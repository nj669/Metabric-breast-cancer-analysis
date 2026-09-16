# =====================================================================
# Breast Cancer Gene Expression Analysis — METABRIC (Kaggle)
# Days 2-5: Load -> Clean -> Statistical Test -> Volcano Plot
# =====================================================================

library(tidyverse)

## ---- Day 2 recap: load the data (skip if `patients` already in memory) ----
patients <- read_csv("METABRIC_RNA_Mutation.csv")


## =====================================================================
## Day 3: Data Cleaning & Pre-processing
## =====================================================================

# --- Step 1: isolate just the gene expression (mRNA z-score) columns ---
# The file mixes clinical fields, gene expression columns, and
# mutation-status columns (which end in "_mut"). We exclude both the known
# numeric clinical fields and the mutation columns, leaving only genes.
# Verified against the real data: this leaves 489 gene columns (not the
# ~331 sometimes quoted for this dataset elsewhere — trust the real data).
clinical_numeric <- c(
  "age_at_diagnosis", "cohort", "neoplasm_histologic_grade",
  "lymph_nodes_examined_positive", "mutation_count",
  "nottingham_prognostic_index", "overall_survival_months",
  "overall_survival", "tumor_size", "tumor_stage",
  "chemotherapy", "hormone_therapy", "radio_therapy"
)

mutation_cols <- names(patients)[str_ends(names(patients), "_mut")]

gene_cols <- patients |>
  select(where(is.numeric)) |>
  select(-any_of(c("patient_id", clinical_numeric))) |>
  select(-any_of(mutation_cols)) |>
  names()

length(gene_cols)   # confirmed: 489

# --- Step 2: eliminate missing values (NAs) ---
expr <- patients |> select(patient_id, er_status, all_of(gene_cols))

na_per_gene <- expr |>
  select(all_of(gene_cols)) |>
  summarise(across(everything(), ~ sum(is.na(.x))))

genes_with_na <- names(na_per_gene)[as.numeric(na_per_gene[1, ]) > 0]
length(genes_with_na)   # confirmed: 0 — this dataset had no missing gene values

expr_clean <- expr |>
  select(-any_of(genes_with_na)) |>   # drop genes with any NA (none, in this case)
  filter(!is.na(er_status))            # drop patients missing our grouping variable (none, in this case)

dim(expr_clean)   # rows = patients kept, columns = genes kept + id + er_status

# --- Step 3: drop low-signal gene profiles (near-zero variance) ---
# Note: this data is already-normalized z-scores (mean ~0, sd ~1 per gene),
# not raw sequencing read counts, so classic "low read-count" filtering
# doesn't apply here. The equivalent QC step for this kind of data is
# dropping genes that barely vary across patients — they carry no
# information for comparing groups either way.
gene_cols_clean <- setdiff(names(expr_clean), c("patient_id", "er_status"))

gene_variance <- expr_clean |>
  select(all_of(gene_cols_clean)) |>
  summarise(across(everything(), ~ var(.x, na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "gene", values_to = "variance")

low_variance_genes <- gene_variance |> filter(variance < 0.01) |> pull(gene)
length(low_variance_genes)   # confirmed: 0 — expected, since z-scored data always has variance = 1

expr_clean <- expr_clean |> select(-any_of(low_variance_genes))
gene_cols_clean <- setdiff(names(expr_clean), c("patient_id", "er_status"))
dim(expr_clean)

# --- Step 4: normalize variance across patients ---
# NOT done as a literal log2 transform. These values are already z-scored
# (mean 0, sd 1 per gene) — a stronger normalization than log2 for this
# purpose — and z-scores are frequently negative, so log2() would return
# NaN for roughly half the data. In the standard pipeline behind this kind
# of dataset, raw microarray intensities are log-transformed FIRST, and
# THEN converted to z-scores — so log2 was already applied upstream, one
# step before the z-scoring visible in this file. Applying it again here
# would double-process the data incorrectly.
#
# Instead, we VERIFY the normalization directly:
gene_means <- expr_clean |> select(all_of(gene_cols_clean)) |> summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
gene_sds   <- expr_clean |> select(all_of(gene_cols_clean)) |> summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

range(as.numeric(gene_means[1, ]))   # should be very close to 0
range(as.numeric(gene_sds[1, ]))     # should be very close to 1


## =====================================================================
## Day 4: Run Statistical Hypothesis Tests
## =====================================================================

# METABRIC has no healthy/normal samples, so we compare two real clinical
# groups instead of "Healthy vs Tumor": ER-positive vs ER-negative tumors.
# Confirmed group sizes: 1459 ER-positive, 445 ER-negative (1904 total).
count(expr_clean, er_status)   # confirm the two group labels and sizes

group_pos <- expr_clean |> filter(er_status == "Positive")
group_neg <- expr_clean |> filter(er_status == "Negative")

run_ttest <- function(gene) {
  test <- t.test(group_pos[[gene]], group_neg[[gene]])
  tibble(
    gene      = gene,
    delta_z   = mean(group_pos[[gene]]) - mean(group_neg[[gene]]),  # ER+ minus ER-
    p_value   = test$p.value
  )
}

results <- map_dfr(gene_cols_clean, run_ttest)

significant_genes <- results |> filter(p_value < 0.05) |> arrange(p_value)

nrow(significant_genes)     # how many genes pass p < 0.05
head(significant_genes, 10) # top 10 — esr1 should be at or near the top


## =====================================================================
## Day 5: Scientific Visualization (Volcano Plot)
## =====================================================================

volcano_data <- results |>
  mutate(
    neg_log10_p = -log10(p_value),
    status = case_when(
      p_value < 0.05 & delta_z >  0.5 ~ "Higher in ER+",
      p_value < 0.05 & delta_z < -0.5 ~ "Higher in ER-",
      TRUE                            ~ "Not significant"
    )
  )

volcano_plot <- ggplot(volcano_data, aes(x = delta_z, y = neg_log10_p, color = status)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(
    "Higher in ER+"   = "#D64550",
    "Higher in ER-"   = "#4472CA",
    "Not significant" = "grey70"
  )) +
  labs(
    title = "Differential Gene Expression: ER-Positive vs ER-Negative Breast Tumors",
    subtitle = "METABRIC dataset (n = 1904 patients, Kaggle)",
    x = "Difference in mean z-score (ER+ minus ER-)",
    y = expression(-log[10]("p-value")),
    color = NULL
  ) +
  theme_minimal(base_size = 13)

volcano_plot

ggsave("volcano_plot.png", volcano_plot, width = 8, height = 6, dpi = 300)
# Check your project folder - volcano_plot.png should now exist there.
