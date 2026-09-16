library(tidyverse)

patients <- read_csv("METABRIC_RNA_Mutation.csv")

class(patients)
dim(patients)
colnames(patients)[1:20]
clinical_numeric <- c("age_at_diagnosis", "cohort", "neoplasm_histologic_grade",
                      "lymph_nodes_examined_positive", "mutation_count",
                      "nottingham_prognostic_index", "overall_survival_months",
                      "overall_survival", "tumor_size", "tumor_stage",
                      "chemotherapy", "hormone_therapy", "radio_therapy")

mutation_cols <- names(patients)[str_ends(names(patients), "_mut")]
gene_cols <- patients |>
  select(where(is.numeric)) |>
  select(-any_of(c("patient_id", clinical_numeric))) |>
  select(-any_of(mutation_cols)) |>
  names()
length(gene_cols)
length(mutation_cols)
sum(clinical_numeric %in% names(patients))
gene_cols[1:20]
tail(gene_cols, 20)
expr <- patients |> select(patient_id, er_status, all_of(gene_cols))
dim(expr)
exists("patients")
exists("gene_cols")
na_per_gene <- expr |>
  select(all_of(gene_cols)) |>
  summarise(across(everything(), ~ sum(is.na(.x))))
genes_with_na <- names(na_per_gene)[as.numeric(na_per_gene[1, ]) > 0]
length(genes_with_na)
expr_clean <- expr |> filter(!is.na(er_status))
dim(expr_clean)
gene_cols_clean <- setdiff(names(expr_clean), c("patient_id", "er_status"))
length(gene_cols_clean)
gene_variance <- expr_clean |>
  select(all_of(gene_cols_clean)) |>
  summarise(across(everything(), ~ var(.x, na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "gene", values_to = "variance")
gene_variance
low_variance_genes <- gene_variance |> filter(variance < 0.01) |> pull(gene)
length(low_variance_genes)
count(expr_clean, er_status)
group_pos <- expr_clean |> filter(er_status == "Positive")
group_neg <- expr_clean |> filter(er_status == "Negative")
nrow(group_pos)
run_ttest <- function(gene) {
  test <- t.test(group_pos[[gene]], group_neg[[gene]])
  tibble(
    gene    = gene,
    delta_z = mean(group_pos[[gene]]) - mean(group_neg[[gene]]),
    p_value = test$p.value
  )
}
results <- map_dfr(gene_cols_clean, run_ttest)
dim(results)
significant_genes <- results |> filter(p_value < 0.05) |> arrange(p_value)

nrow(significant_genes)
head(significant_genes, 10)
volcano_data <- results |>
  mutate(
    neg_log10_p = -log10(p_value),
    status = case_when(
      p_value < 0.05 & delta_z >  0.5 ~ "Higher in ER+",
      p_value < 0.05 & delta_z < -0.5 ~ "Higher in ER-",
      TRUE                            ~ "Not significant"
    )
  )
dim(volcano_data)
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
