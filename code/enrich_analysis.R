library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

result_fash <- enrich_set(genes_selected = genes_highlighted1, background_gene = all_genes)
result_fash_unadj <- enrich_set(genes_selected = genes_highlighted1_before, background_gene = all_genes)
result_strober_linear <- enrich_set(genes_selected = genes_highlighted_strober_linear, background_gene = all_genes)
result_strober_nonlinear <- enrich_set(genes_selected = genes_highlighted_strober_nonlinear, background_gene = all_genes)

result_fash_early <- enrich_set(genes_selected = early_genes, background_gene = all_genes)
result_fash_late <- enrich_set(genes_selected = late_genes, background_gene = all_genes)
result_fash_switch <- enrich_set(genes_selected = switch_genes, background_gene = all_genes)
result_fash_mid <- enrich_set(genes_selected = middle_genes, background_gene = all_genes)

# Extract (Term, PValue) from a clusterProfiler enrichResult
extract_pvals <- function(enrich_res, p_col = c("pvalue", "p.adjust", "qvalue")) {
  p_col <- match.arg(p_col)
  enrich_res@result %>%
    transmute(
      Term   = gsub("^HALLMARK_", "", Description),
      PValue = .data[[p_col]]
    )
}

# Avoid -log10(0) = Inf
safe_neglog10 <- function(p) {
  p <- pmin(pmax(p, .Machine$double.xmin), 1)
  -log10(p)
}

make_pval_scatter <- function(df_a, df_b, name_a, name_b,
                              alpha = 0.6, label_n = 50,
                              sig_cut = 0.05, clip_max = 5) {

  cmp <- full_join(df_a, df_b, by = "Term", suffix = c("_a", "_b")) %>%
    mutate(
      PValue_a = coalesce(PValue_a, 1),
      PValue_b = coalesce(PValue_b, 1),
      negLog_a = safe_neglog10(PValue_a),
      negLog_b = safe_neglog10(PValue_b),
      sig = case_when(
        PValue_a < sig_cut & PValue_b < sig_cut ~ "Both",
        PValue_a < sig_cut & PValue_b >= sig_cut ~ paste0(name_a, " only"),
        PValue_a >= sig_cut & PValue_b < sig_cut ~ paste0(name_b, " only"),
        TRUE ~ "Neither"
      )
    )

  # Label only the strongest (smallest) significant p-values
  to_label <- cmp %>%
    filter(PValue_a < sig_cut | PValue_b < sig_cut) %>%
    arrange(pmin(PValue_a, PValue_b)) %>%
    slice_head(n = label_n)

  ggplot(cmp, aes(x = negLog_a, y = negLog_b)) +
    geom_point(aes(color = sig), alpha = alpha) +
    # geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(sig_cut), linetype = "dashed") +
    geom_vline(xintercept = -log10(sig_cut), linetype = "dashed") +
    coord_equal(xlim = c(0, clip_max), ylim = c(0, clip_max), expand = FALSE) +
    geom_text_repel(
      data = to_label, aes(label = Term),
      size = 3, max.overlaps = Inf
    ) +
    labs(
      title = paste0("GSEA p-value comparison: ", name_a, " vs ", name_b),
      x = paste0("-log10(p) ", name_a),
      y = paste0("-log10(p) ", name_b),
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

# --- Your data frames ---
p_val_df_fash <- extract_pvals(result_fash, p_col = "pvalue")
p_val_df_fash_unadj <- extract_pvals(result_fash_unadj, p_col = "pvalue")
p_val_df_strober_linear <- extract_pvals(result_strober_linear, p_col = "pvalue")
p_val_df_strober_nonlinear <- extract_pvals(result_strober_nonlinear, p_col = "pvalue")

# adjusted fash vs strober
pval_plot_fash_strober_linear <-
  make_pval_scatter(p_val_df_fash, p_val_df_strober_linear, "FASH", "Strober linear")
pval_plot_fash_strober_nonlinear <-
  make_pval_scatter(p_val_df_fash, p_val_df_strober_nonlinear, "FASH", "Strober nonlinear")
pval_plot_fash_strober_linear
pval_plot_fash_strober_nonlinear

# unadjusted fash vs strober
pval_plot_fash_unadj_strober_linear <-
  make_pval_scatter(p_val_df_fash_unadj, p_val_df_strober_linear, "FASH unadj", "Strober linear")
pval_plot_fash_unadj_strober_nonlinear <-
  make_pval_scatter(p_val_df_fash_unadj, p_val_df_strober_nonlinear, "FASH unadj", "Strober nonlinear")
pval_plot_fash_unadj_strober_linear
pval_plot_fash_unadj_strober_nonlinear

# FASH early, late, switch, mid vs strober
pval_plot_fash_early_strober_linear <-
  make_pval_scatter(extract_pvals(result_fash_early, p_col = "pvalue"),
                    p_val_df_strober_linear, "FASH early", "Strober linear")
pval_plot_fash_late_strober_linear <-

