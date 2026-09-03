options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_setting <- Sys.getenv("PHASE6B_OUTPUT_DIR", unset = file.path("submission", "anzjs"))
out <- if (grepl("^(?:[A-Za-z]:|/)", out_setting)) out_setting else file.path(root, out_setting)
fig_dir <- file.path(out, "figures")
tab_dir <- file.path(out, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# The journal-owned class and bibliography style are used only in the clean
# compilation staging directory. They are intentionally not copied into the
# author-supplied submission directory.

read_summary <- function(name) {
  path <- file.path(root, "results", "summaries", name)
  stopifnot(file.exists(path))
  read.csv(path, check.names = FALSE)
}

pair <- read_summary("pairwise_full_summary.csv")
pair_diag <- read_summary("pairwise_full_diagnostics.csv")
boot <- read_summary("bootstrap_focused_summary.csv")
matrix <- read_summary("matrix_full_summary.csv")
matrix_r <- read_summary("matrix_correlation_pair_summary.csv")
matrix_p <- read_summary("matrix_partial_pair_summary.csv")
aec_primary <- read_summary("aec_primary_analysis.csv")
aec_exact <- read_summary("aec_exact_id_sensitivity.csv")
aec_m20 <- read_summary("aec_m20_sensitivity.csv")
aec_small <- read_summary("aec_small_place_sensitivity.csv")
aec_booth <- read_summary("aec_state_bootstrap.csv")
aec_cluster19 <- read_summary("aec_division_cluster_bootstrap_2019.csv")
aec_cluster22 <- read_summary("aec_division_cluster_bootstrap_2022.csv")
thin <- read_summary("aec_thinning_summary.csv")
scenarios <- read.csv(file.path(root, "config", "pairwise_scenarios.csv"))

stopifnot(
  nrow(pair) == 96L,
  nrow(matrix) == 3L,
  nrow(thin) == 5L,
  sum(matrix$theorem_violations) == 0L,
  identical(as.integer(matrix$replicates), rep(5000L, 3L)),
  aec_primary$n == 6874L,
  aec_cluster19$valid_replicates == 4999L,
  aec_cluster22$valid_replicates == 4999L,
  all(thin$label == "SEMI-SYNTHETIC CONTROLLED THINNING")
)

if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("patchwork", quietly = TRUE)) {
  stop("ggplot2 and patchwork are required to regenerate manuscript figures")
}
library(ggplot2)
library(patchwork)

theme_anzjs <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.25),
      strip.background = element_rect(fill = "grey92", colour = "grey35"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.4, "cm"),
      plot.title = element_text(face = "bold", size = rel(1)),
      plot.margin = margin(4, 6, 4, 4)
    )
}

method_labels <- c(
  proposed = "Proposed",
  naive = "Naive",
  oracle_sample = "Oracle sample"
)

pair_plot <- pair[pair$method %in% names(method_labels), ]
pair_plot$method_label <- factor(method_labels[pair_plot$method],
                                 levels = unname(method_labels))
pair_plot$scenario_id <- factor(pair_plot$scenario_id, levels = paste0("S", 1:6))
pair_plot$n <- as.numeric(pair_plot$n)

line_scale <- scale_linetype_manual(values = c("solid", "dashed", "dotted"),
                                    drop = FALSE)
shape_scale <- scale_shape_manual(values = c(16, 1, 15), drop = FALSE)
grey_scale <- scale_colour_manual(values = c("black", "grey40", "grey65"),
                                  drop = FALSE)

p_bias <- ggplot(pair_plot, aes(n, bias, colour = method_label,
                                linetype = method_label, shape = method_label)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_errorbar(aes(ymin = bias - 1.96 * mcse_bias,
                    ymax = bias + 1.96 * mcse_bias),
                width = 8, linewidth = 0.3) +
  geom_line(linewidth = 0.45) + geom_point(size = 1.6) +
  facet_wrap(~ scenario_id, nrow = 2) +
  scale_x_continuous(breaks = c(50, 100, 250, 500)) +
  labs(title = "(a) Bias", x = NULL, y = "Monte Carlo bias",
       colour = NULL, linetype = NULL, shape = NULL) +
  line_scale + shape_scale + grey_scale + theme_anzjs()

p_rmse <- ggplot(pair_plot, aes(n, rmse, colour = method_label,
                                linetype = method_label, shape = method_label)) +
  geom_errorbar(aes(ymin = pmax(0, rmse - 1.96 * mcse_rmse),
                    ymax = rmse + 1.96 * mcse_rmse),
                width = 8, linewidth = 0.3) +
  geom_line(linewidth = 0.45) + geom_point(size = 1.6) +
  facet_wrap(~ scenario_id, nrow = 2) +
  scale_x_continuous(breaks = c(50, 100, 250, 500)) +
  labs(title = "(b) Root mean squared error", x = "Number of units, n",
       y = "RMSE", colour = NULL, linetype = NULL, shape = NULL) +
  line_scale + shape_scale + grey_scale + theme_anzjs()

ggsave(file.path(fig_dir, "figure1_pairwise.pdf"),
       p_bias / p_rmse + plot_layout(guides = "collect"),
       width = 7.1, height = 8.6, units = "in", device = cairo_pdf)

coverage <- pair[pair$method == "proposed", ]
coverage$scenario_id <- factor(coverage$scenario_id, levels = paste0("S", 1:6))
coverage$n <- as.numeric(coverage$n)
boot_long <- rbind(
  data.frame(scenario_id = boot$scenario_id, n = boot$n,
             method = "Analytic (focused)", coverage = boot$analytic_coverage,
             mcse = boot$analytic_mcse_coverage),
  data.frame(scenario_id = boot$scenario_id, n = boot$n,
             method = "Percentile bootstrap", coverage = boot$bootstrap_coverage,
             mcse = boot$bootstrap_mcse_coverage)
)
boot_long$scenario_id <- factor(boot_long$scenario_id, levels = paste0("S", 1:6))

p_cov <- ggplot(coverage, aes(n, coverage, group = scenario_id)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey45") +
  geom_errorbar(aes(ymin = coverage - 1.96 * mcse_coverage,
                    ymax = coverage + 1.96 * mcse_coverage),
                width = 8, linewidth = 0.3) +
  geom_line(linewidth = 0.45, colour = "black") +
  geom_point(size = 1.6, colour = "black") +
  geom_errorbar(data = boot_long,
                aes(x = n, ymin = coverage - 1.96 * mcse,
                    ymax = coverage + 1.96 * mcse),
                inherit.aes = FALSE, position = position_dodge(width = 12),
                width = 5, linewidth = 0.3, colour = "grey35") +
  geom_point(data = boot_long, aes(x = n, y = coverage, shape = method),
             inherit.aes = FALSE, position = position_dodge(width = 12),
             size = 2, colour = "grey25") +
  facet_wrap(~ scenario_id, nrow = 2) +
  scale_x_continuous(breaks = c(50, 100, 250, 500)) +
  scale_y_continuous(limits = c(0.84, 0.99), breaks = c(0.85, 0.90, 0.95)) +
  scale_shape_manual(values = c(4, 17)) +
  labs(title = "(a) Analytic 95% interval coverage",
       subtitle = "Focused analytic and bootstrap points are overlaid for S3 and S6 at n = 50, 100",
       x = NULL, y = "Coverage", shape = NULL) +
  theme_anzjs()

reg <- reshape(pair_diag[, c("scenario_id", "n", "projection_rate", "floor_rate")],
               direction = "long", varying = c("projection_rate", "floor_rate"),
               v.names = "rate", timevar = "repair", times = c("Projection", "Floor"))
reg$scenario_id <- factor(reg$scenario_id, levels = paste0("S", 1:6))
reg$n <- as.numeric(reg$n)
p_reg <- ggplot(reg, aes(n, rate, linetype = repair, shape = repair)) +
  geom_line(linewidth = 0.45, colour = "black") + geom_point(size = 1.6) +
  facet_wrap(~ scenario_id, nrow = 2) +
  scale_x_continuous(breaks = c(50, 100, 250, 500)) +
  scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_shape_manual(values = c(16, 1)) +
  labs(title = "(b) Regularisation activity", x = "Number of units, n",
       y = "Activation rate", linetype = NULL, shape = NULL) +
  theme_anzjs()

ggsave(file.path(fig_dir, "figure2_inference.pdf"),
       p_cov / p_reg + plot_layout(guides = "collect"),
       width = 7.1, height = 8.6, units = "in", device = cairo_pdf)

thin_long <- rbind(
  data.frame(cap = thin$cap, method = "Naive", estimate = thin$mean_naive),
  data.frame(cap = thin$cap, method = "Corrected", estimate = thin$mean_corrected)
)
p_thin <- ggplot(thin_long, aes(cap, estimate, linetype = method, shape = method)) +
  geom_hline(yintercept = unique(thin$reference_full_data_Y_correlation),
             colour = "grey45", linetype = "dotdash", linewidth = 0.5) +
  geom_line(colour = "black", linewidth = 0.55) +
  geom_point(colour = "black", size = 2) +
  scale_x_continuous(breaks = thin$cap) +
  scale_y_continuous(limits = c(0.25, 0.94), breaks = seq(0.3, 0.9, 0.1)) +
  scale_linetype_manual(values = c(Corrected = "solid", Naive = "dashed")) +
  scale_shape_manual(values = c(Corrected = 16, Naive = 1)) +
  annotate("text", x = 58, y = unique(thin$reference_full_data_Y_correlation) + 0.025,
           label = "Full-data observed correlation = 0.897704", size = 3) +
  labs(x = "Denominator cap, m", y = "Mean estimated correlation",
       linetype = NULL, shape = NULL) + theme_anzjs(10)
ggsave(file.path(fig_dir, "figure3_thinning.pdf"), p_thin,
       width = 7.1, height = 4.4, units = "in", device = cairo_pdf)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "--", formatC(x, digits = digits, format = "f"))
}
write_tex <- function(name, lines) {
  writeLines(lines, file.path(tab_dir, name), useBytes = TRUE)
}

scenario_desc <- c(
  S1 = "Positive, moderate non-informative",
  S2 = "Positive, low/imbalanced",
  S3 = "Positive, informative/correlated",
  S4 = "Negative, informative/correlated",
  S5 = "Null, informative/correlated",
  S6 = "Positive, boundary-heavy"
)
truth <- aggregate(truth ~ scenario_id, pair[pair$method == "proposed", ], unique)
scenarios2 <- merge(scenarios, truth, by = "scenario_id", sort = FALSE)
t1 <- c("\\begin{tabular}{llccclc}", "\\toprule",
        "Scenario & Latent construction & $a$ & $b$ & $\\delta$ & Denominator & $\\rho$ \\\\",
        "\\midrule")
for (i in seq_len(nrow(scenarios2))) {
  z <- scenarios2[i, ]
  t1 <- c(t1, sprintf("%s & %s & %.1f & %.1f & %.2f & %s & %s \\\\",
                      z$scenario_id, scenario_desc[z$scenario_id], z$a, z$b,
                      z$delta, z$denominator_id, fmt(z$truth, 3)))
}
write_tex("table1_scenarios.tex", c(t1, "\\bottomrule", "\\end{tabular}"))

selected <- pair[pair$n %in% c(50, 500) & pair$method %in% c("proposed", "naive"), ]
wide <- reshape(selected[, c("scenario_id", "n", "method", "bias", "rmse", "coverage")],
                idvar = c("scenario_id", "n"), timevar = "method", direction = "wide")
wide <- wide[order(wide$scenario_id, wide$n), ]
t2 <- c("\\begin{tabular}{lrrrrrr}", "\\toprule",
        "& & \\multicolumn{2}{c}{Proposed} & \\multicolumn{2}{c}{Naive} & \\\\",
        "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
        "Scenario & $n$ & Bias & RMSE & Bias & RMSE & Coverage \\\\", "\\midrule")
for (i in seq_len(nrow(wide))) {
  z <- wide[i, ]
  t2 <- c(t2, sprintf("%s & %d & %s & %s & %s & %s & %s \\\\", z$scenario_id,
                      z$n, fmt(z$bias.proposed), fmt(z$rmse.proposed),
                      fmt(z$bias.naive), fmt(z$rmse.naive), fmt(z$coverage.proposed)))
}
write_tex("table2_pairwise.tex", c(t2, "\\bottomrule", "\\end{tabular}"))

t3 <- c("\\begin{tabular}{rrrrrr}", "\\toprule",
        "$n$ & PSD activation & Raw Frobenius & PSD Frobenius & PSD/raw ratio & Correlation RMSE \\\\",
        "\\midrule")
for (i in seq_len(nrow(matrix))) {
  z <- matrix[i, ]
  t3 <- c(t3, sprintf("%d & %s & %s & %s & %s & %s \\\\", z$n,
                      fmt(z$projection_rate, 4), fmt(z$raw_fro_mean, 5),
                      fmt(z$psd_fro_mean, 5), fmt(z$loss_ratio_mean, 4),
                      fmt(z$R_rmse_offdiag_mean, 4)))
}
write_tex("table3_matrix.tex", c(t3, "\\bottomrule", "\\end{tabular}"))

aec_cols <- Reduce(intersect, lapply(list(aec_primary, aec_exact, aec_m20, aec_small), names))
aec_rows <- do.call(rbind, lapply(list(aec_primary, aec_exact, aec_m20, aec_small),
                                  function(x) x[, aec_cols, drop = FALSE]))
aec_labels <- c("Primary", "Exact-ID only", "$M_{2019},M_{2022}\\geq20$", "Small polling places")
t4 <- c("\\begin{tabular}{lrrrrl}", "\\toprule",
        "Analysis & $n$ & Naive & Latent & Difference & 95\\% interval \\\\", "\\midrule")
for (i in seq_len(nrow(aec_rows))) {
  z <- aec_rows[i, ]
  t4 <- c(t4, sprintf("%s & %d & %s & %s & %s & [%s, %s] \\\\",
                      aec_labels[i], z$n, fmt(z$naive_correlation, 4),
                      fmt(z$latent_correlation, 4), fmt(z$latent_minus_naive, 4),
                      fmt(z$analytic_ci_lower, 4), fmt(z$analytic_ci_upper, 4)))
}
t4 <- c(t4, "\\midrule",
        sprintf("\\multicolumn{2}{l}{Primary, booth bootstrap} & -- & %s & -- & [%s, %s] \\\\",
                fmt(aec_primary$latent_correlation, 4), fmt(aec_booth$ci_lower, 4), fmt(aec_booth$ci_upper, 4)),
        sprintf("\\multicolumn{2}{l}{Primary, 2019-division bootstrap (151 clusters)} & -- & %s & -- & [%s, %s] \\\\",
                fmt(aec_primary$latent_correlation, 4), fmt(aec_cluster19$ci_lower, 4), fmt(aec_cluster19$ci_upper, 4)),
        sprintf("\\multicolumn{2}{l}{Primary, 2022-division bootstrap (151 clusters)} & -- & %s & -- & [%s, %s] \\\\",
                fmt(aec_primary$latent_correlation, 4), fmt(aec_cluster22$ci_lower, 4), fmt(aec_cluster22$ci_upper, 4)))
write_tex("table4_aec.tex", c(t4, "\\bottomrule", "\\end{tabular}"))

method_display <- c(proposed = "Proposed", naive = "Naive", raw = "Raw", oracle_sample = "Oracle")
sup_pair <- pair[order(pair$scenario_id, pair$n, match(pair$method, names(method_display))), ]
sp <- c("\\begin{longtable}{llrrrrrr}",
        "\\caption{Complete pairwise simulation results.}\\label{tab:sup-pair}\\\\", "\\toprule",
        "Scenario & Method & $n$ & Truth & Mean & Bias & RMSE & Coverage \\\\", "\\midrule", "\\endfirsthead",
        "\\toprule", "Scenario & Method & $n$ & Truth & Mean & Bias & RMSE & Coverage \\\\",
        "\\midrule", "\\endhead")
for (i in seq_len(nrow(sup_pair))) {
  z <- sup_pair[i, ]
  sp <- c(sp, sprintf("%s & %s & %d & %s & %s & %s & %s & %s \\\\", z$scenario_id,
                      method_display[z$method], z$n, fmt(z$truth), fmt(z$mean_estimate),
                      fmt(z$bias), fmt(z$rmse), fmt(z$coverage)))
}
write_tex("supp_pairwise_full.tex", c(sp, "\\bottomrule", "\\end{longtable}"))

sd <- c("\\begin{longtable}{lrrrrrrrr}",
        "\\caption{Pairwise diagnostics.}\\label{tab:sup-pair-diag}\\\\", "\\toprule",
        "Scenario & $n$ & Failures & Warnings & Negative variance & Projection & Floor & $\\overline M_1$ & $\\overline M_2$ \\\\",
        "\\midrule", "\\endfirsthead", "\\toprule",
        "Scenario & $n$ & Failures & Warnings & Negative variance & Projection & Floor & $\\overline M_1$ & $\\overline M_2$ \\\\",
        "\\midrule", "\\endhead")
for (i in seq_len(nrow(pair_diag))) {
  z <- pair_diag[i, ]
  sd <- c(sd, sprintf("%s & %d & %d & %d & %s & %s & %s & %s & %s \\\\",
                      z$scenario_id, z$n, z$failures, z$warnings,
                      fmt(z$raw_negative_variance_rate, 4), fmt(z$projection_rate, 4),
                      fmt(z$floor_rate, 4), fmt(z$mean_M1, 2), fmt(z$mean_M2, 2)))
}
write_tex("supp_pairwise_diagnostics.tex", c(sd, "\\bottomrule", "\\end{longtable}"))

sb <- c("\\begin{tabular}{lrrrrrrrr}", "\\toprule",
        "Scenario & $n$ & Analytic cov. & Bootstrap cov. & Difference & Analytic length & Bootstrap length & Boot. projection & Boot. floor \\\\",
        "\\midrule")
for (i in seq_len(nrow(boot))) {
  z <- boot[i, ]
  sb <- c(sb, sprintf("%s & %d & %s & %s & %s & %s & %s & %s & %s \\\\",
                      z$scenario_id, z$n, fmt(z$analytic_coverage), fmt(z$bootstrap_coverage),
                      fmt(z$coverage_difference_bootstrap_minus_analytic), fmt(z$analytic_mean_length),
                      fmt(z$bootstrap_mean_length), fmt(z$bootstrap_projection_rate), fmt(z$bootstrap_floor_rate)))
}
write_tex("supp_bootstrap.tex", c(sb, "\\bottomrule", "\\end{tabular}"))

sm <- c("\\begin{tabular}{rrrrrrrrr}", "\\toprule",
        "$n$ & Projection & Floor & Raw Frob. & PSD Frob. & $R$ RMSE & Partial RMSE & Partial max $>0.5$ & Violations \\\\",
        "\\midrule")
for (i in seq_len(nrow(matrix))) {
  z <- matrix[i, ]
  sm <- c(sm, sprintf("%d & %s & %s & %s & %s & %s & %s & %s & %d \\\\",
                      z$n, fmt(z$projection_rate, 4), fmt(z$floor_rate, 4),
                      fmt(z$raw_fro_mean, 5), fmt(z$psd_fro_mean, 5),
                      fmt(z$R_rmse_offdiag_mean, 4), fmt(z$partial_rmse_offdiag_mean, 4),
                      fmt(z$partial_max_error_gt_0_5, 4), z$theorem_violations))
}
write_tex("supp_matrix_summary.tex", c(sm, "\\bottomrule", "\\end{tabular}"))

pair_table <- function(dat, caption, label) {
  out_lines <- c("\\begin{longtable}{lllrrrrrr}",
                 sprintf("\\caption{%s}\\label{%s}\\\\", caption, label), "\\toprule",
                 "$n$ & Pair & Truth & Mean & Bias & RMSE & Mean abs. error & $|e|>0.5$ & $|e|>0.75$ \\\\",
                 "\\midrule", "\\endfirsthead", "\\toprule",
                 "$n$ & Pair & Truth & Mean & Bias & RMSE & Mean abs. error & $|e|>0.5$ & $|e|>0.75$ \\\\",
                 "\\midrule", "\\endhead")
  for (i in seq_len(nrow(dat))) {
    z <- dat[i, ]
    out_lines <- c(out_lines, sprintf("%d & %s--%s & %s & %s & %s & %s & %s & %s & %s \\\\",
                                            z$n, z$variable_j, z$variable_k, fmt(z$truth),
                                            fmt(z$mean_estimate), fmt(z$bias), fmt(z$rmse),
                                            fmt(z$mean_absolute_error), fmt(z$abs_error_gt_0_5),
                                            fmt(z$abs_error_gt_0_75)))
  }
  c(out_lines, "\\bottomrule", "\\end{longtable}")
}
write_tex("supp_matrix_correlations.tex",
          pair_table(matrix_r, "Pairwise ordinary-correlation performance in the matrix simulation.",
                     "tab:sup-matrix-r"))
write_tex("supp_matrix_partial.tex",
          pair_table(matrix_p, "Partial-correlation performance in the matrix simulation.",
                     "tab:sup-matrix-partial"))

st <- c("\\begin{tabular}{rrrrrrrr}", "\\toprule",
        "Cap & Reference & Mean naive & Naive bias & Naive RMSE & Mean corrected & Corrected bias & Corrected RMSE \\\\",
        "\\midrule")
for (i in seq_len(nrow(thin))) {
  z <- thin[i, ]
  st <- c(st, sprintf("%d & %s & %s & %s & %s & %s & %s & %s \\\\",
                      z$cap, fmt(z$reference_full_data_Y_correlation, 6),
                      fmt(z$mean_naive, 6), fmt(z$naive_bias, 6), fmt(z$naive_rmse, 6),
                      fmt(z$mean_corrected, 6), fmt(z$corrected_bias, 6), fmt(z$corrected_rmse, 6)))
}
write_tex("supp_thinning.tex", c(st, "\\bottomrule", "\\end{tabular}"))

cat("Phase 6B deterministic artifact build complete\n")
cat("Figures:", length(list.files(fig_dir, pattern = "\\.pdf$")), "\n")
cat("Tables:", length(list.files(tab_dir, pattern = "\\.tex$")), "\n")
