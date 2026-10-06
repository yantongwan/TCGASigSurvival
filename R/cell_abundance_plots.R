#' Forest plot comparing marker and estimated-fraction Cox effects
#' @param result Survival result list or results data frame.
#' @param cohort Cohort label.
#' @param endpoint Clinical endpoint.
#' @param adjustment Adjustment-set label.
#' @param model continuous or median.
#' @return ggplot with successful primary-cell estimates; QC retains failures.
#' @export
plot_cell_abundance_forest <- function(result, cohort = "PAAD", endpoint = "OS",
    adjustment = "unadjusted", model = "continuous") {
  d <- if (is.list(result) && !is.data.frame(result)) result$results else result
  d <- d[d$cohort_label == cohort & d$endpoint == endpoint & d$adjustment == adjustment &
    d$model == model & d$role == "primary", , drop = FALSE]
  d$cell_type <- factor(d$cell_type, levels = rev(names(default_immune_signatures())))
  d$measurement <- factor(d$measurement, levels = c("marker_score", "estimated_fraction"),
    labels = c("Marker score", "Estimated fraction"))
  failed <- sum(d$status != "ok"); d <- d[d$status == "ok", , drop = FALSE]
  dodge <- ggplot2::position_dodge(width = 0.5)
  ggplot2::ggplot(d, ggplot2::aes(x = HR, y = cell_type, colour = measurement, group = measurement)) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey65", linetype = 2, linewidth = 0.4) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = CI_lower, xmax = CI_upper),
      orientation = "y", width = 0.16, linewidth = 0.55, position = dodge) +
    ggplot2::geom_point(size = 2.2, position = dodge) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c("Marker score" = "#50789B", "Estimated fraction" = "#B56357"), drop = FALSE) +
    ggplot2::labs(x = if (model == "continuous") "Hazard ratio per reference-cohort SD (95% CI)" else "Hazard ratio: high versus low (95% CI)",
      y = NULL, colour = NULL, title = paste(cohort, endpoint, "immune-context association"),
      subtitle = paste("Adjustment:", adjustment, "| Not estimable:", failed),
      caption = "Marker scores are not cell percentages. CD4 fraction includes Treg; macrophage fraction is model M1 + M2.") +
    ggplot2::theme_classic(base_size = 10, base_family = "Arial") +
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(size = 12),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0))
}

#' Kaplan-Meier plot of fixed marker/fraction median groups
#' @param result Survival result list with inputs and results.
#' @param cell_type Cell label.
#' @param measurement marker_score or estimated_fraction.
#' @param cohort,endpoint Cohort and clinical endpoint.
#' @return A single-panel ggplot; at-risk counts are attached as risk_table.
#'   Groups and cutoff are fixed before endpoint exclusions, not optimized.
#' @export
plot_cell_abundance_km <- function(result, cell_type, measurement,
    cohort = "PAAD", endpoint = "OS") {
  d <- result$inputs
  d <- d[d$cohort_label == cohort & d$cell_type == cell_type & d$measurement == measurement &
    d$endpoint == endpoint & d$adjustment == "unadjusted" & d$exclusion == "included", , drop = FALSE]
  r <- result$results
  r <- r[r$cohort_label == cohort & r$cell_type == cell_type & r$measurement == measurement &
    r$endpoint == endpoint & r$adjustment == "unadjusted" & r$model == "median", , drop = FALSE]
  if (nrow(r) != 1L || r$status != "ok") stop("Requested median survival contrast is not estimable; consult source QC.", call. = FALSE)
  fit <- survival::survfit(survival::Surv(endpoint_time, endpoint_status) ~ group, data = d)
  s <- summary(fit, censored = TRUE)
  curve <- data.frame(months = s$time / (365.25 / 12), survival_estimate = s$surv,
    lower_survival = s$lower, upper_survival = s$upper, n_censor = s$n.censor,
    curve = sub("^group=", "", as.character(s$strata)))
  curve <- rbind(data.frame(months = 0, survival_estimate = 1, lower_survival = 1,
    upper_survival = 1, n_censor = 0, curve = c("low", "high")), curve)
  curve <- curve[order(curve$curve, curve$months), ]
  times <- seq(0, floor(max(d$endpoint_time) / (365.25 / 12) / 12) * 12, by = 12)
  risk <- do.call(rbind, lapply(c("low", "high"), function(g) data.frame(group = g,
    months = times, n_risk = vapply(times, function(t) sum(d$group == g & d$endpoint_time >= t * (365.25 / 12)), integer(1)))))
  labels <- vapply(c("low", "high"), function(g) paste0(if (g == "low") "Low" else "High",
    " (n=", sum(d$group == g), ", events=", sum(d$group == g & d$endpoint_status == 1), ")"), character(1))
  title <- paste(cohort, cell_type, endpoint)
  subtitle <- paste(if (measurement == "marker_score") "Marker score" else "Estimated cell fraction",
    "| Median cutoff =", signif(r$cutoff, 3), "| Log-rank P =", format.pval(r$logrank_p, digits = 3))
  band <- do.call(rbind, lapply(split(curve, curve$curve), function(z) {
    if (nrow(z) < 2L) return(z)
    old <- z[-nrow(z), , drop = FALSE]; old$months <- z$months[-1L]
    both <- rbind(old, z); both$order <- c(seq_len(nrow(old)) * 2 - 1, seq_len(nrow(z)) * 2 - 2)
    both[order(both$months, both$order), ]
  }))
  p <- ggplot2::ggplot(curve, ggplot2::aes(x = months, y = survival_estimate, colour = curve, fill = curve)) +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(ymin = lower_survival, ymax = upper_survival), alpha = 0.12, colour = NA, na.rm = TRUE) +
    ggplot2::geom_step(linewidth = 0.65) +
    ggplot2::geom_point(data = curve[curve$n_censor > 0, ], shape = 3, size = 1.7, stroke = 0.35) +
    ggplot2::scale_colour_manual(values = c(low = "#50789B", high = "#B56357"), breaks = c("low", "high"), labels = labels) +
    ggplot2::scale_fill_manual(values = c(low = "#50789B", high = "#B56357"), breaks = c("low", "high"), labels = labels) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(x = "Time (months)", y = paste(endpoint, "survival probability"), colour = NULL, fill = NULL,
      title = title, subtitle = subtitle, caption = paste0("95% survival intervals. High > median; ties are low.\n",
        "Descriptive median split; primary inference uses continuous Cox. At-risk counts are in the paired source table.")) +
    ggplot2::theme_classic(base_size = 10, base_family = "Arial") +
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(size = 12),
      plot.subtitle = ggplot2::element_text(size = 9), plot.caption = ggplot2::element_text(size = 8, hjust = 0))
  attr(p, "risk_table") <- risk
  p
}
