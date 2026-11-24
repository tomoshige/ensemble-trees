# Rscript penguins_mob.R [出力先ディレクトリ]
# install.packages("partykit")
# datasets::penguins がない R では install.packages("palmerpenguins") も実行。

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1L]] else "."
stopifnot(dir.exists(output_dir))

penguins <- tryCatch(datasets::penguins, error = function(e) NULL)
if (is.null(penguins)) {
  penguins <- palmerpenguins::penguins[c(
    "body_mass_g", "flipper_length_mm", "species", "sex", "bill_length_mm"
  )]
  names(penguins) <- c("body_mass", "flipper_len", "species", "sex", "bill_len")
}

d <- na.omit(penguins[c("body_mass", "flipper_len", "species", "sex", "bill_len")])
d$flipper_200 <- d$flipper_len - 200

# | の左は各ノードで当てはめるモデル、右は分割変数の候補。
fit <- partykit::lmtree(
  body_mass ~ flipper_200 | species + sex + bill_len,
  data = d, alpha = 0.05, bonferroni = TRUE, trim = 0.10, minsize = 25L
)

cat("Complete observations:", nrow(d), "\n")
print(fit)
cat("\nParameter instability tests at the root\n")
print(strucchange::sctest(fit, node = 1L))

leaf_data <- split(d, predict(fit, type = "node"))
leaf_coef <- coef(fit)
flipper_ranges <- t(vapply(leaf_data, function(group) range(group$flipper_len), numeric(2L)))
colnames(flipper_ranges) <- c("min_mm", "max_mm")
cat("\nCoefficients and observed flipper-length ranges by leaf\n")
print(cbind(n = vapply(leaf_data, nrow, integer(1L)), leaf_coef, flipper_ranges))

draw_tree <- function() {
  ink <- "#1C333D"
  axis <- grid::gpar(col = ink, fontsize = 10, fontface = "bold")
  x_limits <- range(d$flipper_len) + c(-3, 3)
  y_limits <- range(d$body_mass) + c(-150, 150)
  x_ticks <- pretty(x_limits, 3)
  x_ticks <- x_ticks[x_ticks > x_limits[1L] & x_ticks < x_limits[2L]]
  y_ticks <- pretty(y_limits, 4)
  y_ticks <- y_ticks[y_ticks > y_limits[1L] & y_ticks < y_limits[2L]]

  draw_leaf <- function(node) {
    key <- as.character(partykit::id_node(node))
    obs <- leaf_data[[key]]
    beta <- leaf_coef[key, ]
    grid::pushViewport(grid::viewport(
      x = 0.53, y = 0.47, width = 0.87, height = 0.73,
      xscale = x_limits, yscale = y_limits
    ))
    grid::grid.rect(gp = grid::gpar(fill = "#F1F6F7", col = NA))
    grid::grid.segments(
      x0 = 0, x1 = 1, y0 = grid::unit(y_ticks, "native"),
      y1 = grid::unit(y_ticks, "native"), gp = grid::gpar(col = "#D4E2E6")
    )
    grid::grid.points(
      obs$flipper_len, obs$body_mass, default.units = "native", pch = 16,
      size = grid::unit(1.04, "mm"),
      gp = grid::gpar(col = grDevices::adjustcolor("#3F7A91", alpha.f = 0.57))
    )
    observed_range <- range(obs$flipper_len)
    grid::grid.lines(
      observed_range, beta[1L] + beta[2L] * (observed_range - 200),
      default.units = "native", gp = grid::gpar(col = "#005E54", lwd = 2.5)
    )
    grid::grid.rect(gp = grid::gpar(fill = NA, col = "#B0C7CE"))
    grid::grid.xaxis(at = x_ticks, gp = axis)
    if (key == names(leaf_data)[1L]) grid::grid.yaxis(at = y_ticks, gp = axis)
    grid::grid.text(
      paste("n =", nrow(obs)), x = 0.05, y = 0.95, just = c("left", "top"), gp = axis
    )
    grid::popViewport()
  }

  plot(
    fit, terminal_panel = draw_leaf, tnex = 2.05, margins = c(2.5, 4, 0.8, 0.85),
    ip_args = list(id = FALSE, pval = FALSE, fill = "#E7F0F2"),
    ep_args = list(justmin = 9, just = "alternate"),
    gp = grid::gpar(fontface = "bold", fontsize = 12, col = ink)
  )
  grid::grid.text("Flipper length (mm)", y = grid::unit(0.7, "lines"), gp = axis)
  grid::grid.text(
    "Body mass (g)", x = grid::unit(0.9, "lines"), y = 0.25, rot = 90, gp = axis
  )
}

grDevices::pdf(file.path(output_dir, "penguins_mob.pdf"), width = 11, height = 5.65,
               family = "Helvetica", useDingbats = FALSE)
tryCatch(draw_tree(), finally = invisible(grDevices::dev.off()))
