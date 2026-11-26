# Rscript crps_gpu_sgemm.R [UCIのZIPファイルまたはURL] [出力先]
# 出典: Paredes and Ballester-Ripoll (2017), https://doi.org/10.24432/C5MK70
# 追加のRパッケージは不要。

args <- commandArgs(trailingOnly = TRUE)
data_source <- if (length(args)) args[[1L]] else
  "https://archive.ics.uci.edu/static/public/440/sgemm+gpu+kernel+performance.zip"
output_dir <- if (length(args) >= 2L) args[[2L]] else "."
stopifnot(dir.exists(output_dir))

if (grepl("^https?://", data_source)) {
  options(timeout = max(120, getOption("timeout")))
  zip <- tempfile(fileext = ".zip")
  download.file(data_source, zip, mode = "wb", quiet = TRUE)
  raw <- read.csv(unz(zip, "sgemm_product.csv"), check.names = FALSE)
  unlink(zip)
} else {
  raw <- read.csv(unz(data_source, "sgemm_product.csv"), check.names = FALSE)
}

runs <- sprintf("Run%d (ms)", seq_len(4L))
features <- setdiff(names(raw), runs)
stopifnot(length(features) == 14L, all(complete.cases(raw)), !anyDuplicated(raw[features]))

# 小さいタイルの候補に対象を限定する。探索時に見た設定を学習側にまとめる。
scope <- which(raw$MWG <= 64 & raw$NWG <= 64)
set.seed(2026L, kind = "Mersenne-Twister", sample.kind = "Rejection")
explored_full <- sample.int(nrow(raw))[seq_len(30000L)]
set.seed(2026L, kind = "Mersenne-Twister", sample.kind = "Rejection")
shuffled_scope <- scope[sample.int(length(scope))]
train_ids <- union(shuffled_scope[seq_len(30000L)], intersect(explored_full, scope))
test_ids <- setdiff(shuffled_scope, train_ids)

to_runs <- function(ids) {
  settings <- raw[ids, ]
  d <- settings[rep(seq_along(ids), each = length(runs)), features]
  d$runtime_ms <- as.vector(t(as.matrix(settings[runs])))
  d
}

# 各設定の4測定を常に同じ側に置く。1設定当たりの測定数は全設定で等しい。
train <- to_runs(train_ids)
test <- to_runs(test_ids)
min_leaf <- 350L * length(runs)
max_depth <- 2L
deadline <- 100

# 実行時間の経験分布を使った最小平均CRPS。
empirical_entropy <- function(y_sorted) {
  n <- length(y_sorted)
  sum((2 * seq_len(n) - n - 1) * y_sorted) / n^2
}

node_loss <- function(y, criterion) {
  if (criterion == "CRPS") length(y) * empirical_entropy(sort(y)) else sum((y - mean(y))^2)
}

grow_tree <- function(data, criterion, depth = 0L, id = 1L, depth_limit = max_depth) {
  y <- data$runtime_ms
  leaf <- list(id = id, y = sort(y))
  n <- nrow(data)
  if (depth >= depth_limit || n < 2L * min_leaf) return(leaf)

  candidates <- lapply(features, function(variable) {
    order_x <- order(data[[variable]])
    x <- data[[variable]][order_x]
    y_ordered <- y[order_x]
    positions <- which(diff(x) > 0)
    positions <- positions[positions >= min_leaf & n - positions >= min_leaf]
    losses <- vapply(positions, function(k) {
      node_loss(y_ordered[seq_len(k)], criterion) +
        node_loss(y_ordered[seq.int(k + 1L, n)], criterion)
    }, numeric(1L))
    data.frame(variable = rep(variable, length(positions)), loss = losses,
               cutoff = (x[positions] + x[positions + 1L]) / 2)
  })
  candidates <- do.call(rbind, candidates)
  if (!nrow(candidates)) return(leaf)
  best <- candidates[which.min(candidates$loss), ]
  if (best$loss >= node_loss(y, criterion)) return(leaf)

  left <- data[[best$variable]] <= best$cutoff
  list(id = id, variable = best$variable, cutoff = best$cutoff,
       left = grow_tree(data[left, ], criterion, depth + 1L, 2L * id, depth_limit),
       right = grow_tree(data[!left, ], criterion, depth + 1L, 2L * id + 1L, depth_limit))
}

terminal_nodes <- function(tree) {
  if (is.null(tree$variable)) return(setNames(list(tree$y), tree$id))
  c(terminal_nodes(tree$left), terminal_nodes(tree$right))
}

predict_leaf <- function(tree, d) {
  if (is.null(tree$variable)) return(rep.int(tree$id, nrow(d)))
  left <- d[[tree$variable]] <= tree$cutoff
  ids <- integer(nrow(d))
  ids[left] <- predict_leaf(tree$left, d[left, ])
  ids[!left] <- predict_leaf(tree$right, d[!left, ])
  ids
}

predictive_crps <- function(y_sorted, actual) {
  n <- length(y_sorted)
  prefix <- c(0, cumsum(y_sorted))
  k <- findInterval(actual, y_sorted)
  ((2 * k - n) * actual + prefix[n + 1L] - 2 * prefix[k + 1L]) / n -
    empirical_entropy(y_sorted)
}

leaf_stats <- function(y) {
  c(settings = length(y) / length(runs), mean = mean(y), median = median(y),
    q90 = quantile(y, 0.90, type = 1L, names = FALSE), p_late = mean(y > deadline))
}

evaluate <- function(tree) {
  leaves <- terminal_nodes(tree)
  ids <- as.character(predict_leaf(tree, test))
  y <- test$runtime_ms
  point <- vapply(leaves, mean, numeric(1L))[ids]
  prob <- vapply(leaves, function(z) mean(z > deadline), numeric(1L))[ids]
  scores <- numeric(length(y))
  for (key in names(leaves)) {
    rows <- ids == key
    scores[rows] <- predictive_crps(leaves[[key]], y[rows])
  }
  c(CRPS = mean(scores), Brier = mean((prob - (y > deadline))^2),
    RMSE = sqrt(mean((point - y)^2)))
}

print_tree <- function(tree, indent = "") {
  if (is.null(tree$variable)) {
    cat(sprintf("%sLeaf %d: %d settings\n", indent, tree$id, length(tree$y) / length(runs)))
    return(invisible(NULL))
  }
  cat(sprintf("%sNode %d: %s <= %g\n", indent, tree$id, tree$variable, tree$cutoff))
  print_tree(tree$left, paste0(indent, "  L "))
  print_tree(tree$right, paste0(indent, "  R "))
}

trees <- list(Empirical = list(id = 1L, y = sort(train$runtime_ms)),
              CART = grow_tree(train, "SSE"), CRPS_tree = grow_tree(train, "CRPS"))
cat(sprintf("Settings: %d total in scope; %d training; %d test\n",
            length(scope), length(train_ids), length(test_ids)))
for (model in names(trees)) {
  cat("\n", model, "\n", sep = "")
  tree <- trees[[model]]
  print_tree(tree)
  cat("Training leaf summaries\n")
  print(t(vapply(terminal_nodes(tree), leaf_stats, numeric(5L))), digits = 5)
  cat("Test leaf summaries\n")
  test_leaves <- split(test$runtime_ms, predict_leaf(tree, test))
  print(t(vapply(test_leaves, leaf_stats, numeric(5L))), digits = 5)
}
cat("\nTest scores\n")
print(t(vapply(trees, evaluate, c(CRPS = 0, Brier = 0, RMSE = 0))), digits = 6)

test_settings <- raw[test_ids, features]
cart_id <- predict_leaf(trees$CART, test_settings)
crps_id <- predict_leaf(trees$CRPS_tree, test_settings)
cart_leaves <- terminal_nodes(trees$CART)
crps_leaves <- terminal_nodes(trees$CRPS_tree)
late_probability <- function(y) mean(y > deadline)
cart_late <- vapply(cart_leaves, late_probability, numeric(1L))
crps_late <- vapply(crps_leaves, late_probability, numeric(1L))
late_gap <- abs(cart_late[as.character(cart_id)] - crps_late[as.character(crps_id)])
cat(sprintf("\nTest settings with at least a 10 percentage point difference: %.1f%%\n",
            100 * mean(late_gap >= 0.10)))
matched_groups <- table(CRPS = crps_id, CART = cart_id)[c("5", "7"), ]
cart_weights <- prop.table(matched_groups, 1L)
group_cart_late <- drop(cart_weights %*% cart_late)
group_crps_late <- crps_late[rownames(matched_groups)]
cat("\nCART leaf counts for the test settings in CRPS leaves 5 and 7\n")
print(matched_groups)
cat("\nMean predicted exceedance probability over the same test settings\n")
print(cbind(CART = group_cart_late, CRPS = group_crps_late))

# テスト設定から、NWG 以外が同じ組を探す。実測時間は選択に使わない。
shared_features <- setdiff(features, "NWG")
pair_candidates <- merge(
  test_settings[crps_id == 5L & test_settings$NWG == 16, shared_features],
  test_settings[crps_id == 7L & test_settings$NWG == 32, shared_features],
  by = shared_features
)
stopifnot(nrow(pair_candidates) > 0L)
first <- do.call(order, pair_candidates[shared_features])[1L]
paired_settings <- pair_candidates[rep(first, 2L), shared_features]
paired_settings$NWG <- c(16, 32)
pair_cart <- as.character(predict_leaf(trees$CART, paired_settings))
pair_crps <- as.character(predict_leaf(trees$CRPS_tree, paired_settings))
cat("\nMatched test settings (NWG only differs)\n")
print(cbind(paired_settings[features], CART_leaf = pair_cart, CRPS_leaf = pair_crps,
            p_late_CART = unname(cart_late[pair_cart]),
            p_late_CRPS = vapply(crps_leaves[pair_crps], late_probability, numeric(1L))),
      row.names = FALSE)

draw_all_leaves <- function() {
  ink <- "#193340"
  edge <- "#718B97"
  guide <- "#DDE8EC"
  reference <- "#657B86"
  threshold_grid <- seq(0, 250, by = 0.25)
  root_survival <- 1 - ecdf(trees$Empirical$y)(threshold_grid)
  leaf_x <- seq(0.198, 0.875, length.out = 4L)
  x <- c(mean(leaf_x), mean(leaf_x[1:2]), mean(leaf_x[3:4]), leaf_x)
  y <- c(0.919, 0.751, 0.751, rep(0.558, 4L))
  node_height <- 0.043

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(gp = grid::gpar(fontfamily = "sans", col = ink)))

  draw_leaf <- function(key, values, color) {
    center <- x[as.integer(key)]
    probability <- late_probability(values)
    n_settings <- length(values) / length(runs)
    grid::grid.text(sprintf("Leaf %s   |   n = %s", key, format(n_settings, big.mark = ",")),
                    x = center, y = 0.565, gp = grid::gpar(fontsize = 11.25, fontface = "bold"))
    grid::grid.text(sprintf(">100 ms   %.1f%%", 100 * probability),
                    x = center, y = 0.497,
                    gp = grid::gpar(col = color, fontsize = 11.25, fontface = "bold"))

    grid::pushViewport(grid::viewport(
      x = center, y = 0.265, width = 0.200, height = 0.367,
      xscale = c(0, 250), yscale = c(0, 1)
    ))
    grid::grid.rect(gp = grid::gpar(fill = "#FAFCFD", col = NA))
    y_ticks <- grid::unit(c(0, 0.5, 1), "native")
    grid::grid.segments(x0 = 0, x1 = 1, y0 = y_ticks, y1 = y_ticks,
                        gp = grid::gpar(col = guide, lwd = 0.8))
    grid::grid.segments(x0 = grid::unit(deadline, "native"),
                        x1 = grid::unit(deadline, "native"), y0 = 0, y1 = 1,
                        gp = grid::gpar(col = "#9AADBA", lty = 3))
    grid::grid.lines(threshold_grid, root_survival, default.units = "native",
                     gp = grid::gpar(col = reference, lwd = 1.1, lty = 2))
    grid::grid.lines(threshold_grid, 1 - ecdf(values)(threshold_grid),
                     default.units = "native", gp = grid::gpar(col = color, lwd = 2))
    grid::grid.points(deadline, probability, default.units = "native",
                      pch = 19, size = grid::unit(1.45, "mm"), gp = grid::gpar(col = color))
    grid::grid.rect(gp = grid::gpar(fill = NA, col = "#AFC4CC", lwd = 0.75))

    x_ticks <- grid::unit(c(0, 100, 200), "native")
    tick_label <- grid::gpar(col = ink, fontsize = 10, fontface = "bold")
    grid::grid.segments(x_ticks, 0, x_ticks, grid::unit(-1.1, "mm"), gp = grid::gpar(col = edge))
    grid::grid.text(c("0", "100", "200"), x = x_ticks,
                    y = grid::unit(-2.9, "mm"), gp = tick_label)
    if (key == "4") {
      grid::grid.segments(0, y_ticks, grid::unit(-1.1, "mm"), y_ticks,
                          gp = grid::gpar(col = edge))
      grid::grid.text(c("0%", "50%", "100%"), x = grid::unit(-1.9, "mm"),
                      y = y_ticks, just = "right", gp = tick_label)
    }
    grid::popViewport()
  }

  draw_tree_row <- function(tree, name, center, color, node_fill) {
    grid::pushViewport(grid::viewport(x = 0.5, y = center, width = 1, height = 0.447))
    grid::grid.text(name, x = 0.098, y = 0.916, just = "left",
                    gp = grid::gpar(col = color, fontsize = 16, fontface = "bold"))

    for (i in 1:3) {
      children <- 2L * i + 0:1
      top <- y[i] - node_height
      joint <- top - 0.036
      stop_at <- y[children] + if (i == 1L) node_height else 0.077
      grid::grid.segments(c(x[i], x[children[1L]]), c(top, joint),
                          c(x[i], x[children[2L]]), c(joint, joint),
                          gp = grid::gpar(col = edge, lwd = 0.95))
      grid::grid.segments(x[children], joint, x[children], stop_at,
                          gp = grid::gpar(col = edge, lwd = 0.95))
      grid::grid.text(c("Yes", "No"), x = x[children] + c(-0.027, 0.027),
                      y = joint - 0.021, gp = grid::gpar(col = ink, fontsize = 9))
    }

    for (node in list(tree, tree$left, tree$right)) {
      i <- node$id
      grid::grid.rect(x[i], y[i], width = 0.137, height = 2 * node_height,
                      gp = grid::gpar(fill = node_fill, col = "#BBCED4", lwd = 0.95))
      grid::grid.text(bquote(bold(.(node$variable)) <= bold(.(node$cutoff))), x = x[i], y = y[i],
                      gp = grid::gpar(fontsize = 11.5))
    }

    leaves <- terminal_nodes(tree)
    for (key in names(leaves)) draw_leaf(key, leaves[[key]], color)
    grid::popViewport()
  }

  draw_tree_row(trees$CART, "CART", 0.748, "#9E4A24", "#F7EEE8")
  draw_tree_row(trees$CRPS_tree, "CRPS", 0.278, "#087266", "#E8F3F0")
  grid::grid.segments(0.097, 0.511, 0.978, 0.511, gp = grid::gpar(col = "#D4E2E7", lwd = 0.9))
  grid::grid.segments(0.732, 0.936, 0.777, 0.936,
                      gp = grid::gpar(col = reference, lwd = 1.3, lty = 2))
  grid::grid.text("Before splitting", x = 0.788, y = 0.936, just = "left",
                  gp = grid::gpar(fontsize = 10.5))
  grid::grid.text("Runtime threshold (ms)", x = 0.537, y = 0.024,
                  gp = grid::gpar(fontsize = 12, fontface = "bold"))
  grid::grid.text("Probability above threshold", x = 0.022, y = 0.49, rot = 90,
                  gp = grid::gpar(fontsize = 12, fontface = "bold"))
  grid::popViewport()
}

pdf(file.path(output_dir, "crps_gpu_sgemm.pdf"), width = 9.6, height = 6.7,
    family = "Helvetica", useDingbats = FALSE)
tryCatch(draw_all_leaves(), finally = invisible(dev.off()))

draw_paired_forecasts <- function() {
  ink <- "#193340"
  colors <- c(CART = "#9E4A24", CRPS = "#087266")
  line_types <- c(CART = 2, CRPS = 1)
  x_range <- c(0, 220)
  grid <- sort(unique(c(seq(x_range[1L], x_range[2L], by = 0.25), deadline)))
  y_ticks <- seq(0, 1, by = 0.25)
  par(mfrow = c(2L, 2L), mar = c(2.1, 3.35, 2.15, 0.9),
      oma = c(2.15, 2.15, 0.35, 0.38), mgp = c(1.65, 0.4, 0), tcl = -0.24,
      family = "sans", font.axis = 2, fg = ink, col.axis = ink, col.main = ink)

  exceedance_curve <- function(y) 1 - ecdf(y)(grid)
  cart_curves <- vapply(cart_leaves, exceedance_curve, numeric(length(grid)))
  crps_curves <- vapply(crps_leaves, exceedance_curve, numeric(length(grid)))

  panel <- function(heading, cart_curve, crps_curve, probabilities, legend_labels) {
    plot(NA_real_, NA_real_, type = "n", xlim = x_range, ylim = c(0, 1),
         xaxs = "i", yaxs = "i", xlab = "", ylab = "", axes = FALSE)
    abline(h = y_ticks, col = "#DFE9EC", lwd = 0.9)
    abline(v = deadline, col = "#83959F", lty = 3, lwd = 1.15)
    axis(1, at = seq(0, 200, by = 50))
    axis(2, at = y_ticks, labels = c("0", "25%", "50%", "75%", "100%"), las = 1)
    box(col = "#ABC0C9", lwd = 0.85)
    title(heading, line = 0.67, cex.main = 1.14)

    lines(grid, cart_curve, type = "s", lwd = 2.3, lty = line_types[["CART"]], col = colors[["CART"]])
    lines(grid, crps_curve, type = "s", lwd = 2.3, lty = line_types[["CRPS"]], col = colors[["CRPS"]])
    points(rep(deadline, 2L), probabilities, pch = 19, cex = 0.74, col = colors)
    legend("topright", legend = sprintf("%s  (%.1f%%)", legend_labels, 100 * probabilities),
           col = colors, lty = line_types, lwd = 2.3, bg = "white", box.col = "#D8E4E7",
           cex = 0.84, inset = 0.04, seg.len = 2.4)
  }

  for (leaf in rownames(matched_groups)) {
    # 同じテスト設定へのCART予測分布を、葉に送られた設定の比率で平均する。
    mixture <- drop(cart_curves %*% cart_weights[leaf, colnames(cart_curves)])
    panel(sprintf("All test settings in CRPS leaf %s", leaf), mixture, crps_curves[, leaf],
          c(group_cart_late[[leaf]], group_crps_late[[leaf]]), c("CART average", "CRPS"))
  }
  for (i in seq_len(nrow(paired_settings))) {
    panel(sprintf("Matched setting  |  NWG = %d", paired_settings$NWG[i]),
          cart_curves[, pair_cart[i]], crps_curves[, pair_crps[i]],
          c(cart_late[[pair_cart[i]]], group_crps_late[[pair_crps[i]]]),
          c(paste("CART leaf", pair_cart[i]), paste("CRPS leaf", pair_crps[i])))
  }
  mtext("Runtime threshold (ms)", side = 1, outer = TRUE, line = 0.53, font = 2, col = ink)
  mtext("Probability above threshold", side = 2, outer = TRUE, line = 0.63, font = 2, col = ink)
}

pdf(file.path(output_dir, "crps_gpu_sgemm_same_settings.pdf"), width = 8.9, height = 6.45,
    family = "Helvetica", useDingbats = FALSE)
tryCatch(draw_paired_forecasts(), finally = invisible(dev.off()))

depth_three <- list(CART = grow_tree(train, "SSE", depth_limit = 3L),
                    CRPS = grow_tree(train, "CRPS", depth_limit = 3L))
cat("\nTest scores for depth 3\n")
print(t(vapply(depth_three, evaluate, c(CRPS = 0, Brier = 0, RMSE = 0))), digits = 6)

draw_depth_three <- function() {
  ink <- "#193340"
  edge <- "#718B97"
  reference <- "#657B86"
  threshold_grid <- seq(0, 300, by = 0.25)
  root_survival <- 1 - ecdf(trees$Empirical$y)(threshold_grid)
  node_x <- c(0.5, 0.25, 0.75, seq(0.125, 0.875, by = 0.25), (seq_len(8L) - 0.5) / 8)
  node_y <- c(0.925, rep(0.885, 2L), rep(0.845, 4L), rep(0.802, 8L))
  box_half_height <- 0.012

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(gp = grid::gpar(fontfamily = "sans", col = ink)))

  draw_leaf <- function(key, values, index, color, highlight) {
    x <- c(0.25, 0.75)[(index - 1L) %% 2L + 1L]
    heading_y <- 0.747 - 0.170 * ((index - 1L) %/% 2L)
    probability <- late_probability(values)
    count <- format(length(values) / length(runs), big.mark = ",")
    grid::grid.text(sprintf("Leaf %s   |   n = %s", key, count), x, heading_y,
                    gp = grid::gpar(fontsize = 11, fontface = "bold"))
    grid::grid.text(sprintf(">100 ms   %.2f%%", 100 * probability), x, heading_y - 0.019,
                    gp = grid::gpar(col = color, fontsize = 10.7, fontface = "bold"))

    grid::pushViewport(grid::viewport(x = x, y = heading_y - 0.087, width = 0.425,
                                      height = 0.116, xscale = c(0, 300), yscale = c(0, 1)))
    grid::grid.rect(gp = grid::gpar(fill = if (highlight) "#EDF7F3" else "#FAFCFD", col = NA))
    y_ticks <- grid::unit(c(0, 0.5, 1), "native")
    grid::grid.segments(x0 = 0, x1 = 1, y0 = y_ticks, y1 = y_ticks,
                        gp = grid::gpar(col = "#DDE8EC", lwd = 0.8))
    grid::grid.segments(x0 = grid::unit(deadline, "native"),
                        x1 = grid::unit(deadline, "native"), y0 = 0, y1 = 1,
                        gp = grid::gpar(col = "#9AADBA", lty = 3))
    grid::grid.lines(threshold_grid, root_survival, default.units = "native",
                     gp = grid::gpar(col = reference, lwd = 1.1, lty = 2))
    grid::grid.lines(threshold_grid, 1 - ecdf(values)(threshold_grid), default.units = "native",
                     gp = grid::gpar(col = color, lwd = 1.95))
    grid::grid.points(deadline, probability, default.units = "native", pch = 19,
                      size = grid::unit(1.35, "mm"), gp = grid::gpar(col = color))
    grid::grid.rect(gp = grid::gpar(fill = NA, col = if (highlight) color else "#AFC4CC", lwd = 0.8))

    x_ticks <- grid::unit(seq(0, 300, by = 100), "native")
    tick_label <- grid::gpar(fontsize = 9.4, fontface = "bold")
    grid::grid.segments(x_ticks, 0, x_ticks, grid::unit(-1, "mm"), gp = grid::gpar(col = edge))
    grid::grid.text(c("0", "100", "200", "300"), x_ticks, grid::unit(-2.65, "mm"), gp = tick_label)
    if (index %% 2L == 1L) {
      grid::grid.segments(0, y_ticks, grid::unit(-1, "mm"), y_ticks, gp = grid::gpar(col = edge))
      grid::grid.text(c("0%", "50%", "100%"), grid::unit(-1.75, "mm"), y_ticks,
                      just = "right", gp = tick_label)
    }
    grid::popViewport()
  }

  draw_model <- function(tree, name, x_center, color, tint, highlighted = integer()) {
    grid::pushViewport(grid::viewport(x = x_center, width = 0.418))
    grid::grid.text(name, x = 0.5, y = 0.967,
                    gp = grid::gpar(col = color, fontsize = 17, fontface = "bold"))

    for (i in 1:7) {
      kids <- 2L * i + 0:1
      start <- node_y[i] - box_half_height
      joint <- start - 0.009
      end <- node_y[kids] + if (i <= 3L) box_half_height else 0.010
      grid::grid.segments(c(node_x[i], node_x[kids[1L]]), c(start, joint),
                          c(node_x[i], node_x[kids[2L]]), c(joint, joint),
                          gp = grid::gpar(col = edge, lwd = 0.8))
      grid::grid.segments(node_x[kids], joint, node_x[kids], end,
                          gp = grid::gpar(col = edge, lwd = 0.8))
    }

    branches <- list(tree, tree$left, tree$right,
                     tree$left$left, tree$left$right, tree$right$left, tree$right$right)
    for (node in branches) {
      i <- node$id
      emphasized <- all((2L * i + 0:1) %in% highlighted)
      grid::grid.rect(node_x[i], node_y[i], width = 0.228, height = 2 * box_half_height,
                      gp = grid::gpar(fill = tint, col = if (emphasized) color else "#BBCED4"))
      grid::grid.text(bquote(bold(.(node$variable)) <= .(node$cutoff)), node_x[i], node_y[i],
                      gp = grid::gpar(fontsize = 9.6))
    }
    leaves <- terminal_nodes(tree)
    for (i in seq_along(leaves)) {
      key <- names(leaves)[i]
      highlight <- as.integer(key) %in% highlighted
      grid::grid.text(key, node_x[as.integer(key)], node_y[as.integer(key)],
                      gp = grid::gpar(col = if (highlight) color else ink,
                                      fontsize = 10.3, fontface = "bold"))
      draw_leaf(key, leaves[[i]], i, color, highlight)
    }
    grid::popViewport()
  }

  draw_model(depth_three$CART, "CART", 0.277, "#9E4A24", "#F7EEE8")
  draw_model(depth_three$CRPS, "CRPS", 0.755, "#087266", "#E8F3F0", c(12L, 13L))
  grid::grid.segments(0.507, 0.115, 0.507, 0.970, gp = grid::gpar(col = "#D4E2E7", lwd = 0.8))
  grid::grid.text("Probability above threshold", x = 0.020, y = 0.405, rot = 90,
                  gp = grid::gpar(fontsize = 12, fontface = "bold"))
  grid::grid.text("Runtime threshold (ms)", x = 0.516, y = 0.055,
                  gp = grid::gpar(fontsize = 12, fontface = "bold"))
  grid::grid.text("Left branch satisfies the rule", x = 0.28, y = 0.784,
                  gp = grid::gpar(fontsize = 9.4, col = reference))
  grid::grid.text("Left branch satisfies the rule", x = 0.755, y = 0.784,
                  gp = grid::gpar(fontsize = 9.4, col = reference))
  grid::grid.segments(0.399, 0.027, 0.442, 0.027,
                      gp = grid::gpar(col = reference, lwd = 1.3, lty = 2))
  grid::grid.text("Before splitting", x = 0.453, y = 0.027, just = "left",
                  gp = grid::gpar(fontsize = 10.2))
  grid::popViewport()
}

pdf(file.path(output_dir, "crps_gpu_sgemm_depth3.pdf"), width = 10.5, height = 10.4,
    family = "Helvetica", useDingbats = FALSE)
tryCatch(draw_depth_three(), finally = invisible(dev.off()))
