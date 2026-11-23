# UCI Raisin Dataset: Çınar, Koklu, and Tasdemir (2020).
# https://doi.org/10.24432/C5660T
# 実行方法: Rscript raisin_tree.R [ローカルのCSVファイル]
# パッケージが入っていない場合: install.packages(c("rpart", "ragg"))
# macOS 以外では、Cairo 対応の R とシステムの日本語フォントを使用する。

args <- commandArgs(trailingOnly = TRUE)
uci_csv <- "https://archive.ics.uci.edu/static/public/850/data.csv"
data_source <- if (length(args) == 0L) uci_csv else args[[1L]]
raisin <- read.csv(data_source)
raisin$Class <- factor(raisin$Class, levels = c("Kecimen", "Besni"))
stopifnot(nrow(raisin) == 900L, !anyNA(raisin), all(table(raisin$Class) == 450L))

# 品種ごとに450粒のうち360粒を学習用、残り90粒をテスト用とする。
set.seed(42, kind = "Mersenne-Twister", sample.kind = "Rejection")
by_class <- split(seq_len(nrow(raisin)), raisin$Class)
train_index <- sort(unlist(lapply(by_class, sample, size = 360L), use.names = FALSE))
train <- raisin[train_index, ]
test <- raisin[-train_index, ]

# Gini 不純度に基づく分類木。図で追える大きさに制限する。
tree <- rpart::rpart(
  Class ~ ., data = train, method = "class", parms = list(split = "gini"),
  control = rpart::rpart.control(
    maxdepth = 4, minsplit = 40, minbucket = 20, cp = 0, xval = 0,
    maxcompete = 0, maxsurrogate = 0, usesurrogate = 0
  )
)

predicted <- predict(tree, newdata = test, type = "class")
print(tree, digits = 9)
print(table(Actual = test$Class, Predicted = predicted))
cat(sprintf("Test accuracy: %d/%d (%.1f%%)\n",
            sum(predicted == test$Class), nrow(test), 100 * mean(predicted == test$Class)))

# 分割条件と各葉の訓練標本の内訳を図にする。
draw_tree <- function(tree, font_family) {
  frame <- tree$frame
  node_id <- as.integer(rownames(frame))
  is_leaf <- frame$var == "<leaf>"
  parents <- which(!is_leaf)
  children <- cbind(match(2L * node_id[parents], node_id),
                    match(2L * node_id[parents] + 1L, node_id))

  # 葉を走査順に右に並べ、内部ノードを左右の子の高さの中間に置く。
  x <- 1.52 + 2.58 * floor(log2(node_id))
  x[is_leaf] <- 12.62
  y <- numeric(nrow(frame))
  y[is_leaf] <- seq(7.34, 1.31, length.out = sum(is_leaf))
  for (i in rev(seq_along(parents))) {
    y[parents[i]] <- mean(y[children[i, ]])
  }

  ink <- "#182B3A"
  muted <- "#476272"
  edge <- "#8BA0AB"
  border <- "#D2DEE3"
  class_color <- c(Kecimen = "#096B73", Besni = "#A34D22")
  class_count <- frame$yval2[, 2:3, drop = FALSE]
  node_class <- attr(tree, "ylevels")[frame$yval]
  widths <- ifelse(is_leaf, 3.22, 2.04)
  half_height <- ifelse(is_leaf, 0.53, 0.49)

  par(mar = rep(0, 4), family = font_family, fg = ink)
  plot.new()
  plot.window(xlim = c(0, 14.75), ylim = c(0, 9.26), xaxs = "i", yaxs = "i")
  text(0.52, 8.79, "乾燥ブドウの品種分類木", adj = 0,
       font = 2, cex = 1.62, col = ink)
  text(0.52, 8.38, sprintf("UCI Raisin データ  ｜  学習用 %d 粒  ｜  Gini 不純度", frame$n[1]),
       adj = 0, cex = 0.86, col = muted)
  segments(0.52, 8.08, 14.23, 8.08, col = border)

  for (i in seq_along(parents)) {
    parent <- parents[i]
    from <- x[parent] + widths[parent] / 2
    elbow <- from + 0.22
    for (side in 1:2) {
      child <- children[i, side]
      to <- x[child] - widths[child] / 2
      lines(c(from, elbow, elbow, to), c(y[parent], y[parent], y[child], y[child]),
            col = edge, lwd = 1.3)
      text(elbow - 0.055, y[child] + 0.14, c("はい", "いいえ")[side],
           adj = 1, cex = 0.70, col = muted)
    }
  }

  feature_labels <- c(
    Area = "面積", MajorAxisLength = "長軸長",
    MinorAxisLength = "短軸長", Eccentricity = "離心率",
    ConvexArea = "凸包面積", Extent = "外接矩形内の面積比", Perimeter = "周囲長"
  )
  for (i in seq_along(parents)) {
    node <- parents[i]
    feature <- frame$var[node]
    threshold <- tree$splits[i, "index"]
    threshold_label <- if (feature %in% c("Eccentricity", "Extent")) {
      sprintf("%.3f", threshold)
    } else {
      sprintf("%.1f", threshold)
    }
    condition <- if (tree$splits[i, "ncat"] < 0) {
      bquote("" < .(threshold_label))
    } else {
      bquote("" >= .(threshold_label))
    }
    rect(x[node] - widths[node] / 2, y[node] - half_height[node],
         x[node] + widths[node] / 2, y[node] + half_height[node],
         col = "#F0F5F7", border = border)
    text(x[node], y[node] + 0.19, feature_labels[[feature]], font = 2, cex = 0.90)
    text(x[node], y[node] - 0.15, condition, cex = 0.98)
  }

  for (node in which(is_leaf)) {
    left <- x[node] - widths[node] / 2
    right <- x[node] + widths[node] / 2
    label <- node_class[node]
    rect(left, y[node] - half_height[node], right, y[node] + half_height[node],
         col = "#FFFFFF", border = border)
    text(left + 0.16, y[node] + 0.29, label, adj = 0, font = 2,
         cex = 0.98, col = class_color[[label]])
    text(right - 0.16, y[node] + 0.29, sprintf("%d 粒", frame$n[node]),
         adj = 1, cex = 0.83, col = muted)
    bar_left <- left + 0.16
    bar_right <- right - 0.16
    bar_split <- bar_left + (bar_right - bar_left) * class_count[node, 1] / frame$n[node]
    rect(bar_left, y[node] - 0.08, bar_split, y[node] + 0.07,
         col = class_color[["Kecimen"]], border = NA)
    rect(bar_split, y[node] - 0.08, bar_right, y[node] + 0.07,
         col = class_color[["Besni"]], border = NA)
    text(bar_left, y[node] - 0.30, paste("Kecimen", class_count[node, 1]),
         adj = 0, cex = 0.87, col = ink)
    text(bar_right, y[node] - 0.30, paste("Besni", class_count[node, 2]),
         adj = 1, cex = 0.87, col = ink)
  }
  text(0.52, 0.38, "各葉の見出しは予測品種、下段の数値は学習データの品種別の粒数。閾値は丸めて表示。",
       adj = 0, cex = 0.76, col = muted)
}

# 日本語フォントを使用し、TeX に取り込めるベクター PDF を作成する。
if (capabilities("aqua")) {
  grDevices::quartzFonts(raisinJapanese = grDevices::quartzFont(
    c("HiraginoSans-W3", "HiraginoSans-W6", "HiraginoSans-W3", "HiraginoSans-W6")
  ))
  pdf_family <- "raisinJapanese"
  png_family <- "Hiragino Sans"
  grDevices::quartz(type = "pdf", file = "raisin_tree.pdf",
                    width = 12.5, height = 7.85, pointsize = 14, bg = "white",
                    family = pdf_family)
} else {
  pdf_family <- png_family <- "sans"
  grDevices::cairo_pdf("raisin_tree.pdf", width = 12.5, height = 7.85,
                      pointsize = 14, bg = "white", family = pdf_family)
}
draw_tree(tree, pdf_family)
invisible(grDevices::dev.off())

ragg::agg_png("raisin_tree.png", width = 12.5, height = 7.85, units = "in",
              res = 180, pointsize = 14, background = "white")
draw_tree(tree, png_family)
invisible(grDevices::dev.off())
