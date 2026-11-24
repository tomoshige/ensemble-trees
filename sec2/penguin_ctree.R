# Rscript penguins_ctree.R [CSVファイルまたはURL] [出力先ディレクトリ]
# 必要なパッケージ: install.packages(c("partykit", "rpart"))
# データ: https://allisonhorst.github.io/palmerpenguins/reference/penguins.html

args <- commandArgs(trailingOnly = TRUE)
csv <- if (length(args)) args[[1L]] else {
  paste0("https://raw.githubusercontent.com/",
         "allisonhorst/palmerpenguins/main/inst/extdata/penguins.csv")
}
output_dir <- if (length(args) >= 2L) args[[2L]] else "."
stopifnot(dir.exists(output_dir))

penguins <- read.csv(csv)
penguins <- na.omit(penguins[c(
  "species", "island", "sex", "bill_length_mm", "bill_depth_mm",
  "flipper_length_mm", "body_mass_g"
)])
categorical <- c("species", "island", "sex")
penguins[categorical] <- lapply(penguins[categorical], factor)

set.seed(42, kind = "Mersenne-Twister", sample.kind = "Rejection")
by_species <- split(seq_len(nrow(penguins)), penguins$species)
train_id <- sort(unlist(lapply(by_species, function(i) {
  sample(i, size = floor(0.8 * length(i)))
}), use.names = FALSE))
train <- penguins[train_id, ]
test <- penguins[-train_id, ]

min_leaf <- 10L
ct <- partykit::ctree(
  species ~ ., data = train,
  control = partykit::ctree_control(
    teststat = "quadratic", splitstat = "quadratic", testtype = "Bonferroni",
    mincriterion = 0.95, minsplit = 2L * min_leaf, minbucket = min_leaf
  )
)
# 学習標本内で、種ごとの羽数が各分割でなるべく均等になるようにする。
set.seed(2026, kind = "Mersenne-Twister", sample.kind = "Rejection")
fold_id <- integer(nrow(train))
for (species in levels(train$species)) {
  i <- which(train$species == species)
  fold_id[i] <- sample(rep(seq_len(10L), length.out = length(i)))
}

grown_cart <- rpart::rpart(
  species ~ ., data = train, method = "class", parms = list(split = "gini"),
  control = rpart::rpart.control(
    cp = 0, minsplit = 2L * min_leaf, minbucket = min_leaf, xval = fold_id
  )
)

# rpart の交差検証で使う各剪定区間の代表値を求める。
cp_table <- grown_cart$cptable
cp_lower <- cp_table[, "CP"]
cv_cp <- c(10 * cp_lower[1L], sqrt(head(cp_lower, -1L) * tail(cp_lower, -1L)))
root_error <- grown_cart$frame$dev[1L] / nrow(train)
cart_cv <- data.frame(
  leaves = as.integer(cp_table[, "nsplit"]) + 1L,
  cp = cv_cp,
  alpha = root_error * cv_cp,
  error = root_error * cp_table[, "xerror"]
)

# 葉数の昇順なので、同率なら最小の木を選ぶ。
selected <- which.min(cart_cv$error)
cart <- rpart::prune(grown_cart, cp = cart_cv$cp[selected])

node_ids <- partykit::nodeids(ct)
nodes <- partykit::nodeapply(ct, node_ids, identity)
children <- lapply(nodes, function(node) {
  match(vapply(partykit::kids_node(node), partykit::id_node, integer(1)), node_ids)
})
leaves <- which(lengths(children) == 0L)
parents <- which(lengths(children) > 0L)
node_depth <- integer(length(nodes))
for (parent in parents) {
  node_depth[children[[parent]]] <- node_depth[parent] + 1L
}
leaf_counts <- table(predict(ct, type = "node"), train$species)

ct_prediction <- predict(ct, newdata = test)
cart_prediction <- predict(cart, newdata = test, type = "class")
results <- data.frame(
  model = c("CTree", "CART"),
  leaves = c(length(leaves), sum(cart$frame$var == "<leaf>")),
  depth = c(max(node_depth), max(floor(log2(as.integer(rownames(cart$frame)))))),
  correct = c(sum(ct_prediction == test$species), sum(cart_prediction == test$species)),
  test_n = nrow(test)
)
results$accuracy <- results$correct / results$test_n

cat(sprintf("Training: %d; test: %d; partykit %s\n",
            nrow(train), nrow(test), utils::packageVersion("partykit")))
print(ct)
cat("\nCART 10-fold cross-validation\n")
print(cart_cv, row.names = FALSE)
cat(sprintf("Selected alpha = %.9g; validation error = %.2f%%\n",
            cart_cv$alpha[selected], 100 * cart_cv$error[selected]))
print(cart)
print(results, row.names = FALSE)
print(table(Observed = test$species, Predicted = ct_prediction))
print(table(Observed = test$species, Predicted = cart_prediction))
print(partykit::info_node(partykit::node_party(ct))$criterion["p.value", ])
print(table(train$island, train$species))

# 根の候補数。同じくちばし長の標本はまとめ、島は各水準と残りに分ける。
bill_left <- cumsum(table(train$bill_length_mm))
island_left <- table(train$island)
admissible <- function(n_left) {
  sum(n_left >= min_leaf & nrow(train) - n_left >= min_leaf)
}
print(c(island = admissible(island_left), bill_length_mm = admissible(bill_left)))

draw_ctree <- function() {
  width <- 8.85
  x <- numeric(length(nodes))
  x[leaves] <- seq(0.76, width - 0.76, length.out = length(leaves))
  for (parent in rev(parents)) {
    x[parent] <- mean(x[children[[parent]]])
  }
  y <- ifelse(seq_along(nodes) %in% leaves, 0.70, 4.40 - 1.03 * node_depth)
  node_half_height <- 0.285
  leaf_half_height <- 0.525
  ink <- "#182B3A"
  edge <- "#8397A3"
  border <- "#C5D3DA"
  class_col <- c(Adelie = "#096B73", Chinstrap = "#A34D22", Gentoo = "#5D4895")
  feature_label <- c(island = "Island", sex = "Sex", bill_length_mm = "Bill length (mm)",
                     bill_depth_mm = "Bill depth (mm)", flipper_length_mm = "Flipper length (mm)",
                     body_mass_g = "Body mass (g)")
  model_data <- partykit::data_party(ct)

  par(mar = rep(0, 4), family = "sans", font = 2, fg = ink)
  plot.new()
  plot.window(c(0, width), c(0, 4.86), xaxs = "i", yaxs = "i")
  for (parent in parents) {
    start <- y[parent] - node_half_height
    elbow <- start - 0.23
    for (child in children[[parent]]) {
      end <- y[child] + if (child %in% leaves) leaf_half_height else node_half_height
      lines(c(x[parent], x[parent], x[child], x[child]), c(start, elbow, elbow, end),
            lwd = 1.35, col = edge)
    }
    branch_x <- (x[parent] + x[children[[parent]]]) / 2
    text(branch_x, elbow + 0.12, c("Yes", "No"), cex = 0.81, col = "#425E70")
  }

  for (parent in parents) {
    split <- partykit::split_node(nodes[[parent]])
    feature <- names(model_data)[partykit::varid_split(split)]
    half_width <- if (feature == "flipper_length_mm") 0.92 else 0.79
    rect(x[parent] - half_width, y[parent] - node_half_height,
         x[parent] + half_width, y[parent] + node_half_height,
         col = "#EFF5F7", border = border)
    text(x[parent], y[parent] + 0.10, feature_label[[feature]], cex = 0.93)
    if (is.factor(model_data[[feature]])) {
      groups <- levels(model_data[[feature]])[partykit::index_split(split) == 1L]
      text(x[parent], y[parent] - 0.10, paste0(paste(groups, collapse = " / "), "?"), cex = 0.97)
    } else {
      cutoff <- format(partykit::breaks_split(split), trim = TRUE, digits = 8L)
      text(x[parent], y[parent] - 0.10, bquote("" <= bold(.(cutoff))), cex = 0.96)
    }
  }

  for (leaf in leaves) {
    counts <- leaf_counts[as.character(node_ids[leaf]), ]
    prediction <- names(counts)[which.max(counts)]
    half_width <- 0.66
    rect(x[leaf] - half_width, y[leaf] - leaf_half_height,
         x[leaf] + half_width, y[leaf] + leaf_half_height,
         col = "white", border = border)
    text(x[leaf], y[leaf] + 0.375, prediction, col = class_col[[prediction]], cex = 1.07)
    bar <- (x[leaf] - 0.53) + 1.06 * c(0, cumsum(counts) / sum(counts))
    rect(head(bar, -1), y[leaf] + 0.163, tail(bar, -1), y[leaf] + 0.235,
         col = class_col[names(counts)], border = NA)
    for (i in seq_along(counts)) {
      line_y <- y[leaf] - 0.185 * (i - 1L) + 0.035
      text(x[leaf] - 0.53, line_y, names(counts)[i], adj = 0, cex = 0.76)
      text(x[leaf] + 0.53, line_y, counts[[i]], adj = 1, cex = 0.78)
    }
  }
}

grDevices::pdf(file.path(output_dir, "penguins_ctree.pdf"), width = 8.85, height = 4.86,
               family = "Helvetica", pointsize = 12, useDingbats = FALSE)
tryCatch(draw_ctree(), finally = grDevices::dev.off())
