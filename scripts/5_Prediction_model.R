library(svd)
library(nnet)
library(dplyr)
library(tidyr)

set.seed(999)

# ---- Load and align data ----
microbiome_data <- read.csv("bacterial_abundance.csv", row.names = 1)  # with your own data
virome_data     <- read.csv("viral_abundance.csv",    row.names = 1)  # with your own data
microbiome_data <- t(microbiome_data)
virome_data     <- t(virome_data)
common_samples  <- intersect(rownames(microbiome_data), rownames(virome_data))
microbiome_data <- microbiome_data[common_samples, , drop = FALSE]
virome_data     <- virome_data[common_samples, , drop = FALSE]

# ---- Train/test split ----
n <- nrow(microbiome_data)
train_indices <- sample(seq_len(n), size = floor(0.8 * n))
test_indices  <- setdiff(seq_len(n), train_indices)

# ---- Center/scale using training set only ----
center_scale <- function(X, idx_train) {
  mu <- colMeans(X[idx_train, , drop = FALSE])
  sd <- apply(X[idx_train, , drop = FALSE], 2, sd)
  sd[sd == 0] <- 1
  list(
    Xs_train = sweep(sweep(X[idx_train, , drop = FALSE], 2, mu, "-"), 2, sd, "/"),
    Xs_test  = sweep(sweep(X[-idx_train, , drop = FALSE], 2, mu, "-"), 2, sd, "/"),
    mu = mu, sd = sd
  )
}

mic_cs <- center_scale(microbiome_data, train_indices)
vir_cs <- center_scale(virome_data,     train_indices)

# ---- SVD on TRAIN ONLY; then project TEST consistently ----
svd_train <- function(X) {
  # X is centered/scaled train matrix (samples x features)
  s <- svd(as.matrix(X))
  list(u = s$u, d = s$d, v = s$v)
}

# microbiome
mic_svd  <- svd_train(mic_cs$Xs_train)
mic_ev   <- mic_svd$d^2 / sum(mic_svd$d^2)
mic_cum  <- cumsum(mic_ev)
mic_k    <- which(mic_cum >= 0.80)[1]
# scores on train: U * D
mic_scores_train <- mic_svd$u[, 1:mic_k, drop = FALSE] %*% diag(mic_svd$d[1:mic_k])
colnames(mic_scores_train) <- paste0("microbiome_dim", seq_len(mic_k))
# project test onto V: X_test %*% V
mic_scores_test <- as.matrix(mic_cs$Xs_test) %*% mic_svd$v[, 1:mic_k, drop = FALSE]
mic_scores_test <- mic_scores_test %*% diag(mic_svd$d[1:mic_k])
colnames(mic_scores_test) <- colnames(mic_scores_train)

# virome
vir_svd  <- svd_train(vir_cs$Xs_train)
vir_ev   <- vir_svd$d^2 / sum(vir_svd$d^2)
vir_cum  <- cumsum(vir_ev)
vir_k    <- which(vir_cum >= 0.80)[1]
vir_scores_train <- vir_svd$u[, 1:vir_k, drop = FALSE] %*% diag(vir_svd$d[1:vir_k])
colnames(vir_scores_train) <- paste0("viral_dim", seq_len(vir_k))
vir_scores_test  <- as.matrix(vir_cs$Xs_test) %*% vir_svd$v[, 1:vir_k, drop = FALSE]
vir_scores_test  <- vir_scores_test %*% diag(vir_svd$d[1:vir_k])
colnames(vir_scores_test) <- colnames(vir_scores_train)

# ---- Predict microbiome SVD from virome SVD (per-component nnet) ----
run_neural_network <- function(train_x, test_x, y_vec, hidden_units) {
  fit <- nnet(
    x = as.matrix(train_x),
    y = as.numeric(y_vec),
    size = hidden_units,
    linout = TRUE,
    trace = FALSE,
    maxit = 1000,
    decay = 1e-4
  )
  as.numeric(predict(fit, as.matrix(test_x)))
}

calc_err <- function(actual, predicted) {
  mae <- mean(abs(actual - predicted))
  mse <- mean((actual - predicted)^2)
  rmse <- sqrt(mse)
  r2 <- 1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
  c(R2 = r2, MAE = mae, RMSE = rmse)
}

hidden_units_grid <- 1:33
variance_weights  <- mic_ev[1:mic_k]
variance_weights  <- variance_weights / sum(variance_weights)  # <-- normalize

results_list <- vector("list", length(hidden_units_grid))
names(results_list) <- hidden_units_grid

for (hu in hidden_units_grid) {
  model_r2 <- model_mae <- model_rmse <- numeric(mic_k)
  for (i in seq_len(mic_k)) {
    y_tr <- mic_scores_train[, i]
    y_te <- mic_scores_test[,  i]
    yhat <- run_neural_network(vir_scores_train, vir_scores_test, y_tr, hu)
    e    <- calc_err(y_te, yhat)
    model_r2[i]  <- e["R2"]; model_mae[i] <- e["MAE"]; model_rmse[i] <- e["RMSE"]
  }
  results_list[[as.character(hu)]] <- data.frame(
    Hidden_Units = hu,
    Signal       = paste0("Signal_", seq_len(mic_k)),
    R2 = model_r2, MAE = model_mae, RMSE = model_rmse, Weight = variance_weights
  )
}

all_results <- bind_rows(results_list)
weighted_summary <- all_results |>
  group_by(Hidden_Units) |>
  summarise(
    Weighted_R2   = sum(R2  * Weight),
    Weighted_MAE  = sum(MAE * Weight),
    Weighted_RMSE = sum(RMSE* Weight),
    SD_R2 = sd(R2), SD_MAE = sd(MAE), SD_RMSE = sd(RMSE),
    .groups = "drop"
  )
head(weighted_summary)  

z_score_model <- weighted_summary %>%
  mutate(
    z_r2 = scale(Weighted_R2),
    z_mae = scale(-Weighted_MAE),
    z_rmse = scale(-Weighted_RMSE),
    z_total = 0.33 * z_r2 + 0.33 * z_mae + 0.33 * z_rmse
  ) %>%
  arrange(desc(z_total))

best_model_units <- z_score_model$Hidden_Units[1]
message("Best model based on z-score uses ", best_model_units, " hidden units.")

# ---- Train final per-signal models using best units ----
neuralnet_predictions <- vector("list", mic_k)
for (i in seq_len(mic_k)) {
  y_tr <- mic_scores_train[, i]
  neuralnet_predictions[[i]] <- run_neural_network(
    vir_scores_train, vir_scores_test, y_tr, hidden_units = best_model_units
  )
}
predicted_svd_df <- do.call(cbind, neuralnet_predictions)
colnames(predicted_svd_df) <- paste0("Signal_", seq_len(mic_k))

# True test signals
true_svd_df <- mic_scores_test
colnames(true_svd_df) <- colnames(predicted_svd_df)

# Per-signal R2
r2_scores <- sapply(seq_len(mic_k), function(i) {
  calc_err(true_svd_df[, i], predicted_svd_df[, i])["R2"]
})
r2_results <- data.frame(Signal = colnames(true_svd_df), R2 = as.numeric(r2_scores))
print(r2_results)

# ---- Reconstruct abundances from predicted signals ----
reconstruct_abundance_all_signals <- function(pred_signals, train_signals, train_abundance) {
  reconstructed <- vector("list", ncol(train_abundance))
  Xtr <- as.matrix(train_signals)
  Xte <- as.matrix(pred_signals)
  for (i in seq_len(ncol(train_abundance))) {
    fit <- lm(as.numeric(train_abundance[, i]) ~ Xtr)
    coef_i <- coef(fit); coef_i[is.na(coef_i)] <- 0
    reconstructed[[i]] <- cbind(1, Xte) %*% as.numeric(coef_i)
  }
  out <- do.call(cbind, reconstructed)
  colnames(out) <- colnames(train_abundance)
  out
}

mapped_abundance_all_models <- reconstruct_abundance_all_signals(
  predicted_svd_df,
  mic_scores_train,
  microbiome_data[train_indices, , drop = FALSE]
)

actual_abundance <- microbiome_data[test_indices, , drop = FALSE]

comparison_all_models <- data.frame(
  Actual    = as.vector(as.matrix(actual_abundance)),
  Predicted = as.vector(as.matrix(mapped_abundance_all_models))
)

calc_errors_overall <- function(actual, forecast) {
  mae <- mean(abs(actual - forecast))
  mse <- mean((actual - forecast)^2)
  rmse <- sqrt(mse)
  r2 <- 1 - sum((actual - forecast)^2) / sum((actual - mean(actual))^2)
  list(MAE = mae, MSE = mse, RMSE = rmse, R2 = r2)
}
overall <- calc_errors_overall(comparison_all_models$Actual, comparison_all_models$Predicted)
print(overall)

# Per-taxon metrics
per_taxon_metrics <- tibble(
  Taxon = colnames(actual_abundance),
  R2  = apply(actual_abundance, 2, function(a, p) calc_errors_overall(a, p)$R2,  p = mapped_abundance_all_models),
  MAE = apply(actual_abundance, 2, function(a, p) calc_errors_overall(a, p)$MAE, p = mapped_abundance_all_models),
  RMSE= apply(actual_abundance, 2, function(a, p) calc_errors_overall(a, p)$RMSE,p = mapped_abundance_all_models)
)
