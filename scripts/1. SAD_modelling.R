#0. Load required libraries
library(sads)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble) 
set.seed(999)

# 1. Define modified R² function on rank–abundance curves (RAD)
r2m <- function(obs, pred) {
  o <- obs$abund
  p <- pred$abund
  1 - sum((o - p)^2) / sum((o - mean(o))^2)
}

# 2. Fit four SAD models to each sample & compute modified R²
results <- tibble(
  Sample    = colnames(comm),
  r2_bs     = NA_real_,
  r2_ls     = NA_real_,
  r2_poilog = NA_real_,
  r2_zipf   = NA_real_
)

for (j in seq_len(n_samples)) {
  abund <- comm[, j]
  abund <- abund[abund > 0]        # drop zeros
  
  # fit models
  fit_bs     <- fitsad(abund, sad = "bs")
  fit_ls     <- fitsad(abund, sad = "ls")
  fit_poilog <- fitsad(abund, sad = "poilog")
  fit_zipf   <- fitrad(abund, rad = "zipf", N = sum(abund))
  
  # get observed & predicted RADs
  rad_obs    <- rad(abund)
  rad_bs     <- radpred(fit_bs)
  rad_ls     <- radpred(fit_ls)
  rad_poilog <- radpred(fit_poilog)
  rad_zipf   <- radpred(fit_zipf)
  
  # compute and store r2m for each model
  results$r2_bs[j]     <- r2m(rad_obs, rad_bs)
  results$r2_ls[j]     <- r2m(rad_obs, rad_ls)
  results$r2_poilog[j] <- r2m(rad_obs, rad_poilog)
  results$r2_zipf[j]   <- r2m(rad_obs, rad_zipf)
}

# 3. Reshape and summarise R² across all samples

# 3a. Make long table of per-sample fits
r2_long <- results %>%
  pivot_longer(
    cols      = -Sample,
    names_to  = "model",
    values_to = "r2"
  )

write.csv(r2_long,"r2_long_sample_as_mag.csv")

# 3b. Compute mean and SD of R² for each model
r2_stats <- r2_long %>%
  group_by(model) %>%
  summarise(
    mean = mean(r2),
    sd   = sd(r2),
    .groups = "drop"
  )

write.csv(r2_stats, "r2_stats_sample_as_mag.csv")
# 3c. Compute percentage of samples in which each model was best
best_pct <- r2_long %>%
  group_by(Sample) %>%
  slice_max(r2, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  count(model) %>%
  mutate(pct = 100 * n / n_samples) %>%
  select(model, pct)

# 3d. Combine stats and best‐model percentages
r2_summary <- left_join(r2_stats, best_pct, by = "model")
write.csv(r2_summary,"r2_summary_mag_as_mag.csv")

# 4. Aggregate all samples into one meta‐community and fit models
abund_agg <- rowSums(comm)
abund_agg <- abund_agg[abund_agg > 0]

agg_bs     <- fitsad(abund_agg,   sad = "bs")
agg_ls     <- fitsad(abund_agg,   sad = "ls")
agg_poilog <- fitsad(abund_agg,   sad = "poilog")
agg_zipf   <- fitrad(abund_agg,   rad = "zipf", N = sum(abund_agg))

# 5. Get observed & predicted RADs for the aggregated community
rad_obs_agg    <- rad(abund_agg)
rad_bs_agg     <- radpred(agg_bs)
rad_ls_agg     <- radpred(agg_ls)
rad_poilog_agg <- radpred(agg_poilog)
rad_zipf_agg   <- radpred(agg_zipf)

# 6. Build a long data.frame for ggplot
df_all <- bind_rows(
  tibble(rank = rad_obs_agg$rank,    abund = rad_obs_agg$abund,    model = "Observed"),
  tibble(rank = rad_bs_agg$rank,     abund = rad_bs_agg$abund,     model = "Broken-stick"),
  tibble(rank = rad_ls_agg$rank,     abund = rad_ls_agg$abund,     model = "Log-series"),
  tibble(rank = rad_poilog_agg$rank, abund = rad_poilog_agg$abund, model = "Poisson lognormal"),
  tibble(rank = rad_zipf_agg$rank,   abund = rad_zipf_agg$abund,   model = "Zipf")
)


# 7. Plot Fig. 2–style panel with ggplot2
df_wide <- df_all %>%
  pivot_wider(names_from = model, values_from = abund)

ggplot(df_wide, aes(x = rank)) +
 
  geom_line(aes(y = log10(`Poisson lognormal`), color = "Poisson lognormal"),
               size = 1.5) +
  geom_line(aes(y = log10(`Broken-stick`),       color = "Broken-stick"),
              size = 1.5) +
  geom_line(aes(y = log10(`Log-series`),         color = "Log-series"),
              size = 1.5) +
  geom_line(aes(y = log10(Zipf),                 color = "Zipf"),
              size = 1.5) + 
  geom_point(aes(y = log10(Observed),              color = "Observed"),
                                       size = 3, shape = 1,alpha=0.2) +

  labs(
    x = "Abundance rank",
    y = "log10(Abundance)"
  ) +
  theme_classic() +
  scale_color_manual(
    breaks = c("Observed","Broken-stick","Poisson lognormal","Log-series","Zipf"),
    values = c(
      "Observed"           = "black",
      "Broken-stick"       = "#fdae61",
      "Poisson lognormal" = "#d53e4f",
      "Log-series"         = "#5e4fa2",
      "Zipf"               = "#3288bd"
    )
  ) +
  theme_classic() +
  theme(
    legend.title   = element_blank(),
    legend.position= c(0.8, 0.75),
    legend.text    = element_text(size = 10),
    legend.key.size= unit(1.5, "lines"),
    panel.border   = element_rect(fill = NA, size = 1, linetype = "solid")
  )+ theme(axis.text = element_text(size = 10))+theme(axis.title = element_text(size = 15,face="bold"))
