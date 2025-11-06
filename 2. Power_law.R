#0. Load required libraries
library(ecofolio)
library(dplyr)
library(tidyr)
library(ggplot2)


# 1. Fit Taylor’s power law: variance = c * mean^z with comm_t as the community matrix
taylor_res <- fit_taylor(comm_t)
print(taylor_res)
# $c  is the intercept on the log scale
# $z  is the slope (Taylor’s exponent)

# 2. Plot the mean–variance relationship with confidence intervals
plot_mv(comm_t, show = "linear", ci = TRUE)


# 3. Compute per‐taxon mean and variance of abundance
taxon_stats <- comm_t %>%
  pivot_longer(cols = everything(), names_to = "taxon", values_to = "abund") %>%
  group_by(taxon) %>%
  summarise(
    mean_abund = mean(abund, na.rm = TRUE),
    var_abund  = var(abund,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(mean_abund > 0)

# 4. Fit the linear Taylor’s Law model on log–log scale
fit_tl <- lm(log10(var_abund) ~ log10(mean_abund), data = taxon_stats)
tl_coef <- coef(fit_tl)  # intercept = tl_coef[1], slope = tl_coef[2]

# 5. Plot with ggplot2: points + linear fit + 95% CI ribbon
ggplot(taxon_stats, aes(x = mean_abund, y = var_abund)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(
    method  = "lm",
    formula = y ~ x,
    se      = TRUE,
    color   = "#D95F02",
    fill    = "#D95F02",
    size    = 1
  ) +
  labs(
    x = "Mean abundance",
    y = "Variance of abundance",
    subtitle = paste0("slope (z) = ", round(tl_coef[2], 2))
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(fill = NA, color = "black"),
    plot.subtitle = element_text(face = "italic", size = 12)
  )
