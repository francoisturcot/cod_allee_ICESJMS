# Supplementary: GAM production curves and Bmsy proxy
# Requires objects from the main analysis: df, fits_seg, stock_order

library(mgcv)

# Psi from the segmented fits (biomass/Bmax scale)
psi_df <- map_dfr(fits_seg, function(f) {
  tibble(stock = f$stock, psi = f$psi)
})

n_grid <- 500
alpha  <- 0.05

curve_pred  <- data.frame()
curve_class <- data.frame()

for (s in stock_order) {
  
  # Subset one stock (biomass already scaled to Bmax)
  DF <- as.data.frame(df[df$stock == s, ])
  DF <- DF[!is.na(DF$sp) & !is.na(DF$biomass), ]
  
  b_min <- min(DF$biomass)
  b_max <- max(DF$biomass)
  
  # Empirical production curve
  fit <- gam(sp ~ s(biomass), data = DF, family = "gaussian", method = "REML")
  smooth_p   <- summary(fit)$s.table[1, "p-value"]
  smooth_edf <- summary(fit)$s.table[1, "edf"]
  
  # Predictions over observed range
  PF <- data.frame(biomass = seq(b_min, b_max, length.out = n_grid))
  pr <- predict(fit, newdata = PF, se = TRUE)
  PF$fit   <- pr$fit
  PF$lower <- pr$fit - 1.96 * pr$se.fit
  PF$upper <- pr$fit + 1.96 * pr$se.fit
  PF$stock <- s
  curve_pred <- rbind(curve_pred, PF)
  
  # Peak of fitted curve
  i_max  <- which.max(PF$fit)
  b_peak <- PF$biomass[i_max]
  
  # Linear predictor matrix and coefficient covariance
  X    <- predict(fit, newdata = PF, type = "lpmatrix")
  beta <- coef(fit)
  V    <- vcov(fit)
  
  # Rise: P(peak) - P(Bmin)
  d_rise     <- X[i_max, ] - X[1, ]
  rise       <- as.numeric(d_rise %*% beta)
  rise_se    <- sqrt(as.numeric(t(d_rise) %*% V %*% d_rise))
  rise_lower <- rise - 1.96 * rise_se
  
  # Decline: P(peak) - P(Bmax)
  d_fall     <- X[i_max, ] - X[n_grid, ]
  fall       <- as.numeric(d_fall %*% beta)
  fall_se    <- sqrt(as.numeric(t(d_fall) %*% V %*% d_fall))
  fall_lower <- fall - 1.96 * fall_se
  
  # Significant rise / decline (lower 95% bound above zero; zero SE at an edge peak gives FALSE)
  rise_sig <- rise_se > 0 & rise_lower > 0
  fall_sig <- fall_se > 0 & fall_lower > 0
  
  # Classify
  if (smooth_p > alpha) {
    shape  <- "flat"
    interp <- "Bmsy undetermined"
    bmsy   <- NA
  } else if (rise_sig & fall_sig) {
    shape  <- "dome"
    interp <- "Bmsy within observed range"
    bmsy   <- b_peak
  } else if (rise_sig & !fall_sig) {
    shape  <- ifelse(i_max == n_grid, "increasing", "plateau")
    interp <- ifelse(i_max == n_grid, "Bmsy > Bmax", "Bmsy not identified (production levels off)")
    bmsy   <- NA
  } else if (!rise_sig & fall_sig) {
    shape  <- "decreasing"
    interp <- "Bmsy < Bmin"
    bmsy   <- NA
  } else {
    shape  <- "flat"
    interp <- "Bmsy undetermined"
    bmsy   <- NA
  }
  
  curve_class <- rbind(curve_class, data.frame(stock = s, edf = smooth_edf, p_value = smooth_p,
                                               B_peak = b_peak, rise_lower = rise_lower,
                                               fall_lower = fall_lower, shape = shape,
                                               interpretation = interp, Bmsy_proxy = bmsy))
}

# Psi as a fraction of the Bmsy proxy
curve_class <- left_join(curve_class, psi_df, by = "stock")
curve_class$psi_over_Bmsy <- curve_class$psi / curve_class$Bmsy_proxy

# Facet order
curve_pred$stock  <- factor(curve_pred$stock, levels = stock_order)
curve_class$stock <- factor(curve_class$stock, levels = stock_order)
df_supp <- df %>% filter(stock %in% stock_order)
df_supp$stock <- factor(df_supp$stock, levels = stock_order)

bmsy_lines <- curve_class %>% filter(!is.na(Bmsy_proxy))
psi_lines  <- curve_class %>% filter(!is.na(psi))

# Supplementary figure
p_supp <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
  geom_point(data = df_supp, aes(x = biomass, y = sp), size = 1, alpha = 0.7) +
  geom_ribbon(data = curve_pred, aes(x = biomass, ymin = lower, ymax = upper),
              fill = "lightblue", alpha = 0.5) +
  geom_line(data = curve_pred, aes(x = biomass, y = fit), linewidth = 1) +
  geom_vline(data = bmsy_lines, aes(xintercept = Bmsy_proxy), colour = "orange", linewidth = 0.8) +
  geom_vline(data = psi_lines, aes(xintercept = psi), colour = "grey", linetype = "dashed", linewidth = 0.8) +
  geom_text(data = curve_class, aes(x = 0.02, y = Inf, label = shape),
            hjust = 0, vjust = 1.5, size = 2.5) +
  facet_wrap(~ stock, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1),
                     labels = c("0", "0.25", "0.50", "0.75", "1")) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        strip.text = element_text(size = 9),
        panel.grid.minor = element_blank()) +
  labs(x = "Biomass (proportion of maximum)",
       y = "Surplus production")

p_supp
curve_class
