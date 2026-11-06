# Packages: install if missing, then load

pkgs <- c("tidyverse", "segmented", "patchwork", "here", "mgcv")

missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}

library(tidyverse)
library(segmented)
library(patchwork)
library(here)

#read data file
df <- readRDS(here("data", "cod_clean.RDS"))

# Surplus production and production rate

df <- df %>% arrange(stock, year)

df <- df %>%
  group_by(stock) %>%
  mutate(sp = catch + lead(biomass) - biomass) %>%
  ungroup()

df$prod.rate <- df$sp / df$biomass

# drop rows without a production rate (last year of each stock)
df <- df %>% filter(!is.na(prod.rate))

#depletion
df_depletion <- df %>%
  group_by(stock) %>%
  summarise(
    max_biomass = max(biomass, na.rm = TRUE),
    lrp         = 0.4 * max_biomass,
    min_biomass = min(biomass, na.rm = TRUE),
    depletion   = min_biomass / lrp,
    .groups = "drop"
  )

#b/bmax
df <- df %>%
  group_by(stock) %>%
  mutate(biomass_scaled = biomass / max(biomass, na.rm = TRUE)) %>%
  ungroup()
df$biomass = df$biomass_scaled


#stocks by region
nw_stocks <- c("2J3KL", "3NO", "3M", "3Pn4RS", "3Ps", "4TVn", "4VsW",
               "4X5Y", "5Zjm", "GB", "EGOM", "WGOM", "SNE", "1F-XIV", "1IN")
ne_stocks <- c("NEAR", "ICE", "FAPL", "IS", "IIIaW-IV-VIId",
               "BA2532", "BA2224", "VIa", "VIIek")

stock_order <- c(nw_stocks, ne_stocks)
stock_order <- stock_order[stock_order %in% unique(df$stock)]

# fit lm for every stock, and segmented where it converges

fits_all <- map(sort(unique(df$stock)), function(s) {
  temp <- df %>% filter(stock == s)
  m    <- lm(prod.rate ~ biomass, data = temp)
  sm   <- try(segmented(m, seg.Z = ~ biomass), silent = TRUE)
  out  <- list(stock = s, data = temp, lm = m, seg = FALSE)
  if (!(inherits(sm, "try-error") || is.null(sm$psi))) {
    b1 <- as.numeric(coef(sm)[2])
    b2 <- b1 + as.numeric(coef(sm)[3])
    out <- c(out, list(mod = sm, psi = as.numeric(sm$psi[2]), slope2 = b2))
    out$seg <- b1 > 0 & b2 < 0 & AIC(sm) < AIC(m)
  }
  out
})

fits_seg <- fits_all %>% keep(~ .x$seg)
fits_lin <- fits_all %>% keep(~ !.x$seg)

add_region <- function(d) {
  d %>% mutate(
    region = if_else(stock %in% nw_stocks, "Northwest Atlantic", "Northeast Atlantic"),
    stock  = factor(stock, levels = stock_order)
  )
}

# get psi
df_obs <- map_dfr(fits_all, function(f) {
  f$data %>% mutate(psi = if (f$seg) f$psi else NA_real_)
}) %>% add_region()

df_psi <- df_obs %>% filter(!is.na(psi)) %>% distinct(stock, psi, region)

# linear fits (non-segmented stocks)
df_lin <- map_dfr(fits_lin, function(f) {
  b <- seq(min(f$data$biomass), max(f$data$biomass), length.out = 100)
  tibble(stock = f$stock, biomass = b, segment = "linear",
         prod.rate = predict(f$lm, data.frame(biomass = b)))
})

# segmented fits, split at psi for colouring
df_seg <- map_dfr(fits_seg, function(f) {
  b1 <- seq(min(f$data$biomass), f$psi, length.out = 100)
  b2 <- seq(f$psi, max(f$data$biomass), length.out = 100)
  bind_rows(tibble(biomass = b1, segment = "before_break"),
            tibble(biomass = b2, segment = "after_break")) %>%
    mutate(stock = f$stock,
           prod.rate = predict(f$mod, data.frame(biomass = biomass)))
})

# post-break slope extrapolated below psi
df_ext <- map_dfr(fits_seg, function(f) {
  b <- seq(min(f$data$biomass), f$psi, length.out = 100)
  tibble(stock = f$stock, biomass = b,
         prod.rate = predict(f$mod, data.frame(biomass = f$psi)) + f$slope2 * (b - f$psi))
})

df_lines <- bind_rows(df_lin, df_seg) %>% add_region()
df_ext   <- add_region(df_ext)

# significance markers: linear stocks only (slope p < 0.05)
df_sig <- map_dfr(fits_lin, function(f) {
  tibble(stock = f$stock,
         p = summary(f$lm)$coefficients["biomass", "Pr(>|t|)"])
}) %>%
  filter(p < 0.05) %>%
  add_region()

make_panel <- function(reg) {
  ggplot(df_obs %>% filter(region == reg), aes(x = biomass, y = prod.rate)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey", linewidth = 0.5) +
    geom_vline(data = df_psi %>% filter(region == reg), aes(xintercept = psi),
               colour = "grey", linetype = "dashed", linewidth = 0.5) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_line(alpha = 0.4, linewidth = 0.3) +
    geom_line(data = df_ext %>% filter(region == reg),
              linetype = "dashed", linewidth = 0.7, colour = "coral") +
    geom_line(data = df_lines %>% filter(region == reg),
              aes(colour = segment), linewidth = 0.8) +
    geom_text(data = df_sig %>% filter(region == reg),
              x = Inf, y = Inf, label = "*",
              hjust = 1.5, vjust = 1.2, size = 6,
              colour = "black", inherit.aes = FALSE) +
    scale_colour_manual(values = c(before_break = "#00B8E7",
                                   after_break  = "coral",
                                   linear       = "grey30")) +
    facet_wrap(~ stock, ncol = 4, scales = "free_y") +
    scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1),
                       labels = c("0", "0.5", "1"),
                       expand = expansion(mult = 0.08)) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 8),
      axis.title = element_blank(),
      strip.text = element_text(size = 8),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 9, face = "bold")
    ) +
    ggtitle(reg)
}

p_nw <- make_panel("Northwest Atlantic")
p_ne <- make_panel("Northeast Atlantic")

n_rows <- function(reg) ceiling(n_distinct(df_obs$stock[df_obs$region == reg]) / 4)

p_combined <- p_nw / p_ne +
  plot_layout(heights = c(n_rows("Northwest Atlantic"), n_rows("Northeast Atlantic")))

p_final <- wrap_elements(
  grid::textGrob("Relative surplus production", rot = 90, gp = grid::gpar(fontsize = 9))
) +
  wrap_elements(
    p_combined + plot_annotation(
      caption = "Biomass (proportion of maximum)",
      theme = theme(plot.caption = element_text(hjust = 0.5, size = 9)))
  ) +
  plot_layout(widths = c(1, 40))

p_final
