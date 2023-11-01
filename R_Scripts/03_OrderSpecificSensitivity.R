library(dplyr)
library(ggplot2)
library(brms)
library(cmdstanr)
library(tidybayes)
library(stringr)
set_cmdstan_path()

#Script to run continued models
mdf <- read.csv("data/modelData_allSpecies.csv")

mdf_scaled <- mdf %>% 
    mutate(Dev_1 = scale(Dev_1))

traits <- read.csv("data/TraitData.csv")

mdf_scaled <- mdf_scaled %>% left_join(traits, by = c("scientific_name" = "Species"))

Order_abund <- mdf_scaled %>% 
    group_by(Order, Site) %>% 
    summarise(totalAbund = sum(abundance))

atFiveSites <- Order_abund %>% 
    filter(totalAbund > 0) %>% 
    group_by(Order) %>% 
    summarise(count = length(unique(Site))) %>% 
    filter(count >= 5)

mdf_scaled <- mdf_scaled %>% 
    filter(Order %in% atFiveSites$Order)

dev <- read.csv("data/impervious_surface.csv")

mdf2 <- Order_abund %>% 
    filter(Order %in% atFiveSites$Order) %>% 
    left_join(dev)

model1 <- brm(formula = bf(totalAbund ~ Dev_1 +
                               (1+ Dev_1 | Order)),
              family = negbinomial(),
              data = mdf2,
              chains = 4, iter = 2400, warmup = 1000,
              control = list(adapt_delta = 0.99),
              cores = 4, seed = 1234, 
              threads = threading(2),
              backend = "cmdstanr")

# examine model assumptions
plot(model1)
pp_check(model1)
pp_check(model1, type = "stat", stat = "mean")
pp_check(model1, type = "bars")

## examine output
summary(model1)

m_sum <- summary(model1) 
fixed <- m_sum$fixed %>% 
    tibble::rownames_to_column()
random <- ranef(model1)$Order %>% 
    unlist() %>% as.data.frame()

# tab outputs
write.csv(fixed, "Tables/fixedEffects_Orderspecific.csv", row.names = F)

write.csv(random, "Tables/randomEffects_Orderspecific.csv", row.names = F)

# plot conditional effects
ce <- conditional_effects(x = model1, effects = "Dev_1")
ce_df <- ce$Dev_1

ggplot() +
    geom_jitter(mdf2, mapping = aes(x = Dev_1, y = totalAbund, color = Order), alpha = 0.5) +
    geom_line(ce_df, mapping = aes(x = Dev_1, y = estimate__)) +
    geom_ribbon(ce_df, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__), alpha = 0.3) +
    labs(x = "Urban development", y = "Total Abundance") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16))

ggsave(filename = "Figures/OverallResponseToUrb_OrderRE.png", dpi = 450)

# plot Order specific results

draw <- model1 %>%
    spread_draws(b_Dev_1, r_Order[Order,Dev_1]) %>%
    # add the grand mean to the group-specific deviations
    mutate(mu = b_Dev_1 + r_Order) %>%
    ungroup() %>%
    mutate(Order = str_replace_all(Order, "[.]", " ")) 

ggplot(t, aes(x = mu, y = reorder(Order, mu))) +
    geom_vline(xintercept = fixef(model1)[2, 1], color = "#839496", size = 1) +
    geom_vline(xintercept = fixef(model1)[2, 3:4], color = "#839496", linetype = "dotted", size = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    stat_halfeye(.width = c(0.025,0.975), size = 2/3, fill = "dodgerblue", alpha = 0.9) +
    labs(x = expression("Sensitivity to urban development"),
         y = "") +
    theme_classic() 

ggsave(filename = "Figures/OrderResponseToUrb.png", dpi = 450)
