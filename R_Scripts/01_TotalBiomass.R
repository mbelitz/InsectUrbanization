library(tidyverse)
library(brms)
cmdstanr::set_cmdstan_path()

mdf <- read.csv("Data/BycatchIdentification.csv") %>% 
    mutate(Site = case_when(
        Site %in% c("Cofr", "Cofr ") ~ "Cofr",
        Site %in% c("Prcr", "prcr", "PRCR") ~ "Prcr",
        Site %in% c("Demi", "Demi ") ~ "Demi",
        Site %in% c("Rist", "Rist ") ~ "Rist",
        Site %in% c("Biva", "Biva ", "Biar") ~ "Biva",
        Site %in% c("Auca", "Auca ") ~ "Auca",
        Site %in% c("Baca", "Baca ") ~ "Baca",
        Site %in% c("Bowa", "Bowa ") ~ "Bowa",
        Site %in% c("Joma", "Joma ") ~ "Joma"
    ))

biomass <- mdf %>% 
    group_by(Site,Date) %>% 
    summarise(biomass = length(unique(Tube..)))

totBiomass <- biomass %>% 
    group_by(Site) %>% 
    summarise(totBiomass = sum(biomass),
              collectionDays = length(unique(Date)),
              BiomassPerCollection = totBiomass/collectionDays)


urb <- read.csv("Data/impervious_surface.csv")

mdf2 <- totBiomass %>% 
    left_join(urb)

#### Richness Modeling ####
biomass_dev <- brm(formula = bf(BiomassPerCollection ~ Dev_1),
                data = mdf2,
                family = gaussian(),
                chains = 4, iter = 2400, warmup = 1000,
                control = list(adapt_delta = 0.93),
                cores = 4, seed = 1234, 
                threads = threading(2),
                backend = "cmdstanr", 
)

# examine model assumptions
plot(biomass_dev)
pp_check(biomass_dev, ndraws = 50)
pp_check(biomass_dev, type = "stat", stat = "mean")

summary(biomass_dev)
biomass_dev_sum <- summary(biomass_dev)$fixed %>% 
    tibble::rownames_to_column()

write.csv(x = biomass_dev_sum, file = "Tables/totalBiomassModel.csv", row.names = F)

ce <- conditional_effects(x = biomass_dev, effects = "Dev_1")
ce_df <- ce$Dev_1

biomass_plot <- ggplot() +
    geom_jitter(mdf2, mapping = aes(x = Dev_1, y = BiomassPerCollection), alpha = 0.5) +
    geom_line(ce_df, mapping = aes(x = Dev_1, y = estimate__)) +
    geom_ribbon(ce_df, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__), alpha = 0.3) +
    labs(x = "Urban development", y = "Total Biomass per Collection") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16))

biomass_plot  

ggsave(filename = "Figures/TotalBiomass.png", plot = biomass_plot, dpi = 450)
