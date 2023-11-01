library(dplyr)
library(ggplot2)
library(brms)

cmdstanr::set_cmdstan_path()
#### Creating Model Dataframe ####
cm_long_urbanization <- read.csv("data/modelData_allSpecies.csv")

urb <- read.csv("Data/impervious_surface.csv")

rich <- cm_long_urbanization %>% 
    filter(abundance != 0) %>% 
    group_by(Site) %>% 
    summarise(richness = length(unique(scientific_name)))%>% 
    left_join(urb)

#### Richness Modeling ####
rich_dev <- brm(formula = bf(richness ~ Dev_1),
                data = rich,
                family = gaussian(),
                chains = 4, iter = 2400, warmup = 1000,
                control = list(adapt_delta = 0.93),
                cores = 4, seed = 1234, 
                threads = threading(2),
                backend = "cmdstanr", 
)

# examine model assumptions
plot(rich_dev)
pp_check(rich_dev)
pp_check(rich_dev, type = "stat", stat = "mean")

summary(rich_dev)
rich_dev_sum <- summary(rich_dev)$fixed %>% 
    tibble::rownames_to_column()

write.csv(x = rich_dev_sum, file = "Tables/richnessModel.csv", row.names = F)

ce <- conditional_effects(x = rich_dev, effects = "Dev_1")
ce_df <- ce$Dev_1

rich_plot <- ggplot() +
    geom_jitter(rich, mapping = aes(x = Dev_1, y = richness), alpha = 0.5) +
    geom_line(ce_df, mapping = aes(x = Dev_1, y = estimate__)) +
    geom_ribbon(ce_df, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__), alpha = 0.3) +
    labs(x = "Urban development", y = "Richness") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16))

rich_plot  

ggsave(plot = rich_plot, filename = "Figures/Richness.png")

## source in biomass so those can be plotted together
source("R_Scripts/01_TotalBiomass.R")

library(cowplot)

cp <- plot_grid(biomass_plot, rich_plot, labels = c("A","B"),
                nrow = 2, ncol = 1)

ggsave(filename = "Figures/Biomass&Richness.png", plot = cp, width = 3.5, height = 5, dpi = 450)
