library(dplyr)
library(ggplot2)
library(brms)
library(ggpubr)

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
# examine model assumptions
plot(rich_dev)
t <- pp_check(rich_dev, ndraws = 50)
ggsave(filename = "Figures/Supplemental/RichnessCheck.png", plot = t)
pp_check(rich_dev, type = "stat", stat = "mean")
ggsave(filename = "Figures/Supplemental/RichnessCheckHist.png")

# Abundance model
abund_df <- cm_long_urbanization %>% 
    group_by(Site) %>% 
    summarise(totalAbund = sum(abundance))

abund_df <- ungroup(abund_df) %>% 
    left_join(urb)

abund_dev <- brm(formula = bf(totalAbund ~ Dev_1),
                data = abund_df,
                family = negbinomial(),
                chains = 4, iter = 2400, warmup = 1000,
                control = list(adapt_delta = 0.93),
                cores = 4, seed = 1234, 
                threads = threading(2),
                backend = "cmdstanr", 
)

ce_abund <- conditional_effects(x = abund_dev, effects = "Dev_1")
ce_df_abund <- ce$Dev_1

abund_plot <- ggplot() +
    geom_jitter(abund_df, mapping = aes(x = Dev_1, y = totalAbund), alpha = 0.5) +
    geom_line(ce_df_abund, mapping = aes(x = Dev_1, y = estimate__)) +
    geom_ribbon(ce_df_abund, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__), alpha = 0.3) +
    labs(x = "Urban development", y = "Total abundance") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16))

abund_plot  

ga <- ggarrange(rich_plot, abund_plot, labels = c("A", "B"))

ggsave(plot = abund_plot, filename = "Figures/Abundance.png", width = 4.5, height = 3)
