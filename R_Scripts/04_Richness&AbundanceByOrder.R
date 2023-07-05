library(dplyr)
library(ggplot2)
library(cowplot)

#Script to run continued models
mdf <- read.csv("data/modelData_allSpecies.csv")

mdf_scaled <- mdf %>% 
    mutate(Dev_1 = scale(Dev_1),
           Dev_10 = scale(Dev_10),
           temp_niche = scale(temp_niche),
           meanLight = scale(log(meanLight + 1)),
           mean_temp = scale((mean_temp - 1.02) * -1),
           BodySize = scale(BodySize))

traits <- read.csv("data/TraitData.csv")

mdf_scaled <- mdf_scaled %>% left_join(traits, by = c("scientific_name" = "Species"))

# pull out Order specific plots

blattodea_plot <- filter(mdf_scaled, Order == "Blattodea") %>% 
    ggplot() +
    geom_point(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name)) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name),
                method = "lm",
                alpha = 0.1, linewidth = 0.4, se = F) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance), method = "glm", 
                method.args = list(family = "poisson"), linewidth = 1) +
    geom_text(x = 0.9, y = 1.4, label = "Blattella asahinai", color = "black", size = 3) +
    geom_text(x = 0.75, y = -0.1, label = "Periplaneta americana", color = "black", size = 3) +
    scale_color_manual(values = rep("black", 2)) +
    scale_fill_manual(values = rep("black", 2)) +
    theme_classic()+
    theme(legend.position = "none") +
    labs(y = "Abundance", x = "Proportion developed (10-km)")
blattodea_plot

coleoptera_plot <- filter(mdf_scaled, Order == "Coleoptera") %>% 
    ggplot() +
    geom_point(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name)) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name),
                method = "lm",
                alpha = 0.1, linewidth = 0.4, se = F) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance), method = "glm", 
                method.args = list(family = "poisson"), linewidth = 1) +
    geom_text(x = 0.9, y = 140, label = "Hydaticus bimarginatus", color = "black", size = 3) +
    geom_text(x = 0, y = 70, label = "Callistethus marginatus", color = "black", size = 3) +
    scale_color_manual(values = rep("black", 19)) +
    scale_fill_manual(values = rep("black", 19)) +
    theme_classic()+
    theme(legend.position = "none") +
    labs(y = "Abundance", x = "Proportion developed (10-km)")
coleoptera_plot

diptera_plot <- filter(mdf_scaled, Order == "Diptera") %>% 
    ggplot() +
    geom_point(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name)) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name),
                method = "lm",
                alpha = 0.1, linewidth = 0.4, se = F) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance), method = "glm", 
                method.args = list(family = "poisson"), linewidth = 1) +
    geom_text(x = -1, y = 23, label = "Plecia nearctica", color = "black", size = 3) +
    geom_text(x = 0.5, y = 1.8, label = "Tabanus atratus", color = "black", size = 3) +
    scale_color_manual(values = rep("black", 19)) +
    scale_fill_manual(values = rep("black", 19)) +
    theme_classic()+
    theme(legend.position = "none") +
    labs(y = "Abundance", x = "Proportion developed (10-km)")
diptera_plot

diptera_plot <- filter(mdf_scaled, Order == "Diptera") %>% 
    ggplot() +
    geom_point(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name)) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name),
                method = "lm",
                alpha = 0.1, linewidth = 0.4, se = F) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance), method = "glm", 
                method.args = list(family = "poisson"), linewidth = 1) +
    geom_text(x = -1, y = 23, label = "Plecia nearctica", color = "black", size = 3) +
    geom_text(x = 0.5, y = 1.8, label = "Tabanus atratus", color = "black", size = 3) +
    scale_color_manual(values = rep("black", 19)) +
    scale_fill_manual(values = rep("black", 19)) +
    theme_classic()+
    theme(legend.position = "none") +
    labs(y = "Abundance", x = "Proportion developed (10-km)")
diptera_plot

hemimptera_plot <- filter(mdf_scaled, Order == "Hemiptera") %>% 
    ggplot() +
    geom_point(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name)) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name),
                method = "lm",
                alpha = 0.1, linewidth = 0.4, se = F) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance), method = "glm", 
                method.args = list(family = "poisson"), linewidth = 1)  +
    geom_text(x = -0.75, y = 85, label = "Prosapia bicincta", color = "black", size = 3) +
    geom_text(x = -1, y = 50, label = "Draeculacephala inscripta", color = "black", size = 3) +
    scale_color_manual(values = rep("black", 19)) +
    scale_fill_manual(values = rep("black", 19)) +
    theme_classic()+
    theme(legend.position = "none") +
    labs(y = "Abundance", x = "Proportion developed (10-km)")
hemimptera_plot

hymenoptera_plot <- filter(mdf_scaled, Order == "Hymenoptera") %>% 
    ggplot() +
    geom_point(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name)) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name),
                method = "lm",
                alpha = 0.1, linewidth = 0.4, se = F) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance), method = "glm", 
                method.args = list(family = "poisson"), linewidth = 1)  +
    geom_text(x = -0.65, y = 11, label = "Camponotus castaneus", color = "black", size = 3) +
    scale_color_manual(values = rep("black", 19)) +
    scale_fill_manual(values = rep("black", 19)) +
    theme_classic()+
    theme(legend.position = "none") +
    labs(y = "Abundance", x = "Proportion developed (10-km)")
hymenoptera_plot

megaloptera_plot <- filter(mdf_scaled, Order == "Megaloptera") %>% 
    ggplot() +
    geom_point(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name)) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name),
                method = "lm",
                alpha = 0.1, linewidth = 0.4, se = F) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance), method = "glm", 
                method.args = list(family = "poisson"), linewidth = 1)  +
    geom_text(x = 0.75, y = 16, label = "Chauliodes rastricornis", color = "black", size = 3) +
    geom_text(x = -0.5, y = 5, label = "Corydalus cornutus", color = "black", size = 3) +
    scale_color_manual(values = rep("black", 19)) +
    scale_fill_manual(values = rep("black", 19)) +
    theme_classic()+
    theme(legend.position = "none") +
    labs(y = "Abundance", x = "Proportion developed (10-km)")
megaloptera_plot

neuroptera_plot <- filter(mdf_scaled, Order == "Neuroptera") %>% 
    ggplot() +
    geom_point(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name)) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name),
                method = "lm",
                alpha = 0.1, linewidth = 0.4, se = F) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance), method = "glm", 
                method.args = list(family = "poisson"), linewidth = 1, se = F)  +
    scale_color_manual(values = rep("black", 19)) +
    scale_fill_manual(values = rep("black", 19)) +
    theme_classic()+
    theme(legend.position = "none") +
    labs(y = "Abundance", x = "Proportion developed (10-km)")
neuroptera_plot

orthoptera_plot <- filter(mdf_scaled, Order == "Orthoptera") %>% 
    ggplot() +
    geom_point(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name)) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name),
                method = "lm",
                alpha = 0.1, linewidth = 0.4, se = F) +
    geom_smooth(mapping = aes(x = Dev_10, y = abundance), method = "glm", 
                method.args = list(family = "poisson"), linewidth = 1, se = T)  +
    geom_text(x = -1, y = 2.3, label = "Neoscapteriscus borellii", color = "black", size = 3) +
    geom_text(x = -0.85, y = 1.2, label = "Stilpnochlora couloniana", color = "black", size = 3) +
    scale_color_manual(values = rep("black", 19)) +
    scale_fill_manual(values = rep("black", 19)) +
    theme_classic()+
    theme(legend.position = "none") +
    labs(y = "Abundance", x = "Proportion developed (10-km)")
orthoptera_plot

cp <- plot_grid(blattodea_plot, coleoptera_plot, 
                diptera_plot, hemimptera_plot,
                hymenoptera_plot, megaloptera_plot, 
                neuroptera_plot, orthoptera_plot,
                nrow = 4, ncol = 2, 
                labels = c("Blattodea", "Coleoptera",
                           "Diptera", "Hemiptera",
                           "Hymenoptera", "Megaloptera",
                           "Neuroptera", "Orthoptera"), label_x = 0.25)

ggsave(plot = cp, filename = "Figures/AbundanceResponsesByOrder.png",
       width = 9, height = 10, dpi = 400)



ggplot(mdf_scaled, mapping = aes(x = Dev_10, y = abundance, color = scientific_name, fill = scientific_name)) +
    geom_point() +
    geom_smooth(mdf_scaled, 
                mapping = aes(x = Dev_10, y = abundance, 
                              color = scientific_name, fill = scientific_name),
                alpha = 0.1,
                method = "lm",
                size =0.75,
                se = F) +
    theme_bw() +
    scale_color_manual(values = rep("black", 43)) +
    scale_fill_manual(values = rep("black", 43)) +
    labs(x = "Proportion developed (10-km)", y = "Abundance") +
    theme(legend.position = "none") +
    facet_wrap(~NativeStatus, scales = "free") 

ggsave(filename = "Figures/AbundanceResponseNativeInvasive.png", width = 6 , height = 4)

# richness plots by order
rich_zero <- mdf_scaled %>% 
    group_by(Order, Site, .drop = FALSE) %>% 
    summarise(zero = 0)

rich_true <- filter(mdf_scaled, abundance > 0) %>% 
    group_by(Site, Order) %>% 
    summarise(richness = length(unique(scientific_name)))

rich <- left_join(rich_zero, rich_true)
rich <- rich %>% 
    mutate(richness = if_else(condition = is.na(richness), 
                              true = 0, 
                              false = richness))


urb <- read.csv('Data/impervious_surface.csv')
light <- read.csv('Data/lightData.csv')

rich <- left_join(rich, light, by = "Site")

richP <- ggplot(rich, aes(x = log(meanLight + 0.1), y = richness)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_classic() +
    labs(x = "Artificial Light at Night (log10 + 0.1 lux)",
         y = "Richness") +
    facet_wrap(~Order, scales = "free", nrow = 4, ncol = 2)

ggsave(plot = richP, filename = "Figures/RichnessByOrder.png",
       width = 6, height = 8, dpi = 400)
