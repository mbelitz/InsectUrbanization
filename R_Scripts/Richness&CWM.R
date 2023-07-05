library(dplyr)
library(ggplot2)
library(sjPlot)
library(MuMIn)

#### Creating Model Dataframe ####
cm_long_urbanization <- read.csv("data/modelData_allSpecies.csv")

urb <- read.csv("Data/impervious_surface.csv")
tempData <- read.csv("Data/temp_gradient.csv")
light <- read.csv("Data/lightData.csv")

rich <- cm_long_urbanization %>% 
    filter(abundance != 0) %>% 
    group_by(Site) %>% 
    summarise(richness = length(unique(scientific_name)))%>% 
    left_join(urb) %>% 
    left_join(tempData) %>% 
    left_join(light) %>% 
    mutate(meanLight = log(meanLight + 1),
           mean_temp = (mean_temp - 1.02) * -1)

#### Richness Modeling ####
modelRich <- lm(formula = richness ~ Dev_1,
                na.action = "na.fail",
                data = rich)

modelRich2 <- lm(formula = richness ~ Dev_10,
                 na.action = "na.fail",
                 data = rich)

modelRich3 <- lm(formula = richness ~ mean_temp,
                 na.action = "na.fail",
                 data = rich)

modelRich4 <- lm(formula = richness ~ meanLight,
                 na.action = "na.fail",
                 data = rich)

Weights(AICc(modelRich, modelRich2, modelRich3, modelRich4))
table <- MuMIn::model.sel(modelRich, modelRich2, modelRich3, modelRich4)
table
summary(modelRich4)
r.squaredGLMM(modelRich4)

a <- sjPlot::plot_model(modelRich4, terms = "meanLight", type = "pred", title = "") 

a <- a + 
    geom_point(rich, mapping = aes(x = meanLight, y = richness)) +
    labs(x = "Artificial light at night", y = "Species Richness") +
    font_size(title  = 10, axis_title.x = 10, axis_title.y = 10, labels.x = 8, labels.y = 8) + theme_classic()
a


#### Community Weighted Mean Body Size Dataframe ####
cwm <- cm_long_urbanization %>% 
    group_by(Site) %>% 
    summarise(
        BodySize_cwm =weighted.mean(BodySize, abundance, na.rm=T)
    )

cwm_urb <- left_join(cwm, urb) %>% 
    left_join(light) %>% 
    left_join(tempData) %>% 
    mutate(meanLight = log(meanLight + 1),
           mean_temp = (mean_temp - 1.02) * -1)

head(cwm_urb)


## CWM Body Size Modeling##
modelCWM <- lm(formula = BodySize_cwm ~ Dev_1,
               na.action = "na.fail",
               data = cwm_urb)

modelCWM2 <- lm(formula = BodySize_cwm ~ Dev_10,
                na.action = "na.fail",
                data = cwm_urb)

modelCWM3 <- lm(formula = BodySize_cwm ~ mean_temp,
                na.action = "na.fail",
                data = cwm_urb)

modelCWM4 <- lm(formula = BodySize_cwm ~ meanLight,
                na.action = "na.fail",
                data = cwm_urb)

Weights(AICc(modelCWM, modelCWM2, modelCWM3, modelCWM4))
table <- MuMIn::model.sel(modelCWM, modelCWM2, modelCWM3, modelCWM4)
table
summary(modelCWM2)
r.squaredGLMM(modelCWM2)

b <- sjPlot::plot_model(modelCWM2, terms = "Dev_10", type = "pred", title = "")

b <- b + 
    geom_point(cwm_urb, mapping = aes(x = Dev_10, y = BodySize_cwm)) +
    labs(x = "Proportion urbanization (10-km)", y = "Body Size CWM") +
    font_size(title  = 10, axis_title.x = 10, axis_title.y = 10, labels.x = 8, labels.y = 8) + theme_classic()

b

## Community Weight Mean Temperature Niche Dataframe ##
cwm <- cm_long_urbanization %>% 
    group_by(Site) %>% 
    summarise(
        tempNiche_cwm =weighted.mean(temp_niche, abundance, na.rm=T)
    )

cwm_urb <- left_join(cwm, urb) %>% 
    left_join(light) %>% 
    left_join(tempData) %>% 
    mutate(meanLight = log(meanLight + 1),
           mean_temp = (mean_temp - 1.02) * -1)

head(cwm_urb)

## CWM Temperature Niche Modeling ##
modelCWM <- lm(formula = tempNiche_cwm ~ Dev_1,
               na.action = "na.fail",
               data = cwm_urb)

modelCWM2 <- lm(formula = tempNiche_cwm ~ Dev_10,
                na.action = "na.fail",
                data = cwm_urb)

modelCWM3 <- lm(formula = tempNiche_cwm ~ mean_temp,
                na.action = "na.fail",
                data = cwm_urb)

modelCWM4 <- lm(formula = tempNiche_cwm ~ meanLight,
                na.action = "na.fail",
                data = cwm_urb)

Weights(AICc(modelCWM, modelCWM2, modelCWM3, modelCWM4))
table <- MuMIn::model.sel(modelCWM, modelCWM2, modelCWM3, modelCWM4)
table
summary(modelCWM)
r.squaredGLMM(modelCWM)

c <- sjPlot::plot_model(modelCWM, terms = "Dev_1", type = "pred", title = "")

c <- c + 
    geom_point(cwm_urb, mapping = aes(x = Dev_1, y = tempNiche_cwm)) +
    labs(x = "Proportion urbanization (1-km)", y = "Temperature Niche CWM") +
    font_size(title  = 10, axis_title.x = 10, axis_title.y = 10, labels.x = 8, labels.y = 8) + theme_classic()

c

#ggsave("Figures/CWMniche.png", width = 7, height = 7)


Figure4 <- cowplot::plot_grid(a,b,c, labels = c("A", "B", "C"), label_size = 12)
Figure4

ggsave("Figures/Figure4_CWM.png", Figure4,
       width = 6, height = 4.5)
