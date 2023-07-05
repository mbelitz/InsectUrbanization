library(dplyr)
library(lme4)
library(ggplot2)
library(lmerTest)
library(sjPlot)
library(MuMIn)
library(cowplot)

#Script to run continued models
mdf <- read.csv("Data/modelData_allSpecies.csv")

mdf_scaled <- mdf %>% 
    mutate(Dev_1 = scale(Dev_1),
           Dev_10 = scale(Dev_10),
           temp_niche = scale(temp_niche),
           meanLight = scale(log(meanLight + 1)),
           mean_temp = scale((mean_temp - 1.02) * -1),
           BodySize = scale(BodySize))

traits <- read.csv("Data/TraitData.csv") %>% 
    select(Species, Order)

mdf_scaled <- left_join(mdf_scaled, traits, by = c("scientific_name" = "Species"))

# Abundance Modeling #
model1 <- glmer.nb(formula = abundance ~ Dev_1 + meanLight + mean_temp +
                       LHSCategory:Dev_1 + VoltinismCategory:Dev_1 + 
                       LarvalHabitatCategory:Dev_1 + BodySize:Dev_1 +
                       temp_niche:Dev_1 + 
                       (1 | scientific_name),
                   na.action = "na.fail", family = poisson,
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                   data = mdf_scaled)

summary(model1)


model1reduced <- glmer.nb(formula = abundance ~ Dev_1 + 
                              LHSCategory:Dev_1 + 
                              (1 | scientific_name),
                          na.action = "na.fail", family = poisson,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                          data = mdf_scaled)

summary(model1reduced)

Weights(AICc(model1, model1reduced)) #model1 reduced is our first competing model
car::vif(model1reduced)

model2 <- glmer.nb(formula = abundance ~ Dev_1 + meanLight + mean_temp +
                       LHSCategory:meanLight + VoltinismCategory:meanLight + 
                       LarvalHabitatCategory:meanLight + BodySize:meanLight +
                       temp_niche:meanLight + 
                       (1 | scientific_name),
                   na.action = "na.fail", family = poisson,
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                   data = mdf_scaled)

summary(model2)

model2reduced <- glmer.nb(formula = abundance ~ Dev_1 + 
                              (1 | scientific_name),
                          na.action = "na.fail", family = poisson,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                          data = mdf_scaled)
summary(model2reduced)

Weights(AICc(model2, model2reduced)) #reduced model selected

model3 <- glmer.nb(formula = abundance ~ Dev_1 + meanLight + mean_temp +
                       LHSCategory:mean_temp + VoltinismCategory:mean_temp + 
                       LarvalHabitatCategory:mean_temp + BodySize:mean_temp +
                       temp_niche:mean_temp + 
                       (1 | scientific_name),
                   na.action = "na.fail", family = poisson,
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                   data = mdf_scaled)

summary(model3)

model3reduced <- glmer.nb(formula = abundance ~ Dev_1 + mean_temp +
                              LarvalHabitatCategory:mean_temp +
                              (1 | scientific_name),
                          na.action = "na.fail", family = poisson,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                          data = mdf_scaled)
summary(model3reduced)

Weights(AICc(model3, model3reduced)) #Reduced model selected
car::vif(model3reduced)


#10-km scale
model10 <- glmer.nb(formula = abundance ~ Dev_10 + meanLight + mean_temp +
                        LHSCategory:Dev_10 + VoltinismCategory:Dev_10 + 
                        LarvalHabitatCategory:Dev_10 + BodySize:Dev_10 +
                        temp_niche:Dev_10 + 
                        (1 | scientific_name),
                    na.action = "na.fail", family = poisson,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                    data = mdf_scaled)

summary(model10)

model10reduced <- glmer.nb(formula = abundance ~ Dev_10 + 
                               LHSCategory:Dev_10 + 
                               (1 | scientific_name),
                           na.action = "na.fail", family = poisson,
                           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                           data = mdf_scaled)

summary(model10reduced)

Weights(AICc(model10, model10reduced)) #Reduced model selected
car::vif(model10reduced)

model20 <- glmer.nb(formula = abundance ~ Dev_10 + meanLight + mean_temp +
                        LHSCategory:meanLight + VoltinismCategory:meanLight + 
                        LarvalHabitatCategory:meanLight + BodySize:meanLight +
                        temp_niche:meanLight + 
                        (1 | scientific_name),
                    na.action = "na.fail", family = poisson,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                    data = mdf_scaled)

summary(model20)

model20reduced <- glmer.nb(formula = abundance ~ Dev_10 + 
                               (1 | scientific_name),
                           na.action = "na.fail", family = poisson,
                           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                           data = mdf_scaled)
summary(model20reduced)

Weights(AICc(model20, model20reduced)) #reduced model selected
#car::vif(model20reduced)

model30 <- glmer.nb(formula = abundance ~ Dev_10 + meanLight + mean_temp +
                        LHSCategory:mean_temp + VoltinismCategory:mean_temp + 
                        LarvalHabitatCategory:mean_temp + BodySize:mean_temp +
                        temp_niche:mean_temp + 
                        (1 | scientific_name),
                    na.action = "na.fail", family = poisson,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                    data = mdf_scaled)

summary(model30)

model30reduced <- glmer.nb(formula = abundance ~ Dev_10 + mean_temp +
                               LarvalHabitatCategory:mean_temp + 
                               (1 | scientific_name),
                           na.action = "na.fail", family = poisson,
                           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                           data = mdf_scaled)
summary(model30reduced)

Weights(AICc(model30, model30reduced)) #Reduced model selected
car::vif(model30reduced)

Weights(AICc(model1reduced, model2reduced, model3reduced, model10reduced, model20reduced, model30reduced))
summary(model10reduced)
summary(model20reduced)
summary(model30reduced)


# model selection table
table <- model.sel(model1reduced, model2reduced, model3reduced, 
                   model10reduced, model20reduced, model30reduced)
write.csv(table,"Tables/model_selection_table.csv")

# test for overdispersion
library(DHARMa)
sim_m <- simulateResiduals(model20reduced, refit = T)
testDispersion(sim_m)
od_plot <- plot(sim_m) # no large overdispersion detected

# figures of top model (model20reduced)
top_model_output <- plot_model(model20reduced, type = "eff", terms = c("Dev_10")) +
    ggtitle("") +
    labs(x = "Proportion developed (10-km)", y = "Abundance") +
    theme_classic()
top_model_output
ggsave(filename = "Figures/TopModel_proportionDeveloped.png",plot = top_model_output, width = 5, height = 4)

# figure of second best model (model30reduced)
summary(model30reduced)
b <- plot_model(model30reduced, type = "eff", terms = c("mean_temp", "LarvalHabitatCategory")) +
    ggtitle("") +
    labs(x = "Relative temperature", y = "Abundance",
         color = "Larval habitat", fill = "Larval habitat") +
    scale_color_viridis_d(option = "turbo") +
    theme_classic()

b
ggsave(filename = "Figures/secondBestModel.png", dpi = 450, plot = b,
       width = 5, height = 3.5)

# third best model
summary(model10reduced)
a <- plot_model(model10reduced, type = "eff", terms = c("Dev_10", "LHSCategory")) +
    ggtitle("") +
    labs(x = "Proportion developed (10-km)", y = "Abundance", 
         color = "Larval diet", fill = "Larval diet") +
    scale_color_viridis_d() +
    theme_classic()
a

ggsave(filename = "Figures/thirdBestModel.png", dpi = 450, plot = a,
       width = 5, height = 3.5)
