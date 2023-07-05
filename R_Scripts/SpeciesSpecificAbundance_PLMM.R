library(lme4)
library(dplyr)
library(car)
library(MuMIn)
library(sjPlot)
library(ggplot2)
library(ape)
library(phyr)

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
    dplyr::select(Order, Species)

mdf_scaled <- left_join(mdf_scaled, traits, by = c("scientific_name" = "Species"))


# function to capitalize spp names
firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}

tt <- read.tree("Data/insect_tree_wBranches.tre")
tt$tip.label <- stringr::str_replace(tt$tip.label, pattern = "_", " ")
tt$tip.label <- firstup(tt$tip.label)


mdf_phylo <- left_join(mdf_scaled, data.frame(scientific_name = tt$tip.label, Phylo = "Yes"))

sppNotInAnlysis <- data.frame(scientific_name = tt$tip.label) %>% 
    filter(!scientific_name %in% mdf_phylo$scientific_name)

mdf_phylo <- mdf_phylo %>% 
    filter(!is.na(Phylo))

tt <- ape::drop.tip(tt, tip = sppNotInAnlysis$validName)


model20_pglmm <-  pglmm(formula = abundance ~ Dev_10 + 
                            (1 | scientific_name__),
                        data = mdf_phylo, 
                        family = "zeroinflated.poisson",
                        cov_ranef = list(scientific_name = tt), 
                        bayes = TRUE)

summary(model20_pglmm)

model30_pglmm <- pglmm(formula = abundance ~ Dev_10 + mean_temp +
                           LarvalHabitatCategory:mean_temp + 
                           (1 | scientific_name__),
                       data = mdf_phylo, 
                       family = "zeroinflated.poisson",
                       cov_ranef = list(scientific_name = tt), 
                       bayes = TRUE)
summary(model30_pglmm)


model10_pglmm <- pglmm(formula = abundance ~ Dev_10 + 
                           LHSCategory:Dev_10 + 
                           (1 | scientific_name__),
                       data = mdf_phylo, 
                       family = "zeroinflated.poisson",
                       cov_ranef = list(scientific_name = tt), 
                       bayes = TRUE)

summary(model10_pglmm)