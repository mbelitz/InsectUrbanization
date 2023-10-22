library(tidyverse)
library(brms)
library(ape)
library(cmdstanr)
set_cmdstan_path()


#Script to run continued models
mdf <- read.csv("Data/modelData_allSpecies.csv")

mdf_scaled <- mdf %>% 
 #   filter(scientific_name %in% threeSites$scientific_name) %>% 
    mutate(Dev_1 = scale(Dev_1),
           Dev_10 = scale(Dev_10),
           temp_niche = scale(temp_niche),
           mean_temp = scale((mean_temp - 1.02) * -1),
           BodySize = scale(BodySize))

traits <- read.csv("Data/TraitData.csv") %>% 
    select(Species, Order)

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
A <- ape::vcv.phylo(tt)

mdf_phylo <- mdf_phylo %>% 
    mutate(phyloName = scientific_name)

model1 <- brm(formula = bf(abundance ~ Dev_1 +
                      LHSCategory + LarvalHabitatCategory +
                      BodySize + temp_niche +
                      LHSCategory:Dev_1 + LarvalHabitatCategory:Dev_1 + 
                      BodySize:Dev_1 + temp_niche:Dev_1 + 
                      (1 | scientific_name) + (1|gr(phyloName, cov = A))),
                    family = negbinomial(),
                  data = mdf_phylo,
                  data2 = list(A=A),
                  chains = 4, iter = 2400, warmup = 1000,
                  control = list(adapt_delta = 0.99),
                  cores = 4, seed = 1234, 
                  threads = threading(2),
                  backend = "cmdstanr")

summary(model1)

# examine model assumptions
pp_check(model1)
pp_check(model1, type = "stat", stat = "mean")

conditional_effects(x = model1, effects = "Dev_1:LHSCategory")
conditional_effects(x = model1, effects = "Dev_1")

summary(model1)
m_sum <- summary(model1) 
fixed <- m_sum$fixed %>% 
    tibble::rownames_to_column()
random_phy <- m_sum$random$phyloName %>% 
    tibble::rownames_to_column()
random_species <- m_sum$random$scientific_name %>% 
    tibble::rownames_to_column()

m_sum_allSites <- bind_rows(fixed, random_phy, random_species)

# tab outputs
write.csv(m_sum_allSites, "Tables/speciesSpecificResults.csv", row.names = F)
