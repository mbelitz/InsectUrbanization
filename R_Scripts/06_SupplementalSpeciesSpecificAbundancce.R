library(tidyverse)
library(brms)
library(ape)
library(cmdstanr)
library(tidybayes)
set_cmdstan_path()


#Script to run Species specific model without Hydaticus bimarginatus in the dataset
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

# remove most abundant species
mdf_scaled <- mdf_scaled %>% filter(scientific_name != "Hydaticus bimarginatus")

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
    mutate(phyloName = scientific_name)

#tt <- ape::drop.tip(tt, tip = sppNotInAnlysis$validName)
A <- ape::vcv.phylo(tt)
#
#mdf_phylo <- mdf_phylo %>% 
#    mutate(phyloName = scientific_name)

model1 <- brm(formula = bf(abundance ~ Dev_1 +
                               LHSCategory + LarvalHabitatCategory +
                               BodySize + temp_niche +
                               LHSCategory:Dev_1 + LarvalHabitatCategory:Dev_1 + 
                               BodySize:Dev_1 + temp_niche:Dev_1 + 
                               (1|scientific_name) + 
                               (0 + Dev_1|scientific_name)),
              family = negbinomial(),
              data = mdf_phylo,
              data2 = list(A=A),
              chains = 4, iter = 2400, warmup = 1000,
              control = list(adapt_delta = 0.99),
              cores = 4, seed = 1234, 
              threads = threading(2),
              backend = "cmdstanr")

summary(model1)

t <- pp_check(model1, ndraws = 50)
#ggsave(filename = "Figures/Supplemental/SpeciesSpecificCheck.png", plot = t)
pp_check(model1, type = "stat", stat = "mean")
#ggsave(filename = "Figures/Supplemental/SpeciesSpecificCheckHist.png")

## examine output
summary(model1)

m_sum <- summary(model1) 
fixed <- m_sum$fixed %>% 
    tibble::rownames_to_column()
random <- ranef(model1)$Order %>% 
    unlist() %>% as.data.frame()

# plot random and conditional effects
mycols <- viridisLite::turbo(8, begin = 0.05, end = 0.90)
c("#3E378FFF", "#458EFEFF", "#18DAC7FF", "#69FD66FF", "#CDEC34FF", 
  "#FEAA33FF", "#E94D0CFF", "#A31301FF")
colvals <-  c(
    "Coleoptera"  = mycols[1],
    "Diptera"     = mycols[2],
    "Hemiptera"   = mycols[3],
    "Hymenoptera" = mycols[4],
    "Megaloptera" = mycols[5],
    "Blattodea"   = mycols[6],
    "Neuroptera"  = mycols[7],
    "Orthoptera"  = mycols[8]
)

draw <- model1 %>%
    spread_draws(b_Dev_1, r_scientific_name[scientific_name,Dev_1]) %>%
    # add the grand mean to the group-specific deviations
    mutate(mu = b_Dev_1 + r_scientific_name) %>%
    ungroup() %>%
    mutate(Species = str_replace_all(scientific_name, "[.]", " ")) %>% 
    left_join(traits)

ggplot(draw, aes(x = mu, y = reorder(scientific_name, mu), color = Order)) +
    geom_vline(xintercept = fixef(model1)[2, 1], color = "#839496", linewidth = 1) +
    geom_vline(xintercept = fixef(model1)[2, 3:4], color = "#839496", linetype = "dotted", linewidth = 0.75) +
    geom_vline(xintercept = 0) +
    stat_halfeye(.width = c(0.025,0.975), fill = NA, linewidth = 2/3, alpha = 0.9, point_interval = "mean_qi") +
    scale_color_manual(values = colvals) +
    # scale_fill_manual(values = c(rep("transparent", 8))) +
    labs(x = "Sensitivity to urban development", y = "") +
    scale_x_continuous(limits = c(-6,6)) +
    theme_classic() +
    theme(legend.position = c(0.15,0.75)) 

ggsave(filename = "Figures/Supplemental_SpeciesResponseToUrb.png", dpi = 450,
       width = 8.5, height = 11)

# now let's plot predictor posterior distributions

draw1 <- model1 %>% 
    spread_draws(b_Dev_1,                                              
                 b_LHSCategoryDetrivore,                               
                 b_LHSCategoryHerbivore,                                
                 b_LHSCategoryOmnivore,                                
                 b_LarvalHabitatCategoryBelowGround,                   
                 b_LarvalHabitatCategoryFreshwater,                    
                 b_LarvalHabitatCategoryHostOrganisms,                  
                 b_BodySize,                  
                 b_temp_niche,                                          
                 `b_Dev_1:LHSCategoryDetrivore`,                          
                 `b_Dev_1:LHSCategoryHerbivore`,                          
                 `b_Dev_1:LHSCategoryOmnivore`,                          
                 `b_Dev_1:LarvalHabitatCategoryBelowGround`,           
                 `b_Dev_1:LarvalHabitatCategoryFreshwater`,              
                 `b_Dev_1:LarvalHabitatCategoryHostOrganisms`,            
                 `b_Dev_1:BodySize`,            
                 `b_Dev_1:temp_niche`)

draw1_long <- draw1 %>% 
    pivot_longer(cols = ! c(`.draw`, `.chain`, `.iteration`), 
                 names_to = "predictor", values_to = "estimate") %>% 
    mutate(predictor_type = case_when( predictor == "b_Dev_1" ~ "Urban development",                                              
                                       predictor == "b_LHSCategoryDetrivore" ~ "Trait",                               
                                       predictor == "b_LHSCategoryHerbivore"~ "Trait",                                
                                       predictor == "b_LHSCategoryOmnivore" ~ "Trait",                                
                                       predictor == "b_LarvalHabitatCategoryBelowGround"~ "Trait",                   
                                       predictor == "b_LarvalHabitatCategoryFreshwater"~ "Trait",                    
                                       predictor == "b_LarvalHabitatCategoryHostOrganisms"~ "Trait",                  
                                       predictor == "b_BodySize"~ "Trait",                  
                                       predictor == "b_temp_niche"~ "Trait",                                          
                                       predictor == "b_Dev_1:LHSCategoryDetrivore" ~ "Urban development",                          
                                       predictor == "b_Dev_1:LHSCategoryHerbivore"~ "Urban development",                          
                                       predictor == "b_Dev_1:LHSCategoryOmnivore"~ "Urban development",                          
                                       predictor == "b_Dev_1:LarvalHabitatCategoryBelowGround"~ "Urban development",           
                                       predictor == "b_Dev_1:LarvalHabitatCategoryFreshwater"~ "Urban development",              
                                       predictor == "b_Dev_1:LarvalHabitatCategoryHostOrganisms"~ "Urban development",            
                                       predictor == "b_Dev_1:BodySize"~ "Urban development",            
                                       predictor == "b_Dev_1:temp_niche"~ "Urban development")) %>% 
    
    mutate(predictor_rename = case_when( predictor == "b_Dev_1" ~ "Dev",                                              
                                         predictor == "b_LHSCategoryDetrivore" ~ "Detrivore",                               
                                         predictor == "b_LHSCategoryHerbivore"~  "Herbivore",                                
                                         predictor == "b_LHSCategoryOmnivore" ~  "Omnivore",                                
                                         predictor == "b_LarvalHabitatCategoryBelowGround"~ "LH_BelowGround",                   
                                         predictor == "b_LarvalHabitatCategoryFreshwater"~ "LH_Freshwater",                    
                                         predictor == "b_LarvalHabitatCategoryHostOrganisms"~ "LH_HostOrganisms",                  
                                         predictor == "b_BodySize"~ "BodySize",                  
                                         predictor == "b_temp_niche"~ "TempNiche",                                          
                                         predictor == "b_Dev_1:LHSCategoryDetrivore" ~ "Dev:Detrivore",                          
                                         predictor == "b_Dev_1:LHSCategoryHerbivore"~ "Dev:Herbivore",                          
                                         predictor == "b_Dev_1:LHSCategoryOmnivore"~ "Dev:Omnivore",                          
                                         predictor == "b_Dev_1:LarvalHabitatCategoryBelowGround"~ "Dev:LH_BelowGround",           
                                         predictor == "b_Dev_1:LarvalHabitatCategoryFreshwater"~ "Dev:LH_Freshwater",              
                                         predictor == "b_Dev_1:LarvalHabitatCategoryHostOrganisms"~ "Dev:LH_HostOrganisms",            
                                         predictor == "b_Dev_1:BodySize"~ "Dev:BodySize",            
                                         predictor == "b_Dev_1:temp_niche"~ "Dev:TempNiche")) %>% 
    mutate(Sig = case_when(predictor_rename == "Herbivore" | predictor_rename == "LH_Freshwater" ~ "Sig",
                           .default = "Not sig"))

ggplot(draw1_long, aes(x = estimate, y = reorder(predictor_rename, estimate), fill = Sig)) +
    stat_halfeye(.width = c(0.025,0.975), linewidth = 2/3, alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "Coefficient estimate", y = "") +
    scale_fill_manual(values = c("grey40","purple")) +
    facet_wrap(~predictor_type, scales = "free") +
    theme_classic() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"))

ggsave(filename = "Figures/SupplementalSppSpecificResponses.png", dpi = 450,
       width = 8, height = 5)

# Table of results
summary(model1)
m_sum <- summary(model1) 
fixed <- m_sum$fixed %>% 
    tibble::rownames_to_column()
random_species <- m_sum$random$scientific_name %>% 
    tibble::rownames_to_column()

m_sum_allSites <- bind_rows(fixed, random_species)

# tab outputs
write.csv(m_sum_allSites, "Tables/Supplemental_speciesSpecificResults.csv", row.names = F)
