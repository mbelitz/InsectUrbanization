library(tidyverse)

mdf <- read.csv("Data/modelData_allSpecies.csv")

head(mdf)

sppSpecificTable <- mdf %>% 
    select(Site, scientific_name, abundance) %>% 
    pivot_wider(names_from = Site,
                values_from = abundance)

write.csv(x = sppSpecificTable, file = "Tables/speciesAbundancePerSite.csv", row.names = FALSE)

# now order specific
mdf <- read.csv("data/modelData_allSpecies.csv")

mdf_scaled <- mdf %>% 
    mutate(Dev_1 = scale(Dev_1))

traits <- read.csv("data/TraitData.csv")

mdf_scaled <- mdf_scaled %>% left_join(traits, by = c("scientific_name" = "Species"))

Order_abund <- mdf_scaled %>% 
    group_by(Order, Site) %>% 
    summarise(totalAbund = sum(abundance))

orderSpecificTable <- Order_abund %>% 
    pivot_wider(names_from = Site,
                values_from = totalAbund)

write.csv(orderSpecificTable, file = "Tables/orderAbundancePerSite.csv", row.names = FALSE)