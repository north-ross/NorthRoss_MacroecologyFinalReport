library(dplyr)
df <- read_csv("Koster_2023_data_2021_Vascular_plant_biomass.csv")

# Get SR per site and treatment
df$treat_site <- paste(df$Treatment, df$Site, sep="_")
df %>%
  group_by(treat_site) %>%  
  summarise(SR = n_distinct(Species), .groups = 'drop')
# group by plot
plt_SR = df %>%
  group_by(treat_site_plt) %>%  
  summarise(SR = n_distinct(Species), .groups = 'drop')

# get mean, SD, n plots   by treat_site
df %>%
  group_by(treat_site, Plot) %>%
  summarise(SR = n_distinct(Species), .groups = 'drop') %>%
  group_by(treat_site) %>%
  summarise(
    n_plots = n(),
    mean_SR = mean(SR),
    sd_SR = sd(SR),
    .groups = 'drop'
  )

# aggregate FGs
dict = list(
  "Deciduous shrubs" = 'Shrubs', 
  "Evergreen shrubs" = 'Shrubs', 
  "Forbs" = 'Herbs',
  "Graminoids" = "Herbs",
  "Tree seedlings" = "Trees")
dc = df$Fgroup
for (i in 1:length(dict)){dc <- replace(dc, dc == names(dict[i]), dict[[i]])}
df$FG <- dc
table(df$FG)

df$treat_site_FG <- paste(df$treat_site, df$FG, sep="_")
# Report biomass by FGs and category
df %>%
  group_by(treat_site_FG, Plot) %>%
  summarise(BM = mean(Max_BM), .groups = 'drop') %>%
  group_by(treat_site_FG) %>%
  summarise(
    mean_BM = mean(BM),
    n_plots = n(),
    sd_BM = sd(BM),
    .groups = 'drop'
  )
