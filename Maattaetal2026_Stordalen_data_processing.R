###############################################################################################################################################
### Code for: 
### Määttä et al. (2026) Plant belowground traits indicate increased plant-mediated methane transport along a peatland permafrost thaw gradient
### Code created by Tiia Määttä, with parts written with gpt 4, 4o and 5.4
###############################################################################################################################################

### open libraries

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(stringr)
library(lubridate)
library(patchwork)
library(mgcv)
library(tibble)
library(zoo)
library(readr)
library(myClim)


#############################################################
###---------------VEGETATION SURVEY DATA------------------###
#############################################################

### Rubus Chamaemorus area calculations

# load the dataset
photo_data <- read.csv("path/rubus_imagej.csv")

# calculate widthtolength ratios for each diamond
photo_data$R1 <- photo_data$w1 / photo_data$r1
photo_data$R2 <- photo_data$w2 / photo_data$r2
photo_data$R3 <- photo_data$w3 / photo_data$r3

# compute the average ratio across all 50 leaves
avg_ratios <- colMeans(photo_data[, c("R1", "R2", "R3")], na.rm = TRUE)

# print the average ratios
print(avg_ratios)

### import the vga measurements csv

vga <- read.csv("path/VGA_measurements.csv")

# assign pft classes for species

vga <- vga %>%
  mutate(
    aer_class = case_when(
      species %in% c("erivag", "eriang", "equflu", "carrot", "galpal", "carcan") ~ "aerenchymatous",
      TRUE ~ "nonaerenchymatous"
    ),
    PFT_class = case_when(
      species %in% c("erivag", "eriang", "equflu", "drorot", "carrot", "galpal", "compal", "carcan") ~ "herbaceous",
      TRUE ~ "shrub"
    )
  )

# create a new df where leaf and stem areas for each species are shown in separate columns, and then add similar columns for the veg classes

# first, create new df and take out the speciesspecific data 
df_wide <- vga %>%
  distinct(plot_ID, species, .keep_all = TRUE) %>% 
  pivot_wider(
    names_from = species,
    values_from = c(avg.no.of.leaves.per.m2.per.plot, avg.no.of.stems.per.m2.per.plot, 
                    species.avg.leaf.area.per.plot..m2, species.avg.stem.area.per.plot..m2),
    names_glue = "{species}_{.value}"
  )

# selecting columns
final_df <- df_wide %>%
  select(date, plot_ID, chamber_ID, thaw_stage, starts_with("eriang"), starts_with("erivag"), 
         starts_with("carrot"), starts_with("galpal"), starts_with("carcan"), 
         starts_with("equflu"), starts_with("compal"), starts_with("vaculi"), starts_with("vacvit"),
         starts_with("vacoxy"), starts_with("betnan"), starts_with("rubcha"), starts_with("empnig"), 
         starts_with("andpol"), starts_with("drorot"))

final_df <- as.data.frame(final_df)

# remove rows where all speciesspecific columns are na
final_df <- final_df %>%
  filter(if_any(-c(date, plot_ID, chamber_ID, thaw_stage), ~ !is.na(.)))

# replace "." with "_" and ".." with "_" in relevant column names
colnames(final_df) <- gsub("\\.\\.", "_", colnames(final_df))
colnames(final_df) <- gsub("\\.", "_", colnames(final_df))

# convert row.leaf.area and row.stem.area from cm2 to m2 and calculate mean and sd per aer_class group
vga_summary <- vga %>%
  mutate(
    row_leaf_area_m2 = row.leaf.area / 10000,
    row_stem_area_m2 = row.stem.area / 10000
  ) %>%
  group_by(plot_ID, aer_class) %>%
  summarise(
    mean_leaf_area = mean(row_leaf_area_m2, na.rm = TRUE),
    sd_leaf_area = sd(row_leaf_area_m2, na.rm = TRUE),
    mean_stem_area = mean(row_stem_area_m2, na.rm = TRUE),
    sd_stem_area = sd(row_stem_area_m2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = aer_class,
    values_from = c(mean_leaf_area, sd_leaf_area, mean_stem_area, sd_stem_area),
    names_glue = "{aer_class}_{.value}_per_plot_m2"
  )

# merge with final_df
final_df <- final_df %>%
  left_join(vga_summary, by = "plot_ID")


# do the same for pft_class

vga_summary_2 <- vga %>%
  mutate(
    row_leaf_area_m2 = row.leaf.area / 10000,
    row_stem_area_m2 = row.stem.area / 10000
  ) %>%
  group_by(plot_ID, PFT_class) %>%
  summarise(
    mean_leaf_area = mean(row_leaf_area_m2, na.rm = TRUE),
    sd_leaf_area = sd(row_leaf_area_m2, na.rm = TRUE),
    mean_stem_area = mean(row_stem_area_m2, na.rm = TRUE),
    sd_stem_area = sd(row_stem_area_m2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = PFT_class,
    values_from = c(mean_leaf_area, sd_leaf_area, mean_stem_area, sd_stem_area),
    names_glue = "{PFT_class}_{.value}_per_plot_m2"
  )

# merge with final_df
final_df <- final_df %>%
  left_join(vga_summary_2, by = "plot_ID")

# replace specific words in column names
colnames(final_df) <- gsub("nonaerenchymatous", "nonaer", colnames(final_df))
colnames(final_df) <- gsub("aerenchymatous", "aer", colnames(final_df))
colnames(final_df) <- gsub("herbaceous", "herb", colnames(final_df))

# calculate average number of leaves and stems per m2 per plot for the classes
totleaves_stems_summary <- vga %>%
  distinct(plot_ID, RQ_ID, species, aer_class, .keep_all = TRUE) %>% 
  group_by(plot_ID, aer_class) %>%
  summarise(
    mean_totleaves = mean(totleaves, na.rm = TRUE),
    sd_totleaves = sd(totleaves, na.rm = TRUE),
    mean_totstems = mean(totstems, na.rm = TRUE),
    sd_totstems = sd(totstems, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = aer_class,
    values_from = c(mean_totleaves, sd_totleaves, mean_totstems, sd_totstems),
    names_glue = "{aer_class}_{.value}_per_plot"
  )

# merge with final_df
final_df <- final_df %>%
  left_join(totleaves_stems_summary, by = "plot_ID")

# do the same for pft_class
totleaves_stems_summary_2 <- vga %>%
  distinct(plot_ID, RQ_ID, species, PFT_class, .keep_all = TRUE) %>%  
  group_by(plot_ID, PFT_class) %>%
  summarise(
    mean_totleaves = mean(totleaves, na.rm = TRUE),
    sd_totleaves = sd(totleaves, na.rm = TRUE),
    mean_totstems = mean(totstems, na.rm = TRUE),
    sd_totstems = sd(totstems, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = PFT_class,
    values_from = c(mean_totleaves, sd_totleaves, mean_totstems, sd_totstems),
    names_glue = "{PFT_class}_{.value}_per_plot"
  )

# merge with final_df
final_df <- final_df %>%
  left_join(totleaves_stems_summary_2, by = "plot_ID")

# rename some columns
colnames(final_df) <- gsub("nonaerenchymatous", "nonaer", colnames(final_df))
colnames(final_df) <- gsub("aerenchymatous", "aer", colnames(final_df))
colnames(final_df) <- gsub("herbaceous", "herb", colnames(final_df))

# calculate green area (ga)

species_list <- unique(vga$species, na.rm=T)
species_list <- species_list[!is.na(species_list)]

for (species in species_list) {
  leaf_col <- paste0(species, "_species_avg_leaf_area_per_plot_m2")
  leaf_num_col <- paste0(species, "_avg_no_of_leaves_per_m2_per_plot")
  stem_col <- paste0(species, "_species_avg_stem_area_per_plot_m2")
  stem_num_col <- paste0(species, "_avg_no_of_stems_per_m2_per_plot")
  ga_col <- paste0(species, "_GA_m2m2")
  
  # check if all required columns exist in final_df
  existing_cols <- c(leaf_col, leaf_num_col, stem_col, stem_num_col) %in% colnames(final_df)
  
  if (all(existing_cols)) {
    final_df[[ga_col]] <- (final_df[[leaf_col]] * final_df[[leaf_num_col]]) +
      (final_df[[stem_col]] * final_df[[stem_num_col]])
  } else {
    missing_cols <- c(leaf_col, leaf_num_col, stem_col, stem_num_col)[!existing_cols]
    message("Skipping ", species, " due to missing columns: ", paste(missing_cols, collapse = ", "))
  }
}

# remove columns that include "_ga" but not "_ga_m2m2"
final_df <- final_df %>%
  select(-which(grepl("_GA", colnames(final_df)) & !grepl("_GA_m2m2", colnames(final_df))))

final_df_long <- final_df %>%
  pivot_longer(
    cols = matches("_GA_m2m2|_species_avg_leaf_area_per_plot_m2|_species_avg_stem_area_per_plot_m2|_avg_no_of_leaves_per_m2_per_plot|_avg_no_of_stems_per_m2_per_plot"),
    names_to = c("species", ".value"),
    names_pattern = "(.*?)_(GA_m2m2|species_avg_leaf_area_per_plot_m2|species_avg_stem_area_per_plot_m2|avg_no_of_leaves_per_m2_per_plot|avg_no_of_stems_per_m2_per_plot)"
  ) %>%
  rename(
    species_GA_m2m2 = GA_m2m2,
    species_avg_leaf_area_m2 = species_avg_leaf_area_per_plot_m2,
    species_avg_stem_area_m2 = species_avg_stem_area_per_plot_m2,
    species_avg_no_leaves_per_m2 = avg_no_of_leaves_per_m2_per_plot,
    species_avg_no_stems_per_m2 = avg_no_of_stems_per_m2_per_plot
  )

final_df_long <- as.data.frame(final_df_long)

# remove rows where key species-specific columns are na
final_df_long <- final_df_long %>%
  filter(!if_all(c(species_avg_no_stems_per_m2, species_avg_leaf_area_m2, species_avg_stem_area_m2, species_GA_m2m2), is.na))

# update ga calculation based on availability of leaf or stem area
final_df_long <- final_df_long %>%
  mutate(
    species_GA_m2m2 = case_when(
      !is.na(species_avg_leaf_area_m2) & !is.na(species_avg_stem_area_m2) ~
        (species_avg_leaf_area_m2 * species_avg_no_leaves_per_m2) + 
        (species_avg_stem_area_m2 * species_avg_no_stems_per_m2),
      
      !is.na(species_avg_leaf_area_m2) & is.na(species_avg_stem_area_m2) ~
        species_avg_leaf_area_m2 * species_avg_no_leaves_per_m2,
      
      is.na(species_avg_leaf_area_m2) & !is.na(species_avg_stem_area_m2) ~
        species_avg_stem_area_m2 * species_avg_no_stems_per_m2,
      
      TRUE ~ NA_real_
    )
  )

### there are 5 rows where ga is na. there are some missing values which should not be missing. inserting them from excel

# palsa3chamber vaculi leaf area: 0.00002404526685
final_df_long$species_avg_leaf_area_m2[final_df_long$species == "vaculi" & final_df_long$plot_ID == "palsa3chamber" & is.na(final_df_long$species_avg_leaf_area_m2)] <- 0.00002404526685

# palsa3chamber betnan leaf area: 0.00002823506397
final_df_long$species_avg_leaf_area_m2[final_df_long$species == "betnan" & final_df_long$plot_ID == "palsa3chamber" & is.na(final_df_long$species_avg_leaf_area_m2)] <- 0.00002823506397

# bog4core drorot: leaf area: 0.00001822123739 (so not na!)
final_df_long$species_avg_leaf_area_m2[final_df_long$species == "drorot" & final_df_long$plot_ID == "bog4core" & is.na(final_df_long$species_avg_leaf_area_m2)] <- 0.00001822123739

# bog4core empnig: no of leaves: 16217.42036, avg no of stems: 415.8312913
final_df_long$species_avg_no_leaves_per_m2[final_df_long$species == "empnig" & final_df_long$plot_ID == "bog4core" & is.na(final_df_long$species_avg_no_leaves_per_m2)] <- 16217.42036

# palsa3core vaculi: leaf area: 0.000006577709618
final_df_long$species_avg_leaf_area_m2[final_df_long$species == "vaculi" & final_df_long$plot_ID == "palsa3core" & is.na(final_df_long$species_avg_leaf_area_m2)] <- 0.000006577709618

### calculate ga for nonaer, aer, herb and shrub again

calc_ga_m2m2 <- function(leaf_area, stem_area, leaf_n, stem_n) {
  case_when(
    !is.na(leaf_area) & !is.na(stem_area) ~ (leaf_area * leaf_n) + (stem_area * stem_n),
    !is.na(leaf_area) & is.na(stem_area) ~ leaf_area * leaf_n,
    is.na(leaf_area) & !is.na(stem_area) ~ stem_area * stem_n,
    TRUE ~ NA_real_
  )
}

final_df_long <- final_df_long %>%
  mutate(
    nonaer_GA_m2m2 = calc_ga_m2m2(nonaer_mean_leaf_area_per_plot_m2, nonaer_mean_stem_area_per_plot_m2,
                                  nonaer_mean_totleaves_per_plot, nonaer_mean_totstems_per_plot),
    aer_GA_m2m2 = calc_ga_m2m2(aer_mean_leaf_area_per_plot_m2, aer_mean_stem_area_per_plot_m2,
                               aer_mean_totleaves_per_plot, aer_mean_totstems_per_plot),
    herb_GA_m2m2 = calc_ga_m2m2(herb_mean_leaf_area_per_plot_m2, herb_mean_stem_area_per_plot_m2,
                                herb_mean_totleaves_per_plot, herb_mean_totstems_per_plot),
    shrub_GA_m2m2 = calc_ga_m2m2(shrub_mean_leaf_area_per_plot_m2, shrub_mean_stem_area_per_plot_m2,
                                 shrub_mean_totleaves_per_plot, shrub_mean_totstems_per_plot)
  )

# change column names to stem density for shrubs, herbs, nonaer and aer

final_df_long <- final_df_long %>%
  rename(
    aer_stem_density_nom2 = aer_mean_totstems_per_plot,
    nonaer_stem_density_nom2 = nonaer_mean_totstems_per_plot,
    herb_stem_density_nom2 = herb_mean_totstems_per_plot,
    shrub_stem_density_nom2 = shrub_mean_totstems_per_plot,
    aer_stem_density_sd = aer_sd_totstems_per_plot,
    nonaer_stem_density_sd = nonaer_sd_totstems_per_plot,
    herb_stem_density_sd = herb_sd_totstems_per_plot,
    shrub_stem_density_sd = shrub_sd_totstems_per_plot
  )

# make a compact version which can be combined with fch4 data

vga_data <- subset(final_df_long, select=c(date, plot_ID, chamber_ID, thaw_stage, species, 
                                           species_GA_m2m2, species_avg_no_stems_per_m2,
                                           aer_stem_density_nom2, nonaer_stem_density_nom2, 
                                           herb_stem_density_nom2, shrub_stem_density_nom2,
                                           aer_stem_density_sd, nonaer_stem_density_sd,
                                           herb_stem_density_sd, shrub_stem_density_sd,
                                           aer_GA_m2m2, nonaer_GA_m2m2, herb_GA_m2m2,
                                           shrub_GA_m2m2
))

vga_data <- vga_data %>%
  rename(
    species_stem_density_nom2 = species_stem_density
  )

# there are some nan > change to na

vga_data <- vga_data %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))

# add aer_class and pft_class
vga_data <- vga_data %>%
  mutate(
    aer_class = case_when(
      species %in% c("erivag", "eriang", "equflu", "carrot", "galpal", "carcan") ~ "aerenchymatous",
      TRUE ~ "nonaerenchymatous"
    ),
    PFT_class = case_when(
      species %in% c("erivag", "eriang", "equflu", "drorot", "carrot", "galpal", "compal", "carcan") ~ "herbaceous",
      TRUE ~ "shrub"
    )
  )

# move next to species

vga_data <- vga_data %>%
  relocate(aer_class, PFT_class, .after = species)

# save as csv

write.csv(vga_data, "path/VGA_data.csv", row.names = FALSE)

# for Table C2, calculate mean ga and sd per species

species_ga_summary <- vga_data %>%
  group_by(thaw_stage, species) %>%
  summarise(
    mean_ga = mean(species_GA_m2m2, na.rm = TRUE),
    sd_ga = sd(species_GA_m2m2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(thaw_stage, species)

species_ga_summary <- as.data.frame(species_ga_summary)

### point intercept species survey data to get species % cover per plot

survey <- read.csv("path/speciessurvey.csv")

# remove empty rows
survey <- subset(survey, species != "")

survey <- subset(survey, select=-notes)

# create a new version of that df where, instead of number of hits, you only have 1 for one or more hits and 0 for no hits

## first check if there are nas in the hit columns

colSums(is.na(survey[, c("P1_hit", "P2_hit", "P3_hit",
                         "P4_hit", "P5_hit", "P6_hit",
                         "P7_hit", "P8_hit",
                         "P9_hit", "P10_hit", "P11_hit",
                         "P12_hit", "P13_hit", "P14_hit",
                         "P15_hit", "P16_hit", "P17_hit",
                         "P18_hit", "P19_hit", "P20_hit",
                         "P21_hit", "P22_hit", "P23_hit",
                         "P24_hit", "P25_hit")]))
# no NAs

survey_2 <- survey %>%
  mutate(across(ends_with("_hit"), ~ ifelse(. > 0, 1, .)))

survey_2 <- survey_2 %>%
  mutate(
    aer_class = case_when(
      species %in% c("erivag", "eriang", "equflu", "carrot", "galpal") ~ "aerenchymatous",
      species %in% c("sphrip", "sphcap", "polsp", "dicsp", "sphbal", "pticil", "scosp") ~ "moss",
      TRUE ~ "nonaerenchymatous"
    ),
    PFT_class = case_when(
      species %in% c("erivag", "eriang", "equflu", "drorot", "carrot", "galpal", "compal") ~ "herbaceous",
      species %in% c("sphrip", "sphcap", "polsp", "dicsp", "sphbal", "pticil", "scosp") ~ "moss",
      TRUE ~ "shrub"
    )
  )

# calculate species cover (%)
survey_2 <- survey_2 %>%
  rowwise() %>%
  mutate(
    species_hits = sum(c_across(ends_with("_hit"))), 
    species_cover = (species_hits / 25) * 100
  ) %>%
  ungroup() %>%
  group_by(plot_ID) %>%
  mutate(
    aer_cover = sum(ifelse(aer_class == "aerenchymatous", species_hits, 0)) / 25 * 100,
    nonaer_cover = sum(ifelse(aer_class == "nonaerenchymatous", species_hits, 0)) / 25 * 100,
    herb_cover = sum(ifelse(PFT_class == "herbaceous", species_hits, 0)) / 25 * 100,
    shrub_cover = sum(ifelse(PFT_class == "shrub", species_hits, 0)) / 25 * 100,
    moss_cover = sum(ifelse(class == "moss", species_hits, 0)) / 25 * 100,
    vasc_cover = sum(ifelse(class == "vasc", species_hits, 0)) / 25 * 100
  ) %>%
  ungroup()

## remove the extra columns

survey_2 <- survey_2[, !grepl("_hit|_hgt", names(survey_2))]

# save

write.csv(survey_2, "path/species_cover.csv", row.names = FALSE)

## calculate mean species cover per thaw stage

t <- survey_2 %>%
  mutate(plot_type = ifelse(grepl("chamber$", plot_ID), "chamber", 
                            ifelse(grepl("core$", plot_ID), "core", NA))) %>% 
  group_by(thaw_stage, plot_type, species) %>%
  summarise(mean_species_cover = mean(species_cover, na.rm = TRUE), .groups = "drop") %>%
  arrange(thaw_stage, plot_type, desc(mean_species_cover))  

t <- as.data.frame(t)
t


# per herb and shrub
PFT_cover_summary <- survey_2 %>%
  mutate(plot_type = ifelse(grepl("chamber$", plot_ID), "chamber", 
                            ifelse(grepl("core$", plot_ID), "core", NA))) %>%
  distinct(plot_ID, .keep_all = TRUE) %>% 
  group_by(thaw_stage, plot_type) %>%
  summarise(
    mean_herb_cover = mean(herb_cover, na.rm = TRUE),
    mean_shrub_cover = mean(shrub_cover, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(thaw_stage, plot_type) 

PFT_cover_summary <- as.data.frame(PFT_cover_summary)


# calculate mean ga per thaw stage per core and chamber

PFT_ga_summary <- vga_data %>%
  mutate(plot_type = ifelse(grepl("chamber$", plot_ID), "chamber", 
                            ifelse(grepl("core$", plot_ID), "core", NA))) %>%  
  distinct(plot_ID, .keep_all = TRUE) %>% 
  group_by(thaw_stage, plot_type) %>%
  summarise(
    mean_herb_ga = mean(herb_GA_m2m2, na.rm = TRUE),
    mean_shrub_ga = mean(shrub_GA_m2m2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(thaw_stage, plot_type)

PFT_ga_summary <- as.data.frame(PFT_ga_summary)


# calculate across all species ("total")

ga <- vga_data %>%
  mutate(plot_type = ifelse(grepl("chamber$", plot_ID), "chamber", 
                            ifelse(grepl("core$", plot_ID), "core", NA))) %>% 
  group_by(thaw_stage, plot_type) %>%
  summarise(mean_ga = mean(species_GA_m2m2, na.rm = TRUE), .groups = "drop") %>%
  arrange(thaw_stage, plot_type, desc(mean_ga)) 

ga <- as.data.frame(ga)
ga

#############################################################
###-----------------ABIOTIC VARIABLES---------------------###
#############################################################

### variables from manual measurements from both chamber and core plots

# add mean soil temp per plot

temp <- read.csv("path/soiltemp_data.csv")

temp$datetime <- paste0(temp$datetime, ":00")

# convert datetime string to proper datetime format
temp$datetime <- dmy_hms(temp$datetime)

# extract only the date part
temp$date <- date(temp$datetime)

temp <- subset(temp, select = -datetime)

# take an average and sd per plot 

temp_ch <- temp %>%
  group_by(plot_ID) %>%
  summarise(
    mean_TS = mean(temp_10, na.rm = TRUE),
    sd_TS = sd(temp_10, na.rm = TRUE),
    .groups = "drop"
  )

temp_ch <- as.data.frame(temp_ch)
temp_ch

# create thaw_stage

temp <- temp %>%
  mutate(thaw_stage = recode(sub("([a-zA-Z]+).*", "\\1", plot_ID),
                              palsa = "intact",
                              bog = "partly_thawed",
                              fen = "fully_thawed"))

# do this but for chamber vs core

temp_summary <- temp %>%
  mutate(plot_type = ifelse(grepl("chamber$", plot_ID), "chamber", 
                            ifelse(grepl("core$", plot_ID), "core", NA))) %>% 
  group_by(thaw_stage, plot_type) %>%
  summarise(
    mean_TS = mean(temp_10, na.rm = TRUE),
    sd_TS = sd(temp_10, na.rm = TRUE),
    .groups = "drop"
  )

temp_summary <- as.data.frame(temp_summary)


#### combine mean ts, ph, mean species, herb and shrub ga, species, herb and shrub cover, and wtd for comparisons

comparison_data <- vga_data %>%
  left_join(temp %>% select(plot_ID, mean_TS, sd_TS), by = "plot_ID")

# combine with survey data
comparison_data_2 <- comparison_data %>%
  full_join(
    survey_2 %>% select(plot_ID, species, class, aer_class, PFT_class, species_cover, 
                        aer_cover, nonaer_cover, herb_cover, shrub_cover, moss_cover, vasc_cover),
    by = c("plot_ID", "species", "aer_class", "PFT_class")
  )

comparison_data_2 <- subset(comparison_data_2, select= -date)

comparison_data_2 <- comparison_data_2 %>%
  mutate(
    class = case_when(
      species %in% c("erivag", "eriang", "equflu", "carrot", "galpal", "drorot", "compal", "betnan", "vacvit", "carcan", "andpol",
                     "vacoxy", "vaculi", "rubcha", "empnig") ~ "vasc",
      species %in% c("sphrip", "sphcap", "polsp", "dicsp", "sphbal", "pticil", "scosp") ~ "moss"
    )
  )

# create thaw_stage
comparison_data_2 <- comparison_data_2 %>%
  mutate(thaw_stage = case_when(
    str_detect(plot_ID, "palsa") ~ "intact",
    str_detect(plot_ID, "bog") ~ "partly_thawed",
    str_detect(plot_ID, "fen") ~ "fully_thawed"
  ))

# clean plot_id
comparison_data_2 <- comparison_data_2 %>%
  mutate(chamber_ID = case_when(
    str_detect(plot_ID, "1") ~ 1,
    str_detect(plot_ID, "2") ~ 2,
    str_detect(plot_ID, "3") ~ 3,
    str_detect(plot_ID, "4") ~ 4,
    str_detect(plot_ID, "5") ~ 5,
    str_detect(plot_ID, "6") ~ 6,
    str_detect(plot_ID, "7") ~ 7,
    str_detect(plot_ID, "8") ~ 8,
    str_detect(plot_ID, "9") ~ 9
  ))

# arrange based on plot_id

comparison_data_2 <- comparison_data_2 %>%
  arrange(plot_ID)

# add ph, wtd and active layer data manually

comparison_data_2 <- comparison_data_2 %>% mutate(soil_pH =
                                                    case_when(comparison_data_2$plot_ID == "palsa1core" ~ 3.01, 
                                                              comparison_data_2$plot_ID == "palsa3core" ~ 2.98,
                                                              comparison_data_2$plot_ID == "palsa5core" ~ 3.16,
                                                              comparison_data_2$plot_ID == "palsa1chamber" ~ 3.1,
                                                              comparison_data_2$plot_ID == "palsa3chamber" ~ 3.1,
                                                              comparison_data_2$plot_ID == "palsa5chamber" ~ 2.83,
                                                              comparison_data_2$plot_ID == "bog2core" ~ 3.18,
                                                              comparison_data_2$plot_ID == "bog4core" ~ 3.78,
                                                              comparison_data_2$plot_ID == "bog6core" ~ 3.39,
                                                              comparison_data_2$plot_ID == "bog2chamber" ~ 3.3,
                                                              comparison_data_2$plot_ID == "bog4chamber" ~ 3.2,
                                                              comparison_data_2$plot_ID == "bog6chamber" ~ 3.3,
                                                              comparison_data_2$plot_ID == "fen7core" ~ 5.74,
                                                              comparison_data_2$plot_ID == "fen8core" ~ 6.02,
                                                              comparison_data_2$plot_ID == "fen9core" ~ 4.81,
                                                              comparison_data_2$plot_ID == "fen7chamber" ~ 5.1,
                                                              comparison_data_2$plot_ID == "fen8chamber" ~ 5,
                                                              comparison_data_2$plot_ID == "fen9chamber" ~ 4.7)
)


comparison_data_2 <- comparison_data_2 %>% mutate(WTD =
                                                    case_when(comparison_data_2$plot_ID == "palsa1core" ~ NA, 
                                                              comparison_data_2$plot_ID == "palsa3core" ~ NA,
                                                              comparison_data_2$plot_ID == "palsa5core" ~ NA,
                                                              comparison_data_2$plot_ID == "palsa1chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "palsa3chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "palsa5chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "bog2core" ~ -14.67,
                                                              comparison_data_2$plot_ID == "bog4core" ~ -21,
                                                              comparison_data_2$plot_ID == "bog6core" ~ -27,
                                                              comparison_data_2$plot_ID == "bog2chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "bog4chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "bog6chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "fen7core" ~ 0,
                                                              comparison_data_2$plot_ID == "fen8core" ~ 0,
                                                              comparison_data_2$plot_ID == "fen9core" ~ 0,
                                                              comparison_data_2$plot_ID == "fen7chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "fen8chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "fen9chamber" ~ NA)
)


comparison_data_2 <- comparison_data_2 %>% mutate(active_layer_cm =
                                                    case_when(comparison_data_2$plot_ID == "palsa1core" ~ 37.5, 
                                                              comparison_data_2$plot_ID == "palsa3core" ~ 39,
                                                              comparison_data_2$plot_ID == "palsa5core" ~ 85.5,
                                                              comparison_data_2$plot_ID == "palsa1chamber" ~ 47,
                                                              comparison_data_2$plot_ID == "palsa3chamber" ~ 41.5,
                                                              comparison_data_2$plot_ID == "palsa5chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "bog2core" ~ NA,
                                                              comparison_data_2$plot_ID == "bog4core" ~ NA,
                                                              comparison_data_2$plot_ID == "bog6core" ~ NA,
                                                              comparison_data_2$plot_ID == "bog2chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "bog4chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "bog6chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "fen7core" ~ NA,
                                                              comparison_data_2$plot_ID == "fen8core" ~ NA,
                                                              comparison_data_2$plot_ID == "fen9core" ~ NA,
                                                              comparison_data_2$plot_ID == "fen7chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "fen8chamber" ~ NA,
                                                              comparison_data_2$plot_ID == "fen9chamber" ~ NA)
)

# save

write.csv(comparison_data_2, "path/plot_env_comparison_data.csv", row.names = FALSE)


#############################################################
###----------ORG LITTER AND DEAD VEGETATION---------------###
#############################################################

litter_survey <- read.csv("path/litter_survey.csv")

# remove empty rows
litter_survey <- subset(litter_survey, cover_class != "")

litter_survey <- subset(litter_survey, select=-notes)

# check if there are NAs in the hit columns

colSums(is.na(litter_survey[, c("P1", "P2", "P3",
                                "P4", "P5", "P6",
                                "P7", "P8",
                                "P9", "P10", "P11",
                                "P12", "P13", "P14",
                                "P15", "P16", "P17",
                                "P18", "P19", "P20",
                                "P21", "P22", "P23",
                                "P24", "P25")]))

# no NAs

# calculate ground cover class cover (%)

p_cols <- paste0("P", 1:25)


litter_survey_2 <- litter_survey %>%
  rowwise() %>% 
  mutate(
    ground_hits = sum(c_across(all_of(p_cols)), na.rm = TRUE),
    ground_cover = (ground_hits / 25) * 100
  ) %>%
  ungroup()

litter_survey_2 <- as.data.frame(litter_survey_2)

# remove the extra columns

litter_survey_2 <- litter_survey_2 %>% select(-all_of(p_cols))

# remove ground_hits

litter_survey_2 <- subset(litter_survey_2, select=-ground_hits)
litter_survey_2 <- subset(litter_survey_2, select=-date)

# keep only litter and water

litter_survey_2 <- subset(litter_survey_2, cover_class %in% c("litter", "water"))

# make litter_cover into its own column and water as well

litter_survey_2 <- litter_survey_2 %>%
  pivot_wider(
    names_from = cover_class,  
    values_from = ground_cover,  
    names_glue = "{cover_class}_cover", 
    values_fill = list(ground_cover = 0)
  )

litter_survey_2 <- as.data.frame(litter_survey_2)

### deadveg

dead <- read.csv("path/deadveg_survey.csv")

# remove extra stuff

dead <- subset(dead, select=-c(notes, date))

dead <- dead[rowSums(!is.na(dead[, setdiff(names(dead), "plot_ID")])) > 0, ]

# calculate ground cover class cover (%)

dead <- dead %>%
  rowwise() %>% 
  mutate(
    dead_hits = sum(c_across(all_of(p_cols)), na.rm = TRUE),
    dead_cover = (dead_hits / 25) * 100
  ) %>%
  ungroup() 

dead <- as.data.frame(dead)

## remove the extra columns

dead <- dead %>% select(-all_of(p_cols))

# remove dead_hits

dead <- subset(dead, select=-dead_hits)

# combine with comparison_data_2

comparison_data_3 <- comparison_data_2 %>%
  left_join(dead, by = c("plot_ID", "chamber_ID"))

# add litter

comparison_data_4 <- comparison_data_3 %>%
  left_join(litter_survey_2 %>% select(plot_ID, chamber_ID, thaw_stage, litter_cover, water_cover), 
            by = c("plot_ID", "chamber_ID", "thaw_stage"))

# save

write.csv(comparison_data_4, "path/comparison_data.csv", row.names = FALSE)


#############################################################
###--------------ROOT AND RHIZOME TRAITS------------------###
#############################################################

root_rhizo <- read.csv("path/stordalen_root_data_raw_clean.csv")

# remove x from in front of the columns

names(root_rhizo) <- sub("^X", "", names(root_rhizo))

# calculate the percentage of fine roots (< 2mm in each tray and sample)
# need to modify the df so that tray has its own column and d_class_mm is a new column, flipped from the individual columns


root_rhizo_2 <- root_rhizo %>%
  mutate(across(matches("^(\\d|L|SA|V)"), as.character))

# vector of nontrait columns to exclude from pivoting
non_trait_cols <- c("sample_id", "global_tot_len_cm", "global_surface_area_cm2", 
                    "global_avg_d_mm", "global_vol_cm3")

# ensure all trait columns are character to avoid type mismatch
root_rhizo <- root_rhizo %>%
  mutate(across(-all_of(non_trait_cols), as.character))

# proceed with transformation
root_long <- root_rhizo %>%
  # extract 'tray' from sample_id
  mutate(tray = str_extract(sample_id, "tray\\d+")) %>%
  
  # pivot diameterclass trait columns to long format
  pivot_longer(
    cols = -all_of(c(non_trait_cols, "tray")),
    names_to = "raw_dclass",
    values_to = "value"
  ) %>%
  
  # identify trait type
  mutate(
    trait = case_when(
      str_detect(raw_dclass, "L\\d") | str_detect(raw_dclass, "^L\\d") ~ "tot_length_cm_per_d_class",
      str_detect(raw_dclass, "SA\\d") | str_detect(raw_dclass, "^SA\\d") ~ "surface_area_cm2_per_d_class",
      str_detect(raw_dclass, "V\\d") | str_detect(raw_dclass, "^V\\d") ~ "vol_d_class_cm3_per_d_class",
      TRUE ~ NA_character_
    )
  ) %>%
  
  # clean up diameter class label
  mutate(
    d_class_mm = raw_dclass %>%
      str_remove("L") %>%
      str_remove("SA") %>%
      str_remove("V") %>%
      str_replace("^(\\d+\\.\\d{3})(\\d+\\.\\d{3})$", "\\1-\\2") %>%
      str_replace("^(\\d+\\.\\d{3})$", ">\\1")
  ) %>%
  
  dplyr::select(-raw_dclass) %>%
  
  pivot_wider(names_from = trait, values_from = value)

root_long <- as.data.frame(root_long)

# change the tray column
root_long <- root_long %>%
  # remove "tray" from tray values
  mutate(tray = str_remove(tray, "tray")) %>%
  
  # move 'tray' column right after 'sample_id'
  relocate(tray, .after = sample_id)


# add pft
root_long <- root_long %>%
  mutate(PFT = case_when(
    str_detect(sample_id, "H") ~ "herbaceous",
    str_detect(sample_id, "S") ~ "shrub"
  ))


# add organ

root_long <- root_long %>%
  mutate(organ = case_when(
    str_detect(sample_id, "Root") ~ "root",
    str_detect(sample_id, "Rhizome") ~ "rhizome"
  ))

# add plot_id and phys_code

# create the phys_code lookup table
phys_lookup <- tribble(
  ~plot_id,     ~PHYS_code,
  "palsa1core",   "325B23TM1",
  "palsa1core",   "325B23TM2",
  "palsa1core",   "325B23TM3",
  "palsa1core",   "325B23TM4",
  "palsa1core",   "325B23TM5",
  "palsa1core",   "325B23TM6",
  "palsa3core",   "325B23TM9",
  "palsa3core",   "325B23TM10",
  "palsa3core",   "325B23TM11",
  "palsa3core",   "325B23TM12",
  "palsa5core",   "325B23TM13",
  "palsa5core",   "325B23TM14",
  "palsa5core",   "325B23TM15",
  "palsa5core",   "325B23TM16",
  "palsa5core",   "325B23TM17",
  "palsa5core",   "325B23TM18",
  "bog2core",     "325B23TM19",
  "bog2core",     "325B23TM20",
  "bog2core",     "325B23TM21",
  "bog2core",     "325B23TM22",
  "bog4core",     "325B23TM23",
  "bog4core",     "325B23TM24",
  "bog4core",     "325B23TM25",
  "bog4core",     "325B23TM26",
  "bog4core",     "325B23TM27",
  "bog4core",     "325B23TM28",
  "bog6core",     "325B23TM29",
  "bog6core",     "325B23TM30",
  "bog6core",     "325B23TM31",
  "bog6core",     "325B23TM32",
  "bog6core",     "325B23TM33",
  "bog6core",     "325B23TM34",
  "fen7core",     "325B23TM35",
  "fen7core",     "325B23TM36",
  "fen7core",     "325B23TM37",
  "fen7core",     "325B23TM38",
  "fen8core",     "325B23TM39",
  "fen8core",     "325B23TM40",
  "fen8core",     "325B23TM41",
  "fen8core",     "325B23TM42",
  "fen9core",     "325B23TM43",
  "fen9core",     "325B23TM44",
  "fen9core",     "325B23TM45",
  "fen9core",     "325B23TM46"
)

# add phys_code and plot_id to root_long
root_long <- root_long %>%
  # extract number from sample_id 
  mutate(sample_num = str_extract(sample_id, "^\\d+")) %>%
  
  # create phys_code by pasting prefix with the number
  mutate(PHYS_code = paste0("325B23TM", sample_num)) %>%
  
  # join with phys_code to get plot_id
  left_join(phys_lookup, by = "PHYS_code") %>%
  
  # remove sample_num helper column
  dplyr::select(-sample_num)


# move the new cols

root_long <- root_long %>% relocate(c(PHYS_code, plot_id, organ, PFT), .before = sample_id)


## calculate the percentage of <2 mm roots

## per increment

# clean and calculate % of root length < 2mm per phys_code

root_long$tot_length_cm_per_d_class <- as.numeric(root_long$tot_length_cm_per_d_class)
root_long$global_tot_len_cm <- as.numeric(root_long$global_tot_len_cm)

# with <2mm and >=2mm  per depth increment

# summarize per sample_id, size_class, phys_code, and pft
per_sampleid <- root_long %>%
  filter(organ == "root") %>%
  filter(str_detect(d_class_mm, "-")) %>%
  mutate(
    d_max = as.numeric(str_extract(d_class_mm, "(?<=-)\\d+\\.\\d+")),
    size_class = ifelse(d_max < 2, "<2mm", ">=2mm")
  ) %>%
  group_by(sample_id, PHYS_code, PFT, size_class) %>%
  summarise(
    length = sum(tot_length_cm_per_d_class, na.rm = TRUE),
    global_tot_len_cm = first(global_tot_len_cm), 
    .groups = "drop"
  )

# aggregate 
percent_by_size_incr <- per_sampleid %>%
  group_by(PHYS_code, PFT, size_class) %>%
  summarise(
    total_length = sum(length, na.rm = TRUE),
    total_global = sum(global_tot_len_cm, na.rm = TRUE),
    pct_length = (total_length / total_global) * 100,
    .groups = "drop"
  )

percent_by_size_incr <- as.data.frame(percent_by_size_incr)

## add depth info

# create the phys_codetosample_increment lookup
increment_lookup <- tribble(
  ~PHYS_code, ~sample_increment,
  "325B23TM1", "0-10",
  "325B23TM2", "0-10",
  "325B23TM3", "10-20",
  "325B23TM4", "10-20",
  "325B23TM5", "20-35.5",
  "325B23TM6", "20-35.5",
  "325B23TM9", "0-10",
  "325B23TM10", "0-10",
  "325B23TM11", "10-26",
  "325B23TM12", "10-26",
  "325B23TM13", "0-10",
  "325B23TM14", "0-10",
  "325B23TM15", "10-20",
  "325B23TM16", "10-20",
  "325B23TM17", "20-32",
  "325B23TM18", "20-32",
  "325B23TM19", "0-10",
  "325B23TM20", "0-10",
  "325B23TM21", "10-25",
  "325B23TM22", "10-25",
  "325B23TM23", "0-10",
  "325B23TM24", "0-10",
  "325B23TM25", "10-20",
  "325B23TM26", "10-20",
  "325B23TM27", "20-30",
  "325B23TM28", "20-30",
  "325B23TM29", "0-10",
  "325B23TM30", "0-10",
  "325B23TM31", "10-20",
  "325B23TM32", "10-20",
  "325B23TM33", "20-32",
  "325B23TM34", "20-32",
  "325B23TM35", "0-10",
  "325B23TM36", "0-10",
  "325B23TM37", "10-18",
  "325B23TM38", "10-18",
  "325B23TM39", "0-10",
  "325B23TM40", "0-10",
  "325B23TM41", "10-21",
  "325B23TM42", "10-21",
  "325B23TM43", "0-10",
  "325B23TM44", "0-10",
  "325B23TM45", "10-24",
  "325B23TM46", "10-24"
)

# join to your existing percent_by_size dataframe
percent_by_size_incr <- percent_by_size_incr %>%
  left_join(increment_lookup, by = "PHYS_code")

percent_by_size_incr <- percent_by_size_incr %>%
  left_join(phys_lookup, by = "PHYS_code")

percent_by_size_incr <- percent_by_size_incr %>% relocate(c(sample_increment, plot_id), .after = PHYS_code)
percent_by_size_incr <- percent_by_size_incr %>% relocate(c(plot_id), .before = PHYS_code)

percent_by_size_incr <- percent_by_size_incr %>%
  arrange(plot_id)

# save
write.csv(percent_by_size_incr, "path/root_percentages_increment.csv", row.names = FALSE)
pct_incr <- read.csv("path/root_percentages_increment.csv")

## per plot per pft 

per_plotid <- root_long %>%
  filter(organ == "root") %>%
  filter(str_detect(d_class_mm, "-")) %>%
  mutate(
    d_max = as.numeric(str_extract(d_class_mm, "(?<=-)\\d+\\.\\d+")),
    size_class = ifelse(d_max < 2, "<2mm", ">=2mm")
  ) %>%
  group_by(sample_id, plot_id, PFT, size_class) %>%
  summarise(
    length = sum(tot_length_cm_per_d_class, na.rm = TRUE),
    global_tot_len_cm = first(global_tot_len_cm), 
    .groups = "drop"
  )

# aggregate 
percent_by_size_plot <- per_plotid %>%
  group_by(plot_id, PFT, size_class) %>%
  summarise(
    total_length = sum(length, na.rm = TRUE),
    total_global = sum(global_tot_len_cm, na.rm = TRUE),
    pct_length = (total_length / total_global) * 100,
    .groups = "drop"
  )

percent_by_size_plot <- as.data.frame(percent_by_size_plot)

# save
write.csv(percent_by_size_plot, "path/root_percentages_plot.csv", row.names = FALSE)

# plot the plot scale

percent_by_size_plot$plot_id <- factor(
  percent_by_size_plot$plot_id,
  levels = c(
    "palsa1core", "palsa3core", "palsa5core",
    "bog2core", "bog4core", "bog6core",
    "fen7core", "fen8core", "fen9core"
  )
)

ggplot() +
  geom_point(
    data = percent_by_size_plot,
    aes(
      x = plot_id,
      y = pct_length,
      shape = PFT,
      fill = size_class,
      color = size_class
    ),
    size = 10,
    stroke = 1.5
  ) +
  theme_bw() +
  labs(x = "", y = expression(paste("percentage"))) +
  scale_fill_manual(
    values = c("<2mm" = "#165BAA", ">=2mm" = "white"),
    na.value = NA
  ) +
  scale_color_manual(
    values = c("<2mm" = "#165BAA", ">=2mm" = "#165BAA")
  ) +
  scale_shape_manual(values = c("herbaceous" = 21, "shrub" = 22)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme(
    panel.border = element_blank(),
    axis.line.y.left = element_line(),
    axis.line.x.bottom = element_line(),
    axis.line = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(size = 20),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.title = element_text(size = 20)
  )


#### aggregate root_long to traylevel

## first need to calculate the d_class_vol_sum_cm3

root_long$vol_d_class_cm3_per_d_class <- as.numeric(root_long$vol_d_class_cm3_per_d_class)

root_long <- root_long %>%
  group_by(sample_id) %>%
  mutate(vol_d_class_sum_cm3 = sum(vol_d_class_cm3_per_d_class, na.rm = TRUE)) %>%
  ungroup()

root_long <- as.data.frame(root_long)

# create a new df without the d classes and d class specific stuff, leave the repeated rows

root_tray <- subset(root_long, select = -c(d_class_mm, tot_length_cm_per_d_class, surface_area_cm2_per_d_class, vol_d_class_cm3_per_d_class))

# then aggregate to tray level

root_tray <- unique(root_tray)

#### add rhizome vol imagej data to the root_tray df

rhizome_imagej <- read.csv("path/rhizome_vol_len_RTD_imagej.csv")

# for this purpose, can remove the dry wt and tissue density
rhizome_imagej <- subset(rhizome_imagej, select = -c(rhizome_dry_wt_g, rhizome_tissue_density_gcm3_imagej))

# rename 11srhizome with adding _tray1 in root_long

root_long <- root_long %>%
  mutate(sample_id = if_else(sample_id == "11SRhizome", "11SRhizome_tray1", sample_id))

root_long <- root_long %>%
  mutate(tray = if_else(sample_id == "11SRhizome_tray1", "1", tray))


# convert rhizome avg d from cm to mm

rhizome_imagej$rhizome_avg_d_mm_imagej <- rhizome_imagej$rhizome_avg_d_cm_imagej*10

# combine dfs
root_complete <- left_join(root_tray, rhizome_imagej, by = c("sample_id", "PHYS_code"))

# create a new column for vol_cm3 where you combine the vol_d_class_sum_cm3 and imagej volumes

imagej_ids <- c(
  "2SRhizome_tray3", "4SRhizome_tray2", "11SRhizome_tray1", "9SRhizome_tray5",
  "15SRhizome_tray2", "29SRhizome_tray2", "23SRhizome_tray2", "21SRhizome_tray2",
  "14SRhizome_tray2", "14SRhizome_tray1", "15SRhizome_tray1", "21SRhizome_tray1",
  "29SRhizome_tray1", "4SRhizome_tray1", "2HRhizome_tray2", "39HRhizome_tray2",
  "35HRhizome_tray2", "29HRhizome_tray2", "23HRhizome_tray2", "19HRhizome_tray2",
  "19HRhizome_tray1", "35HRhizome_tray1", "2HRhizome_tray1"
)

# create new column based on condition
root_complete_2 <- root_complete %>%
  mutate(vol_cm3 = if_else(
    sample_id %in% imagej_ids,
    rhizome_tot_vol_cm3_imagej,
    vol_d_class_sum_cm3
  ))

# do the same for the other rhizome data as well

root_complete_2 <- root_complete_2 %>%
  mutate(tot_len_cm = if_else(
    sample_id %in% imagej_ids,
    rhizome_tot_length_cm_imagej,
    global_tot_len_cm
  ))

root_complete_2 <- root_complete_2 %>%
  mutate(avg_d_mm = if_else(
    sample_id %in% imagej_ids,
    rhizome_avg_d_mm_imagej,
    global_avg_d_mm
  ))

# remove the imagej columns etc

root_all <- subset(root_complete_2, select = -c(rhizome_avg_d_cm_imagej,
                                                rhizome_tot_length_cm_imagej,
                                                rhizome_tot_vol_cm3_imagej,
                                                rhizome_avg_d_mm_imagej,
                                                global_vol_cm3,
                                                vol_d_class_sum_cm3,
                                                global_tot_len_cm,
                                                global_avg_d_mm))

### aggregate to depthincrement level

root_depth <- root_all %>%
  group_by(PHYS_code, organ, PFT) %>%
  summarize(
    avg_d_mm = mean(avg_d_mm, na.rm = TRUE),
    across(c(tot_len_cm, vol_cm3, global_surface_area_cm2), \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  )

root_depth <- as.data.frame(root_depth)

## import the sample weight df

sample_wts <- read.csv("path/sample_wt_gwc_bd.csv")

# remove x from in front of the columns

names(sample_wts) <- sub("^X", "", names(sample_wts))

# replace the 0 values with na in subsample_volume_m3

sample_wts <- sample_wts %>%
  mutate(subsample_volume_m3 = ifelse(subsample_volume_m3 == 0, NA, subsample_volume_m3))


# create lookup values from sample_75 rows
sample75_values <- sample_wts %>%
  filter(subsample25_sample75 == "sample_75") %>%
  dplyr::select(sample_id, sample_increment, 
                `75sample_volume_cm3`, 
                fresh_peat_no_bag, 
                dry_peat_g) %>%
  rename(
    plot_id = sample_id,
    sample75_volume_cm3 = `75sample_volume_cm3`,
    fresh_peat_no_bag_sample75 = fresh_peat_no_bag,
    dry_peat_g_sample75 = dry_peat_g
  )

# join these values back to the main df by plot_id + sample_increment
sample_wts_updated <- sample_wts %>%
  left_join(sample75_values, by = c("sample_id" = "plot_id", "sample_increment"))

# fill missing 75sample_volume_cm3 from sample_75 into subsample_25
sample_wts_updated <- sample_wts_updated %>%
  mutate(
    `75sample_volume_cm3` = if_else(
      subsample25_sample75 == "subsample_25" & is.na(`75sample_volume_cm3`),
      sample75_volume_cm3,
      `75sample_volume_cm3`
    )
  )

# remove unnecessary cols
sample_wts_updated <- subset(sample_wts_updated, select = -c(dry_wt_5g_subsample_aluminum_bowl_75sample, 
                                                             fresh_wt_5g_subsample_aluminum_bowl_75sample,
                                                             aluminum_bowl_wt_g,
                                                             fresh_wt_g_sample_bag
))

# rename sample_id
colnames(sample_wts_updated)[1] <- "plot_id"

# combine with root_depth

filt_sample_wts_root <- semi_join(sample_wts_updated, root_depth, by = "PHYS_code")

# remove the sample column

filt_sample_wts_root <- subset(filt_sample_wts_root, select = -subsample25_sample75)

# combine with root_depth
root_depth_2 <- left_join(root_depth, filt_sample_wts_root, by = "PHYS_code")

# move some columns

root_depth_2 <- root_depth_2 %>%
  relocate(plot_id, .before = PHYS_code)

root_depth_2 <- root_depth_2 %>%
  relocate(sample_increment, increment_depth_cm, increment_diameter_cm, .after = PHYS_code)

### shrub subsample ###

# make a new column where you give the fresh weight of the shrub subsample

filt_sample_wts_root <- filt_sample_wts_root %>%
  mutate(shrub_subsample_wt_g = case_when(
    PHYS_code == "325B23TM2" ~ 15.6,
    PHYS_code == "325B23TM4" ~ 15.2,
    PHYS_code == "325B23TM9" ~ 15,
    PHYS_code == "325B23TM14" ~ 15,
    PHYS_code == "325B23TM15" ~ 15,
    PHYS_code == "325B23TM17" ~ 15.1,
    TRUE ~ NA_real_  
  ))

root_depth_2 <- root_depth_2 %>%
  mutate(shrub_subsample_wt_g = case_when(
    PHYS_code == "325B23TM2" ~ 15.6,
    PHYS_code == "325B23TM4" ~ 15.2,
    PHYS_code == "325B23TM9" ~ 15,
    PHYS_code == "325B23TM14" ~ 15,
    PHYS_code == "325B23TM15" ~ 15,
    PHYS_code == "325B23TM17" ~ 15.1,
    TRUE ~ NA_real_  
  ))


### import dry weight df

dry_wts <- read.csv("path/root_dry_wt.csv")

## first need to sum the p samples with the normal samples

dry_wts_2 <- dry_wts %>%
  mutate(PHYS_code_clean = sub("-.*", "", PHYS_code)) %>%
  group_by(Plot, PHYS_code_clean, organ, PFT) %>%
  summarize(dry_wt_g = sum(dry_wt_g, na.rm = TRUE), .groups = "drop")

dry_wts_2 <- as.data.frame(dry_wts_2)

# rename
colnames(dry_wts_2)[2] <- "PHYS_code"


dry_wts_2 <- dry_wts_2 %>%
  mutate(PFT = if_else(PFT == "herb", "herbaceous", PFT))

colnames(dry_wts_2)[1] <- "plot_id"


# combine with root_depth

root_depth_3 <- left_join(root_depth_2, dry_wts_2, by = c("plot_id", "PHYS_code", "organ", "PFT"))


## now the shrub subsamples

## scale the shrub root dry wt, vol, sa and tot len to the subsample 25% scale

root_depth_scaled <- root_depth_3 %>%
  mutate(
    # calculate shrub_fraction only when all values are available
    shrub_fraction_of_subsample = case_when(
      PFT == "shrub" & organ == "root" & !is.na(shrub_subsample_wt_g) ~ shrub_subsample_wt_g / fresh_peat_no_bag,
      TRUE ~ NA_real_
    ),
    
    # scale dry_wt_g only for shrub roots with a known subsample weight
    dry_wt_subsample_scaled = case_when(
      PFT == "shrub" & organ == "root" & !is.na(shrub_subsample_wt_g) ~ dry_wt_g / shrub_fraction_of_subsample,
      TRUE ~ dry_wt_g
    ),
    
    tot_len_subsample_scaled = case_when(
      PFT == "shrub" & organ == "root" & !is.na(shrub_subsample_wt_g) ~ tot_len_cm / shrub_fraction_of_subsample,
      TRUE ~ tot_len_cm
    ),
    
    vol_subsample_scaled = case_when(
      PFT == "shrub" & organ == "root" & !is.na(shrub_subsample_wt_g) ~ vol_cm3 / shrub_fraction_of_subsample,
      TRUE ~ vol_cm3
    ),
    
    sa_subsample_scaled = case_when(
      PFT == "shrub" & organ == "root" & !is.na(shrub_subsample_wt_g) ~ global_surface_area_cm2 / shrub_fraction_of_subsample,
      TRUE ~ global_surface_area_cm2
    )
  )


## plot level
# use the root_depth_scaled df

# first fix fen8core core depth
root_depth_scaled <- root_depth_scaled %>%
  mutate(
    increment_depth_cm = if_else(
      plot_id == "fen7core" & increment_depth_cm == 18,
      8,
      increment_depth_cm
    )
  )

# create a table of unique increments per plot_id
core_lengths <- root_depth_scaled %>%
  distinct(plot_id, sample_increment, increment_depth_cm) %>%
  group_by(plot_id) %>%
  summarize(core_length_cm = sum(increment_depth_cm, na.rm = TRUE), .groups = "drop")

# join it back to the main data
root_depth_scaled <- root_depth_scaled %>%
  left_join(core_lengths, by = "plot_id")

root_depth_scaled <- as.data.frame(root_depth_scaled)

root_depth_scaled <- subset(root_depth_scaled, select = -core_length_cm.x)
root_depth_scaled <- subset(root_depth_scaled, select = -increment_depth_unique)
colnames(root_depth_scaled)[33] <- "core_length_cm"

# calculate mean core diameter

root_depth_scaled <- root_depth_scaled %>%
  group_by(plot_id) %>%
  mutate(
    core_d_cm = if (n_distinct(increment_diameter_cm, na.rm = TRUE) > 1) {
      mean(increment_diameter_cm, na.rm = TRUE)
    } else {
      first(increment_diameter_cm)
    }
  ) %>%
  ungroup()

root_depth_scaled <- as.data.frame(root_depth_scaled)


# scale the subsample values to depth increment scale

root_depth_scaled_2 <- root_depth_scaled %>%
  group_by(plot_id, sample_increment, PFT, organ) %>%
  mutate(
    # calculate total dry peat mass of the increment (subsample + sample75)
    dry_peat_g_increment = dry_peat_g + dry_peat_g_sample75,
    
    # calculate fraction of total dry peat mass that the subsample represents
    peat_mass_fraction = dry_peat_g / dry_peat_g_increment,
    
    # scale subsamplelevel root traits to the increment scale
    dry_wt_g_increment_scaled = dry_wt_subsample_scaled / peat_mass_fraction,
    tot_len_cm_increment_scaled = tot_len_subsample_scaled / peat_mass_fraction,
    vol_cm3_increment_scaled = vol_subsample_scaled / peat_mass_fraction,
    sa_cm2_increment_scaled = sa_subsample_scaled / peat_mass_fraction
  ) %>%
  ungroup()

root_depth_scaled_2 <- as.data.frame(root_depth_scaled_2)

### for the depthaggregation, normalize each increment to 0.1 m

# calculate root mass / length / vol / sa per dry peat

root_depth_scaled_2$mass_per_dry_peat_gg <- root_depth_scaled_2$dry_wt_g_increment_scaled / root_depth_scaled_2$dry_peat_g_increment
root_depth_scaled_2$len_per_dry_peat_cmg <- root_depth_scaled_2$tot_len_cm_increment_scaled / root_depth_scaled_2$dry_peat_g_increment
root_depth_scaled_2$vol_per_dry_peat_cm3g <- root_depth_scaled_2$vol_cm3_increment_scaled / root_depth_scaled_2$dry_peat_g_increment
root_depth_scaled_2$sa_per_dry_peat_cm2g <- root_depth_scaled_2$sa_cm2_increment_scaled / root_depth_scaled_2$dry_peat_g_increment

# calculate root trait per soil volume

root_depth_scaled_2$mass_per_soil_vol_gcm3 <- root_depth_scaled_2$mass_per_dry_peat_gg * root_depth_scaled_2$bulk_density_gcm3
root_depth_scaled_2$len_per_soil_vol_cmcm3 <- root_depth_scaled_2$len_per_dry_peat_cmg * root_depth_scaled_2$bulk_density_gcm3
root_depth_scaled_2$vol_per_soil_vol_cm3cm3 <- root_depth_scaled_2$vol_per_dry_peat_cm3g * root_depth_scaled_2$bulk_density_gcm3
root_depth_scaled_2$sa_per_soil_vol_cm2cm3 <- root_depth_scaled_2$sa_per_dry_peat_cm2g * root_depth_scaled_2$bulk_density_gcm3

# convert to per m3

root_depth_scaled_2$mass_per_soil_vol_gm3 <- root_depth_scaled_2$mass_per_soil_vol_gcm3 * 1000000
root_depth_scaled_2$len_per_soil_vol_mm3 <- root_depth_scaled_2$len_per_soil_vol_cmcm3 * 10000
root_depth_scaled_2$vol_per_soil_vol_cm3m3 <- root_depth_scaled_2$vol_per_soil_vol_cm3cm3 * 1000000
root_depth_scaled_2$sa_per_soil_vol_cm2m3 <- root_depth_scaled_2$sa_per_soil_vol_cm2cm3 * 1000000

# calculate root trait per area (m2)
root_depth_scaled_2$mass_per_area_gm2 <- root_depth_scaled_2$mass_per_soil_vol_gm3 * root_depth_scaled_2$increment_depth_cm/100
root_depth_scaled_2$len_per_area_mm2 <- root_depth_scaled_2$len_per_soil_vol_mm3 * root_depth_scaled_2$increment_depth_cm/100
root_depth_scaled_2$vol_per_area_cm3m2 <- root_depth_scaled_2$vol_per_soil_vol_cm3m3 * root_depth_scaled_2$increment_depth_cm/100
root_depth_scaled_2$sa_per_area_cm2m2 <- root_depth_scaled_2$sa_per_soil_vol_cm2m3 * root_depth_scaled_2$increment_depth_cm/100

# standardize to 10 cm
root_depth_scaled_2$mass_gm2_std <- root_depth_scaled_2$mass_per_soil_vol_gm3 * 0.1
root_depth_scaled_2$len_mm2_std <- root_depth_scaled_2$len_per_soil_vol_mm3 * 0.1
root_depth_scaled_2$vol_cm3m2_std <- root_depth_scaled_2$vol_per_soil_vol_cm3m3 * 0.1
root_depth_scaled_2$sa_cm2m2_std <- root_depth_scaled_2$sa_per_soil_vol_cm2m3 * 0.1

# save
write.csv(root_depth_scaled_2, "path/root_depth_scale_extra.csv")

# aggregate a plotlevel version

# first calculate core volume
root_depth_scaled_2$core_volume_cm3 <- pi* ((root_depth_scaled_2$core_d_cm/2)^2) * root_depth_scaled_2$core_length_cm

# aggregate
root_plot <- root_depth_scaled_2 %>%
  group_by(plot_id, PFT, organ) %>%
  summarize(
    across(c(avg_d_mm, GWC), \(x) mean(x, na.rm = TRUE)),
    across(c(dry_wt_g_increment_scaled, tot_len_cm_increment_scaled, vol_cm3_increment_scaled, sa_cm2_increment_scaled,
             tot_len_cm, vol_cm3, global_surface_area_cm2, dry_wt_g, dry_peat_g_increment, total_fresh_peat_mass_g), \(x) sum(x, na.rm = TRUE)),
    dry_wt_g_depth_mean = mean(dry_wt_g_increment_scaled, na.rm = TRUE),
    tot_len_cm_depth_mean = mean(tot_len_cm_increment_scaled, na.rm = TRUE),
    vol_cm3_depth_mean = mean(vol_cm3_increment_scaled, na.rm = TRUE),
    sa_cm2_depth_mean = mean(sa_cm2_increment_scaled, na.rm = TRUE),
    core_d_cm = first(core_d_cm),
    core_volume_cm3 = first(core_volume_cm3),
    core_length_cm = first(core_length_cm),
    .groups = "drop"
  )

root_plot <- as.data.frame(root_plot)

# change some column names
root_plot <- root_plot %>%
  rename(dry_wt_g_scaled = dry_wt_g_increment_scaled,
         tot_len_cm_scaled = tot_len_cm_increment_scaled,
         vol_cm3_scaled = vol_cm3_increment_scaled,
         sa_cm2_scaled = sa_cm2_increment_scaled,
         dry_peat_g = dry_peat_g_increment)

root_plot <- root_plot %>%
  rename(dry_wt_g_raw = dry_wt_g,
         tot_len_cm_raw = tot_len_cm,
         vol_cm3_raw = vol_cm3,
         sa_cm2_raw = global_surface_area_cm2)


#### standardize root len, dry weight, vol and sa

# calculate bulk density

root_plot$bulk_density_gcm3 <- root_plot$dry_peat_g / root_plot$core_volume_cm3

# calculate root mass / length / vol / sa per dry peat

root_plot$mass_per_dry_peat_gg <- root_plot$dry_wt_g_scaled / root_plot$dry_peat_g
root_plot$len_per_dry_peat_cmg <- root_plot$tot_len_cm_scaled / root_plot$dry_peat_g
root_plot$vol_per_dry_peat_cm3g <- root_plot$vol_cm3_scaled / root_plot$dry_peat_g
root_plot$sa_per_dry_peat_cm2g <- root_plot$sa_cm2_scaled / root_plot$dry_peat_g

# calculate root trait per soil volume

root_plot$mass_per_soil_vol_gcm3 <- root_plot$mass_per_dry_peat_gg * root_plot$bulk_density_gcm3
root_plot$len_per_soil_vol_cmcm3 <- root_plot$len_per_dry_peat_cmg * root_plot$bulk_density_gcm3
root_plot$vol_per_soil_vol_cm3cm3 <- root_plot$vol_per_dry_peat_cm3g * root_plot$bulk_density_gcm3
root_plot$sa_per_soil_vol_cm2cm3 <- root_plot$sa_per_dry_peat_cm2g * root_plot$bulk_density_gcm3

# convert to per m3

root_plot$mass_per_soil_vol_gm3 <- root_plot$mass_per_soil_vol_gcm3 * 1000000
root_plot$len_per_soil_vol_mm3 <- root_plot$len_per_soil_vol_cmcm3 * 10000
root_plot$vol_per_soil_vol_cm3m3 <- root_plot$vol_per_soil_vol_cm3cm3 * 1000000
root_plot$sa_per_soil_vol_cm2m3 <- root_plot$sa_per_soil_vol_cm2cm3 * 1000000

# calculate root trait per area (m2)
root_plot$mass_per_area_gm2 <- root_plot$mass_per_soil_vol_gm3 * root_plot$core_length_cm/100
root_plot$len_per_area_mm2 <- root_plot$len_per_soil_vol_mm3 * root_plot$core_length_cm/100
root_plot$vol_per_area_cm3m2 <- root_plot$vol_per_soil_vol_cm3m3 * root_plot$core_length_cm/100
root_plot$sa_per_area_cm2m2 <- root_plot$sa_per_soil_vol_cm2m3 * root_plot$core_length_cm/100

# standardize to 30 cm
root_plot$mass_gm2_std <- root_plot$mass_per_soil_vol_gm3 * 0.3
root_plot$len_mm2_std <- root_plot$len_per_soil_vol_mm3 * 0.3
root_plot$vol_cm3m2_std <- root_plot$vol_per_soil_vol_cm3m3 * 0.3
root_plot$sa_cm2m2_std <- root_plot$sa_per_soil_vol_cm2m3 * 0.3

# standardize also the mean depth traits

root_plot$depth_mean_mass_per_dry_peat_gg <- root_plot$dry_wt_g_depth_mean / root_plot$dry_peat_g
root_plot$depth_mean_len_per_dry_peat_cmg <- root_plot$tot_len_cm_depth_mean / root_plot$dry_peat_g
root_plot$depth_mean_vol_per_dry_peat_cm3g <- root_plot$vol_cm3_depth_mean / root_plot$dry_peat_g
root_plot$depth_mean_sa_per_dry_peat_cm2g <- root_plot$sa_cm2_depth_mean / root_plot$dry_peat_g

# calculate root trait per soil volume

root_plot$depth_mean_mass_per_soil_vol_gcm3 <- root_plot$depth_mean_mass_per_dry_peat_gg * root_plot$bulk_density_gcm3
root_plot$depth_mean_len_per_soil_vol_cmcm3 <- root_plot$depth_mean_len_per_dry_peat_cmg * root_plot$bulk_density_gcm3
root_plot$depth_mean_vol_per_soil_vol_cm3cm3 <- root_plot$depth_mean_vol_per_dry_peat_cm3g * root_plot$bulk_density_gcm3
root_plot$depth_mean_sa_per_soil_vol_cm2cm3 <- root_plot$depth_mean_sa_per_dry_peat_cm2g * root_plot$bulk_density_gcm3

# convert to per m3

root_plot$depth_mean_mass_per_soil_vol_gm3 <- root_plot$depth_mean_mass_per_soil_vol_gcm3 * 1000000
root_plot$depth_mean_len_per_soil_vol_mm3 <- root_plot$depth_mean_len_per_soil_vol_cmcm3 * 10000
root_plot$depth_mean_vol_per_soil_vol_cm3m3 <- root_plot$depth_mean_vol_per_soil_vol_cm3cm3 * 1000000
root_plot$depth_mean_sa_per_soil_vol_cm2m3 <- root_plot$depth_mean_sa_per_soil_vol_cm2cm3 * 1000000

# calculate root trait per area (m2)
root_plot$depth_mean_mass_per_area_gm2 <- root_plot$depth_mean_mass_per_soil_vol_gm3 * root_plot$core_length_cm/100
root_plot$depth_mean_len_per_area_mm2 <- root_plot$depth_mean_len_per_soil_vol_mm3 * root_plot$core_length_cm/100
root_plot$depth_mean_vol_per_area_cm3m2 <- root_plot$depth_mean_vol_per_soil_vol_cm3m3 * root_plot$core_length_cm/100
root_plot$depth_mean_sa_per_area_cm2m2 <- root_plot$depth_mean_sa_per_soil_vol_cm2m3 * root_plot$core_length_cm/100

# standardize to 25 cm
root_plot$depth_mean_mass_gm2_std <- root_plot$depth_mean_mass_per_soil_vol_gm3 * 0.3
root_plot$depth_mean_len_mm2_std <- root_plot$depth_mean_len_per_soil_vol_mm3 * 0.3
root_plot$depth_mean_vol_cm3m2_std <- root_plot$depth_mean_vol_per_soil_vol_cm3m3 * 0.3
root_plot$depth_mean_sa_cm2m2_std <- root_plot$depth_mean_sa_per_soil_vol_cm2m3 * 0.3

# save

write.csv(root_plot, "path/root_plot_scale_extra.csv")

### clean the new dfs

# remove unnecessary columns

root_plot <- subset(root_plot, select = -c(tot_len_cm_raw, vol_cm3_raw, sa_cm2_raw, dry_wt_g_raw,
                                           dry_peat_g, total_fresh_peat_mass_g,
                                           dry_wt_g_scaled, tot_len_cm_scaled, vol_cm3_scaled,
                                           sa_cm2_scaled, core_d_cm, core_volume_cm3,
                                           core_length_cm,mass_per_dry_peat_gg, len_per_dry_peat_cmg,
                                           vol_per_dry_peat_cm3g, sa_per_dry_peat_cm2g,
                                           vol_per_soil_vol_cm3cm3,
                                           vol_per_soil_vol_cm3m3,
                                           depth_mean_mass_per_dry_peat_gg,
                                           depth_mean_len_per_dry_peat_cmg,
                                           depth_mean_vol_per_dry_peat_cm3g,
                                           depth_mean_sa_per_dry_peat_cm2g,
                                           depth_mean_mass_per_soil_vol_gcm3,
                                           depth_mean_len_per_soil_vol_cmcm3,
                                           depth_mean_vol_per_soil_vol_cm3cm3,
                                           depth_mean_sa_per_soil_vol_cm2cm3,
                                           depth_mean_mass_per_soil_vol_gm3,
                                           depth_mean_len_per_soil_vol_mm3,
                                           depth_mean_vol_per_soil_vol_cm3m3,
                                           depth_mean_sa_per_soil_vol_cm2m3,
                                           depth_mean_mass_per_area_gm2,
                                           depth_mean_len_per_area_mm2,
                                           depth_mean_vol_per_area_cm3m2,
                                           depth_mean_sa_per_area_cm2m2))


root_depth <- subset(root_depth_scaled_2, select = -c(tot_len_cm, vol_cm3, global_surface_area_cm2, dry_wt_g,
                                                      increment_diameter_cm, total_fresh_peat_mass_g, 
                                                      total_increment_volume_cm3,
                                                      total_sample_fresh_mass_per_total_increment_volume_gcm3,
                                                      subsample_volume_cm3, subsample_volume_m3,
                                                      `75sample_volume_cm3`,
                                                      fresh_peat_no_bag,
                                                      dry_peat_g, sample75_volume_cm3, fresh_peat_no_bag_sample75,
                                                      dry_peat_g_sample75, shrub_subsample_wt_g,
                                                      shrub_fraction_of_subsample, dry_wt_subsample_scaled, 
                                                      tot_len_subsample_scaled, vol_subsample_scaled, 
                                                      sa_subsample_scaled, core_length_cm, core_d_cm, dry_peat_g_increment, 
                                                      peat_mass_fraction, dry_wt_g_increment_scaled, 
                                                      tot_len_cm_increment_scaled, vol_cm3_increment_scaled,
                                                      sa_cm2_increment_scaled,
                                                      root_mass_per_dry_peat_gg, mass_per_dry_peat_gg, len_per_dry_peat_cmg,
                                                      len_per_dry_peat_gg,
                                                      vol_per_dry_peat_cm3g, sa_per_dry_peat_cm2g,
                                                      core_volume_cm3))

# move some columns

root_depth <- root_depth %>%
  relocate(plot_id, .before = PHYS_code)

root_depth <- root_depth %>%
  relocate(sample_increment, .after = PHYS_code)

root_depth <- subset(root_depth, select = -c(increment_depth_cm))

### calculate rtd and srl based on the "raw" data from root_tray

# depthlevel

root_depth_scaled_2$RTD_gcm3 <- root_depth_scaled_2$dry_wt_g / root_depth_scaled_2$vol_cm3

root_depth_scaled_2$SRL_mg1 <- (root_depth_scaled_2$tot_len_cm/100) / root_depth_scaled_2$dry_wt_g

# add the values to root_depth
root_depth <- root_depth %>%
  left_join(root_depth_scaled_2 %>% dplyr::select(plot_id, PHYS_code, sample_increment, PFT, organ, RTD_gcm3, SRL_mg1),
            by = c("plot_id", "PHYS_code", "sample_increment", "PFT", "organ"))


root_depth <- subset(root_depth, select = -c(RTD_gcm3.y, SRL_mg1.y))
root_depth <- root_depth %>%
  rename(RTD_gcm3 = RTD_gcm3.x,
         SRL_mg1 = SRL_mg1.x)


# set rhizome srl to na (doesn't make sense to calculate with my data)
root_depth <- root_depth %>%
  mutate(SRL_mg1 = case_when(
    organ == "rhizome" ~ NA_real_,
    TRUE ~ SRL_mg1
  ))

root_depth <- subset(root_depth, select = -c(PHYS_code))

# rename columns
root_depth <- root_depth %>%
  rename(len_density_cmcm3 = len_per_soil_vol_cmcm3,
         mass_density_gcm3 = mass_per_soil_vol_gcm3,
         sa_density_cm2cm3 = sa_per_soil_vol_cm2cm3)

root_depth$len_density_mcm3 <- root_depth$len_density_cmcm3 / 100
root_depth <- root_depth %>%
  rename(len_density_mm3 = len_per_soil_vol_mm3)
root_depth <- root_depth %>%
  rename(mass_density_gm3 = mass_per_soil_vol_gm3)
root_depth <- root_depth %>%
  rename(sa_density_cm2m3 = sa_per_soil_vol_cm2m3)


# plotlevel > take the values from root_all

plot_summary <- root_all %>%
  group_by(plot_id, PFT, organ) %>%
  summarize(
    avg_d_mm = mean(avg_d_mm, na.rm = TRUE),
    across(c(tot_len_cm, vol_cm3), \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  )

plot_summary <- as.data.frame(plot_summary)

dry_wts_3 <- dry_wts_2 %>%
  group_by(plot_id, PFT, organ) %>%
  summarize(
    across(dry_wt_g, \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  )

dry_wts_3 <- as.data.frame(dry_wts_3)

# combine with plot_summary

plot_summary <- plot_summary %>%
  left_join(dry_wts_3, by = c("plot_id", "PFT", "organ"))

# calculate rtd

plot_summary$RTD_gcm3 <- plot_summary$dry_wt_g / plot_summary$vol_cm3

# calculate srl
plot_summary$SRL_mg1 <- (plot_summary$tot_len_cm/100) / plot_summary$dry_wt_g

# add rtd and srl to root_plot

root_plot <- root_plot %>%
  left_join(plot_summary %>% dplyr::select(plot_id, PFT, organ, RTD_gcm3, SRL_mg1, avg_d_mm),
            by = c("plot_id", "PFT", "organ"))

# remove the old avg d

root_plot <- subset(root_plot, select = -avg_d_mm.x)

root_plot <- root_plot %>%
  relocate(avg_d_mm.y, .after = organ)

root_plot <- root_plot %>%
  relocate(GWC, .after = SRL_mg1)

root_plot <- root_plot %>%
  rename(avg_d_mm = avg_d_mm.y)

# set rhizome srl to na
root_plot <- root_plot %>%
  mutate(SRL_mg1 = case_when(
    organ == "rhizome" ~ NA_real_,
    TRUE ~ SRL_mg1
  ))

# add columns

root_plot <- root_plot %>%
  mutate(thaw_stage = case_when(
    str_detect(plot_id, "palsa") ~ "intact",
    str_detect(plot_id, "bog") ~ "partly_thawed",
    str_detect(plot_id, "fen") ~ "fully_thawed"
  ))

root_plot <- root_plot %>%
  mutate(chamber_id = case_when(
    str_detect(plot_id, "1") ~ 1,
    str_detect(plot_id, "2") ~ 2,
    str_detect(plot_id, "3") ~ 3,
    str_detect(plot_id, "4") ~ 4,
    str_detect(plot_id, "5") ~ 5,
    str_detect(plot_id, "6") ~ 6,
    str_detect(plot_id, "7") ~ 7,
    str_detect(plot_id, "8") ~ 8,
    str_detect(plot_id, "9") ~ 9
  ))

root_depth <- root_depth %>%
  mutate(thaw_stage = case_when(
    str_detect(plot_id, "palsa") ~ "intact",
    str_detect(plot_id, "bog") ~ "partly_thawed",
    str_detect(plot_id, "fen") ~ "fully_thawed"
  ))

root_depth <- root_depth %>%
  mutate(chamber_id = case_when(
    str_detect(plot_id, "1") ~ 1,
    str_detect(plot_id, "2") ~ 2,
    str_detect(plot_id, "3") ~ 3,
    str_detect(plot_id, "4") ~ 4,
    str_detect(plot_id, "5") ~ 5,
    str_detect(plot_id, "6") ~ 6,
    str_detect(plot_id, "7") ~ 7,
    str_detect(plot_id, "8") ~ 8,
    str_detect(plot_id, "9") ~ 9
  ))

root_plot <- root_plot %>%
  relocate(chamber_id, thaw_stage, .after = plot_id)

root_depth <- root_depth %>%
  relocate(chamber_id, thaw_stage, .after = plot_id)

root_plot <- root_plot %>%
  rename(len_density_cmcm3 = len_per_soil_vol_cmcm3,
         mass_density_gcm3 = mass_per_soil_vol_gcm3,
         sa_density_cm2cm3 = sa_per_soil_vol_cm2cm3)

root_plot$len_density_mcm3 <- root_plot$len_density_cmcm3 / 100
root_plot <- root_plot %>%
  rename(len_density_mm3 = len_per_soil_vol_mm3)
root_plot <- root_plot %>%
  rename(mass_density_gm3 = mass_per_soil_vol_gm3)
root_plot <- root_plot %>%
  rename(sa_density_cm2m3 = sa_per_soil_vol_cm2m3)


# save

write.csv(root_plot, "path/root_plot_scale.csv", row.names = FALSE)

write.csv(root_depth, "path/root_depth_scale.csv", row.names = FALSE)

##############################

# combine vga data to root data

vga_data <- read.csv("path/VGA_data.csv")

# choose only plot_id, chamber_id, thaw_stage, aer_stem etc

vga_data <- subset(vga_data, select = c(plot_ID, chamber_ID, thaw_stage, aer_stem_density_nom2, 
                                        nonaer_stem_density_nom2, herb_stem_density_nom2, 
                                        shrub_stem_density_nom2, 
                                        aer_stem_density_sd, 
                                        nonaer_stem_density_sd, 
                                        herb_stem_density_sd, 
                                        shrub_stem_density_sd,
                                        aer_GA_m2m2, 
                                        nonaer_GA_m2m2,
                                        herb_GA_m2m2, 
                                        shrub_GA_m2m2))


vga_data <- unique(vga_data)

### make two versions: 1. with chamber plots where you take the root data from core and use that as if it was chamber, for ch4

### 2. a core version which has root data and the vga data

# remove core from vga data, change plot_id in both vga data and root_plot and root_depth and take root data and add to vga data

# chamber (vga) + core (root)

vga_data_2 <- vga_data %>%
  filter(!grepl("core", plot_ID))

vga_data_2$plot_ID <- sub("chamber$", "", vga_data_2$plot_ID)

colnames(vga_data_2)[1] <- "plot_id"
colnames(vga_data_2)[2] <- "chamber_id"

root_plot_2 <- root_plot

root_plot_2$plot_id <- sub("core$", "", root_plot_2$plot_id)

root_depth_2 <- root_depth
root_depth_2$plot_id <- sub("core$", "", root_depth_2$plot_id)

root_vga_combined_plot <- root_plot_2 %>%
  left_join(vga_data_2,
            by = c("plot_id", "chamber_id", "thaw_stage"))

root_vga_combined_depth <- root_depth_2 %>%
  left_join(vga_data_2,
            by = c("plot_id", "chamber_id", "thaw_stage"))


# save

write.csv(root_vga_combined_plot, "path/root_vga_plot_scale.csv", row.names = FALSE)

write.csv(root_vga_combined_depth, "path/root_vga_depth_scale.csv", row.names = FALSE)


### 2. remove chamber from vga data, and add to root_plot and root_depth

# core (vga) + core (root)

vga_data_3 <- vga_data %>%
  filter(!grepl("chamber", plot_ID))

colnames(vga_data_3)[1] <- "plot_id"
colnames(vga_data_3)[2] <- "chamber_id"

root_vga_combined_plot_core <- root_plot %>%
  left_join(vga_data_3,
            by = c("plot_id", "chamber_id", "thaw_stage"))

root_vga_combined_depth_core <- root_depth %>%
  left_join(vga_data_3,
            by = c("plot_id", "chamber_id", "thaw_stage"))

# save

write.csv(root_vga_combined_plot_core, "path/root_vga_plot_scale_core.csv", row.names = FALSE)

write.csv(root_vga_combined_depth_core, "path/root_vga_depth_scale_core.csv", row.names = FALSE)

#############################################################
###------------------NEE AND GPP--------------------------###
#############################################################

NEE <- read.csv("path/ICOSETC_SE-Sto_FLUXNET_HH_L2.csv")

# remove scientific notation in the timestamps

NEE$TIMESTAMP_START <- format(NEE$TIMESTAMP_START, scientific = FALSE, trim = TRUE)
NEE$TIMESTAMP_END <- format(NEE$TIMESTAMP_END, scientific = FALSE, trim = TRUE)

NEE$TIMESTAMP_START <- paste0(NEE$TIMESTAMP_START, "00")
NEE$TIMESTAMP_END   <- paste0(NEE$TIMESTAMP_END, "00")
NEE$TIMESTAMP_START <- ymd_hms(NEE$TIMESTAMP_START, tz = "UTC")
NEE$TIMESTAMP_END   <- ymd_hms(NEE$TIMESTAMP_END, tz = "UTC")


# filter for data only from the year 2023
NEE_2023 <- NEE[year(NEE$TIMESTAMP_START) == 2023, ]

# subset to include only relevant columns 

NEE_2023 <- subset(NEE_2023, select = c(TIMESTAMP_START, TIMESTAMP_END, TA_F_MDS, TA_F_MDS_QC, 
                                        TA_ERA, TA_F, TA_F_QC, CO2_F_MDS, CO2_F_MDS_QC, TS_F_MDS_1, TS_F_MDS_2,
                                        TS_F_MDS_3, TS_F_MDS_4, TS_F_MDS_5, TS_F_MDS_6,
                                        TS_F_MDS_1_QC, TS_F_MDS_2_QC, TS_F_MDS_3_QC, TS_F_MDS_4_QC, TS_F_MDS_5_QC, TS_F_MDS_6_QC,
                                        SWC_F_MDS_1, SWC_F_MDS_2, SWC_F_MDS_3, SWC_F_MDS_4, SWC_F_MDS_5, SWC_F_MDS_6, SWC_F_MDS_1_QC,
                                        SWC_F_MDS_2_QC, SWC_F_MDS_3_QC, SWC_F_MDS_4_QC, SWC_F_MDS_5_QC, SWC_F_MDS_6_QC,
                                        NEE_VUT_REF, NEE_VUT_REF_QC
))

# gpp version
GPP_2023 <- subset(NEE_2023, select = c(TIMESTAMP_START, TIMESTAMP_END, TA_F_MDS, TA_F_MDS_QC, 
                                        TA_ERA, TA_F, TA_F_QC, CO2_F_MDS, CO2_F_MDS_QC, TS_F_MDS_1, TS_F_MDS_2,
                                        TS_F_MDS_3, TS_F_MDS_4, TS_F_MDS_5, TS_F_MDS_6,
                                        TS_F_MDS_1_QC, TS_F_MDS_2_QC, TS_F_MDS_3_QC, TS_F_MDS_4_QC, TS_F_MDS_5_QC, TS_F_MDS_6_QC,
                                        SWC_F_MDS_1, SWC_F_MDS_2, SWC_F_MDS_3, SWC_F_MDS_4, SWC_F_MDS_5, SWC_F_MDS_6, SWC_F_MDS_1_QC,
                                        SWC_F_MDS_2_QC, SWC_F_MDS_3_QC, SWC_F_MDS_4_QC, SWC_F_MDS_5_QC, SWC_F_MDS_6_QC,
                                        NEE_VUT_REF, NEE_VUT_REF_QC, GPP_NT_VUT_REF, GPP_DT_VUT_REF
))


# filter out qc=3

NEE_2023_filtered <- subset(NEE_2023, NEE_VUT_REF_QC != 3)

GPP_2023_filtered <- subset(GPP_2023, NEE_VUT_REF_QC != 3)

# aggregate to daily scale
NEE_daily <- NEE_2023_filtered %>%
  mutate(date = date(TIMESTAMP_START)) %>%
  group_by(date) %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)), .groups = "drop")

NEE_daily <- as.data.frame(NEE_daily)


GPP_daily <- GPP_2023_filtered %>%
  mutate(date = date(TIMESTAMP_START)) %>%
  group_by(date) %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)), .groups = "drop")

GPP_daily <- as.data.frame(GPP_daily)

# calculate the productive season start date and end date

# filter to valid days
NEE_daily <- NEE_daily %>%
  filter(!is.na(NEE_VUT_REF))

GPP_daily <- GPP_daily %>%
  filter(!is.na(NEE_VUT_REF))

# identify growing season start and end:
# first day where nee becomes negative (uptake)
start_date <- NEE_daily %>%
  filter(NEE_VUT_REF < 0) %>%
  slice(1) %>% 
  pull(date)

# last day where nee is still negative
end_date <- NEE_daily %>%
  filter(NEE_VUT_REF < 0) %>%
  slice_tail(n = 1) %>%
  pull(date)

# combine as output
productive_season <- tibble(
  start = start_date,
  end = end_date,
  length_days = as.integer(end_date - start_date + 1)
)

#   start:             end:           length_days:
#   2023-05-19    2023-08-30            104

# cut halfhourly gpp df to productive season

GPP_2023_filtered <- GPP_2023_filtered %>%
  mutate(TIMESTAMP_START = force_tz(TIMESTAMP_START, tzone = "UTC"))

GPP_2023_filtered_2 <- GPP_2023_filtered %>%
  filter(
    TIMESTAMP_START >= as.POSIXct("2023-05-19 00:00:00", tz = "UTC") &
      TIMESTAMP_START <= as.POSIXct("2023-08-30 23:59:59", tz = "UTC")
  )

# create a column to distinguish positive vs negative nee
NEE_daily <- NEE_daily %>%
  mutate(NEE_sign = ifelse(NEE_VUT_REF < 0, "uptake", "release"))


# for plotting:
start_date <- as.Date("2023-05-19")
end_date   <- as.Date("2023-08-30")

# monthly ticks covering your data range
month_seq <- seq(
  floor_date(min(NEE_daily$date, na.rm = TRUE), "month"),
  ceiling_date(max(NEE_daily$date, na.rm = TRUE), "month"),
  by = "1 month"
)

month_seq <- seq(
  floor_date(min(GPP_daily$date, na.rm = TRUE), "month"),
  ceiling_date(max(GPP_daily$date, na.rm = TRUE), "month"),
  by = "1 month"
)

# combine month ticks + start/end
x_breaks <- sort(unique(c(month_seq, start_date, end_date)))


# monthly ticks covering your data range
month_seq <- seq(
  floor_date(min(NEE_daily$date, na.rm = TRUE), "month"),
  ceiling_date(max(NEE_daily$date, na.rm = TRUE), "month"),
  by = "1 month"
)

# plot
NEE <- ggplot(NEE_daily, aes(x = date, y = NEE_VUT_REF, fill = NEE_sign)) +
  annotate("rect",
           xmin = as.Date("2023-05-19"),
           xmax = as.Date("2023-08-30"),
           ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.4) +
  geom_col() +
  scale_fill_manual(values = c("uptake" = "darkgreen", "release" = "steelblue")) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  scale_x_date(
    breaks = month_seq,
    labels = scales::date_format("%b"), 
    expand = c(0, 0)
  ) +
  labs(
    x = "",
    y = expression(paste("NEE (", mu, "mol CO"[2], " m"^-2, " s"^-1, ")")),
    fill = "CO2 flux"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text    = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "none"
  )

NEE


GPP <- ggplot(GPP_daily, aes(x = date, y = GPP_NT_VUT_REF)) +
  annotate("rect",
           xmin = as.Date("2023-05-19"),
           xmax = as.Date("2023-08-30"),
           ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.4) +
  geom_col(fill = "grey50") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  scale_x_date(
    breaks = month_seq,
    labels = scales::date_format("%b"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5)) +  
  labs(
    x = "",
    y = expression(paste("GPP (", mu, "mol CO"[2], " m"^-2, " s"^-1, ")"))
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text    = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text  = element_text(size = 18),
    legend.title = element_blank()
  )

GPP

## some stats of gpp

# filter for data only from the year 2023
GPP_daily_filtered <- GPP_daily %>%
  filter(date >= ymd("2023-05-19"),
         date <= ymd("2023-08-30"))

GPP_monthly_summary <- GPP_daily_filtered %>%
  mutate(month = floor_date(date, unit = "month")) %>%
  group_by(month) %>%
  summarise(
    n      = sum(!is.na(GPP_NT_VUT_REF)),
    median = median(GPP_NT_VUT_REF, na.rm = TRUE),
    mean   = mean(GPP_NT_VUT_REF, na.rm = TRUE),
    IQR    = IQR(GPP_NT_VUT_REF, na.rm = TRUE),
    sd     = sd(GPP_NT_VUT_REF, na.rm = TRUE),
    cv     = raster::cv(GPP_NT_VUT_REF, na.rm = TRUE),
    .groups = "drop"
  )

GPP_monthly_summary

#### diurnal scale

# use gpp_2023_filtered

# and filtered_flux (for chamber ch4 flux data)


# chamber data
filtered_flux$datetime <- as_datetime(filtered_flux$datetime)

filtered_flux <- filtered_flux %>%
  mutate(datetime = force_tz(datetime, "UTC"))

# ec data
GPP_2023_filtered_2 <- GPP_2023_filtered_2 %>%
  mutate(
    TIMESTAMP_START = force_tz(TIMESTAMP_START, "UTC"),
    TIMESTAMP_END   = force_tz(TIMESTAMP_END, "UTC")
  )

# make sure timestamps are posixct in utc
filtered_flux[, datetime := force_tz(datetime, "UTC")]

# check histogram of gpp
hist(GPP_2023_filtered_2$GPP_DT_VUT_REF,
     breaks = 50,
     main = "Histogram of GPP_DT_VUT_REF",
     xlab = "GPP_DT_VUT_REF",
     col = "grey")

GPP_2023_filtered_2$hour <- hour(GPP_2023_filtered_2$TIMESTAMP_START)

# aggregate to hourly scale
gpp_summary <- GPP_2023_filtered_2 %>%
  group_by(hour) %>%
  summarise(
    median_GPP = median(GPP_DT_VUT_REF, na.rm = TRUE),
    q25 = quantile(GPP_DT_VUT_REF, 0.25, na.rm = TRUE),
    q75 = quantile(GPP_DT_VUT_REF, 0.75, na.rm = TRUE)
  )

# gpp plot
p_gpp_hr <- ggplot(gpp_summary, aes(x = hour, y = median_GPP)) +
  geom_errorbar(aes(ymin = q25, ymax = q75),
                width = 0.2,
                color = "black") +
  geom_point(color = "grey50", size = 2.5) +
  geom_line(color = "grey50", linewidth = 1) +
  labs(y = expression(paste("GPP (", mu, "mol CO"[2], " m"^-2, " s"^-1, ")")), x = "Hour") +
  theme_bw() +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=18))

p_gpp_hr


# chamber ch4 fluxes

filtered_flux$hour <- lubridate::hour(filtered_flux$datetime)

# create the median and mean summary with uncertainty ranges
ch4_summary <- filtered_flux %>%
  group_by(thaw_stage, hour) %>%
  summarise(
    median_ch4 = median(F_JmgCH4m2d, na.rm = TRUE),
    mean_ch4   = mean(F_JmgCH4m2d, na.rm = TRUE),
    sd_ch4     = sd(F_JmgCH4m2d, na.rm = TRUE),
    q25        = quantile(F_JmgCH4m2d, 0.25, na.rm = TRUE),
    q75        = quantile(F_JmgCH4m2d, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ch4_summary$thaw_stage <- factor(
  ch4_summary$thaw_stage,
  levels = c("intact", "partly_thawed", "fully_thawed")
)

# plot
p_ch4 <- ggplot(ch4_summary, aes(x = hour)) +
  geom_errorbar(aes(ymin = q25,
                    ymax = q75),
                width = 0.2,
                color = "#E30B5C") +
  geom_point(aes(y = median_ch4),
             color = "#E30B5C",
             size = 2.5) +
  geom_line(aes(y = median_ch4),
            color = "#E30B5C",
            linewidth = 1) +
  geom_errorbar(aes(x = hour + 0.2,
                    ymin = mean_ch4 - sd_ch4,
                    ymax = mean_ch4 + sd_ch4),
                width = 0.2,
                color = "#1F77B4") +
  geom_point(aes(x = hour + 0.2,
                 y = mean_ch4),
             color = "#1F77B4",
             size = 2.5) +
  geom_line(aes(x = hour + 0.2,
                y = mean_ch4),
            color = "#1F77B4",
            linewidth = 1) +
  facet_wrap(~thaw_stage) +
  scale_x_continuous(breaks = 0:23, limits = c(0, 23.2)) +
  labs(
    y = expression(paste("CH"[4], " flux (mg CH"[4], " m"^-2, " d"^-1, ")")),
    x = "Hour"
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18)
  )

p_ch4

### combine nee, gpp (on the same row), and
### final_plot
### Fig. B1

# top row: dailyscale nee and gpp
top_row <- NEE + GPP +
  plot_layout(ncol = 2)

# bottom rows: diurnal gpp above diurnal ch4
bottom_row <- p_gpp_hr / p_ch4 +
  plot_layout(heights = c(1, 1.4))

# full figure
final_plot <- top_row / bottom_row +
  plot_layout(heights = c(1, 2))

final_plot

#################################################################
###------------------METEOROLOGICAL DATA----------------------###
#################################################################

# used in manuscript section 2.1

temp <- read.csv("path/smhi-opendata_2_188790_20250709_095940.csv")

precip <- read.csv("path/smhi-opendata_5_188790_20250709_120531.csv")

# mean annual air temperature

temp %>%
  mutate(date = dmy(date),
         year = year(date)) %>%
  filter(year >= 1986, year <= 2023) %>%                 
  group_by(year) %>%
  summarise(mean_temp = mean(temp, na.rm = TRUE)) %>%
  summarise(mean_annual_air_temp = mean(mean_temp, na.rm = TRUE)) %>%
  pull(mean_annual_air_temp)  

# precipitation
precip %>%
  mutate(date = dmy(date),
         year = year(date)) %>%
  filter(year >= 1986, year <= 2023) %>%            
  group_by(year) %>%
  summarise(annual_sum = sum(precip, na.rm = TRUE)) %>%
  summarise(mean_annual_precip_sum = mean(annual_sum, na.rm = TRUE)) %>%
  pull(mean_annual_precip_sum)  

# mean productive season TA

TA_daily$date <- as_date(TA_daily$date)

TA_daily %>%
  filter(date >= as.Date("2023-05-19"), date <= as.Date("2023-08-30")) %>%                 
  summarise(mean_temp = mean(TA_F, na.rm = TRUE),
            sd_temp = sd(TA_F, na.rm = T)) %>%
  pull(sd_temp)  

P_daily$date <- as_date(P_daily$date)

P_daily %>%
  filter(date >= as.Date("2023-05-19"), date <= as.Date("2023-08-30")) %>%                 
  summarise(sum_precip = sum(P_F, na.rm = TRUE)) %>%
  pull(sum_precip)  


#############################################################
###--------------------CH4 FLUX DATA----------------------###
#############################################################

# import data

flux <- read.csv("path/Stordalen_AC_Flagged_2023.csv")

flux <- subset(flux, select = -c(seqday.1, set_num, ppmCH4_prior, ppmCO2_prior, ppmCH4_init, ppmCO2_init,
                                 CO2_r2, avgPAR, avgTair, avgTair_filled,
                                 JumolCO2m2s_unf, JmgCH4m2d_unf,
                                 X90_JumolCO2m2s, X90_JmgCH4m2d, X95_JumolCO2m2s, X95_JmgCH4m2d,
                                 Chamber.flag, Questions.CO2, Questions.CH4, r2.flag.NEE, r2.flag.CH4, Night.Day,
                                 Lo.PAR.Uptake, Hi.CH4.filter, Hi.CO2.filter, prior.CH4, prior.CO2,
                                 Air_Press, T_107, CH4_r2))

flux$sitenum <- as.factor(flux$sitenum)
flux$chmbr <- as.factor(flux$chmbr)

# sitenum: 1 = intact, 2 = partly thawed, 3 = fully thawed

# convert doy to datetime

flux$datetime <- as.POSIXct(
  paste(flux$year, "01-01"),
  format = "%Y %m-%d",
  tz = "UTC"
) + (flux$doy - 1) * 86400

# subset dataset to cover productive season: 
# start: 2023-05-19  end: 2023-08-30

flux <- flux %>%
  filter(as.Date(datetime) >= as.Date("2023-05-19") &
           as.Date(datetime) <= as.Date("2023-08-30"))

# remove the first number from chmbr id
flux$chmbr <- as.character(flux$chmbr)
flux$chmbr <- as.numeric(substr(as.character(flux$chmbr), 
                                nchar(flux$chmbr), 
                                nchar(flux$chmbr)))

# rename sitenum to thaw_stage and change values

flux <- flux %>%
  rename(
    thaw_stage = sitenum,
    chamber_id = chmbr
  )

flux <- flux %>% 
  mutate(thaw_stage = recode(thaw_stage,
                             "1" = "intact",
                             "2" = "partly_thawed",
                             "3" = "fully_thawed"))

# create a date column
flux$date <- as_date(flux$datetime)

# remove columns that are not needed
flux <- subset(flux, select = -c(year, seqday, doy, hour, hour_num, week, F_JumolCO2m2s))

# count the number of days where all soil temp are NA

na_days <- flux %>%
  group_by(chamber_id, date) %>%
  summarise(all_na = all(is.na(avgTgrnd)), .groups = "drop") %>%
  filter(all_na)
nrow(na_days)
# 20 

# check the minimum number of nonna observations within date per chamber

flux %>%
  group_by(chamber_id, date) %>%
  summarise(n_nonNA = sum(!is.na(avgTgrnd)), .groups = "drop") %>%
  filter(n_nonNA > 0) %>% 
  summarise(min_nonNA = min(n_nonNA))
# 2 > too little to rely on the gapfilled estimates

# check which hours of day the temp measurements were taken
flux %>%
  mutate(hour = hour(datetime)) %>%
  group_by(chamber_id, hour) %>%
  summarise(prop_nonNA = mean(!is.na(avgTgrnd)), .groups = "drop") %>%
  ggplot(aes(x = hour, y = prop_nonNA)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ chamber_id) +
  scale_x_continuous(breaks = 0:23) +
  labs(
    x = "Hour of day",
    y = "Proportion of days with temperature data",
    title = "Temporal distribution of avgTgrnd observations",
    subtitle = "Shows when during the day measurements are available"
  ) +
  theme_minimal()

flux_counts <- flux %>%
  group_by(chamber_id, date) %>%
  summarise(n_nonNA = sum(!is.na(avgTgrnd)), .groups = "drop")
days_per_count <- flux_counts %>%
  count(chamber_id, n_nonNA)
days_per_count %>%
  filter(n_nonNA %in% c(4))

ggplot(days_per_count, aes(x = n_nonNA, y = n)) +
  geom_col() +
  facet_wrap(~ chamber_id, scales = "free_y") +
  labs(
    x = "Number of non-NA avgTgrnd observations per day",
    y = "Number of days",
    title = "Data coverage per day per chamber"
  ) +
  theme_minimal()

# max number of temp observations per day: 8, min = 2. 
# > set avgtgrnd to na in dates where there are <3 non-NA osbervations

# soil temp: set gap-filled soil temp to na in dates 
# where non-gap-filled temps are all na for both (because gap-filled temp is uncertain)

flux <- flux %>%
  group_by(chamber_id, date) %>%
  mutate(
    avgTgrnd_filled = ifelse(
      all(is.na(avgTgrnd)),
      NA,
      avgTgrnd_filled
    )
  ) %>%
  ungroup()

# set to NA when there are <3 non-NA observations
flux <- flux %>%
  group_by(chamber_id, date) %>%
  mutate(
    n_nonNA = sum(!is.na(avgTgrnd)),
    avgTgrnd = ifelse(n_nonNA < 3, NA, avgTgrnd),
    avgTgrnd_filled = ifelse(n_nonNA < 3, NA, avgTgrnd_filled)
  ) %>%
  ungroup() %>%
  select(-n_nonNA)


write.csv(flux, "path/flux_data_subset.csv", row.names = FALSE)

# read in the df again
filtered_flux <- read.csv("path/flux_data_subset.csv")

# aggregate per thaw stage to check per thaw stage daily trends
daily_summary_thaw_stage <- filtered_flux %>%
  group_by(thaw_stage, date = as.Date(date)) %>%
  summarise(daily_median_CH4 = median(F_JmgCH4m2d, na.rm = TRUE),
            daily_mean_CH4 = mean(F_JmgCH4m2d, na.rm = TRUE),
            .groups = "drop")

# plot
ggplot(daily_summary_thaw_stage, aes(x = date, y = daily_median_CH4)) +
  geom_point() +
  facet_wrap(~ thaw_stage, scales = "free") +
  theme_bw()

# median fluxes increase almost linearly towards september

### productive season data doesn't fit into gaussian curve

### take: 
### 1. peak ch4 date
### 2. date of 50% peak flux (kind of the mid point fluxes on the rise)
### 3. early third season median flux (baseline flux)
### 4. season median

## calculate these per thaw stage

env_vars <- c("avgTgrnd")

# histogram and qq plot for each variable
par(mfrow = c(2, 1))

for (var in env_vars) {
  hist(filtered_flux[[var]], main = paste("Histogram of", var), xlab = var, col = "lightblue")
  qqnorm(filtered_flux[[var]], main = paste("Q-Q Plot of", var))
  qqline(filtered_flux[[var]], col = "red")
}
par(mfrow = c(1, 1))  

## aggregate dataset to daily scale

# daily medians and means of ch4 and soil temp
flux_daily <- filtered_flux %>%
  group_by(chamber_id, date) %>%
  summarise(
    daily_median_CH4 = median(F_JmgCH4m2d, na.rm = TRUE),
    daily_mean_CH4 = mean(F_JmgCH4m2d, na.rm = TRUE),
    avgTgnd_mean  = mean(avgTgrnd,  na.rm = TRUE),
    thaw_stage = dplyr::first(thaw_stage),
    .groups = "drop"
  )

# ch4 daily medians table
daily_medians <- flux_daily %>% select(chamber_id, date, daily_median_CH4, daily_mean_CH4)

daily_medians$date <- as_date(daily_medians$date)

##### GAMs


daily_medians <- as.data.frame(daily_medians)

# ch4 flux data is not normal > cannot use it as is in gaussian + identity gam
qqnorm(daily_medians$daily_median_CH4, pch = 1, frame = FALSE)
qqline(daily_medians$daily_median_CH4, col = "steelblue", lwd = 2)

# check asinh trasnformation

crit_from_gam <- list()
fit_plots_per_chamber   <- vector("list", length(chambers)); names(fit_plots_per_chamber) <- chambers
resid_plots_per_chamber <- vector("list", length(chambers)); names(resid_plots_per_chamber) <- chambers

for (i in seq_along(chambers)) {
  ch <- chambers[i]
  dch <- daily_medians %>%
    filter(chamber_id == ch) %>%
    arrange(date) %>%
    mutate(daynum = as.numeric(date))
  
  # fit gam on asinh(ch4)
  m <- gam(asinh(daily_median_CH4) ~ s(daynum),
           data = dch, method = "REML")
  
  # prediction grid + backtransform + cis
  pred <- tibble(
    date   = seq(min(dch$date), max(dch$date), by = "1 day"),
    daynum = as.numeric(date)
  )
  pr      <- predict(m, newdata = pred, se.fit = TRUE)
  pred$fit_asinh <- as.numeric(pr$fit)
  pred$se_asinh  <- as.numeric(pr$se.fit)
  pred$fit       <- sinh(pred$fit_asinh)
  pred$fit_low   <- sinh(pred$fit_asinh - 1.96 * pred$se_asinh)
  pred$fit_high  <- sinh(pred$fit_asinh + 1.96 * pred$se_asinh)
  
  # timebased early window (first third of season)
  season_start <- min(dch$date); season_end <- max(dch$date)
  early_end    <- season_start + (season_end - season_start) / 3
  
  # critical values from the fitted curve 
  baseline_fit      <- median(pred$fit[pred$date <= early_end], na.rm = TRUE)
  baseline_fit_low  <- median(pred$fit_low[pred$date <= early_end], na.rm = TRUE)
  baseline_fit_high <- median(pred$fit_high[pred$date <= early_end], na.rm = TRUE)
  
  peak_idx        <- which.max(pred$fit)
  peak_date       <- pred$date[peak_idx]
  peak_value_fit  <- pred$fit[peak_idx]
  peak_value_low  <- pred$fit_low[peak_idx]
  peak_value_high <- pred$fit_high[peak_idx]
  
  half_thr        <- 0.5 * peak_value_fit
  half_idx        <- which(pred$fit >= half_thr)[1]
  half_date       <- pred$date[half_idx]
  half_value_fit  <- pred$fit[half_idx]
  half_value_low  <- pred$fit_low[half_idx]
  half_value_high <- pred$fit_high[half_idx]
  
  # observed medians (+-3 days) centered on nearest nonna to the fitted date 
  choose_center <- function(target_date) {
    on_day <- dch$date == target_date
    if (any(on_day) && is.finite(dch$daily_median_CH4[on_day][1])) return(target_date)
    diffs <- abs(as.numeric(dch$date - target_date))
    cand  <- which(diffs <= 3 & is.finite(dch$daily_median_CH4))
    if (length(cand)) dch$date[cand[which.min(diffs[cand])]] else as.Date(NA)
  }
  
  # centers
  center_peak <- choose_center(peak_date)
  center_half <- choose_center(half_date)
  
  # observed medians around those centers
  peak_value_obs_med <- if (is.na(center_peak)) NA_real_ else {
    median(dch$daily_median_CH4[dch$date >= center_peak - 3 & dch$date <= center_peak + 3], na.rm = TRUE)
  }
  half_value_obs_med <- if (is.na(center_half)) NA_real_ else {
    median(dch$daily_median_CH4[dch$date >= center_half - 3 & dch$date <= center_half + 3], na.rm = TRUE)
  }
  
  baseline_obs <- median(dch$daily_median_CH4[dch$date <= early_end], na.rm = TRUE)
  # build data frames for baseline markers inside the loop, after you compute
  # early_end, baseline_fit, baseline_obs:
  df_baseline <- tibble(
    x_baseline = season_start + 0.5*(early_end - season_start),                
    y_fit      = baseline_fit,
    y_obs      = baseline_obs
  )
  
  # actual observed date-specific peak (single-day max) 
  if (all(!is.finite(dch$daily_median_CH4))) {
    peak_date_obs  <- as.Date(NA)
    peak_value_obs <- NA_real_
  } else {
    obs_peak_idx   <- which.max(dch$daily_median_CH4)  # first max if ties
    peak_date_obs  <- dch$date[obs_peak_idx]
    peak_value_obs <- dch$daily_median_CH4[obs_peak_idx]
  }
  
  # store row for this chamber
  crit_from_gam[[i]] <- tibble(
    chamber = ch,
    season_start, season_end, early_end,
    baseline_fit, baseline_fit_low, baseline_fit_high, baseline_obs,
    peak_date,  peak_value_fit,  peak_value_low,  peak_value_high,  peak_value_obs_med,
    peak_date_obs, peak_value_obs,
    half_date,  half_value_fit,  half_value_low,  half_value_high,  half_value_obs_med,
    center_peak, center_half,                
    gam_r2 = summary(m)$r.sq,
    edf    = summary(m)$s.table[1, "edf"],
    n = nrow(dch)
  )
  
  # plotting data for markers & windows
  fit_markers <- tibble(
    date  = c(half_date, peak_date),
    value = c(half_value_fit, peak_value_fit),
    type  = c("half_fit", "peak_fit")
  )
  obs_markers <- tibble(
    date  = c(center_half, center_peak),
    value = c(half_value_obs_med, peak_value_obs_med),
    type  = c("half_obs", "peak_obs")
  )
  # vertical window lines (+-3 days from centers), drop na centers
  win_lines <- bind_rows(
    tibble(date = c(center_half - 3, center_half + 3), kind = "half") %>% filter(!is.na(date)),
    tibble(date = c(center_peak - 3, center_peak + 3), kind = "peak") %>% filter(!is.na(date))
  )
  
  # plot
  p_fit <-
    ggplot(dch, aes(x = date, y = daily_median_CH4)) +geom_ribbon(data = pred,
                aes(x = date, ymin = fit_low, ymax = fit_high),
                inherit.aes = FALSE, fill = "red", alpha = 0.15) +
    geom_line(data = pred, aes(x = date, y = fit),
              inherit.aes = FALSE, color = "red", linewidth = 1) +
    geom_point(color = "grey40", alpha = 0.8, size = 1.6) +
    geom_vline(xintercept = early_end, linetype = "dashed") +
    geom_point(
      data = df_baseline,
      aes(x = x_baseline, y = y_fit),
      inherit.aes = FALSE,
      shape = 22, fill = "white", color = "black", stroke = 0.8, size = 3.5
    ) +
    geom_point(
      data = df_baseline,
      aes(x = x_baseline, y = y_obs),
      inherit.aes = FALSE,
      shape = 22, fill = "black", color = "black", stroke = 0.8, size = 3.5
    ) +
    geom_vline(data = win_lines, aes(xintercept = date),
               inherit.aes = FALSE, linetype = "dashed", color = "black", linewidth = 0.4, alpha = 0.8) +
    geom_point(data = fit_markers %>% dplyr::filter(type == "half_fit"),
               aes(x = date, y = value),
               inherit.aes = FALSE, shape = 21, fill = "white", color = "black", size = 3, stroke = 0.7) +
    geom_point(data = fit_markers %>% dplyr::filter(type == "peak_fit"),
               aes(x = date, y = value),
               inherit.aes = FALSE, shape = 24, fill = "white", color = "black", size = 3, stroke = 0.7) +
    geom_point(data = obs_markers %>% dplyr::filter(type == "half_obs" & is.finite(value)),
               aes(x = date, y = value),
               inherit.aes = FALSE, shape = 16, color = "black", size = 3) +
    geom_point(data = obs_markers %>% dplyr::filter(type == "peak_obs" & is.finite(value)),
               aes(x = date, y = value),
               inherit.aes = FALSE, shape = 24, fill = "black", color = "black", size = 3, stroke = 0.7) +
    annotate("text",
           x = min(dch$date) + 5, 
           y = max(dch$daily_median_CH4, na.rm = TRUE),
           label = paste0("Adj. R² = ", round(summary(m)$r.sq, 2)),
           hjust = 0, vjust = 1, size = 4.5) +
    labs(title = paste("Chamber", ch),
         x = "Date",
         y = expression(paste("Daily median CH"[4], " flux"))) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  fit_plots_per_chamber[[i]] <- p_fit
  
  # residual diagnostics
  res_df <- tibble(
    fitted   = fitted(m),  # asinh scale
    resid    = residuals(m)
  )
  
  p_resid_fit <- ggplot(res_df, aes(x = fitted, y = resid)) +
    geom_hline(yintercept = 0, color = "grey60") +
    geom_point(alpha = 0.8) +
    labs(title = "Residuals vs Fitted", x = "Fitted (asinh)", y = "Residuals (asinh)") +
    theme_minimal(base_size = 11)
  
  p_hist <- ggplot(res_df, aes(x = resid)) +
    geom_histogram(bins = 30, fill = "grey70", color = "white") +
    labs(title = "Residuals histogram", x = "Residuals (asinh)", y = "Count") +
    theme_minimal(base_size = 11)
  
  p_qq <- ggplot(res_df, aes(sample = resid)) +
    stat_qq(size = 1) +
    stat_qq_line(color = "red") +
    labs(title = "QQ plot of residuals", x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_minimal(base_size = 11)
  
  resid_plots_per_chamber[[i]] <- p_resid_fit + p_hist + p_qq + plot_layout(ncol = 3)
}

# outputs
crit_from_gam <- bind_rows(crit_from_gam)

# table with critical metrics, CIs, diagnostics, and observed markers + centers
crit_from_gam

crit_from_gam <- as.data.frame(crit_from_gam)

# check for each chamber
resid_plots_per_chamber[[ which(chambers == 1) ]]
resid_plots_per_chamber[[ which(chambers == 2) ]]
resid_plots_per_chamber[[ which(chambers == 3) ]]
resid_plots_per_chamber[[ which(chambers == 4) ]]
resid_plots_per_chamber[[ which(chambers == 5) ]]
resid_plots_per_chamber[[ which(chambers == 6) ]]
resid_plots_per_chamber[[ which(chambers == 7) ]]
resid_plots_per_chamber[[ which(chambers == 8) ]]
resid_plots_per_chamber[[ which(chambers == 9) ]]

# asinh seems good

##### checking chamber 3 gam diagnostics

chambers <- sort(unique(daily_medians$chamber_id))

gam_stats <- vector("list", length(chambers))

for (i in seq_along(chambers)) {
  
  ch  <- chambers[i]
  dch <- daily_medians %>% 
    filter(chamber_id == ch) %>% 
    arrange(date)
  
  ## gam fit for all chambers, including chamber 3
  m <- gam(asinh(daily_median_CH4) ~ s(as.numeric(date)),
           data = dch, method = "REML")
  
  ## model summary
  sm       <- summary(m)
  adj_r2   <- sm$r.sq
  smooth_p <- sm$s.table[1, "p-value"]
  
  ## predictions + CI on a daily grid
  pred <- tibble(
    date   = seq(min(dch$date), max(dch$date), by = "1 day"),
    daynum = as.numeric(date)
  )
  
  pr <- predict(m, newdata = pred, se.fit = TRUE)
  
  pred$fit      <- sinh(pr$fit)
  pred$fit_low  <- sinh(pr$fit - 1.96 * pr$se.fit)
  pred$fit_high <- sinh(pr$fit + 1.96 * pr$se.fit)
  
  ## critical points 
  peak_idx       <- which.max(pred$fit)
  peak_date      <- pred$date[peak_idx]
  peak_value_fit <- pred$fit[peak_idx]
  
  left_min_idx <- which.min(pred$fit[pred$date <= peak_date])
  left_min_val <- pred$fit[pred$date <= peak_date][left_min_idx]
  
  thr10 <- left_min_val + 0.10 * (peak_value_fit - left_min_val)
  thr50 <- left_min_val + 0.50 * (peak_value_fit - left_min_val)
  
  baseline_date <- pred$date[which(pred$fit >= thr10)[1]]
  half_date     <- pred$date[which(pred$fit >= thr50 & pred$date < peak_date)[1]]
  
  baseline_fit   <- pred$fit[pred$date == baseline_date]
  half_value_fit <- pred$fit[pred$date == half_date]
  
  ## indices of critical dates
  base_idx <- which(pred$date == baseline_date)
  half_idx <- which(pred$date == half_date)
  peak_idx <- which(pred$date == peak_date)
  
  ## cis at critical points
  baseline_low  <- pred$fit_low[base_idx]
  baseline_high <- pred$fit_high[base_idx]
  
  half_low  <- pred$fit_low[half_idx]
  half_high <- pred$fit_high[half_idx]
  
  peak_low  <- pred$fit_low[peak_idx]
  peak_high <- pred$fit_high[peak_idx]
  
  ## store results
  gam_stats[[i]] <- tibble(
    chamber = ch,
    adj_r2  = adj_r2,
    smooth_p = smooth_p,
    
    baseline_date,
    baseline_fit,
    baseline_low,
    baseline_high,
    
    half_date,
    half_value_fit,
    half_low,
    half_high,
    
    peak_date,
    peak_value_fit,
    peak_low,
    peak_high
  )
}

crit_summary <- bind_rows(gam_stats) %>%
  arrange(chamber)

crit_summary <- as.data.frame(crit_summary)

crit_summary

# chamber 3 R2 < 0.1 so use moving medians


# fit gams again but with the updated handling of chamber 3

daily_medians <- daily_medians %>% arrange(chamber_id, date)
chambers <- sort(unique(daily_medians$chamber_id))

daily_medians <- daily_medians %>%
  mutate(date = as.Date(date))
crit_from_gam <- list()
fit_plots_per_chamber <- vector("list", length(chambers)); names(fit_plots_per_chamber) <- chambers

for (i in seq_along(chambers)) {
  ch <- chambers[i]
  
  dch <- daily_medians %>% filter(chamber_id == ch) %>% arrange(date)
  
  if (ch != 3) {
    m <- gam(asinh(daily_median_CH4) ~ s(as.numeric(date)), data = dch, method = "REML")
    pred <- tibble(date = seq(min(dch$date), max(dch$date), by="1 day"),
                   daynum = as.numeric(date))
    pr <- predict(m, newdata=pred, se.fit=TRUE)
    pred$fit      <- sinh(pr$fit)
    pred$fit_low  <- sinh(pr$fit - 1.96*pr$se.fit)
    pred$fit_high <- sinh(pr$fit + 1.96*pr$se.fit)
    
    peak_idx <- which.max(pred$fit)
    peak_date <- pred$date[peak_idx]; peak_value_fit <- pred$fit[peak_idx]
    left_min_idx <- which.min(pred$fit[pred$date <= peak_date])
    left_min_val <- pred$fit[pred$date <= peak_date][left_min_idx]
    thr10 <- left_min_val + 0.10*(peak_value_fit - left_min_val)
    thr50 <- left_min_val + 0.50*(peak_value_fit - left_min_val)
    baseline_date <- pred$date[which(pred$fit >= thr10)[1]]
    half_date     <- pred$date[which(pred$fit >= thr50 & pred$date < peak_date)[1]]
    baseline_fit  <- pred$fit[pred$date == baseline_date]
    half_value_fit<- pred$fit[pred$date == half_date]
    
    prb <- pred %>% filter(is.finite(fit_low), is.finite(fit_high))
    pred_ribbon <- if(nrow(prb)>0) prb else NULL
    
    r2_label <- paste0("Adj. R\u00B2 = ", round(summary(m)$r.sq, 2))
    
  } else {
    # chamber 3
    pred <- tibble(date = seq(min(dch$date), max(dch$date), by="1 day")) %>%
      left_join(dch %>% select(date, daily_median_CH4), by="date") %>%
      arrange(date)
    # ±7 day rolling median
    pred$fit <- zoo::rollapply(pred$daily_median_CH4, width=7, FUN=median,
                               align="center", partial=TRUE, na.rm=TRUE)
    
    peak_idx <- which.max(pred$fit)
    peak_date <- pred$date[peak_idx]; peak_value_fit <- pred$fit[peak_idx]
    left_min_idx <- which.min(pred$fit[pred$date <= peak_date])
    left_min_val <- pred$fit[pred$date <= peak_date][left_min_idx]
    thr10 <- left_min_val + 0.10*(peak_value_fit - left_min_val)
    thr50 <- left_min_val + 0.50*(peak_value_fit - left_min_val)
    baseline_date <- pred$date[which(pred$fit >= thr10)[1]]
    half_date     <- pred$date[which(pred$fit >= thr50 & pred$date < peak_date)[1]]
    baseline_fit  <- pred$fit[pred$date == baseline_date]
    half_value_fit<- pred$fit[pred$date == half_date]
    
    pred_ribbon <- NULL  # No ribbon for chamber 3
    
    # no r2 for chamber 3
    r2_label <- NULL
  }
  
  # store
  crit_from_gam[[i]] <- tibble(
    chamber = ch,
    baseline_date, baseline_fit,
    half_date,     half_value_fit,
    peak_date,     peak_value_fit
  )
  
  # plot
  p_fit <- ggplot(dch, aes(date, daily_median_CH4)) +
    { if (!is.null(pred_ribbon))
      geom_ribbon(data=pred_ribbon, aes(date, ymin=fit_low, ymax=fit_high),
                  inherit.aes=FALSE, fill="red", alpha=0.15)
      else geom_blank()
    } +
    geom_line(data=pred, aes(date, y=fit), color="red", linewidth=1) +
    geom_point(color="grey40", alpha=0.8, size=1.6) +
    geom_point(data=tibble(date=baseline_date, value=baseline_fit),
               aes(date,value), inherit.aes=FALSE,
               shape=22, fill="white", color="black", size=3.5) +
    geom_point(data=tibble(date=half_date, value=half_value_fit),
               aes(date,value), inherit.aes=FALSE,
               shape=21, fill="white", color="black", size=3) +
    geom_point(data=tibble(date=peak_date, value=peak_value_fit),
               aes(date,value), inherit.aes=FALSE,
               shape=24, fill="white", color="black", size=3) +
  { if (!is.null(r2_label))
    annotate("text",
             x = min(dch$date) + 5,
             y = max(dch$daily_median_CH4, na.rm = TRUE),
             label = r2_label,
             hjust = 0, vjust = 1, size = 4.5)
    else
      geom_blank()
  } +
    labs(title=paste("Chamber", ch),
         x="", y=expression(paste("Daily median CH"[4]," flux"))) +
    theme_bw()+
    theme(axis.text = element_text(size=18),
          axis.title = element_text(size=18))
  
  fit_plots_per_chamber[[i]] <- p_fit
}

# combine results
crit_from_gam <- bind_rows(crit_from_gam)
combined_fit <- wrap_plots(fit_plots_per_chamber, ncol=3)
combined_fit

#### GAM fit stats

gam_stats <- vector("list", length(chambers))

for (i in seq_along(chambers)) {
  ch  <- chambers[i]
  dch <- daily_medians %>% filter(chamber_id == ch) %>% arrange(date)
  
  if (ch != 3) {
    ## gam fit
    m <- gam(asinh(daily_median_CH4) ~ s(as.numeric(date)),
             data = dch, method = "REML")
    
    ## model summary
    sm       <- summary(m)
    adj_r2   <- sm$r.sq                    
    smooth_p <- sm$s.table[1, "p-value"]  
    
    ## predictions + ci on a daily grid
    pred <- tibble(
      date   = seq(min(dch$date), max(dch$date), by = "1 day"),
      daynum = as.numeric(date)
    )
    
    pr <- predict(m, newdata = pred, se.fit = TRUE)
    pred$fit      <- sinh(pr$fit)
    pred$fit_low  <- sinh(pr$fit - 1.96 * pr$se.fit)
    pred$fit_high <- sinh(pr$fit + 1.96 * pr$se.fit)
    
    ## critical points
    peak_idx       <- which.max(pred$fit)
    peak_date      <- pred$date[peak_idx]
    peak_value_fit <- pred$fit[peak_idx]
    
    left_min_idx <- which.min(pred$fit[pred$date <= peak_date])
    left_min_val <- pred$fit[pred$date <= peak_date][left_min_idx]
    
    thr10 <- left_min_val + 0.10 * (peak_value_fit - left_min_val)
    thr50 <- left_min_val + 0.50 * (peak_value_fit - left_min_val)
    
    baseline_date <- pred$date[which(pred$fit >= thr10)[1]]
    half_date     <- pred$date[which(pred$fit >= thr50 & pred$date < peak_date)[1]]
    
    baseline_fit   <- pred$fit[pred$date == baseline_date]
    half_value_fit <- pred$fit[pred$date == half_date]
    
    ## indices of those dates in pred
    base_idx <- which(pred$date == baseline_date)
    half_idx <- which(pred$date == half_date)
    peak_idx <- which(pred$date == peak_date)
    
    ## cis at the critical points
    baseline_low  <- pred$fit_low[base_idx]
    baseline_high <- pred$fit_high[base_idx]
    
    half_low  <- pred$fit_low[half_idx]
    half_high <- pred$fit_high[half_idx]
    
    peak_low  <- pred$fit_low[peak_idx]
    peak_high <- pred$fit_high[peak_idx]
    
  } else {
    ## chamber 3 (rolling median, no gam)
    pred <- tibble(date = seq(min(dch$date), max(dch$date), by = "1 day")) %>%
      left_join(dch %>% select(date, daily_median_CH4), by = "date") %>%
      arrange(date)
    
    pred$fit <- zoo::rollapply(
      pred$daily_median_CH4,
      width   = 7,
      FUN     = median,
      align   = "center",
      partial = TRUE,
      na.rm   = TRUE
    )
    
    peak_idx       <- which.max(pred$fit)
    peak_date      <- pred$date[peak_idx]
    peak_value_fit <- pred$fit[peak_idx]
    
    left_min_idx <- which.min(pred$fit[pred$date <= peak_date])
    left_min_val <- pred$fit[pred$date <= peak_date][left_min_idx]
    
    thr10 <- left_min_val + 0.10 * (peak_value_fit - left_min_val)
    thr50 <- left_min_val + 0.50 * (peak_value_fit - left_min_val)
    
    baseline_date <- pred$date[which(pred$fit >= thr10)[1]]
    half_date     <- pred$date[which(pred$fit >= thr50 & pred$date < peak_date)[1]]
    
    baseline_fit   <- pred$fit[pred$date == baseline_date]
    half_value_fit <- pred$fit[pred$date == half_date]

    adj_r2      <- NA_real_
    smooth_p    <- NA_real_
    baseline_low <- baseline_high <- NA_real_
    half_low     <- half_high     <- NA_real_
    peak_low     <- peak_high     <- NA_real_
  }
  
  gam_stats[[i]] <- tibble(
    chamber = ch,
    adj_r2  = adj_r2,
    smooth_p = smooth_p,
    
    baseline_date,
    baseline_fit,
    baseline_low,
    baseline_high,
    
    half_date,
    half_value_fit,
    half_low,
    half_high,
    
    peak_date,
    peak_value_fit,
    peak_low,
    peak_high
  )
}

## final summary table with all chambers (used in Table C5)
crit_summary <- bind_rows(gam_stats) %>%
  arrange(chamber)

crit_summary <- as.data.frame(crit_summary)
crit_summary


### general ch4 season figure with critical points per thaw stage (Fig. B12)

daily_medians <- daily_medians %>%
  arrange(chamber_id, date) %>%
  mutate(
    daynum  = as.numeric(date),
    chamber = as.factor(chamber_id)
  )

# fit a single gam with chamberspecific smooths (exclude chamber 3)
m <- daily_medians %>%
  filter(chamber != 3) %>%
  gam(asinh(daily_median_CH4) ~ chamber + s(daynum, by = chamber), data = ., method = "REML")

# prediction grid for gam chambers
pred_gam <- daily_medians %>%
  filter(chamber != 3) %>%
  group_by(chamber) %>%
  summarise(date = seq(min(date), max(date), by = "1 day"), .groups = "drop") %>%
  mutate(daynum = as.numeric(date)) %>%
  arrange(chamber, date)

pr <- predict(m, newdata = pred_gam, se.fit = TRUE)
pred_gam <- pred_gam %>%
  mutate(
    fit      = sinh(pr$fit),
    fit_low  = sinh(pr$fit - 1.96 * pr$se.fit),
    fit_high = sinh(pr$fit + 1.96 * pr$se.fit),
    has_ribbon = TRUE
  )

# chamber 3: rolling median
ch3 <- daily_medians %>% filter(chamber == 3)
pred_ch3 <- ch3 %>%
  summarise(date = seq(min(date), max(date), by = "1 day"), .groups = "drop") %>%
  left_join(ch3 %>% select(date, daily_median_CH4), by = "date") %>%
  arrange(date) %>%
  mutate(
    daynum = as.numeric(date),
    fit    = zoo::rollapply(daily_median_CH4, width = 7, FUN = median,
                            align = "center", partial = TRUE, na.rm = TRUE),
    fit_low = NA_real_, fit_high = NA_real_,
    chamber = factor(3, levels = levels(daily_medians$chamber)),
    has_ribbon = FALSE
  ) %>%
  select(chamber, date, daynum, fit, fit_low, fit_high, has_ribbon)

# combine predictions from all chambers
pred_all <- bind_rows(pred_gam, pred_ch3)

# compute critical points per chamber (peak, leftmin, baseline 10%, half 50%)
# peak per chamber
peaks <- pred_all %>%
  group_by(chamber) %>%
  slice_max(fit, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(chamber,
            peak_date = date,
            peak_value_fit = fit)

# leftmin before (or at) peak
left_mins <- pred_all %>%
  inner_join(peaks, by = "chamber") %>%
  filter(date <= peak_date) %>%
  group_by(chamber) %>%
  slice_min(fit, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(chamber, left_min_val = fit, left_min_date = date)

# thresholds
thr <- peaks %>%
  inner_join(left_mins, by = "chamber") %>%
  mutate(
    thr10 = left_min_val + 0.10 * (peak_value_fit - left_min_val),
    thr50 = left_min_val + 0.50 * (peak_value_fit - left_min_val)
  )

# baseline (first date where fit >= 10% threshold, before/at peak)
baseline_pts <- pred_all %>%
  inner_join(thr, by = "chamber") %>%
  filter(date <= peak_date, fit >= thr10) %>%
  group_by(chamber) %>%
  arrange(date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(chamber, baseline_date = date, baseline_fit = fit)

# halfrise (first date where fit >= 50% threshold, strictly before peak)
half_pts <- pred_all %>%
  inner_join(thr, by = "chamber") %>%
  filter(date < peak_date, fit >= thr50) %>%
  group_by(chamber) %>%
  arrange(date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(chamber, half_date = date, half_value_fit = fit)

crit_pts <- peaks %>%
  left_join(baseline_pts, by = "chamber") %>%
  left_join(half_pts,     by = "chamber")

# create thaw_stage column to daily_medians
daily_medians <- daily_medians %>% 
  mutate(thaw_stage = case_when(
    chamber_id %in% c(1, 3, 5) ~ "intact",         
    chamber_id %in% c(2, 4, 6) ~ "partly_thawed",   
    chamber_id %in% c(7, 8, 9) ~ "fully_thawed" 
  ))

# intact
d_stage <- daily_medians %>% filter(thaw_stage == "intact")
pred_stage <- pred_all %>% semi_join(d_stage %>% distinct(chamber), by = "chamber")
crit_stage <- crit_pts  %>% semi_join(d_stage %>% distinct(chamber), by = "chamber")

p_intact <- ggplot(d_stage, aes(date, daily_median_CH4)) +
  geom_point(color = "grey40", alpha = 0.8, size = 1.6) +
  geom_line(data = pred_stage, aes(x=date, y = fit, group = chamber),
            inherit.aes = FALSE, linewidth = 1, color = "#E30B5C", alpha=0.7) +
  geom_point(data = crit_stage %>% select(chamber, date = baseline_date, value = baseline_fit),
             aes(date, value), inherit.aes = FALSE,
             shape = 22, fill = "white", color = "black", size = 3.5, stroke = 0.7) +
  geom_point(data = crit_stage %>% select(chamber, date = half_date, value = half_value_fit),
             aes(date, value), inherit.aes = FALSE,
             shape = 21, fill = "white", color = "black", size = 3, stroke = 0.7) +
  geom_point(data = crit_stage %>% select(chamber, date = peak_date, value = peak_value_fit),
             aes(date, value), inherit.aes = FALSE,
             shape = 24, fill = "white", color = "black", size = 3, stroke = 0.7) +
  labs(title = "Intact", x = "", y = expression(paste("Daily median CH"[4], " flux"))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text  = element_text(size = 18),
        axis.title = element_text(size = 18))

# bog
d_stage <- daily_medians %>% filter(thaw_stage == "partly_thawed")
pred_stage <- pred_all %>% semi_join(d_stage %>% distinct(chamber), by = "chamber")
crit_stage <- crit_pts  %>% semi_join(d_stage %>% distinct(chamber), by = "chamber")

p_partly_thawed <- ggplot(d_stage, aes(date, daily_median_CH4)) +
  geom_point(color = "grey40", alpha = 0.8, size = 1.6) +
  geom_line(data = pred_stage, aes(x = date, y = fit, group = chamber),
            inherit.aes = FALSE, linewidth = 1, color = "#E30B5C", alpha=0.7) +
  geom_point(data = crit_stage %>% select(chamber, date = baseline_date, value = baseline_fit),
             aes(date, value), inherit.aes = FALSE,
             shape = 22, fill = "white", color = "black", size = 3.5, stroke = 0.7) +
  geom_point(data = crit_stage %>% select(chamber, date = half_date, value = half_value_fit),
             aes(date, value), inherit.aes = FALSE,
             shape = 21, fill = "white", color = "black", size = 3, stroke = 0.7) +
  geom_point(data = crit_stage %>% select(chamber, date = peak_date, value = peak_value_fit),
             aes(date, value), inherit.aes = FALSE,
             shape = 24, fill = "white", color = "black", size = 3, stroke = 0.7) +
  labs(title = "Partly thawed", x = "", y = expression(paste("Daily median CH"[4], " flux"))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text  = element_text(size = 18),
        axis.title = element_text(size = 18))

# fully thawed
d_stage <- daily_medians %>% filter(thaw_stage == "fully_thawed")
pred_stage <- pred_all %>% semi_join(d_stage %>% distinct(chamber), by = "chamber")
crit_stage <- crit_pts  %>% semi_join(d_stage %>% distinct(chamber), by = "chamber")

p_fully_thawed <- ggplot(d_stage, aes(date, daily_median_CH4)) +
  geom_point(color = "grey40", alpha = 0.8, size = 1.6) +
  geom_line(data = pred_stage, aes(x = date, y = fit, group = chamber),
            inherit.aes = FALSE, linewidth = 1, color = "#E30B5C", alpha=0.7) +
  geom_point(data = crit_stage %>% select(chamber, date = baseline_date, value = baseline_fit),
             aes(date, value), inherit.aes = FALSE,
             shape = 22, fill = "white", color = "black", size = 3.5, stroke = 0.7) +
  geom_point(data = crit_stage %>% select(chamber, date = half_date, value = half_value_fit),
             aes(date, value), inherit.aes = FALSE,
             shape = 21, fill = "white", color = "black", size = 3, stroke = 0.7) +
  geom_point(data = crit_stage %>% select(chamber, date = peak_date, value = peak_value_fit),
             aes(date, value), inherit.aes = FALSE,
             shape = 24, fill = "white", color = "black", size = 3, stroke = 0.7) +
  labs(title = "Fully thawed", x = "", y = expression(paste("Daily median CH"[4], " flux"))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text  = element_text(size = 18),
        axis.title = element_text(size = 18))

# combine
wrap_plots(p_intact, p_partly_thawed, p_fully_thawed, ncol = 3)


# get season medians

flux_daily <- as.data.frame(flux_daily)

# season median ch4 per chamber
season_medians <- flux_daily %>%
  group_by(chamber_id) %>%
  summarise(season_median_ch4 = median(daily_median_CH4, na.rm = TRUE), .groups = "drop")

season_means_Tgnd <- flux_daily %>%
  group_by(chamber_id) %>%
  summarise(season_mean_avgTgnd = mean(avgTgnd_mean, na.rm = TRUE), 
            .groups = "drop")

colnames(season_means_Tgnd)[1] <- "chamber"
colnames(season_medians)[1] <- "chamber"

# attach season median to critical ch4 table
crit_ch4 <- crit_from_gam %>%
  left_join(season_medians, by = "chamber")

# long format for critical dates and combine with flux_daily
colnames(flux_daily)[1] <- "chamber"

crit_long <- crit_ch4 %>%
  select(chamber, baseline_date, half_date, peak_date) %>%
  pivot_longer(cols = c(baseline_date, half_date, peak_date),
               names_to = "crit_type", values_to = "crit_date") %>%
  left_join(flux_daily %>% select(chamber, thaw_stage) %>% distinct(), by = "chamber")

# soil temp from the exact critical date
temp_crit <- crit_long %>%
  rowwise() %>%
  mutate(
    temp_row = list({
      if (is.na(crit_date)) {
        tibble(avgTgnd_mean=NA_real_)
      } else {
        d <- flux_daily %>% filter(chamber == chamber, date == crit_date)
        if (nrow(d) == 0) {
          tibble(avgTgnd=NA_real_)
        } else {
          summarise(d,
                    avgTgnd  = mean(avgTgnd_mean,  na.rm = TRUE))
        }
      }
    })
  ) %>%
  ungroup() %>%
  unnest(temp_row)

# combine ch4 critical values and temp means
complete_df <- temp_crit %>%
  left_join(
    crit_ch4 %>% select(chamber, baseline_fit, half_value_fit, peak_value_fit, season_median_ch4),
    by = "chamber"
  ) %>%
  select(chamber, thaw_stage, crit_type, crit_date,
         baseline_fit, half_value_fit, peak_value_fit, season_median_ch4,
         avgTgnd)

complete_df

# flip and rename crit_type values and attach the matching ch4 value for that crit
flipped <- complete_df %>%
  mutate(
    crit_label = recode(crit_type,
                        baseline_date = "early",
                        half_date     = "middle",
                        peak_date     = "peak"),
    ch4_value_for_crit = case_when(
      crit_label == "early"  ~ baseline_fit,
      crit_label == "middle" ~ half_value_fit,
      crit_label == "peak"   ~ peak_value_fit,
      TRUE ~ NA_real_
    ),
    ch4_date_for_crit = crit_date
  ) %>%
  # keep only needed columns for widening
  select(chamber, thaw_stage,
         crit_label,
         ch4_value_for_crit, ch4_date_for_crit,
         season_median_ch4,
         avgTgnd)

# pivot wider: ch4 value/date and temp get suffixes by crit_label
wide <- flipped %>%
  pivot_wider(
    id_cols = c(chamber, thaw_stage, season_median_ch4),
    names_from  = crit_label,
    values_from = c(ch4_value_for_crit, ch4_date_for_crit,
                    avgTgnd),
    names_glue = "{.value}_{crit_label}"
  )

# rename columns
final_df <- wide %>%
  rename(
    ch4_early        = ch4_value_for_crit_early,
    ch4_middle       = ch4_value_for_crit_middle,
    ch4_peak         = ch4_value_for_crit_peak,
    ch4_season = season_median_ch4,
    ch4_early_date   = ch4_date_for_crit_early,
    ch4_middle_date  = ch4_date_for_crit_middle,
    ch4_peak_date    = ch4_date_for_crit_peak
  )

final_df <- as.data.frame(final_df)


## combine with the trait data

# read in the root-vga data

plant <- read.csv("path/root_vga_plot_scale.csv")

# rename some columns

final_df <- final_df %>%
  rename(
    chamber_id = chamber
  )

# combine with ch4

plant_ch4 <- plant %>%
  left_join(final_df, by = c("chamber_id", "thaw_stage"))


write.csv(plant_ch4, "path/plant_ch4_neefiltered_ch_critical_points_QAQCd.csv", row.names = FALSE)
plant_ch4 <- read.csv("path/plant_ch4_neefiltered_ch_critical_points_QAQCd.csv")

## calculate weighted pft traits in relation to total root biomass

## weighted means for rtd, avg_d_mm and srl by pft biomass

# which traits to sum/mean when collapsing across pfts 
vars_sum  <- c("mass_gm2_std", "sa_cm2m2_std", "len_mm2_std")
vars_mean <- c("avg_d_mm", "RTD_gcm3", "SRL_mg1")

plant_ch4_noPFT <- plant_ch4 %>%
  group_by(chamber_id, thaw_stage, organ) %>%
  summarise(
    across(all_of(vars_mean), ~ weighted.mean(.x, .data$mass_gm2_std, na.rm = TRUE)),
    across(all_of(vars_sum),  ~ sum(.x, na.rm = TRUE)),
    across(
      c(contains("ch4_"),
        contains("GA"),
        starts_with("avgTgnd_")),
      ~ dplyr::first(.x[!is.na(.x)], default = NA)
    ),
    .groups = "drop"
  )
plant_ch4_noPFT <- as.data.frame(plant_ch4_noPFT)

## save
write.csv(plant_ch4_noPFT, "path/root_rhizome_noPFT_chamber_level_QAQCd.csv", row.names = FALSE)

#############################################################
###--------------------POREWATER DATA---------------------###
#############################################################

pw_ch4 <- read.csv("path/Emerge_2023_porewater_1_standardized_autochamber.csv")

# cleaning the df
pw_ch4 <- pw_ch4 %>%
  dplyr::filter(!DepthMin__ %in% c(30, 40, 50, 60, 70))

pw_ch4 <- pw_ch4 %>%
  dplyr::filter(!DepthMin__ %in% c(80))

pw_ch4 <- pw_ch4 %>%
  mutate(thaw_stage = case_when(
    str_detect(Site__, "Eriophorum") ~ "fully_thawed",
    str_detect(Site__, "Sphagnum") ~ "partly_thawed"
  ))

pw_ch4 <- pw_ch4 %>%
  mutate(depth_std = case_when(
    DepthMin__ == 1 ~ "0-10",
    DepthMin__ == 10 ~ "10-20",
    DepthMin__ == 20 ~ "20-30"
  ))

# with all depths
pw_ch4_2 <- pw_ch4 %>%
  mutate(depth_std = case_when(
    DepthMin__ == 1 ~ "1-5",
    DepthMin__ == 10 ~ "10-14",
    DepthMin__ == 20 ~ "20-24",
    DepthMin__ == 30 ~ "30-34",
    DepthMin__ == 40 ~ "40-44",
    DepthMin__ == 50 ~ "50-54",
    DepthMin__ == 60 ~ "60-64",
    DepthMin__ == 70 ~ "70-74",
    DepthMin__ == 80 ~ "80-84"
  ))


# remove extra columns

pw_ch4 <- subset(pw_ch4, select = -c(Site__, FieldSampling__, Date__, DepthMin__, DepthMax__))
pw_ch4_2 <- subset(pw_ch4_2, select = -c(Site__, FieldSampling__, Date__, DepthMin__, DepthMax__))

# rename

pw_ch4 <- pw_ch4 %>%
  rename(Core_ID = Core__,
         d13C_CH4 = d13C_CH4__,
         CH4_conc_mM = CH4.mM__,
         d13C_CO2 = d13C_CO2__,
         CO2_conc_mM = CO2.mM__)

pw_ch4_2 <- pw_ch4_2 %>%
  rename(Core_ID = Core__,
         d13C_CH4 = d13C_CH4__,
         CH4_conc_mM = CH4.mM__,
         d13C_CO2 = d13C_CO2__,
         CO2_conc_mM = CO2.mM__)

# save

write.csv(pw_ch4, "path/porewater_ch4.csv", row.names = FALSE)
write.csv(pw_ch4_2, "path/porewater_ch4_alldepths.csv", row.names = FALSE)

#############################################################
###-----------------PEAT MOISTURE-------------------------###
#############################################################


#  1) read the semicolon-delimited raw file exactly asis 
# your decimals are "." so we keep default locale
raw_path <- "path/conv_data_95140961_2025_10_09_0_B2.csv"

# many tomst "conv_data_*.csv" rows end with a trailing ";" > an empty last column.
# read_delim handles that fine; we then drop empty trailing columns if present.
df0 <- read_delim(raw_path, delim = ";", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# drop columns that are completely empty (e.g., the trailing ";")
df0 <- df0 %>% select(where(~ !all(is.na(.x) | .x == "")))

# rename columns
nm <- c("idx","datetime_local","tz","T1_C","T2_C","T3_C","moisture_raw","VWC_lolly","shake")
names(df0)[seq_len(min(length(nm), ncol(df0)))] <- nm[seq_len(min(length(nm), ncol(df0)))]

# set datetime with local tz and then convert to utc
df <- df0 %>%
  mutate(
    datetime_local = as.POSIXct(datetime_local, format = "%Y.%m.%d %H:%M", tz = "Europe/Stockholm"),
    datetime = with_tz(datetime_local, "UTC")
  )

# build a longformat table for myclim (one row per sensor reading)
# locality_id = plot id
long_tbl <- bind_rows(
  df %>% transmute(locality_id = "B2", sensor_name = "TMS_T1",    datetime, value = as.numeric(T1_C)),
  df %>% transmute(locality_id = "B2", sensor_name = "TMS_T2",    datetime, value = as.numeric(T2_C)),
  df %>% transmute(locality_id = "B2", sensor_name = "TMS_T3",    datetime, value = as.numeric(T3_C)),
  df %>% transmute(locality_id = "B2", sensor_name = "TMS_moist", datetime, value = as.numeric(moisture_raw))
)

# make myclim version for long table
data_mc <- mc_read_long(
  long_tbl,
  clean  = TRUE,
  silent = TRUE
)

# peat calibration curve parameters: 
mc_data_vwc_parameters %>% dplyr::filter(soiltype == "peat")

# convert raw moisture to vwc using the peat curve + T1 soil temp compensation
data_mc <- mc_calc_vwc(
  data_mc,
  moist_sensor  = "TMS_moist",
  temp_sensor   = "TMS_T1",
  soiltype      = "peat", # peat calibration curve
  output_sensor = "VWC_peat",
  frozen2NA     = TRUE # sets frozen periods to NA
)

# convert to df
vwc_df <- mc_reshape_long(data_mc, sensors = c("VWC_peat")) %>%
  select(datetime, locality_id, sensor_name, value) %>%
  rename(VWC_peat = value)

### subset to may-aug 2024-2025

# work on a copy with a localtime column for filtering clarity
vwc_df <- vwc_df %>%
  mutate(datetime_se = with_tz(datetime, "Europe/Stockholm"))

vwc_may_aug_24_25 <- vwc_df %>%
  filter(year(datetime_se) %in% c(2024, 2025),
         month(datetime_se) %in% 5:8)

# how often are you capped?
vwc_flagged <- vwc_df %>%
  mutate(is_capped = !is.na(VWC_peat) & VWC_peat >= 0.995)

vwc_flagged %>%
  summarise(prop_capped = mean(is_capped, na.rm = TRUE))

# plot raw moisture vs vwc to see the plateau at vwc = 1.0 (100%)
# this shows that the peat calibration curve does not capture some
# soil moisture variation at high moisture levels

raw_moist <- mc_reshape_long(data_mc, sensors = "TMS_moist") %>%
  select(datetime, raw_moist = value)

left_join(raw_moist, vwc_df, by = "datetime") %>%
  ggplot(aes(raw_moist, VWC_peat)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_hline(yintercept = 1.0, linetype = 2) +
  labs(x = "TMS raw moisture", y = "VWC (peat)") +
  theme_minimal()

## flag vwc values that are close to 100%

vwc_df <- vwc_df %>%
  mutate(
    VWC_pct = VWC_peat * 100,
    unc_flag_high = VWC_peat >= 0.995
  )

vwc_df_B2 <- vwc_df
vwc_may_aug_24_25_B2 <- vwc_may_aug_24_25


#  1) read the semicolon-delimited raw file exactly asis 
# your decimals are "." so we keep default locale
raw_path <- "path/conv_data_95140964_2025_10_09_0_B1.csv"

# many tomst "conv_data_*.csv" rows end with a trailing ";" > an empty last column.
# read_delim handles that fine; we then drop empty trailing columns if present.
df0 <- read_delim(raw_path, delim = ";", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# drop columns that are completely empty (e.g., the trailing ";")
df0 <- df0 %>% select(where(~ !all(is.na(.x) | .x == "")))

# expecting 9 useful fields like:
nm <- c("idx","datetime_local","tz","T1_C","T2_C","T3_C","moisture_raw","VWC_lolly","shake")
names(df0)[seq_len(min(length(nm), ncol(df0)))] <- nm[seq_len(min(length(nm), ncol(df0)))]

#  2) parse local time (stockholm) and convert to utc (myclim expects utc) 
df <- df0 %>%
  mutate(
    datetime_local = as.POSIXct(datetime_local, format = "%Y.%m.%d %H:%M", tz = "Europe/Stockholm"),
    datetime = with_tz(datetime_local, "UTC")
  )

#  3) build a long table for myclim (one row per sensor reading) 
# give your site a locality_id (e.g., "b2")
long_tbl <- bind_rows(
  df %>% transmute(locality_id = "B1", sensor_name = "TMS_T1",    datetime, value = as.numeric(T1_C)),
  df %>% transmute(locality_id = "B1", sensor_name = "TMS_T2",    datetime, value = as.numeric(T2_C)),
  df %>% transmute(locality_id = "B1", sensor_name = "TMS_T3",    datetime, value = as.numeric(T3_C)),
  df %>% transmute(locality_id = "B1", sensor_name = "TMS_moist", datetime, value = as.numeric(moisture_raw))
)

# tell myclim which physical units these sensors are in
data_mc <- mc_read_long(
  long_tbl,
  clean  = TRUE,
  silent = TRUE
)

#  5) convert raw moisture to vwc using the peat curve + t1 temp compensation 
data_mc <- mc_calc_vwc(
  data_mc,
  moist_sensor  = "TMS_moist",
  temp_sensor   = "TMS_T1",
  soiltype      = "peat",        # peat-specific curve
  output_sensor = "VWC_peat",
  frozen2NA     = TRUE
)

#  6) back to a data.frame
vwc_df <- mc_reshape_long(data_mc, sensors = c("VWC_peat")) %>%
  select(datetime, locality_id, sensor_name, value) %>%
  rename(VWC_peat = value)

# flag high vwc values close to vwc = 100%
vwc_df <- vwc_df %>%
  mutate(
    VWC_pct = VWC_peat * 100,
    unc_flag = VWC_peat >= 0.995
  )

# work on a copy with a local time column for filtering clarity
vwc_df <- vwc_df %>%
  mutate(datetime_se = with_tz(datetime, "Europe/Stockholm"))

vwc_may_aug_24_25 <- vwc_df %>%
  filter(year(datetime_se) %in% c(2024, 2025),
         month(datetime_se) %in% 5:8)

ggplot(vwc_may_aug_24_25 %>% filter(year(vwc_may_aug_24_25$datetime) == 2025), aes(x = datetime, y = VWC_peat * 100)) +
  geom_point(color = "steelblue", size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)",
    title = "Peat soil volumetric water content (VWC)"
  ) +
  theme_minimal(base_size = 12)


# plot raw moisture vs vwc to see the plateau (needs raw_moist from your myclim object)
raw_moist <- mc_reshape_long(data_mc, sensors = "TMS_moist") %>%
  select(datetime, raw_moist = value)

left_join(raw_moist, vwc_df, by = "datetime") %>%
  ggplot(aes(raw_moist, VWC_peat)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_hline(yintercept = 1.0, linetype = 2) +
  labs(x = "TMS raw moisture", y = "VWC (peat)") +
  theme_minimal()

vwc_df_B1 <- vwc_df
vwc_may_aug_24_25_B1 <- vwc_may_aug_24_25


#  1) read the semicolon-delimited raw file exactly asis 
# your decimals are "." so we keep default locale
raw_path <- "path/conv_data_95140965_2025_10_09_0_B3.csv"

# many tomst "conv_data_*.csv" rows end with a trailing ";" > an empty last column.
# read_delim handles that fine; we then drop empty trailing columns if present.
df0 <- read_delim(raw_path, delim = ";", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# drop columns that are completely empty (e.g., the trailing ";")
df0 <- df0 %>% select(where(~ !all(is.na(.x) | .x == "")))

# expecting 9 useful fields like:
# x1=index, x2=datetime, x3=code(=4), x4=t1, x5=t2, x6=t3, x7=moist_raw, x8=vwc_lolly?, x9=battery?
nm <- c("idx","datetime_local","tz","T1_C","T2_C","T3_C","moisture_raw","VWC_lolly","shake")
names(df0)[seq_len(min(length(nm), ncol(df0)))] <- nm[seq_len(min(length(nm), ncol(df0)))]

#  2) parse local time (stockholm) and convert to utc (myclim expects utc) 
df <- df0 %>%
  mutate(
    datetime_local = as.POSIXct(datetime_local, format = "%Y.%m.%d %H:%M", tz = "Europe/Stockholm"),
    datetime = with_tz(datetime_local, "UTC")
  )

#  3) build a long table for myclim (one row per sensor reading) 
# give your site a locality_id (e.g., "b2")
long_tbl <- bind_rows(
  df %>% transmute(locality_id = "B3", sensor_name = "TMS_T1",    datetime, value = as.numeric(T1_C)),
  df %>% transmute(locality_id = "B3", sensor_name = "TMS_T2",    datetime, value = as.numeric(T2_C)),
  df %>% transmute(locality_id = "B3", sensor_name = "TMS_T3",    datetime, value = as.numeric(T3_C)),
  df %>% transmute(locality_id = "B3", sensor_name = "TMS_moist", datetime, value = as.numeric(moisture_raw))
)

# tell myclim which physical units these sensors are in
data_mc <- mc_read_long(
  long_tbl,
  clean  = TRUE,
  silent = TRUE
)


#  5) convert raw moisture to vwc using the peat curve + t1 temp compensation 
data_mc <- mc_calc_vwc(
  data_mc,
  moist_sensor  = "TMS_moist",
  temp_sensor   = "TMS_T1",
  soiltype      = "peat",        # peat-specific curve
  output_sensor = "VWC_peat",
  frozen2NA     = TRUE
)

#  6) back to a data.frame if you want to join/plot 
vwc_df <- mc_reshape_long(data_mc, sensors = c("VWC_peat")) %>%
  select(datetime, locality_id, sensor_name, value) %>%
  rename(VWC_peat = value)
## flag high vwc
vwc_df <- vwc_df %>%
  mutate(
    VWC_pct = VWC_peat * 100,
    unc_flag = VWC_peat >= 0.995
  )

# work on a copy with a local time column for filtering clarity
vwc_df <- vwc_df %>%
  mutate(datetime_se = with_tz(datetime, "Europe/Stockholm"))

# 2) may–august for years 2024 and 2025 (local)
vwc_may_aug_24_25 <- vwc_df %>%
  filter(year(datetime_se) %in% c(2024, 2025),
         month(datetime_se) %in% 5:8)


ggplot(vwc_may_aug_24_25 %>% filter(year(vwc_may_aug_24_25$datetime) == 2025), aes(x = datetime, y = VWC_peat * 100)) +
  geom_point(color = "steelblue", size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)",
    title = "Peat soil volumetric water content (VWC)"
  ) +
  theme_minimal(base_size = 12)

# plot raw moisture vs vwc to see the plateau (needs raw_moist from your myclim object)
raw_moist <- mc_reshape_long(data_mc, sensors = "TMS_moist") %>%
  select(datetime, raw_moist = value)

left_join(raw_moist, vwc_df, by = "datetime") %>%
  ggplot(aes(raw_moist, VWC_peat)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_hline(yintercept = 1.0, linetype = 2) +
  labs(x = "TMS raw moisture", y = "VWC (peat)") +
  theme_minimal()

vwc_df_B3 <- vwc_df
vwc_may_aug_24_25_B3 <- vwc_may_aug_24_25


#  1) read the semicolon-delimited raw file exactly asis 
# your decimals are "." so we keep default locale
raw_path <- "path/conv_data_95140970_2025_10_09_0_F1.csv"

# many tomst "conv_data_*.csv" rows end with a trailing ";" > an empty last column.
# read_delim handles that fine; we then drop empty trailing columns if present.
df0 <- read_delim(raw_path, delim = ";", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# drop columns that are completely empty (e.g., the trailing ";")
df0 <- df0 %>% select(where(~ !all(is.na(.x) | .x == "")))

# expecting 9 useful fields like:
# x1=index, x2=datetime, x3=code(=4), x4=t1, x5=t2, x6=t3, x7=moist_raw, x8=vwc_lolly?, x9=battery?
nm <- c("idx","datetime_local","tz","T1_C","T2_C","T3_C","moisture_raw","VWC_lolly","shake")
names(df0)[seq_len(min(length(nm), ncol(df0)))] <- nm[seq_len(min(length(nm), ncol(df0)))]

#  2) parse local time (stockholm) and convert to utc (myclim expects utc) 
df <- df0 %>%
  mutate(
    datetime_local = as.POSIXct(datetime_local, format = "%Y.%m.%d %H:%M", tz = "Europe/Stockholm"),
    datetime = with_tz(datetime_local, "UTC")
  )

#  3) build a long table for myclim (one row per sensor reading) 
# give your site a locality_id (e.g., "b2")
long_tbl <- bind_rows(
  df %>% transmute(locality_id = "F1", sensor_name = "TMS_T1",    datetime, value = as.numeric(T1_C)),
  df %>% transmute(locality_id = "F1", sensor_name = "TMS_T2",    datetime, value = as.numeric(T2_C)),
  df %>% transmute(locality_id = "F1", sensor_name = "TMS_T3",    datetime, value = as.numeric(T3_C)),
  df %>% transmute(locality_id = "F1", sensor_name = "TMS_moist", datetime, value = as.numeric(moisture_raw))
)

# tell myclim which physical units these sensors are in
data_mc <- mc_read_long(
  long_tbl,
  clean  = TRUE,
  silent = TRUE
)


#  5) convert raw moisture to vwc using the peat curve + t1 temp compensation 
data_mc <- mc_calc_vwc(
  data_mc,
  moist_sensor  = "TMS_moist",
  temp_sensor   = "TMS_T1",
  soiltype      = "peat",        # peat-specific curve
  output_sensor = "VWC_peat",
  frozen2NA     = TRUE
)

#  6) back to a data.frame if you want to join/plot 
vwc_df <- mc_reshape_long(data_mc, sensors = c("VWC_peat")) %>%
  select(datetime, locality_id, sensor_name, value) %>%
  rename(VWC_peat = value)
vwc_df <- vwc_df %>%
  mutate(
    VWC_pct = VWC_peat * 100,
    unc_flag = VWC_peat >= 0.995
  )

# work on a copy with a local time column for filtering clarity
vwc_df <- vwc_df %>%
  mutate(datetime_se = with_tz(datetime, "Europe/Stockholm"))

vwc_may_aug_24_25 <- vwc_df %>%
  filter(year(datetime_se) %in% c(2024, 2025),
         month(datetime_se) %in% 5:8)


ggplot(vwc_may_aug_24_25 %>% filter(year(vwc_may_aug_24_25$datetime) == 2025), aes(x = datetime, y = VWC_peat * 100)) +
  geom_point(color = "steelblue", size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)",
    title = "Peat soil volumetric water content (VWC)"
  ) +
  theme_minimal(base_size = 12)


# plot raw moisture vs vwc to see the plateau (needs raw_moist from your myclim object)
raw_moist <- mc_reshape_long(data_mc, sensors = "TMS_moist") %>%
  select(datetime, raw_moist = value)

left_join(raw_moist, vwc_df, by = "datetime") %>%
  ggplot(aes(raw_moist, VWC_peat)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_hline(yintercept = 1.0, linetype = 2) +
  labs(x = "TMS raw moisture", y = "VWC (peat)") +
  theme_minimal()

vwc_df_F1 <- vwc_df
vwc_may_aug_24_25_F1 <- vwc_may_aug_24_25


#  1) read the semicolon-delimited raw file exactly asis 
# your decimals are "." so we keep default locale
raw_path <- "path/conv_data_95140969_2025_10_09_0_F2.csv"

# many tomst "conv_data_*.csv" rows end with a trailing ";" > an empty last column.
# read_delim handles that fine; we then drop empty trailing columns if present.
df0 <- read_delim(raw_path, delim = ";", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# drop columns that are completely empty (e.g., the trailing ";")
df0 <- df0 %>% select(where(~ !all(is.na(.x) | .x == "")))

# expecting 9 useful fields like:
# x1=index, x2=datetime, x3=code(=4), x4=t1, x5=t2, x6=t3, x7=moist_raw, x8=vwc_lolly?, x9=battery?
nm <- c("idx","datetime_local","tz","T1_C","T2_C","T3_C","moisture_raw","VWC_lolly","shake")
names(df0)[seq_len(min(length(nm), ncol(df0)))] <- nm[seq_len(min(length(nm), ncol(df0)))]

#  2) parse local time (stockholm) and convert to utc (myclim expects utc) 
df <- df0 %>%
  mutate(
    datetime_local = as.POSIXct(datetime_local, format = "%Y.%m.%d %H:%M", tz = "Europe/Stockholm"),
    datetime = with_tz(datetime_local, "UTC")
  )

#  3) build a long table for myclim (one row per sensor reading) 
# give your site a locality_id (e.g., "b2")
long_tbl <- bind_rows(
  df %>% transmute(locality_id = "F2", sensor_name = "TMS_T1",    datetime, value = as.numeric(T1_C)),
  df %>% transmute(locality_id = "F2", sensor_name = "TMS_T2",    datetime, value = as.numeric(T2_C)),
  df %>% transmute(locality_id = "F2", sensor_name = "TMS_T3",    datetime, value = as.numeric(T3_C)),
  df %>% transmute(locality_id = "F2", sensor_name = "TMS_moist", datetime, value = as.numeric(moisture_raw))
)

# tell myclim which physical units these sensors are in
data_mc <- mc_read_long(
  long_tbl,
  clean  = TRUE,
  silent = TRUE
)


#  5) convert raw moisture to vwc using the peat curve + t1 temp compensation 
data_mc <- mc_calc_vwc(
  data_mc,
  moist_sensor  = "TMS_moist",
  temp_sensor   = "TMS_T1",
  soiltype      = "peat",        # peat-specific curve
  output_sensor = "VWC_peat",
  frozen2NA     = TRUE
)

#  6) back to a data.frame if you want to join/plot 
vwc_df <- mc_reshape_long(data_mc, sensors = c("VWC_peat")) %>%
  select(datetime, locality_id, sensor_name, value) %>%
  rename(VWC_peat = value)
vwc_df <- vwc_df %>%
  mutate(
    VWC_pct = VWC_peat * 100,
    unc_flag = VWC_peat >= 0.995
  )

# work on a copy with a localtime column for filtering clarity
vwc_df <- vwc_df %>%
  mutate(datetime_se = with_tz(datetime, "Europe/Stockholm"))

vwc_may_aug_24_25 <- vwc_df %>%
  filter(year(datetime_se) %in% c(2024, 2025),
         month(datetime_se) %in% 5:8)


ggplot(vwc_may_aug_24_25 %>% filter(year(vwc_may_aug_24_25$datetime) == 2024) 
       , aes(x = datetime, y = VWC_peat * 100)) +
  geom_point(color = "steelblue", size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)",
    title = "Peat soil volumetric water content (VWC)"
  ) +
  theme_minimal(base_size = 12)


# plot raw moisture vs vwc to see the plateau (needs raw_moist from your myclim object)
raw_moist <- mc_reshape_long(data_mc, sensors = "TMS_moist") %>%
  select(datetime, raw_moist = value)

left_join(raw_moist, vwc_df, by = "datetime") %>%
  ggplot(aes(raw_moist, VWC_peat)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_hline(yintercept = 1.0, linetype = 2) +
  labs(x = "TMS raw moisture", y = "VWC (peat)") +
  theme_minimal()


vwc_df_F2 <- vwc_df
vwc_may_aug_24_25_F2 <- vwc_may_aug_24_25


#  1) read the semicolon-delimited raw file exactly asis 
# your decimals are "." so we keep default locale
raw_path <- "path/conv_data_95140968_2025_10_09_0_F3.csv"

# many tomst "conv_data_*.csv" rows end with a trailing ";" > an empty last column.
# read_delim handles that fine; we then drop empty trailing columns if present.
df0 <- read_delim(raw_path, delim = ";", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# drop columns that are completely empty (e.g., the trailing ";")
df0 <- df0 %>% select(where(~ !all(is.na(.x) | .x == "")))

# expecting 9 useful fields like:
# x1=index, x2=datetime, x3=code(=4), x4=t1, x5=t2, x6=t3, x7=moist_raw, x8=vwc_lolly?, x9=battery?
nm <- c("idx","datetime_local","tz","T1_C","T2_C","T3_C","moisture_raw","VWC_lolly","shake")
names(df0)[seq_len(min(length(nm), ncol(df0)))] <- nm[seq_len(min(length(nm), ncol(df0)))]

#  2) parse local time (stockholm) and convert to utc (myclim expects utc) 
df <- df0 %>%
  mutate(
    datetime_local = as.POSIXct(datetime_local, format = "%Y.%m.%d %H:%M", tz = "Europe/Stockholm"),
    datetime = with_tz(datetime_local, "UTC")
  )

#  3) build a long table for myclim (one row per sensor reading) 
# give your site a locality_id (e.g., "b2")
long_tbl <- bind_rows(
  df %>% transmute(locality_id = "F3", sensor_name = "TMS_T1",    datetime, value = as.numeric(T1_C)),
  df %>% transmute(locality_id = "F3", sensor_name = "TMS_T2",    datetime, value = as.numeric(T2_C)),
  df %>% transmute(locality_id = "F3", sensor_name = "TMS_T3",    datetime, value = as.numeric(T3_C)),
  df %>% transmute(locality_id = "F3", sensor_name = "TMS_moist", datetime, value = as.numeric(moisture_raw))
)

# tell myclim which physical units these sensors are in
data_mc <- mc_read_long(
  long_tbl,
  clean  = TRUE,
  silent = TRUE
)

#  5) convert raw moisture to vwc using the peat curve + t1 temp compensation 
data_mc <- mc_calc_vwc(
  data_mc,
  moist_sensor  = "TMS_moist",
  temp_sensor   = "TMS_T1",
  soiltype      = "peat",        # peat-specific curve
  output_sensor = "VWC_peat",
  frozen2NA     = TRUE
)

#  6) back to a data.frame if you want to join/plot 
vwc_df <- mc_reshape_long(data_mc, sensors = c("VWC_peat")) %>%
  select(datetime, locality_id, sensor_name, value) %>%
  rename(VWC_peat = value)
vwc_df <- vwc_df %>%
  mutate(
    VWC_pct = VWC_peat * 100,
    unc_flag = VWC_peat >= 0.995
  )

# work on a copy with a localtime column for filtering clarity
vwc_df <- vwc_df %>%
  mutate(datetime_se = with_tz(datetime, "Europe/Stockholm"))

vwc_may_aug_24_25 <- vwc_df %>%
  filter(year(datetime_se) %in% c(2024, 2025),
         month(datetime_se) %in% 5:8)


ggplot(vwc_may_aug_24_25 %>% filter(year(vwc_may_aug_24_25$datetime) == 2025) 
       , aes(x = datetime, y = VWC_peat * 100)) +
  geom_point(color = "steelblue", size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)",
    title = "Peat soil volumetric water content (VWC)"
  ) +
  theme_minimal(base_size = 12)


# plot raw moisture vs vwc to see the plateau (needs raw_moist from your myclim object)
raw_moist <- mc_reshape_long(data_mc, sensors = "TMS_moist") %>%
  select(datetime, raw_moist = value)

left_join(raw_moist, vwc_df, by = "datetime") %>%
  ggplot(aes(raw_moist, VWC_peat)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_hline(yintercept = 1.0, linetype = 2) +
  labs(x = "TMS raw moisture", y = "VWC (peat)") +
  theme_minimal()

vwc_df_F3 <- vwc_df
vwc_may_aug_24_25_F3 <- vwc_may_aug_24_25


#### intact


#  1) read the semicolon-delimited raw file exactly asis 
# your decimals are "." so we keep default locale
raw_path <- "path/conv_data_95140967_2025_10_09_0_P1.csv"

# many tomst "conv_data_*.csv" rows end with a trailing ";" > an empty last column.
# read_delim handles that fine; we then drop empty trailing columns if present.
df0 <- read_delim(raw_path, delim = ";", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# drop columns that are completely empty (e.g., the trailing ";")
df0 <- df0 %>% select(where(~ !all(is.na(.x) | .x == "")))

# expecting 9 useful fields like:
# x1=index, x2=datetime, x3=code(=4), x4=t1, x5=t2, x6=t3, x7=moist_raw, x8=vwc_lolly?, x9=battery?
nm <- c("idx","datetime_local","tz","T1_C","T2_C","T3_C","moisture_raw","VWC_lolly","shake")
names(df0)[seq_len(min(length(nm), ncol(df0)))] <- nm[seq_len(min(length(nm), ncol(df0)))]

#  2) parse local time (stockholm) and convert to utc (myclim expects utc) 
df <- df0 %>%
  mutate(
    datetime_local = as.POSIXct(datetime_local, format = "%Y.%m.%d %H:%M", tz = "Europe/Stockholm"),
    datetime = with_tz(datetime_local, "UTC")
  )

#  3) build a long table for myclim (one row per sensor reading) 
# give your site a locality_id (e.g., "b2")
long_tbl <- bind_rows(
  df %>% transmute(locality_id = "P1", sensor_name = "TMS_T1",    datetime, value = as.numeric(T1_C)),
  df %>% transmute(locality_id = "P1", sensor_name = "TMS_T2",    datetime, value = as.numeric(T2_C)),
  df %>% transmute(locality_id = "P1", sensor_name = "TMS_T3",    datetime, value = as.numeric(T3_C)),
  df %>% transmute(locality_id = "P1", sensor_name = "TMS_moist", datetime, value = as.numeric(moisture_raw))
)

# tell myclim which physical units these sensors are in
data_mc <- mc_read_long(
  long_tbl,
  clean  = TRUE,
  silent = TRUE
)

#  5) convert raw moisture to vwc using the peat curve + t1 temp compensation 
data_mc <- mc_calc_vwc(
  data_mc,
  moist_sensor  = "TMS_moist",
  temp_sensor   = "TMS_T1",
  soiltype      = "peat",        # peat-specific curve
  output_sensor = "VWC_peat",
  frozen2NA     = TRUE
)

#  6) back to a data.frame if you want to join/plot 
vwc_df <- mc_reshape_long(data_mc, sensors = c("VWC_peat")) %>%
  select(datetime, locality_id, sensor_name, value) %>%
  rename(VWC_peat = value)
vwc_df <- vwc_df %>%
  mutate(
    VWC_pct = VWC_peat * 100,
    unc_flag = VWC_peat >= 0.995
  )

# work on a copy with a localtime column for filtering clarity
vwc_df<- vwc_df %>%
  mutate(datetime_se = with_tz(datetime, "Europe/Stockholm"))

# 2) may–august for years 2024 and 2025 (local)
vwc_may_aug_24_25 <- vwc_df %>%
  filter(year(datetime_se) %in% c(2024, 2025),
         month(datetime_se) %in% 5:8)


ggplot(vwc_may_aug_24_25 %>% filter(year(vwc_may_aug_24_25$datetime) == 2024) 
       , aes(x = datetime, y = VWC_peat * 100)) +
  geom_point(color = "steelblue", size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)",
    title = "Peat soil volumetric water content (VWC)"
  ) +
  theme_minimal(base_size = 12)


# plot raw moisture vs vwc to see the plateau (needs raw_moist from your myclim object)
raw_moist <- mc_reshape_long(data_mc, sensors = "TMS_moist") %>%
  select(datetime, raw_moist = value)

left_join(raw_moist, vwc_df, by = "datetime") %>%
  ggplot(aes(raw_moist, VWC_peat)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_hline(yintercept = 1.0, linetype = 2) +
  labs(x = "TMS raw moisture", y = "VWC (peat)") +
  theme_minimal()


vwc_df_P1 <- vwc_df
vwc_may_aug_24_25_P1 <- vwc_may_aug_24_25


#  1) read the semicolon-delimited raw file exactly asis 
# your decimals are "." so we keep default locale
raw_path <- "path/conv_data_95140962_2025_10_09_0_P2.csv"

# many tomst "conv_data_*.csv" rows end with a trailing ";" > an empty last column.
# read_delim handles that fine; we then drop empty trailing columns if present.
df0 <- read_delim(raw_path, delim = ";", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# drop columns that are completely empty (e.g., the trailing ";")
df0 <- df0 %>% select(where(~ !all(is.na(.x) | .x == "")))

# expecting 9 useful fields like:
# x1=index, x2=datetime, x3=code(=4), x4=t1, x5=t2, x6=t3, x7=moist_raw, x8=vwc_lolly?, x9=battery?
nm <- c("idx","datetime_local","tz","T1_C","T2_C","T3_C","moisture_raw","VWC_lolly","shake")
names(df0)[seq_len(min(length(nm), ncol(df0)))] <- nm[seq_len(min(length(nm), ncol(df0)))]

#  2) parse local time (stockholm) and convert to utc (myclim expects utc) 
df <- df0 %>%
  mutate(
    datetime_local = as.POSIXct(datetime_local, format = "%Y.%m.%d %H:%M", tz = "Europe/Stockholm"),
    datetime = with_tz(datetime_local, "UTC")
  )

#  3) build a long table for myclim (one row per sensor reading) 
# give your site a locality_id (e.g., "b2")
long_tbl <- bind_rows(
  df %>% transmute(locality_id = "P2", sensor_name = "TMS_T1",    datetime, value = as.numeric(T1_C)),
  df %>% transmute(locality_id = "P2", sensor_name = "TMS_T2",    datetime, value = as.numeric(T2_C)),
  df %>% transmute(locality_id = "P2", sensor_name = "TMS_T3",    datetime, value = as.numeric(T3_C)),
  df %>% transmute(locality_id = "P2", sensor_name = "TMS_moist", datetime, value = as.numeric(moisture_raw))
)

# tell myclim which physical units these sensors are in
data_mc <- mc_read_long(
  long_tbl,
  clean  = TRUE,
  silent = TRUE
)

#  5) convert raw moisture to vwc using the peat curve + t1 temp compensation 
data_mc <- mc_calc_vwc(
  data_mc,
  moist_sensor  = "TMS_moist",
  temp_sensor   = "TMS_T1",
  soiltype      = "peat",        # peat-specific curve
  output_sensor = "VWC_peat",
  frozen2NA     = TRUE
)

#  6) back to a data.frame if you want to join/plot 
vwc_df <- mc_reshape_long(data_mc, sensors = c("VWC_peat")) %>%
  select(datetime, locality_id, sensor_name, value) %>%
  rename(VWC_peat = value)
vwc_df <- vwc_df %>%
  mutate(
    VWC_pct = VWC_peat * 100,
    unc_flag = VWC_peat >= 0.995
  )

# work on a copy with a localtime column for filtering clarity
vwc_df <- vwc_df %>%
  mutate(datetime_se = with_tz(datetime, "Europe/Stockholm"))


# 2) may–august for years 2024 and 2025 (local)
vwc_may_aug_24_25 <- vwc_df %>%
  filter(year(datetime_se) %in% c(2024, 2025),
         month(datetime_se) %in% 5:8)


ggplot(vwc_may_aug_24_25 %>% filter(year(vwc_may_aug_24_25$datetime) == 2025) 
       , aes(x = datetime, y = VWC_peat * 100)) +
  geom_point(color = "steelblue", size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)",
    title = "Peat soil volumetric water content (VWC)"
  ) +
  theme_minimal(base_size = 12)


# plot raw moisture vs vwc to see the plateau (needs raw_moist from your myclim object)
raw_moist <- mc_reshape_long(data_mc, sensors = "TMS_moist") %>%
  select(datetime, raw_moist = value)

left_join(raw_moist, vwc_df, by = "datetime") %>%
  ggplot(aes(raw_moist, VWC_peat)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_hline(yintercept = 1.0, linetype = 2) +
  labs(x = "TMS raw moisture", y = "VWC (peat)") +
  theme_minimal()

vwc_df_P2 <- vwc_df
vwc_may_aug_24_25_P2 <- vwc_may_aug_24_25


#  1) read the semicolon-delimited raw file exactly asis 
# your decimals are "." so we keep default locale
raw_path <- "path/conv_data_95140963_2025_10_09_0_P3.csv"

# many tomst "conv_data_*.csv" rows end with a trailing ";" > an empty last column.
# read_delim handles that fine; we then drop empty trailing columns if present.
df0 <- read_delim(raw_path, delim = ";", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# drop columns that are completely empty (e.g., the trailing ";")
df0 <- df0 %>% select(where(~ !all(is.na(.x) | .x == "")))

# expecting 9 useful fields like:
# x1=index, x2=datetime, x3=code(=4), x4=t1, x5=t2, x6=t3, x7=moist_raw, x8=vwc_lolly?, x9=battery?
nm <- c("idx","datetime_local","tz","T1_C","T2_C","T3_C","moisture_raw","VWC_lolly","shake")
names(df0)[seq_len(min(length(nm), ncol(df0)))] <- nm[seq_len(min(length(nm), ncol(df0)))]

#  2) parse local time (stockholm) and convert to utc (myclim expects utc) 
df <- df0 %>%
  mutate(
    datetime_local = as.POSIXct(datetime_local, format = "%Y.%m.%d %H:%M", tz = "Europe/Stockholm"),
    datetime = with_tz(datetime_local, "UTC")
  )

#  3) build a long table for myclim (one row per sensor reading) 
# give your site a locality_id (e.g., "b2")
long_tbl <- bind_rows(
  df %>% transmute(locality_id = "P3", sensor_name = "TMS_T1",    datetime, value = as.numeric(T1_C)),
  df %>% transmute(locality_id = "P3", sensor_name = "TMS_T2",    datetime, value = as.numeric(T2_C)),
  df %>% transmute(locality_id = "P3", sensor_name = "TMS_T3",    datetime, value = as.numeric(T3_C)),
  df %>% transmute(locality_id = "P3", sensor_name = "TMS_moist", datetime, value = as.numeric(moisture_raw))
)

# tell myclim which physical units these sensors are in
data_mc <- mc_read_long(
  long_tbl,
  clean  = TRUE,
  silent = TRUE
)


#  5) convert raw moisture to vwc using the peat curve + t1 temp compensation 
data_mc <- mc_calc_vwc(
  data_mc,
  moist_sensor  = "TMS_moist",
  temp_sensor   = "TMS_T1",
  soiltype      = "peat",        # peat-specific curve
  output_sensor = "VWC_peat",
  frozen2NA     = TRUE
)

#  6) back to a data.frame if you want to join/plot 
vwc_df <- mc_reshape_long(data_mc, sensors = c("VWC_peat")) %>%
  select(datetime, locality_id, sensor_name, value) %>%
  rename(VWC_peat = value)
vwc_df <- vwc_df %>%
  mutate(
    VWC_pct = VWC_peat * 100,
    unc_flag = VWC_peat >= 0.995
  )

# work on a copy with a localtime column for filtering clarity
vwc_df <- vwc_df %>%
  mutate(datetime_se = with_tz(datetime, "Europe/Stockholm"))

# 2) may–august for years 2024 and 2025 (local)
vwc_may_aug_24_25 <- vwc_df %>%
  filter(year(datetime_se) %in% c(2024, 2025),
         month(datetime_se) %in% 5:8)


ggplot(vwc_may_aug_24_25_P3 %>% filter(year(vwc_may_aug_24_25_P3$datetime) == 2024) 
       , aes(x = datetime, y = VWC_peat * 100)) +
  geom_point(color = "steelblue", size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)",
    title = "Peat soil volumetric water content (VWC)"
  ) +
  theme_minimal(base_size = 12)


# plot raw moisture vs vwc to see the plateau (needs raw_moist from your myclim object)
raw_moist <- mc_reshape_long(data_mc, sensors = "TMS_moist") %>%
  select(datetime, raw_moist = value)

left_join(raw_moist, vwc_df, by = "datetime") %>%
  ggplot(aes(raw_moist, VWC_peat)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_hline(yintercept = 1.0, linetype = 2) +
  labs(x = "TMS raw moisture", y = "VWC (peat)") +
  theme_minimal()


vwc_df_P3 <- vwc_df
vwc_may_aug_24_25_P3 <- vwc_may_aug_24_25


##### combine

vwc_combined <-
  bind_rows(
    vwc_df_P1,
    vwc_df_P2,
    vwc_df_P3,
    vwc_df_B1,
    vwc_df_B2,
    vwc_df_B3,   
    vwc_df_F1,
    vwc_df_F2,
    vwc_df_F3
  ) %>%
  arrange(locality_id, datetime)


vwc_may_aug_24_25_combined <-
  bind_rows(
    vwc_may_aug_24_25_P1,
    vwc_may_aug_24_25_P2,
    vwc_may_aug_24_25_P3,
    vwc_may_aug_24_25_B1,
    vwc_may_aug_24_25_B2,
    vwc_may_aug_24_25_B3,
    vwc_may_aug_24_25_F1,
    vwc_may_aug_24_25_F2,
    vwc_may_aug_24_25_F3
  ) %>%
  arrange(locality_id, datetime)


vwc_may_aug_24_25_combined <- vwc_may_aug_24_25_combined %>%
  mutate(
    thaw_stage = case_when(
      locality_id == "P1" ~ "intact",
      locality_id == "P2" ~ "partly_thawed",
      locality_id == "P3" ~ "intact",
      locality_id == "B1" ~ "partly_thawed",
      locality_id == "B2" ~ "intact",
      locality_id == "B3" ~ "partly_thawed",
      locality_id == "F1" ~ "fully_thawed",
      locality_id == "F2" ~ "fully_thawed",
      locality_id == "F3" ~ "fully_thawed",
      TRUE ~ NA_character_
    ))

ggplot(vwc_may_aug_24_25_combined,
       aes(x = datetime, y = VWC_pct, color = as.factor(thaw_stage))) +
  geom_point(size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)"
  ) +
  theme_bw()

vwc_combined <- vwc_combined %>%
  mutate(
    thaw_stage = case_when(
      locality_id == "P1" ~ "intact",
      locality_id == "P2" ~ "partly_thawed",
      locality_id == "P3" ~ "intact",
      locality_id == "B1" ~ "partly_thawed",
      locality_id == "B2" ~ "intact",
      locality_id == "B3" ~ "partly_thawed",
      locality_id == "F1" ~ "fully_thawed",
      locality_id == "F2" ~ "fully_thawed",
      locality_id == "F3" ~ "fully_thawed",
      TRUE ~ NA_character_
    ))

### cut the data from 2025 end because seems like data includes 
### period when the sensors were taken out of soil

# subset to include only data up to and including 20250812 23:45:00 utc
# (based on raw data (lolly vwc = 0.0) and visual inspection)

vwc_trimmed <- vwc_combined %>%
  filter(
    datetime >= ymd_hms("2023-09-30 00:00:00", tz = "UTC"),
    datetime <= ymd_hms("2025-08-12 23:45:00", tz = "UTC")
  )

vwc_may_aug_24_25_combined <- vwc_may_aug_24_25_combined %>%
  filter(
    datetime >= ymd_hms("2023-09-30 00:00:00", tz = "UTC"),
    datetime <= ymd_hms("2025-08-12 23:45:00", tz = "UTC")
  )


ggplot(vwc_trimmed,
       aes(x = datetime, y = VWC_pct, color = as.factor(thaw_stage))) +
  geom_point(size = 1, alpha=0.7) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)"
  ) +
  theme_bw()

# flag vwc < 20% points

vwc_trimmed <- vwc_trimmed %>%
  mutate(
    unc_flag_low = ifelse(!is.na(VWC_peat) & VWC_peat < 0.21, TRUE, FALSE)
  )

vwc_may_aug_24_25_combined <- vwc_may_aug_24_25_combined %>%
  mutate(
    unc_flag_low = ifelse(!is.na(VWC_peat) & VWC_peat < 0.21, TRUE, FALSE)
  )


vwc_trimmed <- vwc_trimmed %>%
  mutate(
    unc_flag_high2 = coalesce(unc_flag_high, unc_flag)
  )

vwc_trimmed <- subset(vwc_trimmed, select = -c(unc_flag, unc_flag_high))

colnames(vwc_trimmed)[9] <- "unc_flag_high"


## may aug is missing it

vwc_may_aug_24_25_combined <- vwc_may_aug_24_25_combined %>%
  mutate(
    VWC_pct = VWC_peat * 100,
    unc_flag_high = VWC_peat >= 0.995
  )

# testing filtering:

ggplot(vwc_may_aug_24_25_combined ,
       aes(x = datetime, y = VWC_pct, color = as.factor(chamber_id))) +
  geom_point(size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)"
  ) +
  theme_bw()


ggplot(
  vwc_may_aug_24_25_combined %>%
    filter(chamber_id == 2, unc_flag_high != TRUE),
  aes(x = as_datetime(datetime), y = VWC_pct)
) +
  geom_point(size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)"
  ) +
  theme_bw()

ggplot(
  vwc_may_aug_24_25_combined %>%
    filter(chamber_id == 2, unc_flag_high != TRUE),
  aes(x = as_datetime(datetime), y = VWC_pct)
) +
  geom_point(size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)"
  ) +
  theme_bw()

ggplot(daily_vwc
       %>% filter(daily_vwc$chamber_id == 2),
       aes(x = date, y = VWC_pct_median)) +
  geom_point(size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)"
  ) +
  theme_bw()

ggplot(vwc_trimmed
       %>% filter(vwc_trimmed$unc_flag_low == FALSE),
       aes(x = datetime, y = VWC_pct, color = as.factor(thaw_stage))) +
  geom_point(size = 1) +
  labs(
    x = "Datetime (UTC)",
    y = "Volumetric water content (%)"
  ) +
  theme_bw()


vwc_trimmed <- vwc_trimmed %>%
  arrange(chamber_id)

vwc_trimmed <- vwc_trimmed %>%
  mutate(
    unc_flag_low = ifelse(
      is.na(VWC_pct) & is.na(VWC_peat),
      NA,                  # set to NA when both moisture values are NA
      unc_flag_low         # otherwise keep the existing flag
    )
  )


## add sensor ids

vwc_trimmed <- vwc_trimmed %>%
  mutate(
    sensor_id = case_when(
      locality_id == "P1" ~ "TMS_95140967",
      locality_id == "P2" ~ "TMS_95140962",
      locality_id == "P3" ~ "TMS_95140963",
      locality_id == "B1" ~ "TMS_95140964",
      locality_id == "B2" ~ "TMS_95140961",
      locality_id == "B3" ~ "TMS_95140965",
      locality_id == "F1" ~ "TMS_95140970",
      locality_id == "F2" ~ "TMS_95140969",
      locality_id == "F3" ~ "TMS_95140968"
    )
  )
vwc_trimmed <- vwc_trimmed %>%
  relocate(sensor_id, .after = locality_id)


# save

write.csv(vwc_trimmed, "path/vwc_combined.csv", row.names = FALSE)

write.csv(vwc_may_aug_24_25_combined, "path/vwc_may_aug_24_25_combined.csv", row.names = FALSE)

vwc_trimmed <- read.csv("path/vwc_combined.csv")

vwc_may_aug_24_25_combined <- read.csv("path/vwc_may_aug_24_25_combined.csv")

### change the plot names

vwc_may_aug_24_25_combined <- vwc_may_aug_24_25_combined %>%
  mutate(
    chamber_id = case_when(
      locality_id == "P1" ~ "3",
      locality_id == "P2" ~ "5",
      locality_id == "P3" ~ "1",
      locality_id == "B1" ~ "4",
      locality_id == "B2" ~ "2",
      locality_id == "B3" ~ "6",
      locality_id == "F1" ~ "8",
      locality_id == "F2" ~ "7",
      locality_id == "F3" ~ "9"
    )
  )

vwc_trimmed <- vwc_trimmed %>%
  mutate(
    chamber_id = case_when(
      locality_id == "P1" ~ "3",
      locality_id == "P2" ~ "5",
      locality_id == "P3" ~ "1",
      locality_id == "B1" ~ "4",
      locality_id == "B2" ~ "2",
      locality_id == "B3" ~ "6",
      locality_id == "F1" ~ "8",
      locality_id == "F2" ~ "7",
      locality_id == "F3" ~ "9"
    )
  )


vwc_trimmed <- vwc_trimmed %>%
  relocate(chamber_id, .after = locality_id)

vwc_may_aug_24_25_combined <- vwc_may_aug_24_25_combined %>%
  relocate(chamber_id, .after = locality_id)

vwc_trimmed <- vwc_trimmed %>%
  relocate(thaw_stage, .after = chamber_id)

vwc_may_aug_24_25_combined <- vwc_may_aug_24_25_combined %>%
  relocate(thaw_stage, .after = chamber_id)


# remove some columns

vwc_trimmed <- subset(vwc_trimmed, select = -c(sensor_name))

vwc_trimmed <- vwc_trimmed %>%
  relocate(datetime_se, .after = datetime)

vwc_may_aug_24_25_combined <- subset(vwc_may_aug_24_25_combined, select = -c(sensor_name, unc_flag))

vwc_may_aug_24_25_combined <- vwc_may_aug_24_25_combined %>%
  relocate(datetime_se, .after = datetime)

### make a daily mean/median aggregation

# set rows with unc_flag_low to NA

# and for chamber_id = 2, set rows to na when unc_flag_high == true
# because those points are not realistic 

vwc_filt <- vwc_may_aug_24_25_combined %>%
  mutate(
    VWC_peat = ifelse(
      unc_flag_low == TRUE | (chamber_id == 2 & unc_flag_high == TRUE),
      NA_real_,
      VWC_peat
    ),
    VWC_pct = ifelse(
      unc_flag_low == TRUE | (chamber_id == 2 & unc_flag_high == TRUE),
      NA_real_,
      VWC_pct
    )
  )


daily_vwc <- vwc_filt %>%
  mutate(date = as.Date(datetime_se)) %>%
  group_by(chamber_id, thaw_stage, date) %>%
  summarise(
    across(
      c(VWC_peat, VWC_pct),
      .fns = list(mean   = ~mean(.x, na.rm = TRUE),
                  median = ~median(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    unc_flag_low  = { any_true <- any(unc_flag_low  %in% TRUE, na.rm = TRUE)
    all_na   <- all(is.na(unc_flag_low))
    if (any_true) TRUE else if (all_na) NA else FALSE },
    unc_flag_high = { any_true <- any(unc_flag_high %in% TRUE, na.rm = TRUE)
    all_na   <- all(is.na(unc_flag_high))
    if (any_true) TRUE else if (all_na) NA else FALSE },
    .groups = "drop"
  ) %>%
  mutate(across(matches("VWC_(peat|pct)_(mean|median)$"), ~na_if(.x, NaN))) %>%
  arrange(chamber_id, date)

daily_vwc <- as.data.frame(daily_vwc)


write.csv(daily_vwc, "path/vwc_may_aug_24_25_daily.csv", row.names = FALSE)
daily_vwc <- read.csv("path/vwc_may_aug_24_25_daily.csv")


# quick stats
daily_vwc %>%
  group_by(chamber_id) %>%
  summarise(
    min_VWC_pct_median   = min(VWC_pct_median, na.rm = TRUE),
    max_VWC_pct_median   = max(VWC_pct_median, na.rm = TRUE),
    mean_VWC_pct_median  = mean(VWC_pct_median, na.rm = TRUE),
    median_VWC_pct_median = median(VWC_pct_median, na.rm = TRUE),
    n_obs = sum(!is.na(VWC_pct_median))
  ) %>%
  arrange(chamber_id)


# dataset for si (table c3): 

# remove unnecessary columns: 

daily_vwc <- subset(daily_vwc, select = -c(VWC_peat_mean, VWC_peat_median, VWC_pct_mean, unc_flag_low, unc_flag_high))

# remove rows with only na
daily_vwc <- daily_vwc %>%
  filter(!is.na(VWC_pct_median))

daily_vwc <- daily_vwc %>%
  select(date, everything())

write.csv(daily_vwc, "path/vwc_may_aug_24_25_daily_SI.csv", row.names = FALSE)


## make a chamber_id soil moisture df to be combined with the ch4 critical values

vwc_chamber_median <- daily_vwc %>%
  group_by(chamber_id) %>%
  summarise(
    VWC_pct_2425 = median(VWC_pct_median, na.rm = TRUE),
    .groups = "drop" )

# combine with plant_ch4

vwc_chamber_median <- vwc_chamber_median %>%
  mutate(chamber_id = as.character(chamber_id))

plant_ch4 <- plant_ch4 %>%
  mutate(chamber_id = as.character(chamber_id))

# join
plant_ch4 <- plant_ch4 %>%
  left_join(vwc_chamber_median, by = "chamber_id")

########## add season mean soil temps

season_temp <- flux_daily %>%
  group_by(chamber) %>%
  summarise(
    avgTgnd_season = mean(avgTgnd_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(chamber_id = chamber)

season_temp <- season_temp %>%
  mutate(chamber_id = as.character(chamber_id))

plant_ch4 <- plant_ch4 %>%
  left_join(season_temp, by = "chamber_id")


write.csv(plant_ch4, "path/plant_ch4_neefiltered_ch_critical_points_QAQCd.csv", row.names = FALSE)

#####################################
###--------ZENODO DATASETS--------###
#####################################

### roots and rhizomes

## combine plot and depth scales

plot <- subset(plant_ch4, select = c(plot_id, chamber_id, thaw_stage, PFT, organ, mass_gm2_std, sa_m2m2_std, 
                                     len_mm2_std, avg_d_mm, 
                                     RTD_gcm3, SRL_mg1, GWC))

plot$len_kmm2_std <- plot$len_mm2_std/1000

root_vga_depth <- as.data.frame(root_vga_depth)
root_vga_depth$sa_m2m2_std <- root_vga_depth$sa_cm2m2_std / 1e4 

depth <- subset(root_vga_depth, select = c(plot_id, chamber_id, thaw_stage, depth_std,
                                           PFT, organ, mass_gm2_std, sa_m2m2_std, 
                                           len_kmm2_std, avg_d_mm, RTD_gcm3, SRL_mg1,
                                           GWC))

depth <- depth %>% drop_na(plot_id)

## combine


traits <- c("mass_gm2_std", "sa_m2m2_std", "len_kmm2_std",
            "avg_d_mm", "RTD_gcm3", "SRL_mg1")

#  plot rows 
plot_scaled <- plot %>%
  select(plot_id, chamber_id, thaw_stage, PFT, organ, GWC, all_of(traits)) %>%
  mutate(scale = "plot")

#  depth rows (with scale from depth_std) 
depth_scaled <- depth %>%
  select(plot_id, chamber_id, thaw_stage, PFT, organ, depth_std, GWC, all_of(traits)) %>%
  mutate(scale = case_when(
    depth_std == "0-10"  ~ "depth_0_10",
    depth_std == "10-20" ~ "depth_10_20",
    depth_std == "20-30" ~ "depth_20_30",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(scale)) %>%  
  select(-depth_std)       

#  combine into final df 
root_df <- bind_rows(plot_scaled, depth_scaled)

root_df <- root_df %>%
  mutate(
    thaw_stage = factor(thaw_stage, levels = c("intact", "partly_thawed", "fully_thawed"))
  ) %>%
  arrange(thaw_stage, chamber_id)

## change thaw stage names

root_df <- root_df %>%
  mutate(thaw_stage = case_when(
    str_detect(plot_id, "palsa") ~ "intact",
    str_detect(plot_id, "bog") ~ "partly_thawed",
    str_detect(plot_id, "fen") ~ "fully_thawed"
  ))

root_df <- subset(root_df, select = -plot_id)

colnames(root_df)[1] <- "plot_id"

# rearrange

root_df <- root_df %>%
  relocate(GWC, .after = last_col())

root_df <- root_df %>%
  relocate(scale, .after = organ)

colnames(root_df)[9] <- "d_mm"

root_df <- root_df %>%
  rename_with(tolower)

## add date of collection
root_df$date <- as_date("2023-07-21")

root_df <- root_df %>%
  relocate(date, .before = plot_id)

# add coordinates


plot_coords <- read.csv("path/24072023plots_stordalen.csv")

plot_coords <- plot_coords %>%
  # rename first 3 columns
  rename(
    plot_ID = 1,
    lat = 2,
    lon = 3
  ) %>%
  # filter out rows where 'plot' contains unwanted terms
  filter(!grepl("EC|test|tesu", plot_ID, ignore.case = TRUE)) %>%
  dplyr::select(plot_ID, lat, lon)


# apply transformation
plot_coords <- plot_coords %>%
  mutate(
    plot_ID = if_else(
      str_detect(plot_ID, "ruth"),
      {
        number <- str_extract(plot_ID, "\\d+")
        type <- class_map[number]
        str_c(type, number, "chamber")
      },
      plot_ID  
    )
  )

plot_coords <- plot_coords %>%
  mutate(
    lat = str_remove(lat, "N$"),     
    lon = str_remove(lon, "E$"),
    lat = as.numeric(lat),        
    lon = as.numeric(lon)        
  )

plot_coords <- plot_coords %>%
  mutate(thaw_stage = case_when(
    str_detect(plot_ID, "palsa") ~ "palsa",
    str_detect(plot_ID, "bog") ~ "bog",
    str_detect(plot_ID, "fen") ~ "fen"
  ))

plot_coords_new <- plot_coords %>%
  filter(plot_type != "chamber") %>%                   
  transmute(
    plot_id = str_to_lower(str_extract(plot_ID, "\\d+")),
    lat,
    lon,
    thaw_stage = recode(thaw_stage,                      
                        palsa = "intact",
                        bog   = "partly_thawed",
                        fen   = "fully_thawed")
  )                                                   

plot_coords_new

root_df$plot_id <- as.character(root_df$plot_id)

root_df <- root_df %>%
  left_join(
    plot_coords_new %>% select(plot_id, lat, lon),
    by = "plot_id"
  )

root_df <- root_df %>%
  relocate(lat, lon, .after = plot_id)

# pct values per plot and depth

pct_incr <- pct_incr %>%
  mutate(depth_std = case_when(
    sample_increment == "10-26" ~ "10-20",
    sample_increment == "20-32" ~ "20-30",
    sample_increment == "10-25" ~ "10-20",
    sample_increment == "10-18" ~ "10-20",
    sample_increment == "10-21" ~ "10-20",
    sample_increment == "10-24" ~ "10-20",
    sample_increment == "20-35.5" ~ "20-30",
    TRUE ~ sample_increment
  ))

pct_incr <- subset(pct_incr, select = -c(PHYS_code, sample_increment, total_length, total_global))

pct <- subset(pct, select = -c(total_length, total_global))

# helper: map depth_std > root_df$scale labels
depth_to_scale <- function(x) {
  dplyr::case_when(
    x == "0-10"  ~ "depth_0_10",
    x == "10-20" ~ "depth_10_20",
    x == "20-30" ~ "depth_20_30",
    TRUE ~ NA_character_
  )
}

# clean + reshape pct (plotscale) to wide (fine/coarse)
pct_wide <- pct %>%
  mutate(
    plot_id = str_extract(plot_id, "\\d+"), 
    pft     = str_to_lower(PFT),
    scale   = "plot",
    size_class = case_when(
      size_class == "<2mm"  ~ "pct_fine_root",
      size_class == ">=2mm" ~ "pct_coarse_root",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(size_class)) %>%
  select(plot_id, pft, scale, size_class, pct_length) %>%
  pivot_wider(names_from = size_class, values_from = pct_length)

# clean + reshape pct_incr (depthscale) to wide (fine/coarse)
pct_incr_wide <- pct_incr %>%
  mutate(
    plot_id = str_extract(plot_id, "\\d+"), 
    pft     = str_to_lower(PFT),
    scale   = depth_to_scale(depth_std),
    size_class = case_when(
      size_class == "<2mm"  ~ "pct_fine_root",
      size_class == ">=2mm" ~ "pct_coarse_root",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(scale), !is.na(size_class)) %>%
  select(plot_id, pft, scale, size_class, pct_length) %>%
  pivot_wider(names_from = size_class, values_from = pct_length)

# bind pct sources + join into root_df (keep all root_df rows)
pct_all <- bind_rows(pct_wide, pct_incr_wide)

root_df <- root_df %>%
  mutate(
    plot_id = as.character(plot_id),
    pft     = tolower(pft)
  ) %>%
  left_join(pct_all, by = c("plot_id", "pft", "scale")) %>%
  mutate(
    pct_fine_root   = if_else(organ == "root", pct_fine_root, NA_real_),
    pct_coarse_root = if_else(organ == "root", pct_coarse_root, NA_real_)
  ) %>%
  relocate(pct_fine_root, pct_coarse_root, .before = gwc)

## add ph and soil temp and active layer

soilvars <- read.csv("path/soilvars_data.csv")

soilvars <- soilvars %>%
  rename(
    plot_id       = plot_ID,
    chamber_id   = chamber_ID,
    soil_ph    = pH
  )

soilvars <- subset(soilvars, select = -c(WTL, notes, chamber_id))

soilvars <- soilvars %>%
  filter(!str_detect(plot_id, "chamber")) %>%      
  mutate(plot_id = str_extract(plot_id, "\\d+")) 

root_df <- root_df %>%
  left_join(
    soilvars %>% select(plot_id, soil_ph, act_layer),
    by = "plot_id"
  ) %>%
  relocate(soil_ph, act_layer, .after = gwc)

root_df <- root_df %>%
  mutate(
    act_layer = case_when(
      thaw_stage == "partly_thawed" ~ ">85",    
      plot_id == 5 & act_layer == 85.5 ~ "85", 
      TRUE ~ as.character(act_layer)           
    )
  )

root_df <- root_df %>%
  mutate(
    soil_ph   = if_else(scale == "plot", soil_ph, NA_real_),
    act_layer = if_else(scale == "plot", act_layer, NA_character_)
  )


root_df <- root_df %>%
  rename(
    mass_gm2       = mass_gm2_std,
    sa_m2m2   = sa_m2m2_std,
    len_kmm2    = len_kmm2_std,
    td_gcm3 = rtd_gcm3,
    ph = soil_ph
  )

write.csv(root_df, "path/belowground_traits.csv", row.names = F)

##### veg survey data

veg <- read.csv("path/comparison_data.csv")

veg <- subset(veg, select = c(plot_ID, thaw_stage, species, PFT_class, species_GA_m2m2, aer_GA_m2m2, nonaer_GA_m2m2, herb_GA_m2m2, 
                              shrub_GA_m2m2, mean_TS, sd_TS, class, species_cover, aer_cover, nonaer_cover, 
                              herb_cover, shrub_cover,
                              moss_cover, vasc_cover, soil_pH, active_layer_cm, dead_cover, litter_cover, water_cover))


veg <- veg %>%
  rename(
    plot_id       = plot_ID,
    pft   = PFT_class,
    ts_mean    = mean_TS,
    ts_sd = sd_TS
  )

veg <- veg %>%
  rename_with(tolower)

veg <- veg %>%
  relocate(plot_type, .after = plot_id)

veg <- veg %>% mutate(plot_id = str_extract(plot_id, "\\d+")) 

plot_veg <- subset(veg, select = c(plot_id, plot_type, thaw_stage, aer_ga_m2m2, nonaer_ga_m2m2, 
                                   herb_ga_m2m2, shrub_ga_m2m2, ts_mean, ts_sd,
                                   soil_ph, active_layer_cm, aer_cover, nonaer_cover, 
                                   herb_cover, shrub_cover, moss_cover, 
                                   vasc_cover, dead_cover, litter_cover, water_cover))

veg_species <- subset(veg, select = c(plot_id, plot_type, thaw_stage, species, pft, species_ga_m2m2, species_cover, soil_ph, active_layer_cm))

veg_species <- veg_species %>%
  mutate(
    active_layer_cm = case_when(
      thaw_stage == "partly_thawed" ~ ">85",
      plot_id == 5 & active_layer_cm == 85.5 ~ "85",   
      TRUE ~ as.character(active_layer_cm)  
    )
  )


plot_veg <- plot_veg %>%
  mutate(
    active_layer_cm = case_when(
      thaw_stage == "partly_thawed" ~ ">85", 
      plot_id == 5 & active_layer_cm == 85.5 ~ "85",   
      TRUE ~ as.character(active_layer_cm)  
    )
  )

plot_veg <- plot_veg %>%
  group_by(plot_id, plot_type) %>%
  slice(1) %>%
  ungroup()

plot_veg <- as.data.frame(plot_veg)

colnames(plot_veg)[11] <- "act_layer"

plot_veg <- plot_veg %>%
  mutate(thaw_stage = case_when(
    thaw_stage == "intact" ~ "intact",
    thaw_stage == "partly_thawed" ~ "partly_thawed",
    thaw_stage == "fully_thawed" ~ "fully_thawed"
  ))

# arrange

plot_veg <- plot_veg %>%
  mutate(
    thaw_stage = factor(thaw_stage, levels = c("intact", "partly_thawed", "fully_thawed"))
  ) %>%
  arrange(thaw_stage, plot_id)

# add date and lat lon

date_lookup <- tibble::tribble(
  ~plot_ID,    ~date,
  "fen7chamber",  "19-07-2023",
  "fen8chamber",  "19-07-2023",
  "fen9chamber",  "19-07-2023",
  "bog4chamber",  "20-07-2023",
  "bog6chamber",  "20-07-2023",
  "bog2chamber",  "20-07-2023",
  "palsa1chamber","20-07-2023",
  "palsa5chamber","21-07-2023",
  "palsa3chamber","22-07-2023",
  "fen9core",  "22-07-2023",
  "bog4core",  "22-07-2023",
  "palsa1core","22-07-2023",
  "palsa3core","22-07-2023",
  "bog2core",  "22-07-2023",
  "bog6core",  "22-07-2023",
  "fen8core",  "23-07-2023",
  "fen7core",  "23-07-2023",
  "palsa5core","23-07-2023"
) %>%
  mutate(
    date = dmy(date),
    plot_id = str_extract(plot_ID, "\\d+"),
    plot_type = case_when(
      str_detect(plot_ID, "chamber") ~ "chamber",
      str_detect(plot_ID, "core") ~ "core"
    )
  ) %>%
  select(plot_id, plot_type, date)

plot_veg <- plot_veg %>%
  mutate(plot_id = as.character(plot_id)) %>%
  left_join(date_lookup, by = c("plot_id", "plot_type")) %>%
  relocate(date, .before = everything())

plot_veg <- plot_veg %>%
  left_join(
    plot_coords_2 %>% select(plot_id, plot_type, lat, lon),
    by = c("plot_id", "plot_type")
  )

plot_veg <- subset(plot_veg, select = -c(lat.x, lon.x))

plot_veg <- plot_veg %>% relocate(lat.y, lon.y, .after = plot_type)

plot_veg <- plot_veg %>% relocate(ts_mean, ts_sd, soil_ph, act_layer, .after = water_cover)

plot_veg <- plot_veg %>%
  rename(
    ph = soil_ph,
    lat = lat.y,
    lon = lon.y
  )

write.csv(plot_veg, "path/plot_veg_survey.csv", row.names = F)


### veg species

colnames(veg_species)[9] <- "act_layer"

veg_species <- veg_species %>%
  mutate(
    thaw_stage = factor(thaw_stage, levels = c("intact", "partly_thawed", "fully_thawed"))
  ) %>%
  arrange(thaw_stage, plot_id)

veg_species <- veg_species %>%
  mutate(plot_id = as.character(plot_id)) %>%
  left_join(date_lookup, by = c("plot_id", "plot_type")) %>%
  relocate(date, .before = everything())

veg_species <- subset(veg_species, select = -c(lat, lon))

veg_species <- veg_species %>%
  left_join(
    plot_coords_2 %>% select(plot_id, plot_type, lat, lon),
    by = c("plot_id", "plot_type")
  )

veg_species <- veg_species %>% relocate(lat, lon, .after = plot_type)

veg_species <- veg_species %>%
  rename(
    ph      = soil_ph
  )

write.csv(veg_species, "path/species_veg_survey.csv", row.names = F)


##### flux data for si (table c4)

flux_daily <- as.data.frame(flux_daily)

flux_daily <- flux_daily %>%
  relocate(date, .before = chamber)

flux_daily_SI <- flux_daily %>%
  rename(
    ch4_flux_mgch4m2d1_median      = daily_median_CH4,
    ch4_flux_mgch4m2d1_mean      = daily_mean_CH4,
    peat_temp   = avgTgnd_mean,
    chamber_id    = chamber
  )

write.csv(flux_daily_SI, "path/flux_data_SI.csv", row.names = F)


