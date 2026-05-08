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
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(vegan)
library(reshape2)
library(broom)
library(scales)
library(nlme)
library(patchwork)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(terra)
library(stringr)
library(ggnewscale)
library(data.table)
library(purrr)
library(grid)

#########################################

# import data sets

### root data: 

# plot-scale
# depth-scale
### root-vga data sets

root_vga_plot <- read.csv("path/root_vga_plot_scale.csv")

root_vga_depth <- read.csv("path/root_vga_depth_scale.csv")

# for more direct VGA - root connections
# chamber-core comparison data

comparison_data <- read.csv("path/comparison_data.csv")

# root trait + ch4 data

# no PFTs
plant_ch4_noPFT <- read.csv("path/root_rhizome_noPFT_chamber_level_QAQCd.csv")

# with PFTs
plant_ch4 <- read.csv("path/plant_ch4_neefiltered_ch_critical_points_QAQCd.csv")

# porewater ch4
########################################

### Descriptive statistics and data exploration ###

### ROOTS AND RHIZOMES ###

## per depth ##

# make a horizontal bar graph where y = 0-10, 10-20, 20-30 cm and x axis is mass_density

# create a new column where depths are standardized

root_vga_depth <- root_vga_depth %>%
  mutate(
    depth_std = case_when(
      sample_increment == "0-10" ~ "0-10",
      grepl("^10-", sample_increment) ~ "10-20",
      grepl("^20-", sample_increment) ~ "20-30",
      TRUE ~ NA_character_
    )
  )

root_vga_depth <- root_vga_depth %>%
  relocate(depth_std, .after = sample_increment)

# convert len_density_cmcm3 to km/m3

root_vga_depth$len_density_kmm3 <- root_vga_depth$len_density_mm3 / 1000

root_vga_depth$depth_std <- factor(root_vga_depth$depth_std, levels = c("0-10", "10-20", "20-30"))
root_vga_depth$thaw_stage <- factor(root_vga_depth$thaw_stage, levels = c("intact", "partly_thawed", "fully_thawed"))

root_vga_depth <- root_vga_depth %>%
  mutate(
    thaw_stage = factor(thaw_stage, levels = c("intact", "partly_thawed", "fully_thawed")),
    depth_std = factor(depth_std, levels = c("0-10", "10-20", "20-30"))
  ) %>%
  arrange(depth_std, thaw_stage) 


#### plot

pft_colors   <- c(herbaceous = "#213448", shrub = "#94B4C1")
pt_fill_herb <- "#ECEFCA"
pt_fill_shrb <- "#EBD3F8"

depth_lv <- c("0-10","10-20","20-30")
thaw_lv  <- c("intact","partly_thawed","fully_thawed")

# top to bottom order (depth first, then thaw)
levels_top_to_bottom <- expand.grid(depth = depth_lv, thaw = thaw_lv, 
                                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) |>
  arrange(match(depth, depth_lv), match(thaw, thaw_lv)) |>
  transmute(lbl = paste(depth, thaw)) |>
  pull(lbl)

# modify the df for plotting
root_vga_depth_complete <- root_vga_depth %>%
  filter(organ == "root") %>%
  complete(depth_std, thaw_stage, PFT, organ, fill = list(mass_gm2_std = 0)) %>%
  filter(organ == "root") %>%
  mutate(
    depth_std  = factor(depth_std, levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    # Make combined y factor and reverse levels so 0–10 rows are at the top
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom))
  )

# plot
mass_depth_root <- ggplot(root_vga_depth_complete,
                          aes(x = mass_gm2_std, y = group_y, fill = PFT)) +
  geom_col(position = "stack") +
  geom_point(
    data = ~ dplyr::filter(.x, PFT == "herbaceous"),
    position = position_stack(vjust = 0.5),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_herb, color = "darkgrey"
  ) +
  geom_point(
    data = ~ dplyr::filter(.x, PFT == "shrub"),
    position = position_stack(vjust = 0.5),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_shrb, color = "darkgrey"
  ) +
  scale_fill_manual(values = pft_colors, name = "PFT") +
  labs(
    x = expression(paste("Root biomass (g m"^-2, ")")),
    y = "Depth (cm)"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)), 
    limits = c(0, NA)
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text  = element_text(size = 18)
  )

mass_depth_root


# rhizome biomass

root_vga_depth_complete <- root_vga_depth %>%
  filter(organ == "rhizome") %>%
  complete(depth_std, thaw_stage, PFT, organ, fill = list(mass_gm2_std = 0)) %>%
  filter(organ == "rhizome") %>%
  mutate(
    depth_std  = factor(depth_std, levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    # Make combined y factor and reverse levels so 0–10 rows are at the top
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom))
  )


# plot
mass_depth_rhiz <- ggplot(root_vga_depth_complete,
                          aes(x = mass_gm2_std, y = group_y, fill = PFT)) +
  geom_col(position = "stack") +
  geom_point(
    data = ~ dplyr::filter(.x, PFT == "herbaceous"),
    position = position_stack(vjust = 0.5),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_herb, color = "#343434cc"
  ) +
  geom_point(
    data = ~ dplyr::filter(.x, PFT == "shrub"),
    position = position_stack(vjust = 0.5),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_shrb, color = "#343434cc"
  ) +
  scale_fill_manual(values = pft_colors, name = "PFT") +
  labs(
    x = expression(paste("Rhizome biomass (g m"^-2, ")")),
    y = "Depth (cm)"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, NA)
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text  = element_text(size = 18)
  )

mass_depth_rhiz


# SRL

pft_colors <- c(herbaceous = "#213448", shrub = "#94B4C1")

depth_lv <- c("0-10","10-20","20-30")
thaw_lv  <- c("intact","partly_thawed","fully_thawed")
group_lv <- as.vector(outer(depth_lv, thaw_lv, paste))

# with weighted means

# summarize per PFT
by_pft <- root_vga_depth %>%
  filter(organ == "root") %>%
  group_by(depth_std, thaw_stage, PFT) %>%
  summarise(
    biomass_sum = sum(mass_gm2_std, na.rm = TRUE),
    SRL_mean    = mean(SRL_mg1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv)
  )

# Compute weighted mean SRL per depth × thaw stage using weighted.mean()
weighted_SRL <- by_pft %>%
  group_by(depth_std, thaw_stage) %>%
  summarise(
    SRL_weighted = weighted.mean(SRL_mean, biomass_sum, na.rm = TRUE),
    .groups = "drop"
  )

# Compute relative contribution of each PFT to the weighted mean
by_pft <- by_pft %>%
  left_join(weighted_SRL, by = c("depth_std", "thaw_stage")) %>%
  group_by(depth_std, thaw_stage) %>%
  mutate(
    biomass_share = biomass_sum / sum(biomass_sum, na.rm = TRUE),
    SRL_contrib = biomass_share * SRL_mean,
    group_y = factor(paste(depth_std, thaw_stage), levels = group_lv)
  ) %>%
  ungroup()



# top to bottom order (depth first, then thaw), then reverse so 0–10 rows appear on top
levels_top_to_bottom <- expand.grid(depth = depth_lv, thaw = thaw_lv,
                                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) |>
  dplyr::arrange(match(depth, depth_lv), match(thaw, thaw_lv)) |>
  dplyr::transmute(lbl = paste(depth, thaw)) |>
  dplyr::pull(lbl)

# make sure factors + group_y use that order
by_pft <- by_pft %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom)) 
  )


# Add 0s for missing combinations (so 20–30 fully_thawed is included)
by_pft_complete <- by_pft %>%
  tidyr::complete(
    depth_std,
    thaw_stage,
    PFT,
    fill = list(SRL_contrib = 0)
  ) %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom))
  )

# Then plot the completed dataset
srl_depth_root <- ggplot(by_pft_complete, aes(x = SRL_contrib, y = group_y, fill = PFT)) +
  geom_col(position = "stack") +
  geom_point(
    data = root_vga_depth %>% filter(organ == "root", PFT == "herbaceous"),
    aes(x = SRL_mg1, y = factor(paste(depth_std, thaw_stage),
                                levels = rev(levels_top_to_bottom))),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_herb, color = "darkgrey"
  ) +
  geom_point(
    data = root_vga_depth %>% filter(organ == "root", PFT == "shrub"),
    aes(x = SRL_mg1, y = factor(paste(depth_std, thaw_stage),
                                levels = rev(levels_top_to_bottom))),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_shrb, color = "darkgrey"
  ) +
  scale_fill_manual(values = pft_colors, name = "PFT") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = expression(paste("Biomass-weighted SRL (m g"^-1, ")")),
    y = "Depth (cm)"
  ) +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text  = element_text(size = 18)
  )

srl_depth_root



# RTD

# with weighted means
# Summarize per PFT
by_pft <- root_vga_depth %>%
  filter(organ == "root") %>%
  group_by(depth_std, thaw_stage, PFT) %>%
  summarise(
    biomass_sum = sum(mass_gm2_std, na.rm = TRUE),
    RTD_mean    = mean(RTD_gcm3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv)
  )

# Compute weighted mean SRL per depth × thaw stage using weighted.mean()
weighted_RTD <- by_pft %>%
  group_by(depth_std, thaw_stage) %>%
  summarise(
    SRL_weighted = weighted.mean(RTD_mean, biomass_sum, na.rm = TRUE),
    .groups = "drop"
  )

# Compute relative contribution of each PFT to the weighted mean
by_pft <- by_pft %>%
  left_join(weighted_RTD, by = c("depth_std", "thaw_stage")) %>%
  group_by(depth_std, thaw_stage) %>%
  mutate(
    biomass_share = biomass_sum / sum(biomass_sum, na.rm = TRUE),
    RTD_contrib = biomass_share * RTD_mean,
    group_y = factor(paste(depth_std, thaw_stage), levels = group_lv)
  ) %>%
  ungroup()

# top to bottom order (depth first, then thaw), then reverse so 0–10 rows appear on top
levels_top_to_bottom <- expand.grid(depth = depth_lv, thaw = thaw_lv,
                                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) |>
  dplyr::arrange(match(depth, depth_lv), match(thaw, thaw_lv)) |>
  dplyr::transmute(lbl = paste(depth, thaw)) |>
  dplyr::pull(lbl)

# make sure factors + group_y use that order
by_pft <- by_pft %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom)) 
  )


# Add 0s for missing combinations (so 20–30 fully_thawed is included)
by_pft_complete <- by_pft %>%
  tidyr::complete(
    depth_std,
    thaw_stage,
    PFT,
    fill = list(RTD_contrib = 0)
  ) %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom))
  )

# Then plot the completed dataset
rtd_depth_root <- ggplot(by_pft_complete, aes(x = RTD_contrib, y = group_y, fill = PFT)) +
  geom_col(position = "stack") +
  geom_point(
    data = root_vga_depth %>% filter(organ == "root", PFT == "herbaceous"),
    aes(x = RTD_gcm3, y = factor(paste(depth_std, thaw_stage),
                                 levels = rev(levels_top_to_bottom))),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_herb, color = "darkgrey"
  ) +
  geom_point(
    data = root_vga_depth %>% filter(organ == "root", PFT == "shrub"),
    aes(x = RTD_gcm3, y = factor(paste(depth_std, thaw_stage),
                                 levels = rev(levels_top_to_bottom))),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_shrb, color = "darkgrey"
  ) +
  scale_fill_manual(values = pft_colors, name = "PFT") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = expression(paste("Tissue density (g cm"^-3, ")")),
    y = "Depth (cm)"
  ) +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text  = element_text(size = 18)
  )

rtd_depth_root

# TD rhizome

# with weighted means
# Summarize per PFT
by_pft <- root_vga_depth %>%
  filter(organ == "rhizome") %>%
  group_by(depth_std, thaw_stage, PFT) %>%
  summarise(
    biomass_sum = sum(mass_gm2_std, na.rm = TRUE),
    RTD_mean    = mean(RTD_gcm3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv)
  )

# Compute weighted mean SRL per depth × thaw stage using weighted.mean()
weighted_RTD <- by_pft %>%
  group_by(depth_std, thaw_stage) %>%
  summarise(
    SRL_weighted = weighted.mean(RTD_mean, biomass_sum, na.rm = TRUE),
    .groups = "drop"
  )

# 3. Compute relative contribution of each PFT to the weighted mean
by_pft <- by_pft %>%
  left_join(weighted_RTD, by = c("depth_std", "thaw_stage")) %>%
  group_by(depth_std, thaw_stage) %>%
  mutate(
    biomass_share = biomass_sum / sum(biomass_sum, na.rm = TRUE),
    RTD_contrib = biomass_share * RTD_mean,
    group_y = factor(paste(depth_std, thaw_stage), levels = group_lv)
  ) %>%
  ungroup()

# top to bottom order (depth first, then thaw), then reverse so 0–10 rows appear on top
levels_top_to_bottom <- expand.grid(depth = depth_lv, thaw = thaw_lv,
                                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) |>
  dplyr::arrange(match(depth, depth_lv), match(thaw, thaw_lv)) |>
  dplyr::transmute(lbl = paste(depth, thaw)) |>
  dplyr::pull(lbl)

# make sure factors + group_y use that order
by_pft <- by_pft %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom))
  )


# Add 0s for missing combinations (so 20–30 fully_thawed is included)
by_pft_complete <- by_pft %>%
  tidyr::complete(
    depth_std,
    thaw_stage,
    PFT,
    fill = list(RTD_contrib = 0)
  ) %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom))
  )

# Then plot the completed dataset
rtd_depth_rhiz <- ggplot(by_pft_complete, aes(x = RTD_contrib, y = group_y, fill = PFT)) +
  geom_col(position = "stack") +
  geom_point(
    data = root_vga_depth %>% filter(organ == "rhizome", PFT == "herbaceous"),
    aes(x = RTD_gcm3, y = factor(paste(depth_std, thaw_stage),
                                 levels = rev(levels_top_to_bottom))),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_herb, color = "#343434cc"
  ) +
  geom_point(
    data = root_vga_depth %>% filter(organ == "rhizome", PFT == "shrub"),
    aes(x = RTD_gcm3, y = factor(paste(depth_std, thaw_stage),
                                 levels = rev(levels_top_to_bottom))),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_shrb, color = "#343434cc"
  ) +
  scale_fill_manual(values = pft_colors, name = "PFT") +
  scale_x_continuous(
    limits = c(0, 3),
    breaks = seq(0, 3, by = 1),
    expand = c(0, 0)            
  ) +
  labs(
    x = expression(paste("Tissue density (g cm"^-3, ")")),
    y = "Depth (cm)"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text  = element_text(size = 18)
  )

rtd_depth_rhiz

## diameter

# with weighted means

# Summarize per PFT
by_pft <- root_vga_depth %>%
  filter(organ == "root") %>%
  group_by(depth_std, thaw_stage, PFT) %>%
  summarise(
    biomass_sum = sum(mass_gm2_std, na.rm = TRUE),
    d_mean    = mean(avg_d_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv)
  )

# Compute weighted mean SRL per depth × thaw stage using weighted.mean()
weighted_d <- by_pft %>%
  group_by(depth_std, thaw_stage) %>%
  summarise(
    SRL_weighted = weighted.mean(d_mean, biomass_sum, na.rm = TRUE),
    .groups = "drop"
  )

# Compute relative contribution of each PFT to the weighted mean
by_pft <- by_pft %>%
  left_join(weighted_d, by = c("depth_std", "thaw_stage")) %>%
  group_by(depth_std, thaw_stage) %>%
  mutate(
    biomass_share = biomass_sum / sum(biomass_sum, na.rm = TRUE),
    d_contrib = biomass_share * d_mean,
    group_y = factor(paste(depth_std, thaw_stage), levels = group_lv)
  ) %>%
  ungroup()


# top to bottom order (depth first, then thaw), then reverse so 0–10 rows appear on top
levels_top_to_bottom <- expand.grid(depth = depth_lv, thaw = thaw_lv,
                                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) |>
  dplyr::arrange(match(depth, depth_lv), match(thaw, thaw_lv)) |>
  dplyr::transmute(lbl = paste(depth, thaw)) |>
  dplyr::pull(lbl)

# make sure factors + group_y use that order
by_pft <- by_pft %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom))
  )


# Add 0s for missing combinations (so 20–30 fully_thawed is included)
by_pft_complete <- by_pft %>%
  tidyr::complete(
    depth_std,
    thaw_stage,
    PFT,
    fill = list(d_contrib = 0)
  ) %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom))
  )

# Then plot the completed dataset
d_depth_root <- ggplot(by_pft_complete, aes(x = d_contrib, y = group_y, fill = PFT)) +
  geom_col(position = "stack") +
  geom_point(
    data = root_vga_depth %>% filter(organ == "root", PFT == "herbaceous"),
    aes(x = avg_d_mm, y = factor(paste(depth_std, thaw_stage),
                                 levels = rev(levels_top_to_bottom))),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_herb, color = "darkgrey"
  ) +
  geom_point(
    data = root_vga_depth %>% filter(organ == "root", PFT == "shrub"),
    aes(x = avg_d_mm, y = factor(paste(depth_std, thaw_stage),
                                 levels = rev(levels_top_to_bottom))),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_shrb, color = "darkgrey"
  ) +
  scale_fill_manual(values = pft_colors, name = "PFT") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = expression(paste("Diameter (mm)")),
    y = "Depth (cm)"
  ) +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text  = element_text(size = 18)
  )

d_depth_root


# rhizome diameter

# with weighted means
# Summarize per PFT
by_pft <- root_vga_depth %>%
  filter(organ == "rhizome") %>%
  group_by(depth_std, thaw_stage, PFT) %>%
  summarise(
    biomass_sum = sum(mass_gm2_std, na.rm = TRUE),
    d_mean    = mean(avg_d_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv)
  )

# Compute weighted mean SRL per depth × thaw stage using weighted.mean()
weighted_d <- by_pft %>%
  group_by(depth_std, thaw_stage) %>%
  summarise(
    SRL_weighted = weighted.mean(d_mean, biomass_sum, na.rm = TRUE),
    .groups = "drop"
  )

# Compute relative contribution of each PFT to the weighted mean
by_pft <- by_pft %>%
  left_join(weighted_d, by = c("depth_std", "thaw_stage")) %>%
  group_by(depth_std, thaw_stage) %>%
  mutate(
    biomass_share = biomass_sum / sum(biomass_sum, na.rm = TRUE),
    d_contrib = biomass_share * d_mean,
    group_y = factor(paste(depth_std, thaw_stage), levels = group_lv)
  ) %>%
  ungroup()

# top to bottom order (depth first, then thaw), then reverse so 0–10 rows appear on top
levels_top_to_bottom <- expand.grid(depth = depth_lv, thaw = thaw_lv,
                                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) |>
  dplyr::arrange(match(depth, depth_lv), match(thaw, thaw_lv)) |>
  dplyr::transmute(lbl = paste(depth, thaw)) |>
  dplyr::pull(lbl)

# make sure factors + group_y use that order
by_pft <- by_pft %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom)) 
  )


# Add 0s for missing combinations (so 20–30 fully_thawed is included)
by_pft_complete <- by_pft %>%
  tidyr::complete(
    depth_std,
    thaw_stage,
    PFT,
    fill = list(d_contrib = 0)
  ) %>%
  mutate(
    depth_std  = factor(depth_std,  levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom))
  )

# Then plot the completed dataset
d_depth_rhiz <- ggplot(by_pft_complete, aes(x = d_contrib, y = group_y, fill = PFT)) +
  geom_col(position = "stack") +
  geom_point(
    data = root_vga_depth %>% filter(organ == "rhizome", PFT == "herbaceous"),
    aes(x = avg_d_mm, y = factor(paste(depth_std, thaw_stage),
                                 levels = rev(levels_top_to_bottom))),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_herb, color = "#343434cc"
  ) +
  geom_point(
    data = root_vga_depth %>% filter(organ == "rhizome", PFT == "shrub"),
    aes(x = avg_d_mm, y = factor(paste(depth_std, thaw_stage),
                                 levels = rev(levels_top_to_bottom))),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_shrb, color = "#343434cc"
  ) +
  scale_fill_manual(values = pft_colors, name = "PFT") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = expression(paste("Diameter (mm)")),
    y = "Depth (cm)"
  ) +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text  = element_text(size = 18)
  )

d_depth_rhiz


## surface area
## root

root_vga_depth_complete <- root_vga_depth %>%
  filter(organ == "root") %>%
  complete(depth_std, thaw_stage, PFT, organ, fill = list(sa_cm2m2_std = 0)) %>%
  filter(organ == "root") %>%
  mutate(
    depth_std  = factor(depth_std, levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    # Make combined y factor and reverse levels so 0–10 rows are at the top
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom))
  )

# plot
sa_depth_root <- ggplot(root_vga_depth_complete,
                        aes(x = sa_cm2m2_std, y = group_y, fill = PFT)) +
  geom_col(position = "stack") +
  geom_point(
    data = ~ dplyr::filter(.x, PFT == "herbaceous"),
    position = position_stack(vjust = 0.5),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_herb, color = "darkgrey"
  ) +
  geom_point(
    data = ~ dplyr::filter(.x, PFT == "shrub"),
    position = position_stack(vjust = 0.5),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_shrb, color = "darkgrey"
  ) +
  scale_fill_manual(values = pft_colors, name = "PFT") +
  labs(
    x = expression(paste("Surface area (cm"^2, " m"^-2, ")")),
    y = "Depth (cm)"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, NA)
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text  = element_text(size = 18)
  )

sa_depth_root


# rhizome surface area
root_vga_depth_complete <- root_vga_depth %>%
  filter(organ == "rhizome") %>%
  complete(depth_std, thaw_stage, PFT, organ, fill = list(sa_cm2m2_std = 0)) %>%
  filter(organ == "rhizome") %>%
  mutate(
    depth_std  = factor(depth_std, levels = depth_lv),
    thaw_stage = factor(thaw_stage, levels = thaw_lv),
    # Make combined y factor and reverse levels so 0–10 rows are at the top
    group_y    = factor(paste(depth_std, thaw_stage),
                        levels = rev(levels_top_to_bottom))
  )

# plot
sa_depth_rhiz <- ggplot(root_vga_depth_complete,
                        aes(x = sa_cm2m2_std, y = group_y, fill = PFT)) +
  geom_col(position = "stack") +
  geom_point(
    data = ~ dplyr::filter(.x, PFT == "herbaceous"),
    position = position_stack(vjust = 0.5),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_herb, color = "#343434cc"
  ) +
  geom_point(
    data = ~ dplyr::filter(.x, PFT == "shrub"),
    position = position_stack(vjust = 0.5),
    shape = 21, size = 3, stroke = 0.3, alpha = 0.8,
    fill = pt_fill_shrb, color = "#343434cc"
  ) +
  scale_fill_manual(values = pft_colors, name = "PFT") +
  labs(
    x = expression(paste("Surface area (cm"^2, " m"^-2, ")")),
    y = "Depth (cm)"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),  # 5% padding on the right
    limits = c(0, NA)
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text  = element_text(size = 18)
  )

sa_depth_rhiz



#### combine into a figure (Fig. B3)
# Remove y-axis text/ticks from the non-leftmost plots
sa_depth_root  <- sa_depth_root  + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank())
d_depth_root   <- d_depth_root   + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank())
srl_depth_root <- srl_depth_root + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank())

# arrange with shared y-axis title
fig_root_depth <- ggarrange(
  mass_depth_root, sa_depth_root, d_depth_root,
  rtd_depth_root, srl_depth_root,
  ncol = 3, nrow = 2,
  align = "hv",
  common.legend = TRUE,
  legend = "top"
)

fig_root_depth <- annotate_figure(
  fig_root_depth,
  left = text_grob("Depth (cm)", rot = 90, vjust = 1, size = 20)
)

fig_root_depth


# rhizome figure (Fig. B5)

# Remove y-axis text/ticks from the non-leftmost plots
sa_depth_rhiz  <- sa_depth_rhiz  + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank())
rtd_depth_rhiz <- rtd_depth_rhiz + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank())

# Now arrange with shared y-axis title
fig_rhiz_depth <- ggarrange(
  mass_depth_rhiz, sa_depth_rhiz, d_depth_rhiz,
  rtd_depth_rhiz,
  ncol = 2, nrow = 2,
  align = "hv",
  common.legend = TRUE,
  legend = "top"
)

fig_rhiz_depth


### depth-resolved trait statistics (descriptive) ###

# biomass
# Calculate within-plot CV% across depths for each PFT and thaw_stage
cv_by_plot <- root_vga_depth %>%
  filter(organ == "root") %>%
  group_by(thaw_stage, PFT, plot_id) %>%
  summarise(
    mean_mass = mean(mass_gm2_std, na.rm = TRUE),
    sd_mass  = sd(mass_gm2_std, na.rm = TRUE),
    cv_percent  = (sd_mass / mean_mass) * 100,
    .groups = "drop"
  )


# Calculate mean CV per thaw_stage × PFT group
mean_cv_per_thaw_stage_PFT <- cv_by_plot %>%
  group_by(thaw_stage, PFT) %>%
  summarise(
    mean_cv_percent = mean(cv_percent, na.rm = TRUE),
    sd_cv_percent   = sd(cv_percent, na.rm = TRUE),
    n               = n(),
    se_cv_percent   = sd_cv_percent / sqrt(n),
    .groups = "drop"
  )


# SRL
cv_by_plot <- root_vga_depth %>%
  filter(organ == "root") %>%
  group_by(thaw_stage, PFT, plot_id) %>%
  summarise(
    mean_length = mean(SRL_mg1, na.rm = TRUE),
    sd_length   = sd(SRL_mg1, na.rm = TRUE),
    cv_percent  = (sd_length / mean_length) * 100,
    .groups = "drop"
  )

# Calculate mean CV per thaw_stage × PFT group
mean_cv_per_thaw_stage_PFT <- cv_by_plot %>%
  group_by(thaw_stage, PFT) %>%
  summarise(
    mean_cv_percent = mean(cv_percent, na.rm = TRUE),
    sd_cv_percent   = sd(cv_percent, na.rm = TRUE),
    n               = n(),
    se_cv_percent   = sd_cv_percent / sqrt(n),
    .groups = "drop"
  )

# Root TD
cv_by_plot <- root_vga_depth %>%
  filter(organ == "root") %>%
  group_by(thaw_stage, PFT, plot_id) %>%
  summarise(
    mean_length = mean(RTD_gcm3, na.rm = TRUE),
    sd_length   = sd(RTD_gcm3, na.rm = TRUE),
    cv_percent  = (sd_length / mean_length) * 100,
    .groups = "drop"
  )


# Calculate mean CV per thaw_stage × PFT group
mean_cv_per_thaw_stage_PFT <- cv_by_plot %>%
  group_by(thaw_stage, PFT) %>%
  summarise(
    mean_cv_percent = mean(cv_percent, na.rm = TRUE),
    sd_cv_percent   = sd(cv_percent, na.rm = TRUE),
    n               = n(),
    se_cv_percent   = sd_cv_percent / sqrt(n),
    .groups = "drop"
  )

# means, min, max per depth per thaw stage and PFT

# total root length
root_length_summary <- root_vga_depth %>%
  filter(organ == "root") %>%
  group_by(depth_std, thaw_stage, PFT) %>%
  summarise(
    mean_length = mean(len_kmm2_std, na.rm = TRUE),
    sd_length   = sd(len_kmm2_std, na.rm = TRUE),
    min_length  = min(len_kmm2_std, na.rm = TRUE),
    max_length  = max(len_kmm2_std, na.rm = TRUE),
    n           = sum(!is.na(len_kmm2_std)),
    .groups = "drop"
  )

# SRL
root_srl_summary <- root_vga_depth %>%
  filter(organ == "root") %>%
  group_by(depth_std, thaw_stage, PFT) %>%
  summarise(
    mean_length = mean(SRL_mg1, na.rm = TRUE),
    sd_length   = sd(SRL_mg1, na.rm = TRUE),
    min_length  = min(SRL_mg1, na.rm = TRUE),
    max_length  = max(SRL_mg1, na.rm = TRUE),
    n           = sum(!is.na(SRL_mg1)),
    .groups = "drop"
  )

# root biomass
root_mass_summary <- root_vga_depth %>%
  filter(organ == "root") %>%
  group_by(depth_std, thaw_stage, PFT) %>%
  summarise(
    mean = mean(mass_gm2_std, na.rm = TRUE),
    sd   = sd(mass_gm2_std, na.rm = TRUE),
    min  = min(mass_gm2_std, na.rm = TRUE),
    max  = max(mass_gm2_std, na.rm = TRUE),
    n           = sum(!is.na(mass_gm2_std)),
    .groups = "drop"
  )

# rhizome biomass
rhizome_mass_summary <- root_vga_depth %>%
  filter(organ == "rhizome") %>%
  group_by(depth_std, thaw_stage, PFT) %>%
  summarise(
    mean_mass = mean(mass_gm2_std, na.rm = TRUE),
    sd_mass   = sd(mass_gm2_std, na.rm = TRUE),
    min_mass  = min(mass_gm2_std, na.rm = TRUE),
    max_mass  = max(mass_gm2_std, na.rm = TRUE),
    n           = sum(!is.na(mass_gm2_std)),
    .groups = "drop"
  )

# rhizome TD
rhizome_rtd_summary <- root_vga_depth %>%
  filter(organ == "rhizome") %>%
  group_by(depth_std, thaw_stage, PFT) %>%
  summarise(
    mean_mass = mean(RTD_gcm3, na.rm = TRUE),
    sd_mass   = sd(RTD_gcm3, na.rm = TRUE),
    min_mass  = min(RTD_gcm3, na.rm = TRUE),
    max_mass  = max(RTD_gcm3, na.rm = TRUE),
    n           = sum(!is.na(RTD_gcm3)),
    .groups = "drop"
  )


# percentage of shrub roots located in the top 0-10 cm

# Calculate per-plot percent of root length in 0–10 cm, grouped by thaw_stage
shrub_root_plot_pct <- root_vga_depth %>%
  filter(organ == "root", PFT == "shrub") %>%
  group_by(thaw_stage, plot_id) %>%
  summarise(
    root_0_10 = sum(mass_gm2_std[depth_std == "0-10"], na.rm = TRUE),
    root_total = sum(mass_gm2_std, na.rm = TRUE),
    percent_0_10 = (root_0_10 / root_total) * 100,
    .groups = "drop"
  )

# Summarize across plots within each thaw_stage
summary_percent_by_thaw <- shrub_root_plot_pct %>%
  group_by(thaw_stage) %>%
  summarise(
    mean_percent = mean(percent_0_10, na.rm = TRUE),
    sd_percent   = sd(percent_0_10, na.rm = TRUE),
    n            = n(),
    se_percent   = sd_percent / sqrt(n),
    .groups = "drop"
  )

# percentage of herbaceous roots
# Step 1: Calculate per-plot percent of root length in 0–10 cm, grouped by thaw_stage
herb_root_plot_pct <- root_vga_depth %>%
  filter(organ == "root", PFT == "herbaceous") %>%
  group_by(thaw_stage, plot_id) %>%
  summarise(
    root_10_20 = sum(mass_gm2_std[depth_std == "0-10"], na.rm = TRUE),
    root_total = sum(mass_gm2_std, na.rm = TRUE),
    percent_10_20 = (root_10_20 / root_total) * 100,
    .groups = "drop"
  )

# Summarize across plots within each thaw_stage
summary_percent_by_thaw <- herb_root_plot_pct %>%
  group_by(thaw_stage) %>%
  summarise(
    mean_percent = mean(percent_10_20, na.rm = TRUE),
    sd_percent   = sd(percent_10_20, na.rm = TRUE),
    n            = n(),
    se_percent   = sd_percent / sqrt(n),
    .groups = "drop"
  )

# percentage of rhizome biomass in the top peat

# shrub
# Calculate per-plot percent of rhizome mass in 0–10 cm, grouped by thaw_stage
shrub_rhi_plot_pct <- root_vga_depth %>%
  filter(organ == "rhizome", PFT == "shrub") %>%
  group_by(thaw_stage, plot_id) %>%
  summarise(
    root_0_10 = sum(mass_gm2_std[depth_std == "0-10"], na.rm = TRUE),
    root_total = sum(mass_gm2_std, na.rm = TRUE),
    percent_0_10 = (root_0_10 / root_total) * 100,
    .groups = "drop"
  )

# Summarize across plots within each thaw_stage
summary_percent_by_thaw <- shrub_rhi_plot_pct %>%
  group_by(thaw_stage) %>%
  summarise(
    mean_percent = mean(percent_0_10, na.rm = TRUE),
    sd_percent   = sd(percent_0_10, na.rm = TRUE),
    n            = n(),
    se_percent   = sd_percent / sqrt(n),
    .groups = "drop"
  )


# herbaceous
# Calculate per-plot percent of rhizome mass in 0–10 cm, grouped by thaw_stage
herb_rhi_plot_pct <- root_vga_depth %>%
  filter(organ == "rhizome", PFT == "herbaceous") %>%
  group_by(thaw_stage, plot_id) %>%
  summarise(
    root_0_10 = sum(mass_gm2_std[depth_std == "0-10"], na.rm = TRUE),
    root_total = sum(mass_gm2_std, na.rm = TRUE),
    percent_0_10 = (root_0_10 / root_total) * 100,
    .groups = "drop"
  )

# Summarize across plots within each thaw_stage
summary_percent_by_thaw <- herb_rhi_plot_pct %>%
  group_by(thaw_stage) %>%
  summarise(
    mean_percent = mean(percent_0_10, na.rm = TRUE),
    sd_percent   = sd(percent_0_10, na.rm = TRUE),
    n            = n(),
    se_percent   = sd_percent / sqrt(n),
    .groups = "drop"
  )


### variation in traits between plots within thaw stage ###
# note: code has been combined with GPT for easier use (the old version of the code created one plot and summary for each individual trait
# that were changed manually across the code, and was not user-friendly. Please contact Tiia Määttä for the old version if interested)

# prepare trait columns
root_vga_plot <- root_vga_plot %>%
  mutate(
    sa_m2m2_std = sa_cm2m2_std / 10000,
    len_kmm2_std = len_mm2_std / 1000,
    thaw_stage = factor(
      thaw_stage,
      levels = c("intact", "partly_thawed", "fully_thawed")
    )
  )

# define traits
trait_info <- tibble::tribble(
  ~trait,          ~label,                                                        ~aggregation,
  "mass_gm2_std",  expression(paste("Biomass (g m"^-2, ")")),                    "sum",
  "sa_m2m2_std",   expression(paste("Surface area (m"^2, " m"^-2, ")")),         "sum",
  "len_kmm2_std",  expression(paste("Total root length (km m"^-2, ")")),         "sum",
  "avg_d_mm",      expression(paste("Diameter (mm)")),                           "weighted",
  "RTD_gcm3",      expression(paste("Tissue density (g cm"^-3, ")")),            "weighted",
  "SRL_mg1",       expression(paste("Specific root length (m g"^-1, ")")),       "weighted"
)

# no SRL for rhizomes
trait_info_by_organ <- list(
  root = trait_info,
  rhizome = trait_info %>% filter(trait != "SRL_mg1")
)

safe_cv <- function(x) {
  v <- raster::cv(x, na.rm = TRUE)
  if (is.finite(v)) v else NA_real_
}

safe_weighted_mean <- function(x, w) {
  if (sum(w, na.rm = TRUE) > 0) {
    weighted.mean(x, w = w, na.rm = TRUE)
  } else {
    NA_real_
  }
}

# calculate variation summaries
calculate_trait_variation <- function(data, organ_sel, traits) {
  
  trait_cols <- traits$trait
  
  pft_variation <- data %>%
    filter(organ == organ_sel) %>%
    group_by(thaw_stage, PFT) %>%
    summarise(
      across(
        all_of(trait_cols),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd   = ~ sd(.x, na.rm = TRUE),
          cv   = ~ safe_cv(.x),
          min  = ~ min(.x, na.rm = TRUE),
          max  = ~ max(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
  
  across_by_chamber <- data %>%
    filter(organ == organ_sel) %>%
    group_by(thaw_stage, chamber_id) %>%
    summarise(
      across(
        all_of(traits$trait[traits$aggregation == "weighted"]),
        ~ safe_weighted_mean(.x, mass_gm2_std),
        .names = "{.col}"
      ),
      across(
        all_of(traits$trait[traits$aggregation == "sum"]),
        ~ sum(.x, na.rm = TRUE),
        .names = "{.col}"
      ),
      .groups = "drop"
    )
  
  across_variation <- across_by_chamber %>%
    group_by(thaw_stage) %>%
    summarise(
      across(
        all_of(trait_cols),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd   = ~ sd(.x, na.rm = TRUE),
          cv   = ~ safe_cv(.x),
          min  = ~ min(.x, na.rm = TRUE),
          max  = ~ max(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    ) %>%
    mutate(PFT = "across") %>%
    relocate(PFT, .after = thaw_stage)
  
  bind_rows(pft_variation, across_variation) %>%
    mutate(organ = organ_sel) %>%
    relocate(organ, .before = thaw_stage)
}

root_trait_variation_all <- calculate_trait_variation(
  root_vga_plot,
  organ_sel = "root",
  traits = trait_info_by_organ$root
)

rhi_trait_variation_all <- calculate_trait_variation(
  root_vga_plot,
  organ_sel = "rhizome",
  traits = trait_info_by_organ$rhizome
)

trait_variation_all <- bind_rows(
  root_trait_variation_all,
  rhi_trait_variation_all
)

trait_variation_all


# visualize

pft_colors <- c(
  "herbaceous" = "#213448",
  "shrub" = "#94B4C1"
)

# across = noPFT
slot_levels <- c("herbaceous", "shrub", "across")

slot_offset <- c(
  herbaceous = -0.28,
  shrub = 0,
  across = 0.28
)

plot_trait_variation <- function(data, organ_sel, trait_col, y_label, aggregation) {
  
  df <- data %>%
    filter(organ == organ_sel)
  
  summary_df <- df %>%
    group_by(thaw_stage, PFT) %>%
    summarise(
      mean_value = mean(.data[[trait_col]], na.rm = TRUE),
      sd_value = sd(.data[[trait_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      slot = factor(PFT, levels = slot_levels),
      x_num = as.numeric(thaw_stage) + slot_offset[as.character(slot)],
      ymin = pmax(mean_value - sd_value, 0),
      ymax = mean_value + sd_value
    )
  
  across_df_chamber <- df %>%
    group_by(thaw_stage, chamber_id) %>%
    summarise(
      value = case_when(
        aggregation == "sum" ~ sum(.data[[trait_col]], na.rm = TRUE),
        aggregation == "weighted" ~ safe_weighted_mean(.data[[trait_col]], mass_gm2_std),
        TRUE ~ NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      slot = factor("across", levels = slot_levels),
      x_num = as.numeric(thaw_stage) + slot_offset[as.character(slot)]
    )
  
  across_stage_summary <- across_df_chamber %>%
    group_by(thaw_stage) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      slot = factor("across", levels = slot_levels),
      x_num = as.numeric(thaw_stage) + slot_offset[as.character(slot)],
      ymin = pmax(mean_value - sd_value, 0),
      ymax = mean_value + sd_value
    )
  
  df_pts <- df %>%
    mutate(
      slot = factor(PFT, levels = slot_levels),
      x_num = as.numeric(thaw_stage) + slot_offset[as.character(slot)]
    )
  
  ggplot() +
    geom_point(
      data = df_pts,
      aes(x = x_num, y = .data[[trait_col]], color = PFT),
      shape = 16, size = 4, alpha = 0.6,
      show.legend = FALSE
    ) +
    geom_errorbar(
      data = summary_df,
      aes(x = x_num, ymin = ymin, ymax = ymax),
      width = 0.06, linewidth = 0.8,
      color = "#636363"
    ) +
    geom_point(
      data = summary_df,
      aes(x = x_num, y = mean_value, fill = PFT),
      shape = 21, size = 8,
      color = "black", stroke = 0.5,
      alpha = 0.85,
      show.legend = FALSE
    ) +
    geom_point(
      data = across_df_chamber,
      aes(x = x_num, y = value),
      shape = 21, fill = "white",
      color = "black", size = 4,
      alpha = 0.6
    ) +
    geom_errorbar(
      data = across_stage_summary,
      aes(x = x_num, ymin = ymin, ymax = ymax),
      width = 0.08, linewidth = 0.8,
      color = "#636363"
    ) +
    geom_point(
      data = across_stage_summary,
      aes(x = x_num, y = mean_value),
      shape = 21, fill = "white",
      color = "black", size = 8,
      stroke = 0.5, alpha = 0.85
    ) +
    scale_fill_manual(values = pft_colors) +
    scale_color_manual(values = pft_colors) +
    scale_x_continuous(
      breaks = seq_along(levels(data$thaw_stage)),
      labels = levels(data$thaw_stage),
      expand = expansion(mult = c(0.06, 0.08))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    coord_cartesian(ylim = c(0, NA)) +
    labs(x = "", y = y_label) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 20),
      plot.title = element_blank()
    )
}

plot_list <- imap(
  trait_info_by_organ,
  ~ pmap(
    .x,
    function(trait, label, aggregation) {
      plot_trait_variation(
        data = root_vga_plot,
        organ_sel = .y,
        trait_col = trait,
        y_label = label,
        aggregation = aggregation
      )
    }
  )
)

names(plot_list$root) <- paste0("root_", trait_info_by_organ$root$trait)
names(plot_list$rhizome) <- paste0("rhizome_", trait_info_by_organ$rhizome$trait)

combined <- ggarrange(
  plotlist = c(
    plot_list$root["root_mass_gm2_std"],
    plot_list$root["root_sa_m2m2_std"],
    plot_list$root["root_avg_d_mm"],
    plot_list$root["root_RTD_gcm3"],
    plot_list$root["root_SRL_mg1"],
    plot_list$rhizome["rhizome_mass_gm2_std"],
    plot_list$rhizome["rhizome_sa_m2m2_std"],
    plot_list$rhizome["rhizome_avg_d_mm"],
    plot_list$rhizome["rhizome_RTD_gcm3"]
  ),
  ncol = 5,
  nrow = 2,
  align = "hv",
  labels = c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)")
)

combined


### trait changes along the thaw gradient ###

# add the missing group for shrub as 0 for fully_thawed

# Define the column names to set to 0
cols_to_zero <- c("avg_d_mm", "dry_wt_g_depth_mean", "tot_len_cm_depth_mean", 
                  "vol_cm3_depth_mean", "sa_cm2_depth_mean", "mass_density_gcm3", 
                  "len_density_cmcm3", "sa_density_cm2cm3", "mass_density_gm3", 
                  "len_density_mm3", "sa_density_cm2m3", "mass_per_area_gm2", 
                  "len_per_area_mm2", "vol_per_area_cm3m2", "sa_per_area_cm2m2", 
                  "mass_gm2_std", "len_mm2_std", "vol_cm3m2_std", "sa_cm2m2_std", 
                  "depth_mean_mass_gm2_std", "depth_mean_len_mm2_std", 
                  "depth_mean_vol_cm3m2_std", "depth_mean_sa_cm2m2_std", 
                  "RTD_gcm3", "SRL_mg1")

# Identify all unique combinations for which shrub rows are missing
fully_thawed_herb_rows <- root_vga_plot %>%
  filter(thaw_stage == "fully_thawed", PFT == "herbaceous", organ %in% c("root", "rhizome")) %>%
  dplyr::select(plot_id, chamber_id, thaw_stage, organ, bulk_density_gcm3, GWC)

# Add PFT and fill zero columns manually
new_shrub_rows <- fully_thawed_herb_rows %>%
  mutate(PFT = "shrub")

# Add the zero columns 
for (col in cols_to_zero) {
  new_shrub_rows[[col]] <- 0
}

# Fill in any other missing columns as NA to match original data
missing_cols <- setdiff(names(root_vga_plot), names(new_shrub_rows))
for (col in missing_cols) {
  new_shrub_rows[[col]] <- 0
}

# Reorder columns to match original dataset
new_shrub_rows <- new_shrub_rows %>% dplyr::select(all_of(names(root_vga_plot)))

# Bind to original data
root_vga_plot <- bind_rows(root_vga_plot, new_shrub_rows)

# fill in the VGA etc data to the shrub group

# List of stem-related variables 
stem_vars <- c(
  "aer_stem_density_nom2", "nonaer_stem_density_nom2", "herb_stem_density_nom2", "shrub_stem_density_nom2",
  "aer_stem_density_sd", "nonaer_stem_density_sd", "herb_stem_density_sd", "shrub_stem_density_sd",
  "aer_GA_m2m2", "nonaer_GA_m2m2", "herb_GA_m2m2", "shrub_GA_m2m2"
)

# Extract the donor data (herbaceous in fully_thawed) — one row per plot_id
donor_stem_data <- root_vga_plot %>%
  filter(thaw_stage == "fully_thawed", PFT == "herbaceous") %>%
  group_by(plot_id) %>%
  summarise(across(all_of(stem_vars), ~ first(.), .names = "{.col}"), .groups = "drop")

# Apply these values to shrub rows in fully_thawed
root_vga_plot <- root_vga_plot %>%
  left_join(donor_stem_data, by = "plot_id", suffix = c("", ".donor")) %>%
  mutate(across(
    all_of(stem_vars),
    ~ ifelse(PFT == "shrub" & thaw_stage == "fully_thawed", get(paste0(cur_column(), ".donor")), .)
  )) %>%
  dplyr::select(-ends_with(".donor"))

# change NA to 0 
root_vga_plot <- root_vga_plot %>%
  mutate(across(all_of(stem_vars), ~ replace_na(., 0)))

# Arrange by plot_id, thaw_stage, organ, and order PFT so herbaceous comes before shrub
root_vga_plot <- root_vga_plot %>%
  mutate(PFT = factor(PFT, levels = c("herbaceous", "shrub"))) %>%  # ensure correct order
  arrange(thaw_stage, plot_id, organ, PFT)


# save the df

write.csv(root_vga_plot, "path/root_plot_scale.csv", row.names = FALSE)

## thaw stage statistics ##

# TD
root_vga_plot %>%
  group_by(thaw_stage, PFT, organ) %>%
  summarise(
    RTD_min = min(RTD_gcm3, na.rm = TRUE),
    RTD_max = max(RTD_gcm3, na.rm = TRUE),
    RTD_mean = mean(RTD_gcm3, na.rm = TRUE),
    .groups = "drop"
  )

# biomass
root_vga_plot %>%
  group_by(thaw_stage, PFT, organ) %>%
  summarise(
    min = min(mass_gm2_std, na.rm = TRUE),
    max = max(mass_gm2_std, na.rm = TRUE),
    mean = mean(mass_gm2_std, na.rm = TRUE),
    .groups = "drop"
  )


#### PLANT GREEN AREA 

# within thaw stages

# Define GA columns of interest
ga_cols <- c("herb_GA_m2m2", "shrub_GA_m2m2")

# Summarize one GA row per plot_id 
ga_per_plot <- root_vga_plot %>%
  dplyr::select(plot_id, thaw_stage, all_of(ga_cols)) %>%
  group_by(plot_id, thaw_stage) %>%
  summarise(across(all_of(ga_cols), ~ first(.), .names = "{.col}"), .groups = "drop")

# descriptive statistics per thaw stage
variation <- ga_per_plot %>%
  group_by(thaw_stage) %>%
  summarise(across(
    all_of(ga_cols),
    list(
      mean = ~mean(.x, na.rm = TRUE),
      sd   = ~sd(.x, na.rm = TRUE),
      cv   = ~sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE) * 100,
      min  = ~min(.x, na.rm = TRUE),
      max  = ~max(.x, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  ungroup()

variation <- as.data.frame(variation)
variation

# visualize so easier to see

cv_long <- variation %>%
  pivot_longer(
    cols = ends_with("_cv"),
    names_to = "trait",
    values_to = "cv_percent"
  ) 

ggplot(cv_long, aes(x = trait, y = cv_percent, fill = thaw_stage)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    x = "Trait",
    y = "Coefficient of Variation (%)",
    fill = "Thaw stage"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

# proportion of shrubs of the total GA in fully thawed stage

# Define GA columns of interest
ga_cols <- c("shrub_GA_m2m2", "herb_GA_m2m2")

# Filter for thaw_stage == "fully_thawed"
ga_fully_thawed_summary <- root_vga_plot %>%
  filter(thaw_stage == "fully_thawed") %>%
  dplyr::select(plot_id, all_of(ga_cols)) %>%
  distinct() 

# Sum shrub and herbaceous GA across unique plots
total_shrub <- sum(ga_fully_thawed_summary$shrub_GA_m2m2, na.rm = TRUE)
total_herb  <- sum(ga_fully_thawed_summary$herb_GA_m2m2, na.rm = TRUE)

# Calculate proportion of shrub GA
total_shrub / (total_shrub + total_herb)


### between thaw stages

# need to first create a new GA column

root_vga_plot <- root_vga_plot %>%
  mutate(GA_m2m2 = case_when(
    PFT == "herbaceous" ~ herb_GA_m2m2,
    PFT == "shrub" ~ shrub_GA_m2m2,
    TRUE ~ NA_real_ 
  ))

# Define PFT color palette
pft_colors <- c("herbaceous" = "#213448", "shrub" = "#94B4C1")

# Factor order for thaw_stage
root_vga_plot$thaw_stage <- factor(root_vga_plot$thaw_stage, levels = c("intact", "partly_thawed", "fully_thawed"))

# settings for ordering within each thaw stage
slot_levels <- c("herbaceous","shrub","across")
slot_offset <- c(herbaceous = -0.28, shrub = 0, across = +0.28)

df <- root_vga_plot

# PFT summary (mean ± SD of GA_m2m2)
summary_df <- df %>%
  group_by(thaw_stage, PFT) %>%
  summarise(
    mean_value = mean(GA_m2m2, na.rm = TRUE),
    sd_value   = sd(GA_m2m2,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    slot  = factor(PFT, levels = slot_levels),
    x_num = as.numeric(thaw_stage) + slot_offset[as.character(slot)],
    ymin  = pmax(mean_value - sd_value, 0), 
    ymax  = mean_value + sd_value
  )

# PFT points with fixed x positions
df_pts <- df %>%
  mutate(
    slot  = factor(PFT, levels = slot_levels),
    x_num = as.numeric(thaw_stage) + slot_offset[as.character(slot)]
  )

# GA_m2m2 sums per chamber (noPFT)
across_df_chamber <- df %>%
  group_by(thaw_stage, chamber_id) %>%
  summarise(value = sum(GA_m2m2, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    slot  = factor("across", levels = slot_levels),
    x_num = as.numeric(thaw_stage) + slot_offset[as.character(slot)]
  )

# noPFT thaw-stage GA mean ± SD across chambers
across_stage_summary <- across_df_chamber %>%
  group_by(thaw_stage) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    slot  = factor("across", levels = slot_levels),
    x_num = as.numeric(thaw_stage) + slot_offset[as.character(slot)],
    ymin  = pmax(mean_value - sd_value, 0),
    ymax  = mean_value + sd_value
  )

# plot (Fig. B9)
p_ga <- ggplot() +
  geom_point(
    data = df_pts,
    aes(x = x_num, y = GA_m2m2, color = PFT),
    shape = 16, size = 4, alpha = 0.6, show.legend = FALSE
  ) +
  geom_point(
    data = summary_df,
    aes(x = x_num, y = mean_value, fill = PFT),
    shape = 21, size = 8, color = "black", stroke = 0.5, alpha = 0.85, show.legend = FALSE
  ) +
  geom_errorbar(
    data = summary_df,
    aes(x = x_num, ymin = ymin, ymax = ymax),
    width = 0.06, linewidth = 0.8, color = "#636363"
  ) +
  geom_point(
    data = across_df_chamber,
    aes(x = x_num, y = value),
    shape = 21, fill = "white", color = "black",
    size = 4, stroke = 0.6, alpha = 0.6
  ) +
  geom_point(
    data = across_stage_summary,
    aes(x = x_num, y = mean_value),
    shape = 21, fill = "white", color = "black",
    size = 9, stroke = 0.5, alpha=0.85
  ) +
  geom_errorbar(
    data = across_stage_summary,
    aes(x = x_num, ymin = ymin, ymax = ymax),
    width = 0.08, linewidth = 0.8, color = "#636363"
  ) +
  scale_fill_manual(values = pft_colors) +
  scale_color_manual(values = pft_colors) +
  labs(
    x = "",
    y = expression(paste("Green area (m"^2, " m"^-2, ")"))
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title  = element_text(size = 20),
    plot.title  = element_blank()
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(root_vga_plot$thaw_stage)),
    labels = levels(root_vga_plot$thaw_stage),
    expand = expansion(mult = c(0.06, 0.08))
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))

p_ga


# Create a de-duplicated version of the data — one GA value per plot_id
ga_unique_shrub <- root_vga_plot %>%
  filter(PFT == "shrub") %>% 
  dplyr::select(plot_id, thaw_stage, GA_m2m2) %>%
  distinct()

ga_unique_herb <- root_vga_plot %>%
  filter(PFT == "herbaceous") %>%
  dplyr::select(plot_id, thaw_stage, GA_m2m2) %>%
  distinct()

df_shrub <- ga_unique_shrub %>%
  group_by(thaw_stage) %>%
  summarise(
    mean_GA = mean(GA_m2m2, na.rm = TRUE),
    min_GA = min(GA_m2m2, na.rm = TRUE),
    max_GA = max(GA_m2m2, na.rm = TRUE),
    cv_GA_percent = (sd(GA_m2m2, na.rm = TRUE) / mean(GA_m2m2, na.rm = TRUE)) * 100
  )

df_shrub <- as.data.frame(df_shrub)
df_shrub


df_herb <- ga_unique_herb %>%
  group_by(thaw_stage) %>%
  summarise(
    mean_GA = mean(GA_m2m2, na.rm = TRUE),
    min_GA = min(GA_m2m2, na.rm = TRUE),
    max_GA = max(GA_m2m2, na.rm = TRUE),
    cv_GA_percent = (sd(GA_m2m2, na.rm = TRUE) / mean(GA_m2m2, na.rm = TRUE)) * 100
  )

df_herb <- as.data.frame(df_herb)
df_herb


#### check plot similarity ####

### simple statistics based on vegetation surveys

survey_2 <- read.csv("")

                              
## Table 1:
## calculate mean species cover per thaw stage
                              
t <- survey_2 %>%
     mutate(plot_type = ifelse(grepl("rv$", plot_ID), "rv", 
                               ifelse(grepl("tm$", plot_ID), "tm", NA))) %>% 
     group_by(thaw_stage, plot_type, species) %>%
     summarise(mean_species_cover = mean(species_cover, na.rm = TRUE), .groups = "drop") %>%
     arrange(thaw_stage, plot_type, desc(mean_species_cover))
                              
t <- as.data.frame(t)
t
                              
# herb and shrub
                              
PFT_cover_summary <- survey_2 %>%
                     mutate(plot_type = ifelse(grepl("rv$", plot_ID), "rv", 
                                               ifelse(grepl("tm$", plot_ID), "tm", NA))) %>% 
                     distinct(plot_ID, .keep_all = TRUE) %>% 
                     group_by(thaw_stage, plot_type) %>%
                     summarise(
                     mean_herb_cover = mean(herb_cover, na.rm = TRUE),
                     mean_shrub_cover = mean(shrub_cover, na.rm = TRUE),
                           .groups = "drop"
                           ) %>%
                     arrange(thaw_stage, plot_type)
                              
PFT_cover_summary <- as.data.frame(PFT_cover_summary)
                              
# for similarity comparisons:
summary_df <- comparison_data %>% filter(str_detect(plot_ID, "tm$|rv$")) %>% distinct(plot_ID, .keep_all = TRUE) %>% dplyr::select( plot_ID, thaw_stage, herb_cover, shrub_cover, herb_GA_m2m2, shrub_GA_m2m2, soil_pH, mean_TS ) %>% mutate( type = ifelse(str_detect(plot_ID, "tm$"), "core", "chamber") )

summary_df <- comparison_data %>%
  filter(str_detect(plot_ID, "core$|chamber$")) %>%
  distinct(plot_ID, .keep_all = TRUE) %>%
  dplyr::select(
    plot_ID,
    thaw_stage,
    herb_cover,
    shrub_cover,
    herb_GA_m2m2,
    shrub_GA_m2m2,
    soil_pH,
    mean_TS
  ) %>%
  mutate(
    type = ifelse(str_detect(plot_ID, "core$"), "core", "chamber")
  )

summary_df <- as.data.frame(summary_df)

summary_df <- summary_df %>%
  mutate(shrub_GA_m2m2 = ifelse(is.na(shrub_GA_m2m2), 0, shrub_GA_m2m2))

ggplot(summary_df, aes(x = type, y = herb_cover)) +
  geom_boxplot()

ggplot(summary_df, aes(x = type, y = shrub_cover)) +
  geom_boxplot() 
### run a PCA

summary_df$type <- as.factor(summary_df$type)

# Create the input data for PCA
pca_input <- summary_df %>%
  dplyr::select(plot_ID, type, thaw_stage, herb_cover, shrub_cover, herb_GA_m2m2, shrub_GA_m2m2, soil_pH, mean_TS)

# Set rownames to match plot_IDs
rownames(pca_input) <- pca_input$plot_ID

# Extract only numeric data for PCA
pca_data <- pca_input %>%
  dplyr::select(herb_cover, shrub_cover, herb_GA_m2m2, shrub_GA_m2m2, soil_pH, mean_TS)

# Run PCA
pca_result <- PCA(pca_data, scale.unit = TRUE, graph = FALSE)

# Extract PCA scores
pca_coords <- as.data.frame(pca_result$ind$coord)
pca_coords$plot_ID <- rownames(pca_coords)

# Merge back metadata (type and thaw_stage) using plot_ID
pca_coords <- left_join(pca_coords, pca_input %>% dplyr::select(plot_ID, type, thaw_stage), by = "plot_ID")

# Make sure thaw_stage is a factor 
pca_coords$thaw_stage <- as.factor(pca_coords$thaw_stage)

# Extract and scale loadings
loadings <- as.data.frame(pca_result$var$coord)
loadings$variable <- rownames(loadings)
arrow_scale <- 3
loadings <- loadings %>%
  mutate(xend = Dim.1 * arrow_scale,
         yend = Dim.2 * arrow_scale)
type_colors <- c(
  chamber = "#161E54",
  core    = "#F16D34"
)

# plot (Fig. B2- modified further in Inkscape for the manuscript)
ggplot(pca_coords, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = type, shape = thaw_stage), size = 4) +
  stat_ellipse(aes(group = type, color = type), level = 0.95) +
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black"
  ) +
  geom_text(
    data = loadings,
    aes(x = xend, y = yend, label = variable),
    color = "black", vjust = 1.2
  ) +
  scale_color_manual(
    values = type_colors,
    name = "Type"
  ) +
  labs(,
       x = paste0("PC1 (", round(pca_result$eig[1, 2], 1), "%)"),
       y = paste0("PC2 (", round(pca_result$eig[2, 2], 1), "%)")
  ) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
# Euclidean distances

# Extract and scale the numeric predictor variables
scaled_data <- pca_data %>%
  scale() %>% 
  as.data.frame()

# Add group info ("type") from pca_input
scaled_data$type <- pca_input$type

# Compute Euclidean distance on scaled variables
dist_matrix <- dist(scaled_data %>% dplyr::select(-type), method = "euclidean")

# Run adonis2 with the scaled distance matrix and grouping
adonis_result <- adonis2(
  dist_matrix ~ type,
  data = scaled_data,
  permutations = 999,
  method = "euclidean" 
)
# quick visualization of the core and chamber differences in the environmental variable values

# Melt summary_df to long format
summary_dt <- as.data.table(summary_df)
long_df <- melt(summary_dt,
                id.vars = c("plot_ID", "type", "thaw_stage"),
                measure.vars = c("herb_cover", "shrub_cover", "herb_GA_m2m2", "shrub_GA_m2m2", "soil_pH"))

ggplot(long_df, aes(x = variable, y = value, group = plot_ID, color = type)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ thaw_stage) +
  theme_minimal()

# similarity of plots within thaw stages

# fully thawed
# Subset to one thaw stage at a time
fully_thawed_data <- subset(summary_df, thaw_stage == "fully_thawed")

# Select numeric variables
fully_thawed_matrix <- fully_thawed_data %>%
  dplyr::select(herb_cover, shrub_cover, herb_GA_m2m2, shrub_GA_m2m2, soil_pH, mean_TS)

# Standardize the variables
fully_thawed_matrix_scaled <- scale(fully_thawed_matrix)

# Then calculate Euclidean distances
dist_matrix <- dist(fully_thawed_matrix_scaled)
as.matrix(dist_matrix)

# partly thawed

partly_thawed_data <- subset(summary_df, thaw_stage == "partly_thawed")

# Select numeric variables
partly_thawed_matrix <- partly_thawed_data %>%
  dplyr::select(herb_cover, shrub_cover, herb_GA_m2m2, shrub_GA_m2m2, soil_pH, mean_TS)

# Standardize the variables
partly_thawed_matrix_scaled <- scale(partly_thawed_matrix)

# Then calculate Euclidean distances
dist_matrix <- dist(partly_thawed_matrix_scaled)
as.matrix(dist_matrix)

# intact

intact_data <- subset(summary_df, thaw_stage == "intact")

# Select numeric variables
intact_matrix <- intact_data %>%
  dplyr::select(herb_cover, shrub_cover, herb_GA_m2m2, shrub_GA_m2m2, soil_pH, mean_TS)

# Standardize the variables
intact_matrix_scaled <- scale(intact_matrix)

# Then calculate Euclidean distances
dist_matrix <- dist(intact_matrix_scaled)
as.matrix(dist_matrix)


#####################################################
#### RELATIONSHIPS BETWEEN TRAITS AND CH4 FLUXES ####
#####################################################

#### create a grid of trait-CH4 relationships where ch4_ values are on the y axis and 
#### trait values are on the x axis
#### the x axis and y axis scales should be the same per row (y) and per column (x)

#### traits in order of: mass, SA, avg_d_mm, RTD, SRL (root)

#### mass, SA, avg_d_mm, RTD (rhizome)

pft_colors <- c("herbaceous" = "#213448", "shrub" = "#94B4C1")
shape_vals <- c(partly_thawed = 21, fully_thawed = 24, intact = 22)

ch4_order    <- c("ch4_early", "ch4_middle","ch4_peak","ch4_season")
traits_root  <- c("mass_gm2_std","sa_cm2m2_std","avg_d_mm","RTD_gcm3","SRL_mg1")
traits_rhiz  <- c("mass_gm2_std","sa_cm2m2_std","avg_d_mm","RTD_gcm3") 

# Helper to build long data for a given organ and trait set

## ROOT GRID
df_root <- plant_ch4 %>%
  filter(organ == "root", PFT %in% c("herbaceous","shrub")) %>%
  select(chamber_id, thaw_stage, PFT, organ, all_of(traits_root), all_of(ch4_order)) %>%
  pivot_longer(cols = all_of(traits_root), names_to = "trait", values_to = "x") %>%
  pivot_longer(cols = all_of(ch4_order),  names_to = "ch4",   values_to = "y") %>%
  filter(is.finite(x), is.finite(y)) %>%
  mutate(trait = factor(trait, levels = traits_root),
         ch4   = factor(ch4,   levels = ch4_order))

# plot 
p_root_grid <-
  ggplot(df_root, aes(x = x, y = y, shape = thaw_stage)) +
  geom_point(aes(fill = PFT), color = "black", size = 3.5, alpha = 0.9, stroke = 0.5) +
  geom_smooth(aes(color = PFT, group = PFT), method = "lm", se = FALSE, linewidth = 0.5) +
  scale_shape_manual(values = shape_vals) +
  scale_fill_manual(values = pft_colors, name = "PFT") +
  scale_color_manual(values = pft_colors, guide = "none") +
  facet_grid(ch4 ~ trait, scales = "free") +
  scale_x_continuous(labels = label_number()) +
  scale_y_continuous(labels = label_number(),
                     expand = expansion(mult = c(0.1, 0.15))) +
  labs(x = "", y = expression(paste("CH"[4], " flux (mg CH"[4], " m"^-2, " d"^-1, ")"))) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size=16)
  )

p_root_grid


## Rhizome grid
df_rhiz <- plant_ch4 %>%
  filter(organ == "rhizome", PFT %in% c("herbaceous","shrub")) %>%
  select(chamber_id, thaw_stage, PFT, organ, all_of(traits_rhiz), all_of(ch4_order)) %>%
  pivot_longer(cols = all_of(traits_rhiz), names_to = "trait", values_to = "x") %>%
  pivot_longer(cols = all_of(ch4_order),  names_to = "ch4",   values_to = "y") %>%
  filter(is.finite(x), is.finite(y)) %>%
  mutate(trait = factor(trait, levels = traits_rhiz),
         ch4   = factor(ch4,   levels = ch4_order))

df_rhiz <- df_rhiz %>%
  mutate(PFT = factor(PFT, levels = c("shrub", "herbaceous"))) %>%
  arrange(PFT)

df_rhiz <- df_rhiz %>%
  mutate(ch4 = factor(ch4, levels = c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season"))) %>%
  arrange(ch4)

# plot
p_rhiz_grid <-
  ggplot(df_rhiz, aes(x = x, y = y, shape = thaw_stage)) +
  geom_point(aes(fill = PFT), color = "black", size = 3.5, alpha = 0.9, stroke = 0.5) +
  geom_smooth(aes(color = PFT, group = PFT), method = "lm", se = FALSE, linewidth = 0.7) +
  scale_shape_manual(values = shape_vals) +
  scale_fill_manual(values = pft_colors, name = "PFT") +
  scale_color_manual(values = pft_colors, guide = "none") +
  facet_grid(ch4 ~ trait, scales = "free") +
  scale_x_continuous(labels = label_number()) +
  scale_y_continuous(labels = label_number(),
                     expand = expansion(mult = c(0.1, 0.15))) +
  labs(x = "", y = expression(paste("FCH"[4], " (nmol CH"[4], " m"^-2, " s"^-1, ")"))) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    legend.text = element_text(size=16)
  )

p_rhiz_grid


## calculate simple linear regression R2s and slopes
## based on centered and scaled traits and ch4

## root

## per PFT (herbaceous / shrub)
stats_pft <- df_root %>%
  group_by(trait, ch4, PFT) %>%
  reframe({
    d <- dplyr::filter(cur_data_all(), is.finite(x), is.finite(y))
    if (nrow(d) >= 2 && sd(d$x) > 0 && sd(d$y) > 0) {
      d <- mutate(d, xs = as.numeric(scale(x)), ys = as.numeric(scale(y)))
      m <- lm(ys ~ xs, d)
      s <- summary(m)
      tibble(
        n          = nrow(d),
        slope_std  = unname(coef(m)[2]),             # standardized slope (β)
        r2         = s$r.squared,                    # model R²
        p_slope    = coef(s)[2, 4],                  # p-value for slope
        p_model    = pf(s$fstatistic[1],             # p-value for overall model
                        s$fstatistic[2],
                        s$fstatistic[3],
                        lower.tail = FALSE)
      )
    } else {
      tibble(n = nrow(d),
             slope_std = NA_real_,
             r2 = NA_real_,
             p_slope = NA_real_,
             p_model = NA_real_)
    }
  }) %>%
  ungroup() %>%
  mutate(group = paste0("PFT:", PFT)) %>%
  select(trait, ch4, group, n, slope_std, r2, p_slope, p_model)

stats_pft <- as.data.frame(stats_pft)
stats_pft

stats_pft %>% dplyr::filter(r2 > 0.2, p_model < 0.1)

# this was used to check the models that passed the criteria per PFT and trait
# by changing the PFT and trait each time. group=PFT:shrub and trait=SRL_mg1 are left as an example
stats_pft %>%
  dplyr::filter(
    r2 > 0.2,
    p_model < 0.1,
    group == "PFT:shrub",
    trait == "SRL_mg1"
  )

#### check model assumptions ####

# Build a per-group residuals data frame
diag_residuals_pft <- df_root %>%
  group_by(trait, ch4, PFT) %>%
  reframe({
    d <- dplyr::filter(cur_data_all(), is.finite(x), is.finite(y))
    if (nrow(d) >= 3 && sd(d$x) > 0 && sd(d$y) > 0) {
      d <- mutate(d, xs = as.numeric(scale(x)), ys = as.numeric(scale(y)))
      m <- lm(ys ~ xs, d)
      tibble(
        fitted = fitted(m),
        resid  = resid(m),
        stdres = rstandard(m),
        xs     = d$xs,
        ys     = d$ys
      )
    } else {
      tibble(fitted = numeric(0), resid = numeric(0), stdres = numeric(0), xs = numeric(0), ys = numeric(0))
    }
  }) %>%
  ungroup()

# Residuals vs Fitted
ggplot(diag_residuals_pft, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, span = 0.8) +
  facet_grid(PFT ~ trait + ch4, scales = "free_x") +
  labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
  theme_bw()

# QQ-plot of standardized residuals 
ggplot(diag_residuals_pft, aes(sample = stdres)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(PFT ~ trait + ch4) +
  labs(title = "QQ-plots of standardized residuals") +
  theme_bw()

## rhizome models

## per PFT
stats_pft_rhiz <- df_rhiz %>%
  group_by(trait, ch4, PFT) %>%
  reframe({
    d <- dplyr::filter(cur_data_all(), is.finite(x), is.finite(y))
    if (nrow(d) >= 2 && sd(d$x) > 0 && sd(d$y) > 0) {
      d <- mutate(d, xs = as.numeric(scale(x)), ys = as.numeric(scale(y)))
      m <- lm(ys ~ xs, d)
      s <- summary(m)
      tibble(
        n          = nrow(d),
        slope_std  = unname(coef(m)[2]),             # standardized slope (β)
        r2         = s$r.squared,                    # model R²
        p_slope    = coef(s)[2, 4],                  # p-value for slope
        p_model    = pf(s$fstatistic[1],             # p-value for overall model
                        s$fstatistic[2],
                        s$fstatistic[3],
                        lower.tail = FALSE)
      )
    } else {
      tibble(n = nrow(d),
             slope_std = NA_real_,
             r2 = NA_real_,
             p_slope = NA_real_,
             p_model = NA_real_)
    }
  }) %>%
  ungroup() %>%
  mutate(group = paste0("PFT:", PFT)) %>%
  select(trait, ch4, group, n, slope_std, r2, p_slope, p_model)


stats_pft_rhiz <- as.data.frame(stats_pft_rhiz)
stats_pft_rhiz

# this was used to check the models that passed the criteria per PFT and trait
# by changing the PFT and trait each time. group=PFT:shrub and trait=SRL_mg1 are left as an example
stats_pft_rhiz %>%
  dplyr::filter(
    r2 > 0.2,
    p_model < 0.1,
    group == "PFT:herbaceous",
    trait == "mass_gm2_std"
  )

#### check model diagnostics ####

# Build a per-group residuals data frame
diag_residuals_pft <- df_rhiz %>%
  group_by(trait, ch4, PFT) %>%
  reframe({
    d <- dplyr::filter(cur_data_all(), is.finite(x), is.finite(y))
    if (nrow(d) >= 3 && sd(d$x) > 0 && sd(d$y) > 0) {
      d <- mutate(d, xs = as.numeric(scale(x)), ys = as.numeric(scale(y)))
      m <- lm(ys ~ xs, d)
      tibble(
        fitted = fitted(m),
        resid  = resid(m),
        stdres = rstandard(m),
        xs     = d$xs,
        ys     = d$ys
      )
    } else {
      tibble(fitted = numeric(0), resid = numeric(0), stdres = numeric(0), xs = numeric(0), ys = numeric(0))
    }
  }) %>%
  ungroup()

# Residuals vs Fitted 
ggplot(diag_residuals_pft, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, span = 0.8) +
  facet_grid(PFT ~ trait + ch4, scales = "free_x") +
  labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
  theme_bw()

# QQ-plot of standardized residuals
ggplot(diag_residuals_pft, aes(sample = stdres)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(PFT ~ trait + ch4) +
  labs(title = "QQ-plots of standardized residuals") +
  theme_bw()


### no PFT models and plots

thaw_stage_colors <- c(
  "intact" = "#71C0BB",
  "partly_thawed"   = "#4E6688",
  "fully_thawed"   = "#332D56"
)

# ROOT

# make a long df version and choose only root
root_noPFT <- plant_ch4_noPFT %>%
  filter(organ == "root") %>%
  select(chamber_id, thaw_stage, organ, c("mass_gm2_std", "sa_cm2m2_std", "RTD_gcm3", "SRL_mg1", "avg_d_mm"), 
         c("ch4_season", "ch4_peak", "ch4_middle", "ch4_early")) %>%
  pivot_longer(cols = c("mass_gm2_std", "sa_cm2m2_std", "RTD_gcm3", "SRL_mg1", "avg_d_mm"), 
               names_to = "trait", values_to = "value_trait") %>%
  pivot_longer(cols = c("ch4_season", "ch4_peak", "ch4_middle", "ch4_early"),     
               names_to = "ch4",   values_to = "value_ch4") %>%
  filter(is.finite(value_trait), is.finite(value_ch4))

root_noPFT <- root_noPFT %>%
  mutate(
    trait = factor(
      trait,
      levels = c("mass_gm2_std", "sa_cm2m2_std", "avg_d_mm", "RTD_gcm3", "SRL_mg1")
    ),
    ch4 = factor(
      ch4,
      levels = c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")
    )
  )

# plot
p_root_noPFT <-
  ggplot(root_noPFT, aes(x = value_trait, y = value_ch4, shape = thaw_stage)) +
  geom_point(aes(fill = thaw_stage), size = 3.8, color = "darkgrey") +
  scale_shape_manual(values = shape_vals) +
  scale_fill_manual(values = thaw_stage_colors, name = "Thaw stage") +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              color = "darkgrey", linewidth = 0.9) +
  facet_grid(ch4 ~ trait, scales = "free") +
  scale_x_continuous(labels = label_number(),
                     expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(labels = label_number(),
                     expand = expansion(mult = c(0.1, 0.2))) +
  labs(
    x = "",
    y = expression(paste("CH"[4], " flux (mg CH"[4], " m"^-2, " d"^-1, ")")),
    title = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 16)
  )

p_root_noPFT


#### RHIZOME

# make a long df version and choose only rhizome
rhi_noPFT <- plant_ch4_noPFT %>%
  filter(organ == "rhizome") %>%
  select(chamber_id, thaw_stage, organ, c("mass_gm2_std", "sa_cm2m2_std", "RTD_gcm3", "avg_d_mm"), 
         c("ch4_season", "ch4_peak", "ch4_middle", "ch4_early")) %>%
  pivot_longer(cols = c("mass_gm2_std", "sa_cm2m2_std", "RTD_gcm3", "avg_d_mm"), 
               names_to = "trait", values_to = "value_trait") %>%
  pivot_longer(cols = c("ch4_season", "ch4_peak", "ch4_middle", "ch4_early"),     
               names_to = "ch4",   values_to = "value_ch4") %>%
  filter(is.finite(value_trait), is.finite(value_ch4))

rhi_noPFT <- rhi_noPFT %>%
  mutate(
    trait = factor(
      trait,
      levels = c("mass_gm2_std", "sa_cm2m2_std", "avg_d_mm", "RTD_gcm3")
    ),
    ch4 = factor(
      ch4,
      levels = c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")
    )
  )

# plot
p_rhi_noPFT <-
  ggplot(rhi_noPFT, aes(x = value_trait, y = value_ch4, shape = thaw_stage)) +
  geom_point(aes(fill = thaw_stage), size = 3.8, color = "darkgrey") +
  scale_shape_manual(values = shape_vals) +
  scale_fill_manual(values = thaw_stage_colors, name = "Thaw stage") +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              color = "darkgrey", linewidth = 0.9) +
  facet_grid(ch4 ~ trait, scales = "free") +
  scale_x_continuous(labels = label_number(),
                     expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(labels = label_number(),
                     expand = expansion(mult = c(0.1, 0.2))) +
  labs(
    x = "",
    y = expression(paste("CH"[4], " flux (mg CH"[4], " m"^-2, " d"^-1, ")")),
    title = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 16)
  )


p_rhi_noPFT
# linear regressions based on centered and scaled ch4 and trait values

# ROOT
stats_root_noPFT <- root_noPFT %>%
  group_by(trait, ch4) %>%
  filter(is.finite(value_ch4), is.finite(value_trait)) %>%
  mutate(
    x_scaled = as.numeric(scale(value_trait, center = TRUE, scale = TRUE)),
    y_scaled = as.numeric(scale(value_ch4, center = TRUE, scale = TRUE))
  ) %>%
  do({
    fit <- lm(y_scaled ~ x_scaled, data = .)
    tibble(
      slope = coef(fit)[["x_scaled"]],
      r2 = summary(fit)$r.squared,
      p_model = glance(fit)$p.value,
      p_x = tidy(fit) %>% filter(term == "x_scaled") %>% pull(p.value)
    )
  }) %>%
  ungroup()

stats_root_noPFT <- as.data.frame(stats_root_noPFT)
stats_root_noPFT

stats_root_noPFT %>% dplyr::filter(r2 > 0.2, p_model < 0.1)

# this was used to check the models that passed the criteria per PFT and trait
# by changing the PFT and trait each time. trait=SRL_mg1 is left as an example
stats_root_noPFT %>%
  dplyr::filter(
    r2 > 0.2,
    p_model < 0.1,
    trait == "SRL_mg1"
  )

#### check model diagnostics ####

# Build a per-group residuals data frame
diag_residuals_pft <- root_noPFT %>%
  group_by(trait, ch4) %>%
  reframe({
    d <- dplyr::filter(cur_data_all(), is.finite(value_ch4), is.finite(value_trait))
    if (nrow(d) >= 3 && sd(d$value_ch4) > 0 && sd(d$value_trait) > 0) {
      d <- mutate(d, xs = as.numeric(scale(value_trait)), ys = as.numeric(scale(value_ch4)))
      m <- lm(ys ~ xs, d)
      tibble(
        fitted = fitted(m),
        resid  = resid(m),
        stdres = rstandard(m),
        xs     = d$xs,
        ys     = d$ys
      )
    } else {
      tibble(fitted = numeric(0), resid = numeric(0), stdres = numeric(0), xs = numeric(0), ys = numeric(0))
    }
  }) %>%
  ungroup()

# Residuals vs Fitted
ggplot(diag_residuals_pft, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.6) +
  facet_grid(trait ~ ch4, scales = "free_x") +
  labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
  theme_bw()

# QQ-plot of standardized residuals
ggplot(diag_residuals_pft, aes(sample = stdres)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(trait ~ ch4, scales = "free_x") +
  labs(title = "QQ-plots of standardized residuals") +
  theme_bw()


### updated plot based on model diagnostics:

# make a long df version and choose only root
root_noPFT <- plant_ch4_noPFT %>%
  filter(organ == "root") %>%
  select(chamber_id, thaw_stage, organ, c("RTD_gcm3", "SRL_mg1"), 
         c("ch4_season", "ch4_peak", "ch4_middle", "ch4_early")) %>%
  pivot_longer(cols = c("RTD_gcm3", "SRL_mg1"), 
               names_to = "trait", values_to = "value_trait") %>%
  pivot_longer(cols = c("ch4_season", "ch4_peak", "ch4_middle", "ch4_early"),     
               names_to = "ch4",   values_to = "value_ch4") %>%
  filter(is.finite(value_trait), is.finite(value_ch4))

root_noPFT <- root_noPFT %>%
  mutate(
    trait = factor(
      trait,
      levels = c("RTD_gcm3", "SRL_mg1")
    ),
    ch4 = factor(
      ch4,
      levels = c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")
    )
  )

# plot
p_root_noPFT_RTD_SRL <-
  ggplot(root_noPFT, aes(x = value_trait, y = value_ch4, shape = thaw_stage)) +
  geom_point(aes(fill = thaw_stage), size = 3.8, color = "black", alpha=0.9) +
  scale_shape_manual(values = shape_vals) +
  scale_fill_manual(values = thaw_stage_colors, name = "Thaw stage") +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              color = "grey30", linewidth = 0.9) +
  facet_grid(ch4 ~ trait, scales = "free") +
  scale_x_continuous(labels = label_number(),
                     expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(labels = label_number(),
                     expand = expansion(mult = c(0.1, 0.2))) +
  labs(
    x = "",
    y = expression(paste("CH"[4], " flux (mg CH"[4], " m"^-2, " d"^-1, ")")),
    title = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 16)
  )

p_root_noPFT_RTD_SRL

### now contains only RTD and SRL


### version for SI with the rest of the traits

# make a long df version and choose only root
root_noPFT <- plant_ch4_noPFT %>%
  filter(organ == "root") %>%
  select(chamber_id, thaw_stage, organ, c("mass_gm2_std", "sa_cm2m2_std", "avg_d_mm", "SRL_mg1"), 
         c("ch4_season", "ch4_peak", "ch4_middle", "ch4_early")) %>%
  pivot_longer(cols = c("mass_gm2_std", "sa_cm2m2_std", "avg_d_mm", "SRL_mg1"), 
               names_to = "trait", values_to = "value_trait") %>%
  pivot_longer(cols = c("ch4_season", "ch4_peak", "ch4_middle", "ch4_early"),     
               names_to = "ch4",   values_to = "value_ch4") %>%
  filter(is.finite(value_trait), is.finite(value_ch4))

root_noPFT <- root_noPFT %>%
  mutate(
    trait = factor(
      trait,
      levels = c("mass_gm2_std", "sa_cm2m2_std", "avg_d_mm", "SRL_mg1")
    ),
    ch4 = factor(
      ch4,
      levels = c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")
    )
  )

# plot
p_root_noPFT_rest <-
  ggplot(root_noPFT, aes(x = value_trait, y = value_ch4, shape = thaw_stage)) +
  geom_point(aes(fill = thaw_stage), size = 3.8, color = "black", alpha = 0.9) +
  scale_shape_manual(values = shape_vals) +
  scale_fill_manual(values = thaw_stage_colors, name = "Thaw stage") +
  facet_grid(ch4 ~ trait, scales = "free") +
  scale_x_continuous(labels = label_number(),
                     expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(labels = label_number(),
                     expand = expansion(mult = c(0.1, 0.2))) +
  labs(
    x = "",
    y = expression(paste("CH"[4], " flux (mg CH"[4], " m"^-2, " d"^-1, ")")),
    title = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 16)
  )

p_root_noPFT_rest


# models for rhizome noPFT

stats_rhi_noPFT <- rhi_noPFT %>%
  group_by(trait, ch4) %>%
  filter(is.finite(value_ch4), is.finite(value_trait)) %>%
  mutate(
    x_scaled = as.numeric(scale(value_ch4, center = TRUE, scale = TRUE)),
    y_scaled = as.numeric(scale(value_trait, center = TRUE, scale = TRUE))
  ) %>%
  do({
    fit <- lm(y_scaled ~ x_scaled, data = .)
    tibble(
      slope = coef(fit)[["x_scaled"]],
      r2 = summary(fit)$r.squared,
      p_model = glance(fit)$p.value,
      p_x = tidy(fit) %>% filter(term == "x_scaled") %>% pull(p.value)
    )
  }) %>%
  ungroup()

stats_rhi_noPFT <- as.data.frame(stats_rhi_noPFT)
stats_rhi_noPFT

stats_rhi_noPFT %>% dplyr::filter(r2 > 0.2, p_model < 0.1)

# this was used to check the models that passed the criteria per PFT and trait
# by changing the PFT and trait each time. trait=mass_gm2_std is left as an example
stats_rhi_noPFT %>%
  dplyr::filter(
    r2 > 0.2,
    p_model < 0.1,
    trait == "mass_gm2_std"
  )

#### check model assumptions ####

# Build a per-group residuals data frame
diag_residuals_pft <- rhi_noPFT %>%
  group_by(trait, ch4) %>%
  reframe({
    d <- dplyr::filter(cur_data_all(), is.finite(value_ch4), is.finite(value_trait))
    if (nrow(d) >= 3 && sd(d$value_ch4) > 0 && sd(d$value_trait) > 0) {
      d <- mutate(d, xs = as.numeric(scale(value_trait)), ys = as.numeric(scale(value_ch4)))
      m <- lm(ys ~ xs, d)
      tibble(
        fitted = fitted(m),
        resid  = resid(m),
        stdres = rstandard(m),
        xs     = d$xs,
        ys     = d$ys
      )
    } else {
      tibble(fitted = numeric(0), resid = numeric(0), stdres = numeric(0), xs = numeric(0), ys = numeric(0))
    }
  }) %>%
  ungroup()

# Residuals vs Fitted
ggplot(diag_residuals_pft, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.6) +
  facet_grid(trait ~ ch4, scales = "free_x") +
  labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
  theme_bw()

# QQ-plot of standardized residuals
ggplot(diag_residuals_pft, aes(sample = stdres)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(trait ~ ch4, scales = "free_x") +
  labs(title = "QQ-plots of standardized residuals") +
  theme_bw()


##### create a combined trait-ch4 flux grid with both PFT and noPFT trends
##### one figure for roots and one figure for rhizomes
##### show only the chosen traits (based on model diagnostics, R2s and slopes)

### ROOT

## traits to keep + desired column order
traits_keep <- c("avg_d_mm", "RTD_gcm3", "SRL_mg1")
ch4_keep    <- c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")

## PFT-grouped data
root_PFT <- plant_ch4 %>%
  filter(organ == "root", PFT %in% c("herbaceous", "shrub")) %>%
  select(chamber_id, thaw_stage, PFT, organ,
         all_of(traits_keep), all_of(ch4_keep)) %>%
  pivot_longer(
    cols = all_of(traits_keep),
    names_to = "trait",
    values_to = "value_trait"
  ) %>%
  pivot_longer(
    cols = all_of(ch4_keep),
    names_to = "ch4",
    values_to = "value_ch4"
  ) %>%
  filter(is.finite(value_trait), is.finite(value_ch4)) %>%
  mutate(
    grouping = "PFT",
    trait_panel = paste0(trait, "_PFT")
  )

## no-PFT data
root_noPFT <- plant_ch4_noPFT %>%
  filter(organ == "root") %>%
  select(chamber_id, thaw_stage, organ,
         all_of(traits_keep), all_of(ch4_keep)) %>%
  pivot_longer(
    cols = all_of(traits_keep),
    names_to = "trait",
    values_to = "value_trait"
  ) %>%
  pivot_longer(
    cols = all_of(ch4_keep),
    names_to = "ch4",
    values_to = "value_ch4"
  ) %>%
  filter(is.finite(value_trait), is.finite(value_ch4)) %>%
  mutate(
    grouping = "noPFT",
    PFT = NA_character_,
    trait_panel = paste0(trait, "_noPFT")
  )

## combine
root_combined <- bind_rows(root_PFT, root_noPFT) %>%
  mutate(
    trait_panel = factor(
      trait_panel,
      levels = c(
        "avg_d_mm_PFT", "avg_d_mm_noPFT",
        "RTD_gcm3_PFT", "RTD_gcm3_noPFT",
        "SRL_mg1_PFT", "SRL_mg1_noPFT"
      ),
      labels = c(
        "avg_d_mm\nPFT",
        "avg_d_mm\nnoPFT",
        "RTD_gcm3\nPFT",
        "RTD_gcm3\nnoPFT",
        "SRL_mg1\nPFT",
        "SRL_mg1\nnoPFT"
      )
    ),
    ch4 = factor(
      ch4,
      levels = c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")
    )
  )

## plot the combined grid (Fig. 3- improved further in Inkscape)
p_root_combined <-
  ggplot(root_combined, aes(x = value_trait, y = value_ch4, shape = thaw_stage)) +
  
  ## PFT-grouped points and lines
  geom_point(
    data = root_combined %>% filter(grouping == "PFT"),
    aes(fill = PFT),
    color = "black",
    size = 3.5,
    alpha = 0.9,
    stroke = 0.5
  ) +
  geom_smooth(
    data = root_combined %>% filter(grouping == "PFT"),
    aes(color = PFT, group = PFT),
    method = "lm",
    se = FALSE,
    linewidth = 0.5
  ) +
  
  ## noPFT points and lines
  geom_point(
    data = root_combined %>% filter(grouping == "noPFT"),
    aes(fill = thaw_stage),
    color = "black",
    size = 3.8,
    alpha = 0.9,
    stroke = 0.5
  ) +
  geom_smooth(
    data = root_combined %>% filter(grouping == "noPFT"),
    aes(group = 1),
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  
  scale_shape_manual(values = shape_vals) +
  scale_fill_manual(
    values = c(pft_colors, thaw_stage_colors),
    name = ""
  ) +
  scale_color_manual(values = pft_colors, guide = "none") +
  
  facet_grid(ch4 ~ trait_panel, scales = "free") +
  
  scale_x_continuous(
    labels = label_number(),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  scale_y_continuous(
    labels = label_number(),
    expand = expansion(mult = c(0.1, 0.2))
  ) +
  
  labs(
    x = "",
    y = expression(paste("CH"[4], " flux (mg CH"[4], " m"^-2, " d"^-1, ")"))
  ) +
  
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 16)
  )

p_root_combined


# mass and SA for SI (Fig. B6)

## traits to keep (note renamed SA variable)
traits_keep2 <- c("mass_gm2_std", "sa_m2m2_std")
ch4_keep     <- c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")

## PFT data
root_PFT_2 <- plant_ch4 %>%
  filter(organ == "root", PFT %in% c("herbaceous", "shrub")) %>%
  select(chamber_id, thaw_stage, PFT, organ,
         mass_gm2_std, sa_cm2m2_std, all_of(ch4_keep)) %>%
  mutate(
    sa_m2m2_std = sa_cm2m2_std / 10000 
  ) %>%
  select(-sa_cm2m2_std) %>%
  pivot_longer(
    cols = all_of(traits_keep2),
    names_to = "trait",
    values_to = "value_trait"
  ) %>%
  pivot_longer(
    cols = all_of(ch4_keep),
    names_to = "ch4",
    values_to = "value_ch4"
  ) %>%
  filter(is.finite(value_trait), is.finite(value_ch4)) %>%
  mutate(
    grouping = "PFT",
    trait_panel = paste0(trait, "_PFT")
  )

## noPFT data
root_noPFT_2 <- plant_ch4_noPFT %>%
  filter(organ == "root") %>%
  select(chamber_id, thaw_stage, organ,
         mass_gm2_std, sa_cm2m2_std, all_of(ch4_keep)) %>%
  mutate(
    sa_m2m2_std = sa_cm2m2_std / 10000 
  ) %>%
  select(-sa_cm2m2_std) %>%
  pivot_longer(
    cols = all_of(traits_keep2),
    names_to = "trait",
    values_to = "value_trait"
  ) %>%
  pivot_longer(
    cols = all_of(ch4_keep),
    names_to = "ch4",
    values_to = "value_ch4"
  ) %>%
  filter(is.finite(value_trait), is.finite(value_ch4)) %>%
  mutate(
    grouping = "noPFT",
    PFT = NA_character_,
    trait_panel = paste0(trait, "_noPFT")
  )

## combine 
root_combined_2 <- bind_rows(root_PFT_2, root_noPFT_2) %>%
  mutate(
    trait_panel = factor(
      trait_panel,
      levels = c(
        "mass_gm2_std_PFT", "mass_gm2_std_noPFT",
        "sa_m2m2_std_PFT",  "sa_m2m2_std_noPFT"
      )
    ),
    ch4 = factor(
      ch4,
      levels = c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")
    )
  )

## plot (Fig. B6)
p_root_mass_sa <-
  ggplot(root_combined_2, aes(x = value_trait, y = value_ch4, shape = thaw_stage)) +
  
  ## PFT points
  geom_point(
    data = root_combined_2 %>% filter(grouping == "PFT"),
    aes(fill = PFT),
    color = "black",
    size = 3.5,
    alpha = 0.9,
    stroke = 0.5
  ) +
  
  ## noPFT points
  geom_point(
    data = root_combined_2 %>% filter(grouping == "noPFT"),
    aes(fill = thaw_stage),
    color = "black",
    size = 3.8,
    alpha = 0.9,
    stroke = 0.5
  ) +
  
  scale_shape_manual(values = shape_vals) +
  scale_fill_manual(values = c(pft_colors, thaw_stage_colors), name = "") +
  
  facet_grid(
    ch4 ~ trait_panel,
    scales = "free",
    labeller = labeller(trait_panel = c(
      "mass_gm2_std_PFT"   = "mass (g m^-2)\nPFT",
      "mass_gm2_std_noPFT" = "mass (g m^-2)\nnoPFT",
      "sa_m2m2_std_PFT"    = expression(SA~(m^2~m^-2)*"\nPFT"),
      "sa_m2m2_std_noPFT"  = expression(SA~(m^2~m^-2)*"\nnoPFT")
    ))
  ) +
  
  scale_x_continuous(
    labels = label_number(),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  scale_y_continuous(
    labels = label_number(),
    expand = expansion(mult = c(0.1, 0.2))
  ) +
  
  labs(
    x = "",
    y = expression(paste("CH"[4], " flux (mg CH"[4], " m"^-2, " d"^-1, ")")),
    title = "Root"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 16)
  )

p_root_mass_sa


### RHIZOME

traits_keep_rhiz <- c("sa_m2m2_std", "avg_d_mm")
ch4_keep         <- c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")

## PFT data 
rhiz_PFT <- plant_ch4 %>%
  filter(organ == "rhizome", PFT %in% c("herbaceous", "shrub")) %>%
  select(
    chamber_id, thaw_stage, PFT, organ,
    sa_cm2m2_std, avg_d_mm,
    all_of(ch4_keep)
  ) %>%
  mutate(
    sa_m2m2_std = sa_cm2m2_std / 10000
  ) %>%
  select(-sa_cm2m2_std) %>%
  pivot_longer(
    cols = all_of(traits_keep_rhiz),
    names_to = "trait",
    values_to = "value_trait"
  ) %>%
  pivot_longer(
    cols = all_of(ch4_keep),
    names_to = "ch4",
    values_to = "value_ch4"
  ) %>%
  filter(is.finite(value_trait), is.finite(value_ch4)) %>%
  mutate(
    grouping = "PFT",
    PFT = factor(PFT, levels = c("shrub", "herbaceous")),
    trait_panel = paste0(trait, "_PFT")
  )

## noPFT data 
rhiz_noPFT <- plant_ch4_noPFT %>%
  filter(organ == "rhizome") %>%
  select(
    chamber_id, thaw_stage, organ,
    sa_cm2m2_std, avg_d_mm,
    all_of(ch4_keep)
  ) %>%
  mutate(
    sa_m2m2_std = sa_cm2m2_std / 10000
  ) %>%
  select(-sa_cm2m2_std) %>%
  pivot_longer(
    cols = all_of(traits_keep_rhiz),
    names_to = "trait",
    values_to = "value_trait"
  ) %>%
  pivot_longer(
    cols = all_of(ch4_keep),
    names_to = "ch4",
    values_to = "value_ch4"
  ) %>%
  filter(is.finite(value_trait), is.finite(value_ch4)) %>%
  mutate(
    grouping = "noPFT",
    PFT = NA,
    trait_panel = paste0(trait, "_noPFT")
  )

## combine 
rhiz_combined <- bind_rows(rhiz_PFT, rhiz_noPFT) %>%
  mutate(
    trait_panel = factor(
      trait_panel,
      levels = c(
        "sa_m2m2_std_PFT",  "sa_m2m2_std_noPFT",
        "avg_d_mm_PFT",     "avg_d_mm_noPFT"
      )
    ),
    ch4 = factor(
      ch4,
      levels = c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")
    )
  )

## plot (Fig. 4)
p_rhiz_combined <-
  ggplot(rhiz_combined, aes(x = value_trait, y = value_ch4, shape = thaw_stage)) +
  
  ## PFT points
  geom_point(
    data = rhiz_combined %>% filter(grouping == "PFT"),
    aes(fill = PFT),
    color = "black",
    size = 3.5,
    alpha = 0.9,
    stroke = 0.5
  ) +
  
  ## PFT regression lines, except avg_d_mm
  geom_smooth(
    data = rhiz_combined %>%
      filter(grouping == "PFT", trait != "avg_d_mm"),
    aes(color = PFT, group = PFT),
    method = "lm",
    se = FALSE,
    linewidth = 0.7
  ) +
  
  ## noPFT points
  geom_point(
    data = rhiz_combined %>% filter(grouping == "noPFT"),
    aes(fill = thaw_stage),
    color = "black",
    size = 3.8,
    alpha = 0.9,
    stroke = 0.5
  ) +
  
  ## noPFT regression lines, except avg_d_mm
  geom_smooth(
    data = rhiz_combined %>%
      filter(grouping == "noPFT", trait != "avg_d_mm"),
    aes(group = 1),
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  
  scale_shape_manual(values = shape_vals) +
  scale_fill_manual(
    values = c(pft_colors, thaw_stage_colors),
    name = ""
  ) +
  scale_color_manual(values = pft_colors, guide = "none") +
  
  facet_grid(
    ch4 ~ trait_panel,
    scales = "free",
    labeller = labeller(
      trait_panel = as_labeller(c(
        "sa_m2m2_std_PFT"   = "SA~(m^2~m^{-2})*'\nPFT'",
        "sa_m2m2_std_noPFT" = "SA~(m^2~m^{-2})*'\nnoPFT'",
        "avg_d_mm_PFT"      = "avg~diameter~(mm)*'\nPFT'",
        "avg_d_mm_noPFT"    = "avg~diameter~(mm)*'\nnoPFT'"
      ), label_parsed)
    )
  ) +
  
  scale_x_continuous(
    labels = label_number(),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  scale_y_continuous(
    labels = label_number(),
    expand = expansion(mult = c(0.1, 0.2))
  ) +
  
  labs(
    x = "",
    y = expression(paste("CH"[4], " flux (mg CH"[4], " m"^-2, " d"^-1, ")")),
    title = "Rhizome"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16)
  )

p_rhiz_combined


# only mass and RTD for SI (Fig. B7)

## traits to keep
traits_keep_rhiz2 <- c("mass_gm2_std", "RTD_gcm3")
ch4_keep          <- c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")

## PFT data
rhiz_PFT2 <- plant_ch4 %>%
  filter(organ == "rhizome", PFT %in% c("herbaceous", "shrub")) %>%
  select(
    chamber_id, thaw_stage, PFT, organ,
    all_of(traits_keep_rhiz2),
    all_of(ch4_keep)
  ) %>%
  pivot_longer(
    cols = all_of(traits_keep_rhiz2),
    names_to = "trait",
    values_to = "value_trait"
  ) %>%
  pivot_longer(
    cols = all_of(ch4_keep),
    names_to = "ch4",
    values_to = "value_ch4"
  ) %>%
  filter(is.finite(value_trait), is.finite(value_ch4)) %>%
  mutate(
    grouping = "PFT",
    PFT = factor(PFT, levels = c("shrub", "herbaceous")),
    trait_panel = paste0(trait, "_PFT")
  )

## noPFT data
rhiz_noPFT2 <- plant_ch4_noPFT %>%
  filter(organ == "rhizome") %>%
  select(
    chamber_id, thaw_stage, organ,
    all_of(traits_keep_rhiz2),
    all_of(ch4_keep)
  ) %>%
  pivot_longer(
    cols = all_of(traits_keep_rhiz2),
    names_to = "trait",
    values_to = "value_trait"
  ) %>%
  pivot_longer(
    cols = all_of(ch4_keep),
    names_to = "ch4",
    values_to = "value_ch4"
  ) %>%
  filter(is.finite(value_trait), is.finite(value_ch4)) %>%
  mutate(
    grouping = "noPFT",
    PFT = NA_character_,
    trait_panel = paste0(trait, "_noPFT")
  )

## combine
rhiz_combined2 <- bind_rows(rhiz_PFT2, rhiz_noPFT2) %>%
  mutate(
    trait_panel = factor(
      trait_panel,
      levels = c(
        "mass_gm2_std_PFT", "mass_gm2_std_noPFT",
        "RTD_gcm3_PFT",     "RTD_gcm3_noPFT"
      )
    ),
    ch4 = factor(
      ch4,
      levels = c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")
    )
  )

## plot (Fig. B7)
p_rhiz_mass_rtd <-
  ggplot(rhiz_combined2, aes(x = value_trait, y = value_ch4, shape = thaw_stage)) +
  
  ## PFT points
  geom_point(
    data = rhiz_combined2 %>% filter(grouping == "PFT"),
    aes(fill = PFT),
    color = "black",
    size = 3.5,
    alpha = 0.9,
    stroke = 0.5
  ) +
  
  ## PFT regression lines, except RTD_gcm3
  geom_smooth(
    data = rhiz_combined2 %>%
      filter(grouping == "PFT", trait != "RTD_gcm3"),
    aes(color = PFT, group = PFT),
    method = "lm",
    se = FALSE,
    linewidth = 0.7
  ) +
  
  ## noPFT points
  geom_point(
    data = rhiz_combined2 %>% filter(grouping == "noPFT"),
    aes(fill = thaw_stage),
    color = "black",
    size = 3.8,
    alpha = 0.9,
    stroke = 0.5
  ) +
  
  ## noPFT regression lines, except RTD_gcm3
  geom_smooth(
    data = rhiz_combined2 %>%
      filter(grouping == "noPFT", trait != "RTD_gcm3"),
    aes(group = 1),
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  
  scale_shape_manual(values = shape_vals) +
  scale_fill_manual(
    values = c(pft_colors, thaw_stage_colors),
    name = ""
  ) +
  scale_color_manual(values = pft_colors, guide = "none") +
  
  facet_grid(
    ch4 ~ trait_panel,
    scales = "free",
    labeller = labeller(
      trait_panel = as_labeller(c(
        "mass_gm2_std_PFT"   = "mass~(g~m^{-2})*'\nPFT'",
        "mass_gm2_std_noPFT" = "mass~(g~m^{-2})*'\nnoPFT'",
        "RTD_gcm3_PFT"       = "RTD~(g~cm^{-3})*'\nPFT'",
        "RTD_gcm3_noPFT"     = "RTD~(g~cm^{-3})*'\nnoPFT'"
      ), label_parsed)
    )
  ) +
  
  scale_x_continuous(
    labels = label_number(),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  scale_y_continuous(
    labels = label_number(),
    expand = expansion(mult = c(0.1, 0.2))
  ) +
  
  labs(
    x = "",
    y = expression(paste("CH"[4], " flux (mg CH"[4], " m"^-2, " d"^-1, ")")),
    title = "Rhizome"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16)
  )

p_rhiz_mass_rtd


## check Kendall correlations for relationships where linear regressions are not valid
# the following code was used for each organ and trait at a time by manually changing the PFT, trait and organ name
# PFT = shrub, organ = rhizome, trait = RTD_gcm3 are left here as an example

ch4_cols <- c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")

plant_long <- plant_ch4 %>%
  filter(PFT == "shrub", organ == "rhizome") %>%
  select(chamber_id, value_trait = RTD_gcm3, all_of(ch4_cols)) %>%
  pivot_longer(
    cols = all_of(ch4_cols),
    names_to = "ch4",
    values_to = "value_ch4"
  ) %>%
  filter(is.finite(value_trait), is.finite(value_ch4))

kendall_by_ch4 <- plant_long %>%
  group_by(ch4) %>%
  nest() %>%
  mutate(
    n_pairs = map_int(data, ~ nrow(.x)),
    test = map(data, ~ {
      if (nrow(.x) >= 2) {
        cor.test(
          .x$value_trait,
          .x$value_ch4,
          method = "kendall",
          exact = FALSE
        )
      } else {
        NULL
      }
    }),
    tau = map_dbl(test, ~ {
      if (is.null(.x)) NA_real_ else unname(.x$estimate)
    }),
    p_value = map_dbl(test, ~ {
      if (is.null(.x)) NA_real_ else .x$p.value
    })
  ) %>%
  select(ch4, n_pairs, tau, p_value) %>%
  arrange(ch4)

kendall_by_ch4

# no PFTs

ch4_cols <- c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")

plant_long <- plant_ch4_noPFT %>%
  filter(organ == "rhizome") %>%
  select(chamber_id, value_trait = RTD_gcm3, all_of(ch4_cols)) %>%
  pivot_longer(
    cols = all_of(ch4_cols),
    names_to = "ch4",
    values_to = "value_ch4"
  ) %>%
  filter(is.finite(value_trait), is.finite(value_ch4))

kendall_by_ch4 <- plant_long %>%
  group_by(ch4) %>%
  nest() %>%
  mutate(
    n_pairs = map_int(data, ~ nrow(.x)),
    test = map(data, ~ {
      if (nrow(.x) >= 2) {
        cor.test(
          .x$value_trait,
          .x$value_ch4,
          method = "kendall",
          exact = FALSE
        )
      } else {
        NULL
      }
    }),
    tau = map_dbl(test, ~ {
      if (is.null(.x)) NA_real_ else unname(.x$estimate)
    }),
    p_value = map_dbl(test, ~ {
      if (is.null(.x)) NA_real_ else .x$p.value
    })
  ) %>%
  select(ch4, n_pairs, tau, p_value) %>%
  arrange(ch4)

kendall_by_ch4

###### changes in traits along the thaw gradient ######

#### calculate % increase or decrease from intact to fully_thawed 

trait_vars <- c("mass_gm2_std", 
                "sa_cm2m2_std",   
                "avg_d_mm",            
                "RTD_gcm3",            
                "SRL_mg1")            

# group means per trait per thaw stage
stage_means <-
  plant_ch4 %>%
  filter(thaw_stage %in% c("intact","partly_thawed","fully_thawed")) %>%
  mutate(thaw_stage = factor(thaw_stage, levels = c("intact","partly_thawed","fully_thawed"))) %>%
  group_by(PFT, organ, thaw_stage) %>%
  summarise(across(all_of(trait_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# wide df version
wide_means <-
  stage_means %>%
  pivot_longer(all_of(trait_vars), names_to = "trait", values_to = "mean_value") %>%
  pivot_wider(names_from = thaw_stage, values_from = mean_value)

# percentage change per PFT per organ per trait
df <- wide_means %>%
  mutate(
    start_val = intact,
    end_val   = if_else(PFT == "shrub", partly_thawed, fully_thawed),
    abs_change = end_val - start_val,
    pct_change = 100 * (end_val - start_val) / start_val
  ) %>%
  select(PFT, organ, trait, intact, partly_thawed, fully_thawed, start_val, end_val, abs_change, pct_change) %>%
  arrange(trait, PFT, organ) %>%
  mutate(pct_change = ifelse(is.finite(pct_change), pct_change, NA_real_))

df <- df %>%  
  mutate(
    fold_change = end_val / start_val,
    log_pct_change = 100 * (log(end_val) - log(start_val)),
    arc_pct_change = 100 * (end_val - start_val) / ((end_val + start_val)/2)
  )

df <- as.data.frame(df)

### across PFTs

trait_vars <- c("mass_gm2_std", "sa_m2m2_std", "avg_d_mm", "RTD_gcm3", "SRL_mg1")

plant_ch4$sa_m2m2_std <- plant_ch4$sa_cm2m2_std / 10000

plant_ch4_2 <- plant_ch4 %>%
  filter(thaw_stage %in% c("intact","partly_thawed","fully_thawed")) %>%
  mutate(thaw_stage = factor(thaw_stage, levels = c("intact","partly_thawed","fully_thawed")))

# Aggregate across PFTs per chamber: sums (mass, sa) + weighted means (avg_d, RTD, SRL)
across_by_chamber <- plant_ch4_2 %>%
  group_by(organ, thaw_stage, chamber_id) %>%
  summarise(
    # weighted means first
    avg_d_mm = if (sum(mass_gm2_std, na.rm = TRUE) > 0)
      weighted.mean(avg_d_mm, w = mass_gm2_std, na.rm = TRUE) else NA_real_,
    RTD_gcm3 = if (sum(mass_gm2_std, na.rm = TRUE) > 0)
      weighted.mean(RTD_gcm3, w = mass_gm2_std, na.rm = TRUE) else NA_real_,
    SRL_mg1  = if (sum(mass_gm2_std, na.rm = TRUE) > 0)
      weighted.mean(SRL_mg1,  w = mass_gm2_std, na.rm = TRUE) else NA_real_,
    # sums
    mass_gm2_std = sum(mass_gm2_std, na.rm = TRUE),
    sa_m2m2_std = sum(sa_m2m2_std, na.rm = TRUE),
    .groups = "drop"
  )

# Mean per thaw_stage × organ
stage_means_across <- across_by_chamber %>%
  group_by(organ, thaw_stage) %>%
  summarise(across(all_of(trait_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Wide table with intact/partly_thawed/fully_thawed columns
wide_means_across <- stage_means_across %>%
  pivot_longer(all_of(trait_vars), names_to = "trait", values_to = "mean_value") %>%
  pivot_wider(names_from = thaw_stage, values_from = mean_value)

# % change intact to fully_thawed (across PFTs). If fully_thawed is NA or start is 0, return NA
df_across <- wide_means_across %>%
  mutate(
    start_val   = intact,
    end_val     = fully_thawed,
    abs_change  = end_val - start_val,
    pct_change  = 100 * (end_val - start_val) / start_val,
    fold_change = end_val / start_val,
    log_pct_change = 100 * (log(end_val) - log(start_val)),
    arc_pct_change = 100 * (end_val - start_val) / ((end_val + start_val) / 2)
  ) %>%
  mutate(
    pct_change     = ifelse(is.finite(pct_change), pct_change, NA_real_),
    fold_change    = ifelse(is.finite(fold_change), fold_change, NA_real_),
    log_pct_change = ifelse(is.finite(log_pct_change), log_pct_change, NA_real_),
    arc_pct_change = ifelse(is.finite(arc_pct_change), arc_pct_change, NA_real_)
  ) %>%
  select(organ, trait, intact, partly_thawed, fully_thawed,
         start_val, end_val, abs_change, pct_change,
         fold_change, log_pct_change, arc_pct_change) %>%
  arrange(trait, organ) %>%
  as.data.frame()

df_across


# root biomass
partly_thawed_vs_intact_fully_thawed_by_PFT_biomass <- 
  plant_ch4 %>%
  filter(organ == "root", thaw_stage %in% c("intact","partly_thawed","fully_thawed")) %>%
  mutate(thaw_stage = factor(thaw_stage, levels = c("intact","partly_thawed","fully_thawed"))) %>%
  group_by(PFT, thaw_stage) %>%
  summarise(root_biomass_mean_gm2 = mean(mass_gm2_std, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = thaw_stage, values_from = root_biomass_mean_gm2) %>%
  mutate(
    # partly_thawed vs intact
    abs_change_partly_thawed_vs_intact = partly_thawed - intact,
    pct_change_partly_thawed_vs_intact = 100 * (partly_thawed - intact) / intact,
    fold_partly_thawed_over_intact     = partly_thawed / intact,
    
    # partly_thawed vs fully_thawed (will be NA for PFTs absent in fully_thawed, i.e., shrubs)
    abs_change_partly_thawed_vs_fully_thawed = partly_thawed - fully_thawed,
    pct_change_partly_thawed_vs_fully_thawed = 100 * (partly_thawed - fully_thawed) / fully_thawed,
    fold_partly_thawed_over_fully_thawed     = partly_thawed / fully_thawed
  ) %>%
  mutate(
    pct_change_partly_thawed_vs_intact = ifelse(is.finite(pct_change_partly_thawed_vs_intact), pct_change_partly_thawed_vs_intact, NA_real_),
    pct_change_partly_thawed_vs_fully_thawed   = ifelse(is.finite(pct_change_partly_thawed_vs_fully_thawed),   pct_change_partly_thawed_vs_fully_thawed,   NA_real_),
    fold_partly_thawed_over_intact     = ifelse(is.finite(fold_partly_thawed_over_intact),     fold_partly_thawed_over_intact,     NA_real_),
    fold_partly_thawed_over_fully_thawed       = ifelse(is.finite(fold_partly_thawed_over_fully_thawed),       fold_partly_thawed_over_fully_thawed,       NA_real_)
  )

partly_thawed_vs_intact_fully_thawed_by_PFT_biomass <- as.data.frame(partly_thawed_vs_intact_fully_thawed_by_PFT_biomass)


# root SA
partly_thawed_vs_intact_fully_thawed_by_PFT_sa <- 
  plant_ch4 %>%
  filter(organ == "root", thaw_stage %in% c("intact","partly_thawed","fully_thawed")) %>%
  mutate(thaw_stage = factor(thaw_stage, levels = c("intact","partly_thawed","fully_thawed"))) %>%
  group_by(PFT, thaw_stage) %>%
  summarise(root_sa_mean = mean(sa_cm2m2_std, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = thaw_stage, values_from = root_sa_mean) %>%
  mutate(
    # partly_thawed vs intact
    abs_change_partly_thawed_vs_intact = partly_thawed - intact,
    pct_change_partly_thawed_vs_intact = 100 * (partly_thawed - intact) / intact,
    fold_partly_thawed_over_intact     = partly_thawed / intact,
    
    # partly_thawed vs fully_thawed (will be NA for PFTs absent in fully_thawed, e.g., shrubs)
    abs_change_partly_thawed_vs_fully_thawed = partly_thawed - fully_thawed,
    pct_change_partly_thawed_vs_fully_thawed = 100 * (partly_thawed - fully_thawed) / fully_thawed,
    fold_partly_thawed_over_fully_thawed     = partly_thawed / fully_thawed
  ) %>%
  mutate(
    pct_change_partly_thawed_vs_intact = ifelse(is.finite(pct_change_partly_thawed_vs_intact), pct_change_partly_thawed_vs_intact, NA_real_),
    pct_change_partly_thawed_vs_fully_thawed   = ifelse(is.finite(pct_change_partly_thawed_vs_fully_thawed),   pct_change_partly_thawed_vs_fully_thawed,   NA_real_),
    fold_partly_thawed_over_intact     = ifelse(is.finite(fold_partly_thawed_over_intact),     fold_partly_thawed_over_intact,     NA_real_),
    fold_partly_thawed_over_fully_thawed       = ifelse(is.finite(fold_partly_thawed_over_fully_thawed),       fold_partly_thawed_over_fully_thawed,       NA_real_)
  )

partly_thawed_vs_intact_fully_thawed_by_PFT_sa <- as.data.frame(partly_thawed_vs_intact_fully_thawed_by_PFT_sa)

# SRL
partly_thawed_vs_intact_fully_thawed_by_PFT_srl <- 
  plant_ch4 %>%
  filter(organ == "root", thaw_stage %in% c("intact","partly_thawed","fully_thawed")) %>%
  mutate(thaw_stage = factor(thaw_stage, levels = c("intact","partly_thawed","fully_thawed"))) %>%
  group_by(PFT, thaw_stage) %>%
  summarise(root_srl_mean = mean(SRL_mg1, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = thaw_stage, values_from = root_srl_mean) %>%
  mutate(
    # partly_thawed vs intact
    abs_change_partly_thawed_vs_intact = partly_thawed - intact,
    pct_change_partly_thawed_vs_intact = 100 * (partly_thawed - intact) / intact,
    fold_partly_thawed_over_intact     = partly_thawed / intact,
    # partly_thawed vs fully_thawed (will be NA for PFTs absent in fully_thawed, i.e., shrubs)
    abs_change_partly_thawed_vs_fully_thawed = partly_thawed - fully_thawed,
    pct_change_partly_thawed_vs_fully_thawed = 100 * (partly_thawed - fully_thawed) / fully_thawed,
    fold_partly_thawed_over_fully_thawed     = partly_thawed / fully_thawed
  ) %>%
  mutate(
    pct_change_partly_thawed_vs_intact = ifelse(is.finite(pct_change_partly_thawed_vs_intact), pct_change_partly_thawed_vs_intact, NA_real_),
    pct_change_partly_thawed_vs_fully_thawed   = ifelse(is.finite(pct_change_partly_thawed_vs_fully_thawed),   pct_change_partly_thawed_vs_fully_thawed,   NA_real_),
    fold_partly_thawed_over_intact     = ifelse(is.finite(fold_partly_thawed_over_intact),     fold_partly_thawed_over_intact,     NA_real_),
    fold_partly_thawed_over_fully_thawed       = ifelse(is.finite(fold_partly_thawed_over_fully_thawed),       fold_partly_thawed_over_fully_thawed,       NA_real_)
  )

partly_thawed_vs_intact_fully_thawed_by_PFT_srl <- as.data.frame(partly_thawed_vs_intact_fully_thawed_by_PFT_srl)



#### RELATIONSHIPS BETWEEN GA, PEAT TEMP, MOISTURE AND CH4 FLUXES ####

## create a new df without duplicates (which come from the multi-level root trait data within each chamber)
plant_ch4_noroot <- plant_ch4 %>%
  group_by(thaw_stage, chamber_id) %>%
  summarise(
    across(
      c(herb_GA_m2m2, shrub_GA_m2m2,
        ch4_season, ch4_early, ch4_middle, ch4_peak,
        ch4_early_date, ch4_middle_date, ch4_peak_date,
        avgTgnd_early, avgTgnd_middle, avgTgnd_peak, avgTgnd_season, VWC_pct_2425
      ),
      ~ unique(.x)[1]  
    ) 
  ) %>%
  ungroup()

plant_ch4_noroot <- as.data.frame(plant_ch4_noroot)

# for plotting:

pft_colors <- c("herbaceous" = "#213448", "shrub" = "#94B4C1")
shape_vals <- c(partly_thawed = 21, fully_thawed = 24, intact = 22)

# Order of CH4 rows
ch4_order <- c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")

# Helper to keep one CH4 per chamber-stage if needed
unique1 <- function(v){
  vv <- unique(stats::na.omit(v))
  if (length(vv) == 1) return(vv)
  warning("Non-unique values; taking the first.")
  vv[1]
}
# CH4 long
ch4_long <- plant_ch4_noroot %>%
  pivot_longer(cols = all_of(ch4_order), names_to = "ch4", values_to = "y") %>%
  mutate(ch4 = factor(ch4, levels = ch4_order)) %>%
  as.data.frame()

# GA by PFT 
df_ga_pft <- plant_ch4_noroot %>%
  select(chamber_id, thaw_stage, herb_GA_m2m2, shrub_GA_m2m2, all_of(ch4_order)) %>%
  pivot_longer(
    cols = c(herb_GA_m2m2, shrub_GA_m2m2),
    names_to = "PFT_raw",
    values_to = "GA_m2m2"
  ) %>%
  mutate(
    PFT = if_else(PFT_raw == "herb_GA_m2m2", "herbaceous", "shrub"),
    var = "GA_PFT"
  ) %>%
  pivot_longer(cols = all_of(ch4_order), names_to = "ch4", values_to = "y") %>%
  transmute(
    chamber_id, thaw_stage, PFT, trait,
    ch4 = factor(ch4, levels = ch4_order),
    x = GA_m2m2,
    y = y
  )

# GA sum (across PFTs)
df_ga_sum <- plant_ch4_noroot %>%
  mutate(
    GA_sum_m2m2 = if_else(
      is.na(herb_GA_m2m2) & is.na(shrub_GA_m2m2),
      NA_real_,
      coalesce(herb_GA_m2m2, 0) + coalesce(shrub_GA_m2m2, 0)
    )
  ) %>%
  select(chamber_id, thaw_stage, GA_sum_m2m2, all_of(ch4_order)) %>%
  pivot_longer(cols = all_of(ch4_order), names_to = "ch4", values_to = "y") %>%
  transmute(
    chamber_id, thaw_stage,
    var = "GA_sum",
    ch4 = factor(ch4, levels = ch4_order),
    x = GA_sum_m2m2,
    y = y
  )

# peat temp (avgTgnd)

df_tgnd <- plant_ch4_noroot %>%
  select(chamber_id, thaw_stage,
         avgTgnd_early, avgTgnd_middle, avgTgnd_peak, avgTgnd_season,
         all_of(ch4_order)) %>%
  pivot_longer(cols = all_of(ch4_order), names_to = "ch4", values_to = "y") %>%
  mutate(
    x = case_when(
      ch4 == "ch4_early"  ~ avgTgnd_early,
      ch4 == "ch4_middle" ~ avgTgnd_middle,
      ch4 == "ch4_peak"   ~ avgTgnd_peak,
      ch4 == "ch4_season" ~ avgTgnd_season
    ),
    var = "avgTgnd"
  ) %>%
  select(chamber_id, thaw_stage, trait, ch4, x, y) %>%
  mutate(ch4 = factor(ch4, levels = ch4_order))

# VWC
df_vwc <- plant_ch4_noroot %>%
  select(chamber_id, thaw_stage, VWC_pct_2425, all_of(ch4_order)) %>%
  distinct() %>%
  pivot_longer(cols = all_of(ch4_order), names_to = "ch4", values_to = "y") %>%
  transmute(
    chamber_id, thaw_stage,
    var = "VWC",
    ch4 = factor(ch4, levels = ch4_order),
    x = VWC_pct_2425,
    y = y
  ) %>%
  filter(is.finite(x), is.finite(y))

# Ensure facet column order: GA_PFT (col 1), GA_sum (col 2), then others
df_ga_pft$var <- factor(df_ga_pft$var, levels = c("GA_PFT", "GA_sum", "avgTgnd", "VWC"))
df_ga_sum$var <- factor(df_ga_sum$var, levels = c("GA_PFT", "GA_sum", "avgTgnd", "VWC"))
df_tgnd$var   <- factor(df_tgnd$var,   levels = c("GA_PFT", "GA_sum", "avgTgnd", "VWC"))
df_vwc$var    <- factor(df_vwc$var,    levels = c("GA_PFT", "GA_sum", "avgTgnd", "VWC"))

# Plot (Fig. B8)
p_grid <-
  ggplot() +
  # GA_PFT: herb + shrub
  geom_point(
    data = df_ga_pft,
    aes(x = x, y = y, fill = PFT, shape = thaw_stage),
    color = "black",
    size = 3.5, stroke = 0.7, alpha = 0.85
  ) +
  # GA_sum: across-PFT, white fill
  geom_point(
    data = df_ga_sum,
    aes(x = x, y = y, shape = thaw_stage),
    fill = "lightgrey", color = "black",
    size = 3.5, stroke = 0.7, alpha = 0.85
  ) +
  # avgTgnd 
  geom_point(
    data = df_tgnd,
    aes(x = x, y = y, shape = thaw_stage),
    fill = "lightgrey", color = "black",
    size = 3.5, stroke = 0.7, alpha = 0.85
  ) +
  # VWC
  geom_point(
    data = df_vwc,
    aes(x = x, y = y, shape = thaw_stage),
    fill = "lightgrey", color = "black",
    size = 3.5, stroke = 0.7, alpha = 0.85
  ) +
  facet_grid(ch4 ~ var, scales = "free_x") +
  scale_shape_manual(values = shape_vals, name = "Thaw stage") +
  scale_fill_manual(values = pft_colors, name = "PFT") +
  scale_x_continuous(labels = label_number(),
                     expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(labels = label_number(),
                     expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = NULL, y = expression(paste("CH"[4], " flux (mg CH"[4], " m"^-2, " d"^-1, ")"))) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),  # <- added
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 16),
    legend.position = "right",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

p_grid


# Kendall correlations

# One row per chamber with only needed columns
dat <- plant_ch4 %>%
  distinct(
    chamber_id,
    herb_GA_m2m2, shrub_GA_m2m2, VWC_pct_2425,
    ch4_early, ch4_middle, ch4_peak, ch4_season,
    avgTgnd_early, avgTgnd_middle, avgTgnd_peak, avgTgnd_season
  )

# Robust Kendall helper
kendall_corr <- function(x, y) {
  ok <- stats::complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  
  # Too few pairs or a constant variable -> return NA with note
  if (n < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(tibble(n = n, tau = NA_real_, p_value = NA_real_, note = "too_few_or_constant"))
  }
  
  ct <- suppressWarnings(stats::cor.test(x, y, method = "kendall", exact = FALSE))
  tibble(n = n, tau = unname(ct$estimate), p_value = ct$p.value, note = NA_character_)
}

# All combos: x in {herb_GA, shrub_GA, VWC} × y in {early, middle, peak, season}
x_vars <- c("herb_GA_m2m2", "shrub_GA_m2m2", "VWC_pct_2425")
y_vars <- c("ch4_early", "ch4_middle", "ch4_peak", "ch4_season")

grid_res <- tidyr::crossing(x = x_vars, y = y_vars) %>%
  mutate(res = purrr::map2(x, y, ~ kendall_corr(dat[[.x]], dat[[.y]]))) %>%
  tidyr::unnest(res) %>%
  mutate(type = "X-by-all-Y")

# Special paired cases (avgTgnd_* with matching CH4 period)
paired_x <- c("avgTgnd_early", "avgTgnd_middle", "avgTgnd_peak", "avgTgnd_season")
paired_y <- c("ch4_early",      "ch4_middle",      "ch4_peak",      "ch4_season")

paired_res <- tibble(x = paired_x, y = paired_y) %>%
  mutate(res = purrr::map2(x, y, ~ kendall_corr(dat[[.x]], dat[[.y]]))) %>%
  tidyr::unnest(res) %>%
  mutate(type = "paired_Tgnd_vs_CH4")

# Combine
all_results <- bind_rows(grid_res, paired_res) %>%
  arrange(type, x, y) %>%
  mutate(
    tau = round(tau, 3),
    p_value = signif(p_value, 3)
  )


# PFTs summed for GA

# helper
kend_fun <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 3 || dplyr::n_distinct(x) < 2 || dplyr::n_distinct(y) < 2) {
    return(tibble(n = n, tau = NA_real_, p_value = NA_real_, note = "too_few_or_constant"))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "kendall", exact = FALSE))
  tibble(n = n, tau = unname(ct$estimate), p_value = ct$p.value, note = NA_character_)
}

# Kendall per CH4 period 
kendall_by_ch4 <- df_ga_sum %>%
  group_by(ch4) %>%
  summarise(res = list(kend_fun(x, y)), .groups = "drop") %>%
  unnest(res) %>%
  mutate(tau = round(tau, 3), p_value = signif(p_value, 3)) %>%
  arrange(ch4)

#### GENERAL CH4 FLUX TRENDS (Supplementary Methods A2) ####

#### general CH4 flux stats

daily_medians <- daily_medians %>%
  mutate(
    thaw_stage = case_when(
      chamber == 1 ~ "intact",
      chamber == 2 ~ "partly_thawed",
      chamber == 3 ~ "intact",
      chamber == 4 ~ "partly_thawed",
      chamber == 5 ~ "intact",
      chamber == 6 ~ "partly_thawed",
      chamber == 7 ~ "fully_thawed",
      chamber == 8 ~ "fully_thawed",
      chamber == 9 ~ "fully_thawed",
      TRUE ~ NA_character_
    ))

# Make sure date is in Date format
daily_medians <- daily_medians %>%
  mutate(date = as.Date(date),
         month = month(date))

### Overall summary statistics (chamber-level, across whole period)
summary_overall <- daily_medians %>%
  group_by(thaw_stage) %>%
  summarise(
    median_flux = median(daily_median_CH4, na.rm = TRUE),
    IQR_flux = IQR(daily_median_CH4, na.rm = TRUE),
    mean_flux = mean(daily_median_CH4, na.rm = TRUE),
    mean_flux_mean = mean(daily_mean_CH4, na.rm = T),
    sd_flux = sd(daily_median_CH4, na.rm = TRUE),
    max_flux = max(daily_median_CH4, na.rm=T),
    min_flux = min(daily_median_CH4, na.rm = T)
  )

### Summary statistics per month (chamber-level)
summary_monthly <- daily_medians %>%
  group_by(month, thaw_stage) %>%
  summarise(
    median_flux = median(daily_median_CH4, na.rm = TRUE),
    IQR_flux = IQR(daily_median_CH4, na.rm = TRUE),
    mean_flux = mean(daily_median_CH4, na.rm = TRUE),
    sd_flux = sd(daily_median_CH4, na.rm = TRUE),
    max_flux = max(daily_median_CH4, na.rm = T),
    min_flux = min(daily_median_CH4, na.rm = T),
    .groups = "drop"
  )

# Across whole period: compute daily CV per thaw_stage, then average over time
cv_within_stage_daily <- daily_medians %>%
  group_by(thaw_stage, date) %>%
  summarise(
    CV_percent = raster::cv(daily_median_CH4, na.rm = TRUE, aszero = TRUE),
    .groups = "drop"
  )

cv_within_stage_overall <- cv_within_stage_daily %>%
  group_by(thaw_stage) %>%
  summarise(
    CV_percent = mean(CV_percent, na.rm = TRUE),
    .groups = "drop"
  )

# Per month: average the daily CVs within each month
cv_within_stage_monthly <- cv_within_stage_daily %>%
  mutate(month = month(date)) %>%
  group_by(thaw_stage, month) %>%
  summarise(
    CV_percent = mean(CV_percent, na.rm = TRUE),
    .groups = "drop"
  )

### CV% between thaw stages
# For between-stage CV you need one value per thaw_stage per day, then cv across stages.
# Here we aggregate chambers to thaw-stage daily means first
agg_stage_daily <- daily_medians %>%
  group_by(thaw_stage, date) %>%
  summarise(
    stage_flux = mean(daily_median_CH4, na.rm = TRUE),
    .groups = "drop"
  )

# Across whole period: compute daily CV across thaw stages, then average
cv_between_stages_daily <- agg_stage_daily %>%
  group_by(date) %>%
  summarise(
    CV_percent = raster::cv(stage_flux, na.rm = TRUE, aszero = TRUE),
    .groups = "drop"
  )

cv_between_stages_overall <- cv_between_stages_daily %>%
  summarise(
    CV_percent = mean(CV_percent, na.rm = TRUE)
  )

# Per month: average the daily between-stage CVs within each month
cv_between_stages_monthly <- cv_between_stages_daily %>%
  mutate(month = month(date)) %>%
  group_by(month) %>%
  summarise(
    CV_percent = mean(CV_percent, na.rm = TRUE),
    .groups = "drop"
  )

summary_overall
summary_monthly <- as.data.frame(summary_monthly)
summary_monthly

cv_within_stage_overall
cv_within_stage_monthly
cv_between_stages_overall
cv_between_stages_monthly


## Linear mixed models for differences in CH4 fluxes between thaw stages
model <- lme(asinh(daily_median_CH4) ~ thaw_stage,
             random = ~ 1 | chamber,
             data = daily_medians,
             method = "REML")

# check residuals
qqnorm(residuals(model, type = "normalized"))
qqline(residuals(model, type = "normalized"))

summary(model)

# Extract normalized residuals
res <- residuals(model, type = "normalized")

# temporal autocorrelation:
# Simple ACF of all residuals
acf(res, main = "ACF of normalized residuals (all chambers)")

daily_medians %>%
  mutate(resid = residuals(model, type = "normalized")) %>%
  group_by(chamber) %>%
  group_walk(~ {
    acf(.x$resid, main = paste("ACF residuals - Chamber", unique(.x$chamber)))
  })

## there is clear temporal autocorrelation, in general and chamber-dependent
## --> try corAR1
model_ar1 <- lme(asinh(daily_median_CH4) ~ thaw_stage,
                 random = ~ 1 | chamber,
                 correlation = corAR1(form = ~ date | chamber),
                 data = daily_medians,
                 method = "REML")

anova(model, model_ar1)

# normalized residuals from AR(1) model
res_ar1 <- residuals(model_ar1, type = "normalized")

# ACF for all chambers combined
acf(res_ar1, main = "ACF of normalized residuals (AR1 model)")

# ACF per chamber
daily_medians %>%
  mutate(resid_ar1 = residuals(model_ar1, type = "normalized")) %>%
  group_by(chamber) %>%
  group_walk(~ {
    acf(.x$resid_ar1, main = paste("ACF residuals - Chamber", unique(.x$chamber)))
  })

### it does look better with corAR1
### --> date / chamber tells the model that each chamber is temporally autocorrelated with date

### differences between thaw stages:

# Set fully_thawed as reference
daily_medians$thaw_stage <- as.factor(daily_medians$thaw_stage)
daily_medians$thaw_stage <- relevel(daily_medians$thaw_stage, ref = "fully_thawed")

model_fully_thawed <- update(model_ar1, . ~ thaw_stage)
summary(model_fully_thawed)
# intact significantly different from fully thawed
# partly thawed not significantly different from fully thawed

# intact as reference level:
daily_medians$thaw_stage <- relevel(daily_medians$thaw_stage, ref = "intact")

model_intact_ref <- update(model_ar1, . ~ thaw_stage)
summary(model_intact_ref)

# fully and partly thawed differ significantly from intact

# partly thawed as reference level:

daily_medians$thaw_stage <- relevel(daily_medians$thaw_stage, ref = "partly_thawed")

model_partly_thawed <- update(model_ar1, . ~ thaw_stage)
summary(model_partly_thawed)

# intact differs significantly from partly thawed
# fully thawed does not differ significantly from partly thawed

#############################
### POREWATER CH4 and CO2 ###
#############################

pw_ch4 <- read.csv("path/porewater_ch4_alldepths.csv")

# colors
stage_cols <- c(
  "partly_thawed"   = "#7F96B5",
  "fully_thawed"   = "#332D56"
)

# order depth bins top to bottom 
depth_order <- c("1-5","10-14","20-24", "30-34", "40-44", "50-54", "60-64", "70-74", "80-84") 

pw_ch4 <- pw_ch4 %>%
  dplyr::mutate(depth_std = factor(depth_std, levels = depth_order))

# ranges + linear transform so secondary axis can be used
conc_rng <- range(pw_ch4$CH4_conc_mM, na.rm = TRUE)
iso_rng  <- range(pw_ch4$d13C_CH4,   na.rm = TRUE)
iso_to_conc <- function(z) (z - iso_rng[1]) / diff(iso_rng) * diff(conc_rng) + conc_rng[1]
conc_to_iso <- function(x) (x - conc_rng[1]) / diff(conc_rng) * diff(iso_rng) + iso_rng[1]

conc_rng_co2 <- range(pw_ch4$CO2_conc_mM, na.rm = TRUE)
iso_rng_co2  <- range(pw_ch4$d13C_CO2,   na.rm = TRUE)
iso_to_conc_co2 <- function(z) (z - iso_rng[1]) / diff(iso_rng) * diff(conc_rng) + conc_rng[1]
conc_to_iso_co2 <- function(x) (x - conc_rng[1]) / diff(conc_rng) * diff(iso_rng) + iso_rng[1]

# get the means and SDs for CH4 concentrations and 13C
# (and precomputed offsets & x-values for error bars)
pw_means_ch4 <- pw_ch4 %>%
  dplyr::group_by(thaw_stage, depth_std) %>%
  dplyr::summarise(
    CH4_mean = mean(CH4_conc_mM, na.rm = TRUE),
    CH4_sd   = sd(  CH4_conc_mM, na.rm = TRUE),
    iso_mean = mean(d13C_CH4,    na.rm = TRUE),
    iso_sd   = sd(  d13C_CH4,    na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(
    y_conc = as.numeric(depth_std) + 0.15,
    y_iso  = as.numeric(depth_std) - 0.15,
    # map isotope mean and ±SD to the primary (concentration) x-axis
    iso_x  = iso_to_conc(iso_mean),
    iso_xmin = iso_to_conc(iso_mean - iso_sd),
    iso_xmax = iso_to_conc(iso_mean + iso_sd),
    conc_xmin = pmax(0, CH4_mean - CH4_sd), 
    conc_xmax = CH4_mean + CH4_sd
  )

pw_means_ch4 <- as.data.frame(pw_means_ch4)
pw_means_ch4

# get the means and SDs for CO2 concentrations and 13C
# (and precomputed offsets & x-values for error bars)

pw_means_co2 <- pw_ch4 %>%
  dplyr::group_by(thaw_stage, depth_std) %>%
  dplyr::summarise(
    CO2_mean = mean(CO2_conc_mM, na.rm = TRUE),
    CO2_sd   = sd(  CO2_conc_mM, na.rm = TRUE),
    iso_mean = mean(d13C_CO2,    na.rm = TRUE),
    iso_sd   = sd(  d13C_CO2,    na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(
    y_conc = as.numeric(depth_std) + 0.15,
    y_iso  = as.numeric(depth_std) - 0.15,
    # map isotope mean and ±SD to the primary (concentration) x-axis
    iso_x  = iso_to_conc(iso_mean),
    iso_xmin = iso_to_conc(iso_mean - iso_sd),
    iso_xmax = iso_to_conc(iso_mean + iso_sd),
    conc_xmin = pmax(0, CO2_mean - CO2_sd),    # cap at 0 if you prefer
    conc_xmax = CO2_mean + CO2_sd
  )

pw_means_co2 <- as.data.frame(pw_means_co2)
pw_means_co2

### plot

# compute alphaC and CH4:CO2 
pw_iso <- pw_ch4 %>%
  mutate(
    alphaC       = (d13C_CO2 + 1000) / (d13C_CH4 + 1000),
    ch4_co2_rat = CH4_conc_mM / CO2_conc_mM
  )

# summarize means + SD 
pw_sum <- pw_iso %>%
  group_by(thaw_stage, depth_std) %>%
  summarise(
    CH4_mean_mM   = mean(CH4_conc_mM, na.rm = TRUE),
    CH4_sd_mM     = sd(CH4_conc_mM, na.rm = TRUE),
    
    CO2_mean_mM   = mean(CO2_conc_mM, na.rm = TRUE),
    CO2_sd_mM     = sd(CO2_conc_mM, na.rm = TRUE),
    
    d13C_CH4_mean = mean(d13C_CH4, na.rm = TRUE),
    d13C_CH4_sd   = sd(d13C_CH4, na.rm = TRUE),
    
    d13C_CO2_mean = mean(d13C_CO2, na.rm = TRUE),
    d13C_CO2_sd   = sd(d13C_CO2, na.rm = TRUE),
    
    alphaC_mean   = mean(alphaC, na.rm = TRUE),
    alphaC_sd     = sd(alphaC, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    CH4_sd_mM   = replace_na(CH4_sd_mM, 0),
    CO2_sd_mM   = replace_na(CO2_sd_mM, 0),
    d13C_CH4_sd = replace_na(d13C_CH4_sd, 0),
    d13C_CO2_sd = replace_na(d13C_CO2_sd, 0),
    alphaC_sd   = replace_na(alphaC_sd, 0)
  )

depth_lvls <- levels(pw_sum$depth_std)

# ch4 concentration plot
p_ch4_conc <- ggplot(pw_sum, aes(y = as.numeric(depth_std), group = thaw_stage)) +
  geom_path(aes(x = CH4_mean_mM, color = thaw_stage), linewidth = 0.9) +
  geom_errorbarh(aes(
    xmin = CH4_mean_mM - CH4_sd_mM,
    xmax = CH4_mean_mM + CH4_sd_mM,
    color = thaw_stage
  ), height = 0.18, linewidth = 0.8) +
  geom_point(aes(x = CH4_mean_mM, fill = thaw_stage),
             shape = 21, size = 4, color = "black", stroke = 0.4) +
  scale_color_manual(values = stage_cols, name = "Thaw stage") +
  scale_fill_manual(values = stage_cols, guide = "none") +
  scale_y_continuous(
    breaks = seq_along(depth_lvls),
    labels = depth_lvls,
    trans = "reverse",
    expand = expansion(add = 0.4)
  ) +
  scale_x_continuous(
    name = expression(CH[4]*" concentration (mM)"),
    expand = expansion(mult = 0.02)
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_text(size = 18),
    axis.title.x.bottom = element_text(size = 18),
    legend.title = element_blank()
  )

# 13C CH4 plot

p_d13C_ch4 <- ggplot(pw_sum, aes(y = as.numeric(depth_std), group = thaw_stage)) +
  geom_path(aes(x = d13C_CH4_mean, color = thaw_stage), linewidth = 0.9) +
  geom_errorbarh(aes(
    xmin = d13C_CH4_mean - d13C_CH4_sd,
    xmax = d13C_CH4_mean + d13C_CH4_sd,
    color = thaw_stage
  ), height = 0.18, linewidth = 0.8) +
  geom_point(aes(x = d13C_CH4_mean, fill = thaw_stage),
             shape = 21, size = 4.2, color = "black", stroke = 0.4) +
  scale_color_manual(values = stage_cols, name = "Thaw stage") +
  scale_fill_manual(values = stage_cols, guide = "none") +
  scale_y_continuous(
    breaks = seq_along(depth_lvls),
    labels = depth_lvls,
    trans = "reverse",
    expand = expansion(add = 0.4)
  ) +
  scale_x_continuous(
    name = expression(delta^{13}*C*"-CH"[4]*" ("*"\u2030"*")"),
    expand = expansion(mult = 0.02)
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_text(size = 18),
    axis.text.y = element_blank(),
    axis.title.x.bottom = element_text(size = 18),
    legend.title = element_blank()
  )


# CO2 concentration plot

p_co2_conc <- ggplot(pw_sum, aes(y = as.numeric(depth_std), group = thaw_stage)) +
  geom_path(aes(x = CO2_mean_mM, color = thaw_stage), linewidth = 0.9) +
  geom_errorbarh(aes(
    xmin = CO2_mean_mM - CO2_sd_mM,
    xmax = CO2_mean_mM + CO2_sd_mM,
    color = thaw_stage
  ), height = 0.18, linewidth = 0.8) +
  geom_point(aes(x = CO2_mean_mM, fill = thaw_stage),
             shape = 21, size = 4, color = "black", stroke = 0.4) +
  scale_color_manual(values = stage_cols, name = "Thaw stage") +
  scale_fill_manual(values = stage_cols, guide = "none") +
  scale_y_continuous(
    breaks = seq_along(depth_lvls),
    labels = depth_lvls,
    trans = "reverse",
    expand = expansion(add = 0.4)
  ) +
  scale_x_continuous(
    name = expression(CO[2]*" concentration (mM)"),
    breaks = seq(0, 12, by = 2),   
    expand = expansion(mult = 0.02)
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_text(size = 18),
    axis.title.x.bottom = element_text(size = 18),
    legend.title = element_blank()
  )


# 13C CO2 plot

p_d13C_co2 <- ggplot(pw_sum, aes(y = as.numeric(depth_std), group = thaw_stage)) +
  geom_path(aes(x = d13C_CO2_mean, color = thaw_stage), linewidth = 0.9) +
  geom_errorbarh(aes(
    xmin = d13C_CO2_mean - d13C_CO2_sd,
    xmax = d13C_CO2_mean + d13C_CO2_sd,
    color = thaw_stage
  ), height = 0.18, linewidth = 0.8) +
  geom_point(aes(x = d13C_CO2_mean, fill = thaw_stage),
             shape = 21, size = 4.2, color = "black", stroke = 0.4) +
  scale_color_manual(values = stage_cols, name = "Thaw stage") +
  scale_fill_manual(values = stage_cols, guide = "none") +
  scale_y_continuous(
    breaks = seq_along(depth_lvls),
    labels = depth_lvls,
    trans = "reverse",
    expand = expansion(add = 0.4)
  ) +
  scale_x_continuous(
    name = expression(delta^{13}*C*"-CO"[2]*" ("*"\u2030"*")"),
    expand = expansion(mult = 0.02)
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_text(size = 18),
    axis.text.y = element_blank(),
    axis.title.x.bottom = element_text(size = 18),
    legend.title = element_blank()
  )

# aplhaC plot

p_alpha <- ggplot(pw_sum, aes(y = as.numeric(depth_std), group = thaw_stage)) +
  geom_path(aes(x = alphaC_mean, color = thaw_stage), linewidth = 0.9) +
  geom_errorbarh(aes(
    xmin = alphaC_mean - alphaC_sd,
    xmax = alphaC_mean + alphaC_sd,
    color = thaw_stage
  ), height = 0.18, linewidth = 0.8) +
  geom_point(aes(x = alphaC_mean, fill = thaw_stage),
             shape = 21, size = 4, color = "black", stroke = 0.4) +
  scale_color_manual(values = stage_cols, name = "Thaw stage") +
  scale_fill_manual(values = stage_cols, guide = "none") +
  scale_y_continuous(
    breaks = seq_along(depth_lvls),
    labels = depth_lvls,
    trans = "reverse",
    expand = expansion(add = 0.4)
  ) +
  scale_x_continuous(
    name = expression(alpha[C]*" ("*CH[4]*"/"*CO[2]*")"),
    breaks = sort(unique(c(pretty(pw_sum$alphaC_mean, n = 5), 1.07))),
    expand = expansion(mult = 0.02)
  ) +
  coord_cartesian(
    xlim = range(c(
      pw_sum$alphaC_mean - pw_sum$alphaC_sd,
      pw_sum$alphaC_mean + pw_sum$alphaC_sd,
      1.07
    ), na.rm = TRUE)
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_text(size = 18),
    axis.title.x.bottom = element_text(size = 18),
    legend.title = element_blank()
  )

# combine into one figure (Fig. B11)
final_plot <- 
  (p_ch4_conc + p_d13C_ch4) /
  (p_co2_conc + p_d13C_co2) /
  (p_alpha + plot_spacer()) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

final_plot


#################

### check an average percentage of how much root was <2 mm ###

pct <- read.csv("path/root_percentages_plot.csv")

pct_2 <- pct %>%
  mutate(thaw_stage = str_extract(plot_id, "^[a-zA-Z]+")) %>% 
  group_by(thaw_stage, PFT, size_class) %>%
  summarise(mean_pct_length = mean(pct_length, na.rm = TRUE), .groups = "drop")


#################

### peat moisture ###

#  Mean and SD per thaw_stage 
plant_ch4 %>%
  distinct(thaw_stage, chamber_id, VWC_pct_2425, .keep_all = TRUE) %>%
  group_by(thaw_stage) %>%
  summarise(
    mean_VWC = mean(VWC_pct_2425, na.rm = TRUE),
    sd_VWC   = sd(VWC_pct_2425, na.rm = TRUE),
    n        = n()
  )

# Mean and SD per chamber_id 
plant_ch4 %>%
  distinct(thaw_stage, chamber_id, VWC_pct_2425, .keep_all = TRUE) %>%
  group_by(chamber_id) %>%
  summarise(
    thaw_stage = first(thaw_stage),
    mean_VWC = mean(VWC_pct_2425, na.rm = TRUE),
    sd_VWC   = sd(VWC_pct_2425, na.rm = TRUE),
    n        = n()
  )


######################
#### MAP (Fig. 1) ####
######################
# Get country outlines (scale = "medium" or "large")
countries <- ne_countries(scale = "medium", returnclass = "sf")

# Filter for Norway, Sweden, Finland
scandinavia <- countries %>%
  dplyr::filter(admin %in% c("Sweden", "Norway", "Finland", "Russia", "Denmark",
                             "Estonia", "Latvia", "Lithuania", "Belarus",
                             "Poland", "Kaliningrad", "Germany"))

# Set CRS for curved polar projection (you can adjust lat_0/lon_0)
crs_arctic <- "+proj=laea +lat_0=65 +lon_0=20 +ellps=GRS80 +units=m +no_defs"

# Transform for the projection
scandinavia_laea <- st_transform(scandinavia, crs = crs_arctic)

# Define your custom bounds (units in meters since you're in projected CRS)
x_bounds <- c(-900000, 600000)     # good full width for Scandinavia
y_bounds <- c(-1100000, 800000)     # cut off Svalbard by reducing ymax

# stordalen coord point
point_coords <- st_sfc(
  st_point(c(19.05, 68.3333)),
  crs = 4326
)
# Project point
point_proj <- st_transform(point_coords, crs = crs_arctic)


# Create long sequence of longitude points at fixed latitude 66.56
arctic_longitudes <- seq(-180, 180, by = 1)
arctic_latitude <- rep(66.56, length(arctic_longitudes))

# Create sf LINESTRING in WGS84
arctic_circle <- st_sfc(
  st_linestring(cbind(arctic_longitudes, arctic_latitude)),
  crs = 4326
)

# Project it into your LAEA Arctic CRS
arctic_circle_proj <- st_transform(arctic_circle, crs = crs_arctic)

inset_map <- ggplot() +
  geom_sf(data = scandinavia_laea, fill = "white", color = "black") +  # land stays white
  geom_sf(data = point_proj, fill = "yellow", color = "black", shape = 21, size = 3) +
  coord_sf(crs = crs_arctic, xlim = x_bounds, ylim = y_bounds, expand = FALSE) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = "white", size = 2),
    plot.margin = margin(0, 0, 0, 0)
  )

# import plot coordinates
plot_coords <- read.csv("path/24072023plots_stordalen.csv")


plot_coords <- plot_coords %>%
  # Rename first 3 columns
  rename(
    plot_ID = 1,
    lat = 2,
    lon = 3
  ) %>%
  # Filter out rows where 'plot' contains unwanted terms
  filter(!grepl("EC|test|tesu", plot_ID, ignore.case = TRUE)) %>%
  # Keep only the first 3 columns
  dplyr::select(plot_ID, lat, lon)

# rename plots
# Create a lookup table for the numeric class mapping
class_map <- c("1" = "intact", "2" = "partly_thawed", "3" = "intact",
               "4" = "partly_thawed", "5" = "intact", "6" = "partly_thawed",
               "7" = "fully_thawed", "8" = "fully_thawed", "9" = "fully_thawed")

# Apply transformation
plot_coords <- plot_coords %>%
  mutate(
    plot_ID = if_else(
      str_detect(plot_ID, "ruth"),
      {
        number <- str_extract(plot_ID, "\\d+")
        type <- class_map[number]
        str_c(type, number, "rv")
      },
      plot_ID  # Leave others unchanged
    )
  )

plot_coords <- plot_coords %>%
  mutate(
    lat = str_remove(lat, "N$"),     # remove 'N' from end
    lon = str_remove(lon, "E$"),
    lat = as.numeric(lat),           # convert to numeric
    lon = as.numeric(lon)            # ensure lon is numeric too
  )

plot_coords <- plot_coords %>%
  mutate(thaw_stage = case_when(
    str_detect(plot_ID, "intact") ~ "intact",
    str_detect(plot_ID, "partly_thawed") ~ "partly_thawed",
    str_detect(plot_ID, "fully_thawed") ~ "fully_thawed"
  ))


plot_coords <- plot_coords %>%
  mutate(plot_type = case_when(
    str_detect(plot_ID, "rv") ~ "chamber",
    str_detect(plot_ID, "tm") ~ "core"
  ))


# thaw stage classes from EMERGE (https://emerge-db.asc.ohio-state.edu/datasources/0145_Wv2-2014_GroundCoverClassifications)

class_map <- rast("path/mire_classification.tif")

# Convert to data frame with coordinates
class_map_df <- as.data.frame(class_map, xy = TRUE, na.rm = TRUE)

# Define bounds
xmin <- 419450
xmax <- 420000
ymin <- 7583450
ymax <- 7584000

# Crop the dataframe
class_map_crop_df <- class_map_df %>%
  dplyr::filter(
    x >= xmin, x <= xmax,
    y >= ymin, y <= ymax
  )

# Plot with manual fill scale
ggplot(class_map_crop_df, aes(x = x, y = y, fill = as.factor(mire_classification))) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual(
    values = class_colors,
    name = "Vegetation Class"
  ) +
  theme_minimal()

# the numbers correspond to: 

class_map_crop_df <- class_map_crop_df %>%
  mutate(
    veg_label = case_when(
      mire_classification == 1 ~ "open_water",
      mire_classification == 2 ~ "hummock",
      mire_classification == 3 ~ "other",
      mire_classification == 4 ~ "rock",
      mire_classification == 5 ~ "semiwet",
      mire_classification == 6 ~ "tall_graminoid",
      mire_classification == 7 ~ "tall_shrub",
      mire_classification == 8 ~ "wet",
      TRUE ~ NA_character_
    )
  )

# create intact partly_thawed fully_thawed, open water rock and other classes

class_map_crop_df <- class_map_crop_df %>%
  mutate(
    thaw_stage = case_when(
      veg_label %in% c("tall_shrub", "hummock")       ~ "intact",
      veg_label %in% c("semiwet", "wet")              ~ "partly_thawed",
      veg_label == "tall_graminoid"                   ~ "fully_thawed",
      veg_label == "open_water"                       ~ "open_water",
      veg_label == "rock"                             ~ "rock",
      veg_label == "other"                            ~ "other",
      TRUE ~ NA_character_
    )
  )
# Convert to sf object with WGS84 CRS
plot_coords_sf <- st_as_sf(plot_coords, coords = c("lon", "lat"), crs = 4326)

crs(class_map)

# Transform to match the raster CRS (UTM zone 34N)
plot_coords_proj <- st_transform(plot_coords_sf, crs = 32634)
class_colors <- c(
  "open_water" = "#5C8CE0",  
  "intact" = "#9BD9D5", 
  "partly_thawed" = "#4E6688", 
  "fully_thawed" = "#332D56", 
  "rock" = "white", 
  "other" = "white" 
)

class_map_crop_df$thaw_stage <- factor(
  class_map_crop_df$thaw_stage,
  levels = c("intact", "partly_thawed", "fully_thawed", "open_water", "rock", "other")
)

# Create a one-row data frame for EC tower
EC_df <- data.frame(
  plot_ID = "ECreal",
  lon = 19.04518,
  lat = 68.35594061
)

# Convert to sf object in WGS84 (lat/lon)
EC_sf <- st_as_sf(EC_df, coords = c("lon", "lat"), crs = 4326)

# Transform to match your map's CRS (EPSG:32634)
EC_proj <- st_transform(EC_sf, crs = 32634)


# plot the map
raster_map <- ggplot() +
  # Raster
  geom_raster(data = class_map_crop_df, aes(x = x, y = y, fill = thaw_stage)) +
  scale_fill_manual(
    values = class_colors,
    name = ""
  ) +
  
  new_scale_fill() +
  
  # Points
  geom_sf(
    data = plot_coords_proj,
    aes(shape = plot_type, fill = thaw_stage), 
    color = "yellow", 
    size = 4.5,
    stroke = 1   
  ) +
  scale_shape_manual(
    values = c("core" = 21, "chamber" = 22), 
    name = ""
  ) +
  scale_fill_manual(
    values = c(
      "intact" = "#71C0BB",  
      "partly_thawed" = "#4E6688",  
      "fully_thawed" = "#332D56"  
    ),
    name = ""
  ) +
  geom_sf(data = EC_proj, shape = 4, color = "yellow", size = 4, stroke = 2.5) +
  annotation_scale(location = "br", width_hint = 0.2, style = "bar") + 
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.3, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_minimal) +
  theme_minimal() +
  coord_sf(crs = 32634, label_graticule = "SE", expand = FALSE) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 10), 
    axis.ticks = element_line(),
    panel.grid.major = element_blank()
  )

raster_map

# Combine inset and raster map (Fig. 1- cleaned in Inkscape)
combined_map <- raster_map +
  inset_element(
    inset_map,
    left = 0.02, bottom = 0.75,
    right = 0.4, top = 0.98,
    align_to = "full" 
  )

combined_map


