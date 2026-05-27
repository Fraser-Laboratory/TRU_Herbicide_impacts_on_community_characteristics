# CRAN packages
ls <- c("rstatix", "tidyverse", "ggpubr", "ggtext", "emmeans", "vegan", "betapart",
        "janitor", "lme4", "lmerTest", "easystats", "devtools", "BiocManager",
        "patchwork", "picante", "DHARMa", "GUniFrac", "glmmTMB", "sf", "bcdata",
        "bcmaps")

new_packages <- ls[!(ls %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# pairwiseAdonis, V.Phylomaker2, ggtree, and ggtreeExtra packages are not on CRAN, 
# install using devtools/BiocManager if needed
if(!"pairwiseAdonis" %in% installed.packages()[,"Package"]) {
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}
if(!"V.PhyloMaker2" %in% installed.packages()[,"Package"]) {
  devtools::install_github("jinyizju/V.PhyloMaker2")
}
if(!"ggtree" %in% installed.packages()[,"Package"]) {
  BiocManager::install("ggtree", update = FALSE)
}
if(!"ggtreeExtra" %in% installed.packages()[,"Package"]) {
  BiocManager::install("ggtreeExtra", update = FALSE)
}

# Load packages, don't need devtools or BiocManager loaded
ls <- c(ls[!ls %in% c("devtools", "BiocManager")], 
        "pairwiseAdonis", "V.PhyloMaker2", "ggtree", "ggtreeExtra")
invisible(suppressPackageStartupMessages(
  lapply(ls, function(x) library(x, character.only = TRUE))))
rm(ls, new_packages)

# Import data -------------------------------------------------------------

veg <- read.csv("Data/All_Diversity_Common.csv") %>% 
  slice(-137) %>% # Remove an empty row of data
  mutate(time = case_when(
    year == "2018" ~ "baseline",
    year == "2019" ~ "year_1",
    .default = "year_2"
  ),
  row_num = paste0("a", 1:nrow(.))) %>% 
  column_to_rownames(., var = "row_num")

# Metadata
meta <- veg %>% 
  select(!spotted_knapweed:meadow_brome) %>% 
  rownames_to_column(., var = "id")

# Phylogenetic information
sp_plnt <- readxl::read_xlsx("Data/list_spp.xlsx") %>% 
  janitor::clean_names() %>% 
  dplyr::select(species = latin_name, family) %>% 
  mutate(genus = str_remove(species, pattern = " .*$")) %>% 
  relocate(genus, .before = family)

# Create output directories for data and plots
invisible(lapply(c("Plots", "Output"), dir.create, showWarnings = FALSE))

# Alpha diversity ---------------------------------------------------------

df_plant <- veg %>% 
  mutate(richness = specnumber(select(.,spotted_knapweed:meadow_brome)),
         shann = diversity(select(., spotted_knapweed:meadow_brome), index = "shannon"),
         simp = diversity(select(., spotted_knapweed:meadow_brome), index = "invsimpson"),
         evenn = shann/log(richness))

veg_tree <- read.csv("Data/All_Diversity_Latin.csv") %>% 
  janitor::clean_names(case = "sentence") %>% 
  slice(-137) %>% 
  mutate(time = case_when(
    Year == "2018" ~ "baseline",
    Year == "2019" ~ "year_1",
    .default = "year_2"
  ),
  row_num = paste0("a", 1:nrow(.))) %>% 
  column_to_rownames(., var = "row_num")


meta_tree <- veg_tree %>% 
  select(!`Centaurea stoebe`:`Bromus commutatus`) %>% 
  rownames_to_column(., var = "id")

df_plant_tree <- veg_tree %>% 
  select(`Centaurea stoebe`:`Bromus commutatus`)


sp_name <- colnames(df_plant_tree) %>% 
  gsub(., pattern = " ", replacement = "_")

names(df_plant_tree) <- sp_name

# Alpha diversity analysis ----------------------------------------------

# Model with richness data subsetted to 2019 and 2020 to prevent dropping of 
# spraying as a factor for all downstream analysis

# Seeding had no impact so dropped it and dropped 2018 year
m1 <- lmer(richness ~ sprayed*ash*time + (1|site),  subset(df_plant, year != 2018)) 
car::Anova(m1, test.statistic = "F") %>% 
  broom::tidy() %>% 
  write.csv("Output/Richness_lmer.csv", row.names = FALSE)
simulateResiduals(m1, plot = TRUE)

# Posthoc plot for species richness differences when sprayed
ph_s_rich <- emmeans(m1, ~ sprayed, type  = "response") %>% 
  multcomp::cld(Letters = letters) %>% 
  as.data.frame() %>% 
  mutate(.group = str_trim(.group, side = "both")) %>% 
  ggplot(., aes(x = sprayed, y = emmean, label = .group)) +
  geom_pointrange(aes(ymin = emmean - SE, ymax = emmean + SE),
                 pch = 23, fill = "steelblue", size =  1.1, linewidth = 1.2 ) +
  geom_text(aes( y = emmean + SE + 0.5)) +
  scale_x_discrete(breaks = c("FALSE", "TRUE"),
                   labels = c("Unsprayed", "Sprayed")) +
  coord_cartesian(
    ylim = c(0, 10), 
    expand = c(top = FALSE, left = TRUE, bottom = FALSE, right = TRUE)) +
  theme_bw(base_size = 16) +
  theme(
    text = element_text(family = "Times New Roman"),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black")
  ) +
  labs(y = "Observed species richness")
ph_s_rich

# Posthoc for ash*sprayed
ph_as_rich <- emmeans::emmeans(m1, ~ ash*sprayed) %>% 
  multcomp::cld(Letters = letters) %>% 
  as.data.frame() %>% 
  mutate(.group = str_trim(.group, side = "both"),
         ash = factor(ash, levels = c("none", "low", "high")),
         sprayed = ifelse(sprayed, "Sprayed", "Unsprayed")) %>% 
  ggplot(., aes(x = ash, y = emmean, label = .group)) +
  geom_pointrange(aes(ymin = emmean - SE, ymax = emmean + SE),
                  pch = 23, fill = "steelblue", size =  1.1, linewidth = 1.2 ) +
  facet_wrap(~ sprayed) +
  geom_text(aes( y = emmean + SE + 0.5)) +
  scale_x_discrete(breaks = c("none", "low", "high"),
                   labels = c("None", "Low", "High")) +
  coord_cartesian(
    ylim = c(0, 11), 
    expand = c(top = FALSE, left = TRUE, bottom = FALSE, right = TRUE)) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(y = "Observed species richness", x = NULL)
ph_as_rich

# Posthoc for time*sprayed
ph_ts_rich <- emmeans(m1, ~ time*sprayed) %>% 
  multcomp::cld(Letters = letters) %>% 
  as.data.frame() %>% 
  mutate(.group = str_trim(.group, side = "both"),
         time = factor(time, levels = c("year_1", "year_2")),
         sprayed = ifelse(sprayed, "Sprayed", "Unsprayed")) %>% 
  ggplot(., aes(x = time, y = emmean, label = .group)) +
  geom_pointrange(aes(ymin = emmean - SE, ymax = emmean + SE),
                  pch = 23, fill = "steelblue", size =  1.1, linewidth = 1.2 ) +
  facet_wrap(~ sprayed) +
  geom_text(aes( y = emmean + SE + 0.5)) +
  scale_x_discrete(breaks = c("year_1", "year_2"),
                   labels = c("2019", "2020")) +
  coord_cartesian(
    ylim = c(0, 11), 
    expand = c(top = FALSE, left = TRUE, bottom = FALSE, right = TRUE)) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(y = "Observed species richness", x = NULL)
ph_ts_rich


# Phylogenetic alpha diversity analysis ----------------------------------------

tree <- phylo.maker(sp_plnt, output.tree = TRUE)
matched <- match.phylo.comm(tree$scenario.3, df_plant_tree)
mytree <- multi2di.phylo(matched$phy)
pd_veg <- pd(matched$comm, tree = mytree, include.root = TRUE) %>%
  rownames_to_column(., var = "id") %>%
  inner_join(., meta, by = "id")

phy_m1 <- lmer(PD ~ sprayed*ash*time + (1|site),  subset(pd_veg, year != 2018)) # seeding had no impact so dropped it and dropped 2018 year
car::Anova(phy_m1, test.statistic = "F") %>% 
  broom::tidy() %>% 
  write.csv("Output/PD_lmer.csv", row.names = FALSE)
DHARMa::simulateResiduals(phy_m1, plot = T)
effectsize::eta_squared(phy_m1 , partial = TRUE)

ph_s_pd <- emmeans(phy_m1, ~ sprayed, type  = "response") %>% 
  multcomp::cld(Letters = letters) %>% 
  as.data.frame() %>% 
  mutate(.group = str_trim(.group, side = "both")) %>% 
  ggplot(., aes(x = sprayed, y = emmean, label = .group)) +
  geom_pointrange(aes(ymin = emmean - SE, ymax = emmean + SE),
                  pch = 23, fill = "steelblue", size =  1.1, linewidth = 1.2 ) +
  geom_text(aes( y = emmean + SE + 30)) +
  scale_x_discrete(breaks = c("FALSE", "TRUE"),
                   labels = c("Unsprayed", "Sprayed")) +
   theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black")
  ) +
  labs(y = "Faith's phylogenetic diversity")
ph_s_pd

# Posthoc for ash*sprayed
ph_as_pd <- emmeans::emmeans(phy_m1, ~ ash*sprayed) %>% 
  multcomp::cld(Letters = letters) %>% 
  as.data.frame() %>% 
  mutate(.group = str_trim(.group, side = "both"),
         ash = factor(ash, levels = c("none", "low", "high")),
         sprayed = ifelse(sprayed, "Sprayed", "Unsprayed")) %>% 
  ggplot(., aes(x = ash, y = emmean, label = .group)) +
  geom_pointrange(aes(ymin = emmean - SE, ymax = emmean + SE),
                  pch = 23, fill = "steelblue", size =  1.1, linewidth = 1.2 ) +
  facet_wrap(~ sprayed) +
  geom_text(aes( y = emmean + SE + 30)) +
  scale_x_discrete(breaks = c("none", "low", "high"),
                   labels = c("None", "Low", "High")) +
  # coord_cartesian(
  #   ylim = c(0, 11), 
  #   expand = c(top = FALSE, left = TRUE, bottom = FALSE, right = TRUE)) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    strip.text.x = element_blank(),
    strip.background = element_blank()
  ) +
  labs(y = "Faith's phylogenetic diversity", x = "Ash treatment")

# Posthoc for time*sprayed
ph_ts_pd <- emmeans(phy_m1, ~ time*sprayed) %>% 
  multcomp::cld(Letters = letters) %>% 
  as.data.frame() %>% 
  mutate(.group = str_trim(.group, side = "both"),
         time = factor(time, levels = c("year_1", "year_2")),
         sprayed = ifelse(sprayed, "Sprayed", "Unsprayed")) %>% 
  ggplot(., aes(x = time, y = emmean, label = .group)) +
  geom_pointrange(aes(ymin = emmean - SE, ymax = emmean + SE),
                  pch = 23, fill = "steelblue", size =  1.1, linewidth = 1.2 ) +
  facet_wrap(~ sprayed) +
  geom_text(aes( y = emmean + SE + 30)) +
  scale_x_discrete(breaks = c("year_1", "year_2"),
                   labels = c("2019", "2020")) +
  # coord_cartesian(
  #   ylim = c(0, 11), 
  #   expand = c(top = FALSE, left = TRUE, bottom = FALSE, right = TRUE)) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    strip.text.x = element_blank(),
    strip.background = element_blank()
  ) +
  labs(y = "Faith's phylogenetic diversity", x = NULL)

# Combine alpha diversity and phylogenetic diversity plots
alpha_div_s <- ph_s_rich + ph_s_pd
alpha_div_as <- ph_as_rich / ph_as_pd
alpha_div_ts <- ph_ts_rich / ph_ts_pd

# Save plots
ggsave("Plots/alpha_diversity_sprayed.png", alpha_div_s, height = 4, width = 8, dpi = 800)
ggsave("Plots/alpha_diversity_sprayed_ash.png", alpha_div_as, height = 8, width = 8, dpi = 800)
ggsave("Plots/alpha_diversity_sprayed_time.png", alpha_div_ts, height = 8, width = 8, dpi = 800)

# Beta diversity analysis ------------------------------------------------

# Bray-Curtis dissimilarity calculation
bray_dist <- veg %>% 
  select(spotted_knapweed:meadow_brome) %>% 
  decostand(., method = "hellinger") %>% 
  vegdist(., method = "bray")

# NMDS
set.seed(111) 
nmds <- metaMDS(bray_dist, k = 3, trymax = 999)  
veg_stress <- round(nmds$stress, 2)

# PERMANOVA - to put in supplementary info as a table
set.seed(111111)
mod_ado <- adonis2(bray_dist ~ sprayed*ash*time, by = "term",  meta)
mod_ado_df <- data.frame(mod_ado)
write.csv(mod_ado_df, "Output/veg_abundance_PERMANOVA.csv", row.names = TRUE)

# Prepare data for plotting
plot_df <- data.frame(scores(nmds), meta) %>% 
  mutate(ash = factor(ash, levels = c("none", "low", "high")))

# Plots for each year, showing separation due to herbicide (and ash) treatments
veg_abun_nmds <- ggplot(
  plot_df, aes(x = NMDS1,
               y = NMDS2,
               fill = sprayed
  )) +
  geom_point(aes(shape = ash), size  = 2) +
  facet_grid(~ year, scales = "free") +
  stat_ellipse() +
  theme_bw(base_size = 18) +
  scale_shape_manual(name = "Ash treatment", 
                     values = c(21, 22, 25),
                     label = c("None", "Low", "High")) +
  scale_fill_viridis_d(name = "           Herbicide",
                       labels = c("FALSE" = "Unsprayed",
                                  "TRUE" = "Sprayed"),
                       begin = 0.2, end = 0.8) +
  guides(fill = guide_legend(override.aes = list(shape =21, size = 4)),
         shape = guide_legend(override.aes = list( size = 4))) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = 0),
        legend.key.size =  unit(5, "pt"),
        legend.key.spacing =  unit(5, "pt"),
        legend.spacing = unit(30, "pt"),
        #legend.box = "vertical",
        legend.box.just = "left",
        legend.title = element_text(size =15, face = "bold"),
        text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid = element_blank()) +
  labs(title = "NMDS of Bray-Curtis dissimilarity")

# ggsave("Plots/veg_abun_nmds.png", veg_abun_nmds, width = 10, height = 6, dpi = 800)

# Phylogenetic beta diversity ---------------------------------------------

unifracs  <- GUniFrac(matched$comm, matched$phy, alpha = 0.5)$unifracs
weighted_unifrac <- unifracs[, , "d_0.5"]

set.seed(213)
phy_ado2 <- adonis2( weighted_unifrac ~ sprayed*ash*time, by = "term", meta)
phy_ado2_df <- data.frame(phy_ado2)
write.csv(phy_ado2_df, "Output/Weighed_Unifrac_PERMANOVA.csv", row.names = TRUE)

set.seed(2112)
w_nmds <- metaMDS(weighted_unifrac, k =3, trymax = 999)
stressplot(w_nmds)
w_nmds$stress

phy_plt1 <- data.frame(scores(w_nmds), meta) %>% 
  mutate(ash = factor(ash, levels = c("none", "low", "high"), labels = c("None", "Low", "High")))

time <- c("baseline", "year_1", "year_2")
new_labs <- c("2018", "2019", "2020")
names(new_labs) <-  time

phy_nmds <- ggplot(
  phy_plt1, aes(x = NMDS1,
                y = NMDS2,
                fill = sprayed
  )) +
  geom_point(aes(shape = ash), size  = 2) +
  facet_grid(~ time, labeller = labeller(time = new_labs), scales = "free") +
  stat_ellipse() +
  theme_bw(base_size = 18) +
  scale_shape_manual(name = "Ash treatment", 
                     values = c(21,22, 25)) +
  scale_fill_viridis_d(name= "           Herbicide",
                       labels = c("FALSE"= "Unsprayed",
                                  "TRUE"="Sprayed"),
                       begin = 0.2, end = 0.8) +
  guides(fill = guide_legend(override.aes = list(shape =21, size = 4)),
         shape = guide_legend(override.aes = list( size = 4))) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = 0),
        legend.key.size =  unit(5, "pt"),
        legend.key.spacing =  unit(5, "pt"),
        legend.spacing = unit(30, "pt"),
        legend.box.just = "left",
        legend.title = element_text(size =15, face = "bold"),
        text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid = element_blank()) +
  labs(title = "NMDS of weighed Unifrac distance")
phy_nmds

# Combine plots
nmds_comb <- veg_abun_nmds / phy_nmds + plot_layout(
  guides = "collect", axes = "collect") & theme(legend.position = "bottom")

ggsave("Plots/Beta_diversity_NMDS.png", nmds_comb, width = 14, height = 10, dpi = 800)

# Betapart analysis -------------------------------------------------------

soren <- veg %>% 
  select(spotted_knapweed:meadow_brome) %>% 
  decostand(., method = "pa") %>% 
  beta.pair(., index.family = "sorensen")

# Species turnover component

set.seed(1111111)
mod_ado1 <- adonis2(soren$beta.sim ~ sprayed*ash*time, by = "term",  meta)
mod_ado1_df <- data.frame(mod_ado1)
write.csv(mod_ado1_df, "Output/Species_turnover_PERMANOVA.csv", row.names = TRUE)

#Pairwise adonis 
rough_meta <- meta %>% 
  mutate(trt = paste(sprayed, time, sep = "_"))


set.seed(20260525)
pair_ado_turn <- pairwise.adonis2(soren$beta.sim ~ trt, by = "term",  rough_meta)
write.csv(pair_ado_turn$FALSE_year_1_vs_TRUE_year_1,
          "Output/Species_turnover_PERMANOVA_year2019.csv",
          row.names = TRUE)
write.csv(pair_ado_turn$FALSE_year_2_vs_TRUE_year_2,
          "Output/Species_turnover_PERMANOVA_year2020.csv",
          row.names = TRUE)

set.seed(202605252)
pair_ado_nest <- pairwise.adonis2(soren$beta.sne ~ trt, by = "term",  rough_meta)
write.csv(pair_ado_nest$FALSE_year_1_vs_TRUE_year_1,
          "Output/Species_nest_PERMANOVA_year2019.csv",
          row.names = TRUE)
write.csv(pair_ado_nest$FALSE_year_2_vs_TRUE_year_2,
          "Output/Species_nest_PERMANOVA_year2020.csv",
          row.names = TRUE)


set.seed(11)
nmds_sim <- metaMDS(soren$beta.sim, k =3,  trymax = 999)

plot_df_sim <- data.frame(scores(nmds_sim), meta)

turnover_plot <- ggplot(
  plot_df_sim, aes(x = NMDS1,
                   y = NMDS2,
                   fill = sprayed
  )) +
  geom_point(aes(shape = ash), size  = 2) +
  facet_grid(~ year) +
  stat_ellipse() +
  theme_bw(base_size = 18) +
  scale_shape_manual(name = "Ash treatment", 
                     values = c(21, 22, 25),
                     label = c("None", "Low", "High")) +
  scale_fill_viridis_d(name = "           Herbicide",
                       labels = c("FALSE" = "Unsprayed",
                                  "TRUE" = "Sprayed"),
                       begin = 0.2, end = 0.8) +
  guides(fill = guide_legend(override.aes = list(shape =21, size = 4)),
         shape = guide_legend(override.aes = list(size = 4))) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = 0),
        legend.key.size =  unit(5, "pt"),
        legend.key.spacing =  unit(5, "pt"),
        legend.spacing = unit(30, "pt"),
        legend.box.just = "left",
        legend.title = element_text(size =15, face = "bold"),
        text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0)) +
  labs(title = "Turnover")

# Species nestedness component
set.seed(121)
nmds_sne <- metaMDS(soren$beta.sne, k =3, trymax = 999, autotransform =TRUE)

set.seed(221)
mod_ado2 <- adonis2(soren$beta.sne ~ sprayed*ash*time, by = "term",  meta)
mod_ado2_df <- data.frame(mod_ado2)
write.csv(mod_ado2_df, "Output/Species_nestedness_PERMANOVA.csv", row.names = TRUE)

plot_df_sne <- data.frame(scores(nmds_sne), meta)

nestedness_plot <- ggplot(
  plot_df_sne, aes(x = NMDS1,
                   y =NMDS2,
                   fill = sprayed
  )) +
  geom_point(aes(shape = ash), size  = 2) +
  facet_grid(~ year) +
  stat_ellipse() +
  theme_bw(base_size = 18) +
  scale_shape_manual(name = "Ash treatment", 
                     values = c(21, 22, 25),
                     label = c("None", "Low", "High")) +
  scale_fill_viridis_d(name = "           Herbicide",
                       labels = c("FALSE" = "Unsprayed",
                                  "TRUE" = "Sprayed"),
                       begin = 0.2, end = 0.8) +
  guides(fill = guide_legend(override.aes = list(shape =21, size = 4)),
         shape = guide_legend(override.aes = list( size = 4))) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = 0),
        legend.key.size =  unit(5, "pt"),
        legend.key.spacing =  unit(5, "pt"),
        legend.spacing = unit(30, "pt"),
        legend.box.just = "left",
        legend.title = element_text(size =15, face = "bold"),
        text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0)) +
  labs(title = "Nestedness")

# Combine turnover and nestedness figures
comb_nest_turn <- turnover_plot/nestedness_plot + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom", panel.spacing = unit(1.5, "lines"))
ggsave("Plots/turnover_nestedness.png", comb_nest_turn, height = 10, width = 14, dpi = 800)

# Analysis of plant coverage ----------------------------------------------

df_plant1 <- df_plant %>% 
  mutate(spot_knap = (spotted_knapweed - min(spotted_knapweed))/(max(spotted_knapweed)- min(spotted_knapweed)),
         spot_total = spotted_knapweed/ rowSums(select(.,spotted_knapweed:meadow_brome)),
         cheat_total = cheatgrass/ rowSums(select(.,spotted_knapweed:meadow_brome)))

plnt_df <- df_plant1 %>% 
  select(time,sprayed, ash,site, cheat_total, spot_total, year) %>% 
  mutate(sprayed = as.factor(sprayed),
         ash = as.factor(ash),
         inv_tot = cheat_total + spot_total)

# Model for relative proportion of dominant invasive plants
inv_res <- df_plant1 %>% 
  select(year, site, cheat_total, spot_total) %>% 
  pivot_longer(-c(year,site)) %>% 
  mutate(year = as.factor(year))

mod_inv <- glmmTMB::glmmTMB(value ~ name * year + (1|site), family = "ordbeta", inv_res)
car::Anova(mod_inv) %>% 
  broom::tidy() %>% 
  write.csv("Output/Relative_proportion_glmmTMB.csv", row.names = FALSE)
DHARMa::simulateResiduals(mod_inv, plot = T)

# Plot
inv_prop_plot <- emmeans(mod_inv, ~ name | year, type = "response") %>%
  multcomp::cld(Letters = letters) %>%
  as.data.frame() %>%
  mutate(.group = str_trim(.group, side = "both")) %>%
  ggplot(
    .,
    aes(
      x = year,
      y = response,
      color = name,
      ymin = response - SE,
      ymax = response + SE,
      label = .group,
      group = name,
    )
  ) +
  geom_point(position = position_dodge(width = 0.05),
             size = 5,
             show.legend = FALSE)+
  geom_errorbar(width = 0,
                linewidth = 2,
                position = position_dodge(width = 0.05)) +
  geom_text(nudge_y = -0.005,
            nudge_x = -0.075,
            size = 5, 
            show.legend = FALSE,
            fontface = "bold") +
  geom_vline(xintercept = 1.25, linetype = "dotdash", color = "grey40") +
  annotate(geom = "text", label = "Herbicide\napplied", x = 1.45, y = 0.17, size = 5, color = "grey40") +
  geom_line() +
  scale_x_discrete(expand = expansion(mult = 0, add = 0.2)) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0.15,0.6),
                     breaks = seq(0.2,0.6, 0.1)) +
  theme_bw(base_size = 18) +
  scale_color_viridis_d(
    name = NULL,
    label = c("*Bromus tectorum*", "*Centaurea stoebe*"),
    begin = 0.2,
    end = 0.8
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.9),
    legend.text = element_markdown(face = "bold"),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold")
  ) +
  labs(y = "Relative proportion")
inv_prop_plot

ggsave("Plots/inv_prop_plot.png", inv_prop_plot, height = 6, width = 8, dpi = 800)


# Map data ----------------------------------------------

# Get the BC Boundary layer
bc_boundary <- bc_bound(ask = FALSE, force = TRUE) %>% 
  st_union() %>% 
  st_transform(3153) %>% 
  st_sf() %>% 
  mutate(name = "British Columbia")

# Get the municipal boundary for the City of Merritt
mer <- municipalities() %>% 
  filter(ADMIN_AREA_ABBREVIATION == "Merritt") %>% 
  select(name = ADMIN_AREA_ABBREVIATION) %>% 
  st_transform(3153)

# Laurie Guichon Memorial Grasslands Interpretive Site
lgm <- bcdc_query_geodata("263338a7-93ee-49c1-83e8-13f0bde70833", crs = 3153) %>% 
  filter(PROJECT_NAME == "GRASSLANDS INTERPRETIVE AREA") %>% 
  collect() %>% 
  mutate(name = "Laurie Guichon Memorial Grasslands Interpretive Site") %>% 
  select(name)

# Write to geopackage for importing into QGIS
sites <- rbind(mer, lgm)
st_write(sites, "Output/BC_Data_Layers.gpkg", layer = "general_sites", delete_layer = TRUE)
st_write(bc_boundary, "Output/BC_Data_Layers.gpkg", layer = "bc", delete_layer = TRUE)
