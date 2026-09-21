source("Master_Code_Publication.R")
library(rstatix)

##no significant different in PD or SR for controls (i.e. absence of any treatment) across different years
## We can use this figure to show that the SR and PD in baseline and untreated sites are not significantly different
pd_veg %>% 
  filter(sprayed == FALSE & ash == "none") %>% 
  group_by(year) %>% 
  # summarise(mean_Pd = mean(PD),
  #           sd_Pd = sd(PD),
  #           mean_rich = mean(SR),
  #           sd_rich = sd(SR)) %>% 
  ggplot(., aes(x = as.factor(year), y = SR)) +
  geom_boxplot(fill = "steelblue") +
  stat_summary(geom = "point", fun = "mean", shape = 23, size = 5, fill = "black")+
  theme_bw() +
 scale_y_continuous(breaks = seq(0, 15, 1),
                    label =  seq(0, 15, 1)) +
  ggpubr::geom_pwc(ref.group = "2018",
                   label = "p.signif")

pd_veg %>% 
  filter(sprayed == FALSE & ash == "none") %>% 
  group_by(year) %>% 
  # summarise(mean_Pd = mean(PD),
  #           sd_Pd = sd(PD),
  #           mean_rich = mean(SR),
  #           sd_rich = sd(SR)) %>% 
  ggplot(., aes(x = as.factor(year), y = PD)) +
  geom_boxplot(fill = "steelblue") +
  stat_summary(geom = "point", fun = "mean", shape = 23, size = 5, fill = "black")+
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 1200, 100),
                     label =  seq(0, 1200, 100)) +
  ggpubr::geom_pwc(ref.group = "2018",
                   label = "p.signif")



##calculating pairwise dissimilarity between baseline and rest of the samples across timepoints
## we can use the same code to calculate pairwise distance for turnover, nestedness and unifrac
## that gives a more compelling effidence that betadiversity and partitioning is impacted by herbicide application
dist_analysis <- function(dat_mat){
dist_from_baseline <- dat_mat %>% 
 pull_lower_triangle() %>% 
  as.data.frame() %>% 
  pivot_longer(-rowname) %>% 
  filter(value != "") %>% 
  mutate(name_time = factor(name, levels = meta$id, labels = meta$year),
         rowname_time = factor(rowname, levels = meta$id, labels = meta$year),
         sprayed = factor(rowname, levels = meta$id, labels = meta$sprayed),
         value = as.numeric(value)) %>% 
  filter(rowname_time != "2018",
         name_time  == "2018") #%>% 
  # group_by(name,rowname_time,sprayed ) %>% 
  # summarise(mean_val = mean(value)) %>% ungroup()

str(dist_from_baseline)

# dist_baseline <- dat_mat %>% 
#   pull_lower_triangle() %>% 
#   as.data.frame() %>% 
#   pivot_longer(-rowname) %>% 
#   filter(value != "") %>% 
#   mutate(name_time = factor(name, levels = meta$id, labels = meta$year),
#          rowname_time = factor(rowname, levels = meta$id, labels = meta$year),
#          sprayed = factor(rowname, levels = meta$id, labels = meta$sprayed),
#          value = as.numeric(value)) %>% 
#   filter(rowname_time == "2018",
#          name_time  == "2018")
# 
# comb_dist <- bind_rows(dist_from_baseline,dist_baseline)

# m_dist <- lmer((mean_val) ~ rowname_time * sprayed + (1|name), 
#                    dist_from_baseline) 
# DHARMa::simulateResiduals(m_dist, plot = T)
# car::Anova(m_dist, test.statistic = "F")
  plt <- ggplot(dist_from_baseline, aes(x = sprayed, y = value)) +
  # stat_summary(geom = "pointrange", fun.data = "mean_sdl",
  #              linewidth = 1.5, size = 1.5, color = "black",
  #              shape = 23, fill = "steelblue", fun.args = list(mult = 1)) +
  geom_boxplot()+
  facet_grid( ~ rowname_time) +
    ylim(c(0, NA)) +
    ggpubr::geom_pwc(label = "p.adj.signif") +
    scale_x_discrete(labels = c("FALSE" = "Unsprayed",
                                "TRUE" = "Sprayed")) +
    theme_bw() +
    theme(text = element_text(color = "black", face = "bold", size = 14),
          axis.text = element_textbox(color = "black", face = "bold", size = 14)) #+
    # labs(x = NULL,
    #      y = "Mean Bray-Curtis dissimilarity from baseline")
    # 
    # 
  
  return(plt)
}
dist_bray_pair <- dist_analysis(bray_dist) +
  labs(y = "Pairwise Bray-Curtis\ndissimilarity from baseline", x = NULL)

unifrac_pair <- dist_analysis(weighted_unifrac) +
  labs(y = "Pairwise weighted-unifrac\ndistance from baseline", x = NULL)

turn_pair <- dist_analysis(soren$beta.sim) +
  labs(y = "Pairwise Sorensen distance \n(Turnover) from baseline", x = NULL) 

nest_pair <- dist_analysis(soren$beta.sne) +
  labs(y = "Pairwise Sorensen distance \n(Nestedness) from baseline", x = NULL) 


comb_pair_resp <- (dist_bray_pair + unifrac_pair) / (turn_pair + nest_pair)
  
ggsave("Plots/comb_pair_resp.png", comb_pair_resp, 
       height = 8, width = 12, dpi = 800, units = "in")

# Dispersion analysis -----------------------------------------------------

get_disp <- function(dist_mat) {
  meta_2019 <- meta %>% 
    filter(year == 2019)
  meta_2020 <-  meta %>% 
    filter(year == 2020)
  
  dist_2019 <-dist_mat %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    filter(rownames(.) %in% meta_2019$id) %>% 
    select(meta_2019$id) %>% 
    as.dist()
  
  set.seed(111)
  disp_2019 <- betadisper(dist_2019, meta_2019$sprayed, type = "centroid")
  anova(disp_2019)
  
  
  dist_2020 <- dist_mat %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    filter(rownames(.) %in% meta_2020$id) %>% 
    select(meta_2020$id) %>% 
    as.dist()
  
  set.seed(1111)
  disp_2020 <- betadisper(dist_2020, meta_2020$sprayed, type = "centroid")
  anova(disp_2020)
  
  return(list(anova(disp_2019), anova(disp_2020)))
}

get_disp(bray_dist)
get_disp(weighted_unifrac)
get_disp(soren$beta.sim)
get_disp(soren$beta.sne)



