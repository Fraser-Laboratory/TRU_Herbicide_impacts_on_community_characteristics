library(rstatix)




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
         name_time  == "2018") %>% 
  group_by(name,rowname_time,sprayed ) %>% 
  summarise(mean_val = mean(value)) %>% ungroup()

str(dist_from_baseline)


m_dist <- lmer((mean_val) ~ rowname_time * sprayed + (1|name), 
                   dist_from_baseline) 
DHARMa::simulateResiduals(m_dist, plot = T)
car::Anova(m_dist, test.statistic = "F")
  plt <- ggplot(dist_from_baseline, aes(x = sprayed, y = mean_val)) +
  # stat_summary(geom = "pointrange", fun.data = "mean_sdl",
  #              linewidth = 1.5, size = 1.5, color = "black",
  #              shape = 23, fill = "steelblue", fun.args = list(mult = 1)) +
  geom_boxplot()+
  facet_grid(~ rowname_time) +
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
dist_analysis(bray_dist)
dist_analysis(weighted_unifrac)
dist_analysis(soren$beta.sor)
dist_analysis(soren$beta.sne)

