library(tidyverse)
library(tidybayes)


theme_set(theme_bw())


source("R/misc.R")
source("R/data_processing.R")
## Vis qol, Figure 3


pars_tte <- tar_read(pars_tte_ph)
pars_qol <- tar_read(pars_qol)

dat_qol <- read_csv(here::here("data", "processed", "qol_uk_jitter.csv")) %>% 
  mutate(SID = "Set A", PID = 1:n())



label_cluster <- function(dat_qol) {
  dat_qol <- dat_qol %>% filter(!is.na(Age))
  ds_nz <- dat_qol %>% filter(EQ5D < 1)
  ds_z <- dat_qol %>% filter(EQ5D >= 1)
  clu <- kmeans(ds_nz$EQ5D, 2)
  
    
    if (clu$centers[1] > clu$centers[2]) {
      ds_nz <- ds_nz %>% mutate(cluster = clu$cluster)
    } else {
      ds_nz <- ds_nz %>% mutate(cluster = 3 - clu$cluster)
    }
  bind_rows(
    ds_nz, ds_z %>% mutate(cluster = 0)
  ) %>% 
    arrange(SID, PID, ti) %>% 
    mutate(visit = 1:n(), .by = PID)
}




fix_eff <- lapply(1:100, \(i) {
  post_qol <- tar_read(post_qol, i)
  
  post_qol %>% 
    filter(Model %in% c("C1", "C2")) %>% 
    group_by(Model)  %>% 
    summarise(
      across(b0, amlu)
    )
}) %>% 
bind_rows() %>% 
  group_by(Model) %>% 
  summarise_all(mean)

  
raw_eff <- lapply(1:100, \(i) {
  
  data_gen <- tar_read(data_gen_pn, i) %>% format_qol() %>%  label_cluster()

  data_gen %>% 
    filter(visit == 1) %>% 
    group_by(cluster) %>% 
    summarise(
      across(EQ5D, amlu)
    )
}) %>% 
  bind_rows() %>% 
  group_by(cluster) %>% 
  summarise_all(mean)


raw_prop_eff <- lapply(1:100, \(i) {
  
  data_gen <- tar_read(data_gen_pn, i) %>% format_qol() %>%  label_cluster()
  
  d90 <- data_gen %>% 
    mutate(Intv = ifelse(ti < 90 / 365.25, "< 90", ">= 90")) %>% 
    summarise(
      N = n(), .by = c(cluster, Intv)
    ) %>% 
    mutate(
      Prop = N / sum(N), .by = c(Intv)
    )
  
  d30 <- data_gen %>% 
    mutate(
      Intv = case_when(
        ti < 30 / 365.25 ~ "0 mo.",
        ti < 90 / 365.25 ~ "1~3 mo.",
        ti < 180 / 365.25 ~ "4~6 mo.",
        T ~ "7 mo.~"
      )
    ) %>% 
    summarise(
      N = n(), .by = c(cluster, Intv)
    ) %>% 
    mutate(
      Prop = N / sum(N), .by = c(Intv)
    )
  
  bind_rows(d90, d30)
}) %>% 
  bind_rows() %>% 
  group_by(cluster, Intv) %>% 
  summarise_all(mean)




data_gen_ph %>% 
  filter(ti < disuti) %>% 
  label_cluster() %>% 
  mutate(cluster = as.factor(cluster)) %>% 
  filter(EQ5D < 1) %>% 
  filter(EQ5D >= -0.5) %>% 
  ggplot() +
  stat_halfeye(aes(x = EQ5D, y = SID, fill = cluster), density = "unbounded", alpha = 0.2) +
  scale_fill_discrete("", labels = c("1" = "mild discomfort", "2" = "severe discomfort")) +
  scale_colour_discrete("", labels = c("1" = "mild discomfort", "2" = "severe discomfort")) +
  scale_y_discrete("", labels = with(list_studies, {
    set_names(paste0(Cite, "\n", Setting), SID)
  }))





colour_sch <- c("0" =  "cadetblue1", "1" = "chartreuse2", "2" = "coral")


dat_qol %>% label_cluster()
dat_fix <- fix_eff %>% mutate(cluster = ifelse(Model == "C1", "1", "2"))


g_qol_dat <- dat_qol %>% label_cluster() %>% 
  mutate(cluster = factor(cluster, c("0", "1", "2"))) %>% 
  ggplot() +
  geom_point(aes(x = ti, y = EQ5D, colour = cluster)) + 
  scale_colour_manual("", labels = c("0" = "temporally well", "1" = "mild discomfort", "2" = "severe discomfort"), values = colour_sch) +
  scale_y_continuous("EQ-5D score", breaks = seq(-0.5, 1, 0.5)) +
  scale_x_continuous("Time since the onset of rash", breaks = seq(0, 1.25, 0.25), 
                     labels = scales::number_format(scale = 12, suffix = " mo."),
                     limits = c(0, 1.25)) 

g_qol_sim <- data_gen_ph %>% 
  filter(ti < disuti) %>% 
  label_cluster() %>% 
  mutate(cluster = factor(cluster, c("0", "1", "2"))) %>% 
  filter(EQ5D < 1) %>% 
  filter(EQ5D >= -0.5) %>% 
  ggplot() +
  geom_density(aes(y = EQ5D, fill = cluster, colour = cluster)) +
  geom_hline(data = fix_eff, aes(yintercept = b0_M)) +
  geom_hline(data = fix_eff, aes(yintercept = b0_L), linetype = 2) +
  geom_hline(data = fix_eff, aes(yintercept = b0_U), linetype = 2) +
  geom_text(data = fix_eff, aes(y = b0_L, x = 0, label = sprintf("Fixed effect: %.3f (%.3f, %.3f)", b0_M, b0_L, b0_U)), hjust = -0.2, vjust = 1.3) +
  scale_fill_manual("", labels = c("1" = "mild discomfort", "2" = "severe discomfort"), values = colour_sch) +
  scale_colour_manual("", labels = c("1" = "mild discomfort", "2" = "severe discomfort"), values = colour_sch) +
  scale_y_continuous("EQ-5D score", breaks = seq(-0.5, 1, 0.5)) +
  scale_x_continuous("Density") +
  theme(legend.position  = c(1, 0.2), legend.justification = c(1.2, 0))


g_qol_prop <- raw_prop_eff %>% 
  filter(Intv %in% c("0 mo.", "1~3 mo.", "4~6 mo.", "7 mo.~")) %>% 
  mutate(cluster = factor(cluster, c("0", "1", "2"))) %>% 
  ggplot() +
  geom_bar(aes(x = Intv, y = Prop, fill = cluster), stat = "identity", position = "dodge", width = 0.8) +
  scale_x_discrete("Time since the onset of rash") +
  scale_y_continuous("Proportion in cluster, %", labels = scales::percent) +
  scale_fill_manual("", labels = c("0" = "temporally well", "1" = "mild discomfort", "2" = "severe discomfort"), values = colour_sch) +
  theme(legend.position = c(0, 1), legend.justification = c(-0.2, 1.2))


g_qol_sim
g_qol_prop

g_qol <- ggpubr::ggarrange(g_qol_dat + labs(subtitle = "(A)"), 
                  g_qol_sim + labs(subtitle = "(B)"), 
                  g_qol_prop + labs(subtitle = "(C)"), 
                  nrow = 2, ncol = 2, common.legend = T, legend = "bottom", widths = c(5, 3) , heights = c(5, 3))




ggsave(g_qol, filename = here::here("docs", "gen", "figs", "g_qol.png"), width = 12, height = 11)

raw_eff  
raw_prop_eff

fix_eff




