library(tidyverse)
library(cowplot)
library(magick)
library(ggpubr)

theme_set(theme_bw())


q_0 <- read_csv(here::here("subset_analysis", "docs", "tabs", "summary_qloss_0_uk.csv"))
q_pn <- read_csv(here::here("docs", "gen", "tabs", "tab_age_pn.csv"))

q_ph <- read_csv(here::here("docs", "gen", "tabs", "tab_age.csv"))



tab_qloss <- bind_rows(
  q_ph %>% mutate(Base = "Perfect health"),
  q_pn %>% mutate(Base = "Population norm"),
  q_0 %>% 
    select(Age, disuti_M = dur_M, disuti_L = dur_L, disuti_U = dur_U, starts_with("QL")) %>% 
    select(- ends_with("_A")) %>% 
    mutate(Base = "Baseline")
) %>% 
  filter(Age %in% seq(50, 100, 10)) %>% 
  relocate(Base, Age) %>% 
  mutate(
    Disuti = sprintf("%.0f (%.0f - %.0f)", disuti_M * 365, disuti_L * 365, disuti_U * 365),
    QL00 = sprintf("%.3f (%.3f - %.3f)", QL00_M, QL00_L, QL00_U),
    QL15 = sprintf("%.3f (%.3f - %.3f)", QL15_M, QL15_L, QL15_U),
    QL35 = sprintf("%.3f (%.3f - %.3f)", QL35_M, QL35_L, QL35_U)
  ) %>% 
  select(Base, Age, Disuti:QL35)

tab_qloss

write_csv(tab_qloss, here::here("docs", "gen", "tabs", "summary_qloss_uk.csv"))




load(here::here("subset_analysis", "docs", "figs", "gs_tte.rdata"))



sims_dur <-lapply(1:100, \(i) {
  tar_read(post_tte, i) %>% 
    crossing(Age = seq(20, 95, 5)) %>% 
    mutate(
      dur = 1 / r0 * exp(-Age * ba1)
    ) %>% 
    summarise(
      M = median(dur),
      L = quantile(dur, 0.025),
      U = quantile(dur, 0.975), .by = Age
    )
}) %>% 
  bind_rows() %>% 
  summarise(across(everything(), mean), .by = Age)


gs$g_tte_bind2


g_tte <- ggplot(sims_dur, aes(x = Age)) +
  geom_ribbon(aes(ymin = L, ymax = U), alpha = 0.3) +
  geom_line(aes(y = M)) +
  scale_y_continuous("Disutility period, months", breaks = seq(0, 2.5, 0.5), 
                     labels = scales::number_format(scale = 12, suffix = " mo."),
                     limits = c(0, 1.5)) +
  scale_x_continuous("Age", limits = c(0, 100))



g_bind <- ggpubr::ggarrange(
  gs$g_imp + labs(subtitle = "(A) Survival curve, examplar data") + theme(legend.position = c(1, 1), legend.justification = c(1.1, 1.1)),
  gs$g_tte_data + labs(subtitle = "(B) Simulated values") + 
    scale_y_continuous("Disutility period, months", breaks = seq(0, 2.5, 0.5), 
                       labels = scales::number_format(scale = 12, suffix = " mo."),
                       limits = c(0, 1.5)) +
    theme(legend.position = c(1, 1), legend.justification = c(1.1, 1.1)),
  g_tte + labs(subtitle = "(C) Predicted values"),
  ncol = 1
)


ggsave(g_bind, file = here::here("docs", "gen", "figs", "g_tte_bind.png"), width = 9, height = 12)



