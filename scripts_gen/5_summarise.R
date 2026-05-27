library(tidyverse)


source(here::here("R", "gen_analyses.R"))


## Set A runs
pars_tte <- read_csv(here::here("posteriors", "post_ph_tte.csv")) %>% 
  mutate(Key = 1:n()) %>% 
  select(- lp__)


pars_qol <- read_csv(here::here("posteriors", "post_qol_b_uk.csv"))


summ_r1 <- pars_qol %>% 
  select(- lp__) %>% 
  mutate(Key = 1:n(), .by = Model) %>% 
  pivot_longer(-c(Model, Key)) %>% 
  filter(!is.na(value)) %>% 
  mutate(name = paste0(Model, "_", name)) %>% 
  select(- Model) %>% 
  pivot_wider() %>% 
  left_join(pars_tte) %>% summary_proj(
  ages = seq(50, 90, 5),
  times = c(0, 14, 30, 60, 90)
)


root_src <- here::here("data", "inputs_for_gen")
pars_qol <- read_csv(here::here(root_src, "pars_qol_km.csv")) %>% 
  filter(SID == "ESP") %>% 
  select(- c(SID, Model))

pars_tte <- read_csv(here::here(root_src, "pars_tte_2gp.csv")) %>% 
  filter(SID == "ESP") %>% 
  select(- SID)

pars_tte_slopes <- read_csv(here::here(root_src, "pars_tte_slopes_2gp.csv")) %>% 
  filter(SID == "ESP") %>% 
  select(- SID)



boot_pars <- revive_pars(n_mc = 200, pars_qol, pars_tte, pars_tte_slopes)

summary_proj_gsk <- boot_pars %>% summary_proj(
  ages = seq(50, 90, 5),
  times = c(0, 14, 30, 60, 90)
)

summary_proj_gsk



pars_qol <- read_csv(here::here(root_src, "pars_qol.csv")) %>% 
  filter(SID == "ScottUK") %>% 
  select(- c(SID, Model))

pars_tte <- read_csv(here::here(root_src, "pars_tte_2gp.csv")) %>% 
  filter(SID == "ScottUK") %>% 
  select(- SID)

pars_tte_slopes <- read_csv(here::here(root_src, "pars_tte_slopes_2gp.csv")) %>% 
  filter(SID == "ScottUK") %>% 
  select(- SID)



boot_pars <- revive_pars(n_mc = 200, pars_qol, pars_tte, pars_tte_slopes)

summary_proj_r1 <- boot_pars %>% summary_proj(
  ages = seq(50, 90, 5),
  times = c(0, 14, 30, 60, 90)
)

summary_proj_r1





bind_rows(
  summary_proj_gsk$tab_age %>% mutate(Set = "GSK"),
  summary_proj_r1$tab_age %>% mutate(Set = "R1"),
  summary_proj_v1$tab_age %>% mutate(Set = "V1")
) %>% 
  ggplot(aes(x = Age)) + 
  geom_ribbon(aes(ymin = disuti_L, ymax = disuti_U, fill = Set), alpha = 0.3) + 
  geom_line(aes(y = disuti_M, colour = Set)) + 
  expand_limits(y = 0)


bind_rows(
  summary_proj_gsk$tab_age %>% mutate(Set = "GSK"),
  summary_proj_r1$tab_age %>% mutate(Set = "R1"),
  summary_proj_v1$tab_age %>% mutate(Set = "V1")
) %>% 
  ggplot(aes(x = Age)) + 
  geom_ribbon(aes(ymin = QL35_L, ymax = QL35_U, fill = Set), alpha = 0.3) + 
  geom_line(aes(y = QL35_M, colour = Set)) + 
  expand_limits(y = 0)



summary_proj_gsk$tab_time
summary_proj_r1$tab_time


