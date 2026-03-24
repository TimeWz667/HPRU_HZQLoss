library(tidyverse)
library(rstan)




## Background population

prof_a <- read_csv(here::here("data", "inputs_studies", "tab_agp_by_studies.csv"))
prof_t <- read_csv(here::here("data", "inputs_studies", "tab_t_by_studies.csv"))


prof_a


pars_tte_slopes <- read_csv(here::here("data", "inputs_r1_post", "slopes_tte.csv"))
pars_tte_slopes


studies <- prof_a %>% pull(SID) %>% unique()


## Simulation
set.seed(11667)

data_gen <- lapply(studies, \(sid) {
  ## Parameters
  
  
  if (T) {  # (sid == "Rampakakis et al. 2017 Costa Rica") {
    stats_tte <- read_csv(here::here("data", "inputs_r1_post", "stats_tte_ph.csv"))
    slopes <- pars_tte_slopes %>% filter(SID == "Combined")
    stats_qol <- read_csv(here::here("data", "inputs_r1_post", "stats_qol_b_orig.csv")) %>% 
      mutate(Var = paste0(Model, "_", Var))
  } else {
    stats_tte <- read_csv(here::here("data", "inputs_r1_post", paste0("stats_tte_sub_", sid, ".csv")))
    slopes <- pars_tte_slopes %>% filter(SID == sid)
    stats_qol <- read_csv(here::here("data", "inputs_r1_post", paste0("stats_qol_b_orig_sub_", sid, ".csv"))) %>% 
      mutate(Var = paste0(Model, "_", Var))
  }
  
  sims_pars <- prof_a %>% 
    filter(SID == sid) %>% 
    filter(!is.na(Agp)) %>% 
    select(Agp, n) %>% 
    uncount(n) %>% 
    tidyr::extract(Agp, c("A0", "A1"), "\\[(\\d+),(\\d+)\\)", convert = T, remove = F)  %>% 
    mutate(
      PID = 1:n(),
      A0 = pmax(A0, 18),
      Age = floor(runif(n(), A0, A1))
    ) %>% 
    select(-c("A0", "A1")) %>% 
    crossing(
      bind_rows(
        stats_tte %>% 
          select(Var, mean, sd) %>% 
          filter(Var != "lp__"),
        stats_qol %>% 
          select(Var, mean, sd) %>% 
          filter(Var != "lp__")
      )
    ) %>% 
    mutate(value = rnorm(n(), mean, sd)) %>% 
    select(-c(mean, sd)) %>% 
    pivot_wider(names_from = Var, values_from = value) %>% 
    bind_cols(
      slopes %>%
        select(Var, m) %>% 
        pivot_wider(names_from = Var, values_from = m)
    ) %>% 
    mutate(
      r0 = pmax(r0, 0.0001),
      ba1_sim = b0 + log(r0) * b_lr0 + rnorm(n(), 0, rsd)
    )
  
  
  sims <- sims_pars %>%
    mutate(
      rate = r0 * exp(Age * ba1_sim),
      disuti = rexp(n(), rate = rate)
    ) %>% 
    crossing(
      prof_t %>% filter(SID == sid) %>% 
        select(visit, t_l, t_u, Prop)
    ) %>% 
    mutate(
      SID = sid,
      time_points = floor(runif(n(), t_l, t_u + 1)),
      ti = time_points / 365.25,
      mis = 1 - rbinom(n(), 1, Prop),
      Q_0 = 1,
      Q_1 = pmin(rnorm(n(), C1_b0, C1_sigma), 1),
      Q_2 = pmin(rnorm(n(), C2_b0, C2_sigma), 1),
      PZ = PZ_b0 + PZ_bd15 * (ti <= 15 / 365.25) + PZ_bd30 * (ti <= 30 / 365.25),
      PZ = 1 / (1 + exp(-PZ)),
      PC1 = PC1_b0 + PC1_bd15 * (ti <= 15 / 365.25) + PC1_bd30 * (ti <= 30 / 365.25),
      PC1 = 1 / (1 + exp(-PC1))
    ) %>% 
    group_by(SID, PID) %>%  
    mutate(
      last = ti == max(ti < disuti),
      PZ = ifelse(last, 0, PZ),
      Prop_0 = PZ,
      Prop_1 = PC1 * (1 - PZ),
      Prop_2 = 1 - Prop_0 - Prop_1,
      seed = runif(n()),
      EQ5D = case_when(
        ti > disuti ~ Q_0, 
        seed < Prop_0 ~ Q_0,
        seed < Prop_0 + Prop_1 ~ Q_1,
        T ~ Q_2
      )
    ) %>% 
    ungroup() %>% 
    filter(!mis) %>% 
    select(
      SID, PID, Agp, Age, visit, ti, time_points, EQ5D, disuti, ba1_sim, r0, rate
    )
  
  sims
}) %>% 
  bind_rows()


data_gen %>% 
  summarise(disuti = mean(1 / rate, na.rm = T), .by = "Age") %>% 
  ggplot() +
  geom_line(aes(x = Age, y = disuti))



post_tte <- lapply(studies, \(sid) {
  po <- read_csv(here::here("posteriors", paste0("post_tte_sub_", sid, ".csv"))) %>% 
    mutate(SID = sid)
}) %>% 
  bind_rows()


post_tte  %>% 
  filter(SID != "Rampakakis et al. 2017 Costa Rica") %>% 
  ggplot() +
  geom_point(aes(x = ba1, y = r0)) +
  geom_point(data = data_gen %>% filter(visit == 1), 
             aes(x = ba1_sim, y = r0, colour = ifelse(ba1_sim > 0, "u", "l"))) +
  facet_wrap(.~SID)


post_tte %>% 
  filter(SID != "Rampakakis et al. 2017 Costa Rica") %>% 
  ggplot() +
  geom_density(aes(x = ba1)) +
  facet_wrap(.~SID)

data_gen %>% 
  filter(visit == 1) %>% 
  filter(SID != "Rampakakis et al. 2017 Costa Rica") %>% 
  ggplot() +
  geom_point(aes(x = ba1_sim, y = r0)) +
  facet_wrap(.~SID)
  

data_gen %>% 
  filter(visit == 1) %>% 
  ggplot() +
  geom_density(aes(x = disuti)) +
  facet_wrap(.~SID, scale = "free_x")


write_csv(data_gen, file = here::here("data", "generated", "data_generated_11667.csv"))






