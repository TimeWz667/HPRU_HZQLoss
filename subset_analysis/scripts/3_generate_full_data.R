library(tidyverse)
library(rstan)


## Function

gen_patients <- function(prof, pars_qol, pars_tte, pars_tte_slopes) {
  data_gen <- prof %>% 
    uncount(n) %>% 
    mutate(
      PID = paste0(SID, "_", 1:n()),
      A0 = pmax(A0, 18),
      Age = floor(runif(n(), A0, A1))
    ) %>% 
    select(-c("Agp", "A0", "A1")) %>% 
    full_join(bind_rows(pars_tte, pars_qol) %>% select(- Model), by = "SID", relationship = "many-to-many")  %>% 
    mutate(value = rnorm(n(), mean, sd)) %>% 
    select(-c(mean, sd)) %>% 
    pivot_wider(names_from = Var, values_from = value) %>% 
    left_join(pars_tte_slopes) %>% 
    mutate(
      r0 = pmax(r0, 0.0001),
      ba1_sim = b0 + log(r0) * b_lr0 + rnorm(n(), 0, rsd),
      rate = r0 * exp(Age * ba1_sim),
      disuti = rexp(n(), rate = rate)
    ) %>% 
    crossing(ti = time_steps) %>% 
    arrange(SID, PID, ti) %>% 
    mutate(
      time_points = ti * 365.25,
      Q_0 = 1,
      Q_1 = pmin(rnorm(n(), C1_b0 + C1_bg * ti, C1_sigma), 1),
      Q_2 = pmin(rnorm(n(), C2_b0 + C2_bg * ti, C2_sigma), 1),
      PZ = PZ_b0 + PZ_bd15 * (ti <= 15 / 365.25) + PZ_bd30 * (ti <= 30 / 365.25) + PZ_bg,
      PZ = 1 / (1 + exp(-PZ)),
      PC1 = PC1_b0 + PC1_bd15 * (ti <= 15 / 365.25) + PC1_bd30 * (ti <= 30 / 365.25) + PC1_bg,
      PC1 = 1 / (1 + exp(-PC1))
    ) %>% 
    group_by(SID, PID) %>%  
    mutate(
      visit = 1:n(),
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
    select(
      SID, PID, Age, visit, ti, time_points, EQ5D, disuti, ba1_sim, r0, rate
    )
  
  return(data_gen)
}


## Background population
list_studies <- read_csv(here::here("data", "list_studies.csv")) %>% 
  select(study, SID, SIK)
list_studies


root_src <- here::here("data", "inputs_for_gen")



## Time steps 
time_steps <- (seq(0, by = 14, length.out = 20) + 1) / 365.25


## load Shared inputs

prof_a <- read_csv(here::here(root_src, "prof_a.csv"))

pars_qol <- read_csv(here::here(root_src, "pars_qol.csv"))

set.seed(11667)
seeds <- round(runif(100, 0, 1e6))

names(seeds) <- paste0("s", 1:length(seeds))



## Generate data
for (type in c("ind", "2gp")) {
  type <- glue::as_glue(type)
  print(type)
  pars_tte <- read_csv(here::here(root_src, "pars_tte_" + type + ".csv"))
  pars_tte_slopes <- read_csv(here::here(root_src, "pars_tte_slopes_" + type + ".csv"))
  
  root_tar <- here::here("out", "gen_" + type)
  
  dir.create(root_tar, showWarnings = F )
  
  
  for (v in names(seeds)) {
    set.seed(seeds[v])
    print(seeds[v])
    
    data_gen <- gen_patients(prof_a, pars_qol, pars_tte, pars_tte_slopes) %>% 
      mutate(
        SID = factor(SID, list_studies$SID),
        Agp = cut(Age, c(18, seq(30, 100, by = 10), right = F))
      ) %>% 
      arrange(SID, PID, visit)
    
    
    write_csv(data_gen, file = here::here(root_tar, "generated_" + glue::as_glue(v) + ".csv"))
  }
}
