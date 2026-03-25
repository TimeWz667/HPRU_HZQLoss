library(tidyverse)
library(rstan)




## Background population

list_studies <- read_csv(here::here("data", "list_studies.csv")) %>% 
  select(study, SID, SIK)
list_studies


prof_a <- read_csv(here::here("data", "inputs_studies", "tab_agp_by_studies.csv"))
prof_a <- prof_a %>% rename(study = SID) %>% left_join(list_studies) %>% select(-c(study, SIK)) %>% 
  mutate(SID = factor(SID, list_studies$SID))

prof_t <- read_csv(here::here("data", "inputs_studies", "tab_t_by_studies.csv"))
prof_t <- prof_t %>% rename(study = SID) %>% left_join(list_studies) %>% select(-c(study, SIK)) %>% 
  mutate(SID = factor(SID, list_studies$SID))


prof_a
prof_t


pars_tte_slopes <- read_csv(here::here("data", "inputs_r1_post", "stats_tte_slopes.csv"))
pars_tte_slopes %>% 
  pivot_wider(names_from = Var, values_from = c(m, l, u)) %>%
  filter(SID != "RamCRI") %>% 
  ggplot() +
  geom_pointrange(aes(x = m_b0, y = m_b_lr0, ymin = l_b_lr0, ymax = u_b_lr0, colour = SID)) +
  geom_pointrange(aes(x = m_b0, xmin = l_b0, xmax = u_b0, y = m_b_lr0, colour = SID)) 


studies <- list_studies$SID


stats_tte0 <- read_csv(here::here("data", "inputs_r1_post", "stats_ph_tte.csv"))

stats_qol0 <- read_csv(here::here("data", "inputs_r1_post", "stats_qol_b_orig.csv")) %>% 
  mutate(Var = paste0(Model, "_", Var))

stats_qol_bg0 <- read_csv(here::here("data", "inputs_r1_post", "stats_qol_b_bg_orig.csv")) %>% 
  mutate(Var = paste0(Model, "_", Var)) %>% 
  tidyr::extract(Var, "SIK", "\\[(\\d+)\\]", convert = T, remove = F) %>% 
  mutate(Var = gsub("\\[(\\d+)\\]", "", Var))

## Simulation

seeds <- c(
  s1 = 1828, 
  s2 = 78197, 
  s3 = 764, 
  s4 = 83877, 
  s5 = 3045,
  s6 = 11667
)


for (v in names(seeds)) {
  set.seed(seeds[v])
  print(seeds[v])
  
  
  
  data_gen <- lapply(studies, \(sid) {
    ## Parameters
    sik <- list_studies %>% filter(SID == sid) %>% pull(SIK)
    print(sid)
    print(sik)
    if (T) {  # (sid == "Rampakakis et al. 2017 Costa Rica") {
      stats_tte <- stats_tte0
      slopes <- pars_tte_slopes %>% filter(SID == "Combined")
      #stats_qol <- stats_qol0 %>% bind_rows(stats_qol_bg0 %>% filter(SIK == sik) %>% select(- SIK))
      
      stats_qol <- bind_rows(
        read_csv(here::here("data", "inputs_r1_post", paste0("stats_qol_b_orig_sub_", sid, ".csv"))) %>% 
          mutate(Var = paste0(Model, "_", Var)),
        read_csv(here::here("data", "inputs_r1_post", paste0("stats_qol_b_bg_orig_sub_", sid, ".csv"))) %>% 
          mutate(Var = gsub("\\[(\\d+)\\]", "", paste0(Model, "_", Var)))
      )
    } else {
      stats_tte <- read_csv(here::here("data", "inputs_r1_post", paste0("stats_tte_sub_", sid, ".csv")))
      slopes <- pars_tte_slopes %>% filter(SID == sid)
      stats_qol <- bind_rows(
        read_csv(here::here("data", "inputs_r1_post", paste0("stats_qol_b_orig_sub_", sid, ".csv"))) %>% 
          mutate(Var = paste0(Model, "_", Var)),
        read_csv(here::here("data", "inputs_r1_post", paste0("stats_qol_b_bg_orig_sub_", sid, ".csv"))) %>% 
          mutate(Var = gsub("\\[(\\d+)\\]", "", paste0(Model, "_", Var)))
      )
    }
    
    sims_pars <- prof_a %>% 
      filter(SID == sid) %>% 
      filter(!is.na(Agp)) %>% 
      select(Agp, n) %>% 
      uncount(n) %>% 
      tidyr::extract(Agp, c("A0", "A1"), "\\[(\\d+),(\\d+)\\)", convert = T, remove = F)  %>% 
      mutate(
        PID = paste0(sid, "_", 1:n()),
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
        Q_1 = pmin(rnorm(n(), C1_b0 + C1_bg * ti, C1_sigma), 1),
        Q_2 = pmin(rnorm(n(), C2_b0 + C2_bg * ti, C2_sigma), 1),
        PZ = PZ_b0 + PZ_bd15 * (ti <= 15 / 365.25) + PZ_bd30 * (ti <= 30 / 365.25) + PZ_bg,
        PZ = 1 / (1 + exp(-PZ)),
        PC1 = PC1_b0 + PC1_bd15 * (ti <= 15 / 365.25) + PC1_bd30 * (ti <= 30 / 365.25) + PC1_bg,
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
  
  
  write_csv(data_gen, file = here::here("data", "generated", "data_generated_" + glue::as_glue(v) + ".csv"))
  
}







