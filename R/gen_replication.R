
gen_patients <- function(prof_a, time_steps, pars_qol, pars_tte, pars_tte_slopes, list_studies, seed = 11667) {
  set.seed(seed)
  
  dat_gen <- prof_a %>% 
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
    ) %>% 
    mutate(
      SID = factor(SID, list_studies$SID),
      Agp = cut(Age, c(18, seq(30, 100, by = 10), right = F))
    ) %>% 
    arrange(SID, PID, visit)
  
  return(dat_gen)
}

