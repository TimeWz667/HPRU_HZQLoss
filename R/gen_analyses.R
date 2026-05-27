

mlu <- list(
  M = median,
  L = \(x) quantile(x, 0.025),
  U = \(x) quantile(x, 0.975)
)


revive_pars <- function(n_mc = 2000, pars_qol, pars_tte, pars_tte_slopes) {
  tibble(Key = 1:n_mc) %>% 
    cross_join(bind_rows(pars_qol, pars_tte)) %>% 
    mutate(value = rnorm(n(), mean, sd)) %>% 
    select(-c(mean, sd)) %>% 
    pivot_wider(names_from = Var, values_from = value) %>% 
    bind_cols(pars_tte_slopes) %>% 
    mutate(
      r0 = pmax(r0, 0.000001),
      ba1 = b0 + log(r0) * b_lr0 + rnorm(n(), 0, rsd)
    )
}

summary_proj <- function(boot_pars, 
                         ages = seq(50, 90, 5),times = c(0, 14, 30, 60, 90),
                         dt = 0.01, ti_until = 5) {
  
  tab_tte <- tibble(Age = ages) %>% 
    cross_join(boot_pars) %>% 
    mutate(
      rate = r0 * exp(Age * ba1),
      disuti = 1 / rate
    ) %>% 
    summarise(
      across(c(disuti, disuti), mlu),
      .by = c(Age)
    )
  
  tab_qol <- tibble(Time = times) %>% 
    cross_join(boot_pars) %>% 
    mutate(
      ti = Time / 365.25,
      Q_0 = 1,
      Q_1 = pmin(C1_b0, 1),
      Q_2 = pmin(C2_b0, 1),
      PZ = PZ_b0 + PZ_bd15 * (ti <= 15 / 365.25) + PZ_bd30 * (ti <= 30 / 365.25),
      PZ = 1 / (1 + exp(-PZ)),
      PC1 = PC1_b0 + PC1_bd15 * (ti <= 15 / 365.25) + PC1_bd30 * (ti <= 30 / 365.25),
      PC1 = 1 / (1 + exp(-PC1)),
      QE = (1 - PZ) * (PC1 * Q_1 + (1 - PC1) * Q_2),
      P_C1 = (1 - PZ) * PC1,
      P_C2 = 1 - PZ - P_C1,
      QL = 1 - QE
    ) %>% 
    summarise(
      across(c(Q_1, Q_2, PZ, P_C1, P_C2, QE, QL), mlu),
      .by = c(Time)
    )
  
  
  tab_ql <- tibble(Age = ages) %>% 
    cross_join(boot_pars) %>% 
    mutate(
      rate = r0 * exp(Age * ba1)
    ) %>% 
    crossing(ti = seq(0, ti_until, dt)) %>% 
    mutate(
      AgeT = floor(Age + ti),
      p_hz = (exp(-rate * ti) - exp(-rate * (ti + dt))) / rate,
      p_health = dt - p_hz,
      Q_0 = 1,
      Q_1 = C1_b0,
      Q_2 = C2_b0,
      PZ = PZ_b0 + PZ_bd15 * (ti <= 15 / 365.25) + PZ_bd30 * (ti <= 30 / 365.25),
      PZ = 1 / (1 + exp(-PZ)),
      PC1 = PC1_b0 + PC1_bd15 * (ti <= 15 / 365.25) + PC1_bd30 * (ti <= 30 / 365.25),
      PC1 = 1 / (1 + exp(-PC1)),
      Prop_0 = PZ,
      Prop_1 = PC1 * (1 - PZ),
      Prop_2 = 1 - Prop_0 - Prop_1,
      ql = Q_0 * Prop_0 + Q_1 * Prop_1 + Q_2 * Prop_2,
      ql = (1 - ql) * p_hz,
      ql15 = ql * exp(- 0.015 * (AgeT - Age)),
      ql35 = ql * exp(- 0.035 * (AgeT - Age))
    ) %>% 
    summarise(
      QL00 = sum(ql),
      QL15 = sum(ql15),
      QL35 = sum(ql35), .by = c("Age", "Key")
    ) %>% 
    summarise(
      across(starts_with("QL"), mlu), .by = "Age"
    )
  
  return(list(
    tab_age = tab_tte %>% left_join(tab_ql, by = "Age"),
    tab_time = tab_qol
  ))
}


summarise_qls <- function(pars_qol, pars_tte, ages = seq(50, 90, 5), times = c(0, 14, 30, 60, 90), batch) {
  summ <- pars_qol %>% 
    select(- lp__) %>% 
    mutate(Key = 1:n(), .by = Model) %>%  
    pivot_longer(-c(Model, Key)) %>% 
    filter(!is.na(value)) %>% 
    mutate(name = paste0(Model, "_", name)) %>%
    select(- Model) %>% 
    pivot_wider() %>% 
    left_join(
      pars_tte %>% 
        mutate(Key = 1:n()), by = "Key"
    ) %>% summary_proj(
      ages = ages,
      times = times
    )
  
  summ <- lapply(summ, \(x) x %>% mutate(Batch = batch))
  
  return(summ)
}






