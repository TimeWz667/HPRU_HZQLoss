
sim_pre_ql_k <- function(age0 = 50, pars, pars_demo) {
  n_sim <- max(pars$Key)
  
  crossing(Key = 1:n_sim, Age = age0, ti = 0:(100 - age0)) %>% 
    mutate(
      AgeT = floor(Age + ti)
    ) %>% 
    left_join(pars, by = "Key") %>% 
    left_join(pars_demo %>% rename(AgeT = Age), by = "AgeT") %>% 
    mutate(
      rate = r0 * exp(Age * ba1),
      p_health = 1 + (exp(-rate * (ti + 1)) - exp(-rate * ti)) / rate,
      p_hz = 1 - p_health
    )
      
}


sim_pre_ql_f <- function(age0 = 50, pars, dt = 0.01) {
  n_sim <- max(pars$Key)
  #n_sim = 100
  df <- crossing(Key = 1:n_sim, Age = age0, ti = seq(0, 100-age0, dt)) %>% 
    mutate(
      AgeT = floor(Age + ti)
    ) %>% 
    left_join(pars, by = "Key")  %>% 
    # left_join(pars_demo %>% rename(AgeT = Age), by = "AgeT") %>% 
    mutate(
      rate = r0 * exp(Age * ba1),
      p_hz = (exp(-rate * ti) - exp(-rate * (ti + dt))) / rate,
      p_health = dt - p_hz,
      Q_0 = 1,
      Q_1 = C1_b0,
      Q_2 = C2_b0,
      #PZ = PZ_ba2 * Age ^2 + PZ_ba1 * Age + PZ_b0 + PZ_bd15 * (ti <= 15 / 365.25) + PZ_bd30 * (ti <= 30 / 365.25),
      PZ = PZ_b0 + PZ_bd15 * (ti <= 15 / 365.25) + PZ_bd30 * (ti <= 30 / 365.25),
      PZ = 1 / (1 + exp(-PZ)),
      PC1 = PC1_b0 + PC1_bd15 * (ti <= 15 / 365.25) + PC1_bd30 * (ti <= 30 / 365.25),
      PC1 = 1 / (1 + exp(-PC1)),
      Prop_0 = PZ,
      Prop_1 = PC1 * (1 - PZ),
      Prop_2 = 1 - Prop_0 - Prop_1,
      ql = Q_0 * Prop_0 + Q_1 * Prop_1 + Q_2 * Prop_2,
      ql = (1 - ql) * p_hz,
      # ql_pn = pmin(Q_0, norm) * Prop_0 + pmin(Q_1, norm) * Prop_1 + pmin(Q_2, norm) * Prop_2,
      # ql_pn = (norm - ql_pn) * p_hz
    ) %>% 
    group_by(Key, Age, AgeT) %>% 
    #select(Key, Age, AgeT, ti, ql_ph, ql_pn, rate, p_health, p_hz)
    summarise(
      across(c(p_hz, p_health, ql), sum)
      # across(c(mr, norm), mean)
    )
  
  return(df)
  
}


sim_pre_ql_b <- function(age0 = 50, pars, dt = 0.05) {
  n_sim <- max(pars$Key)
  #n_sim = 100
  df <- crossing(Key = 1:n_sim, Age = age0, ti = seq(0, 100-age0, dt)) %>% 
    mutate(
      AgeT = floor(Age + ti)
    ) %>% 
    left_join(pars, by = "Key")  %>% 
    mutate(
      rate = r0 * exp(Age * ba1),
      p_hz = (exp(-rate * ti) - exp(-rate * (ti + dt))) / rate,
      p_health = dt - p_hz,
      Q_0 = 1,
      Q_1 = C1_b0,
      Q_2 = C2_b0,
      # PZ = PZ_ba2 * Age ^2 + PZ_ba1 * Age + PZ_b0 + PZ_bd15 * (ti <= 15 / 365.25) + PZ_bd30 * (ti <= 30 / 365.25),
      PZ = PZ_b0 + PZ_bd15 * (ti <= 15 / 365.25) + PZ_bd30 * (ti <= 30 / 365.25),
      PZ = 1 / (1 + exp(-PZ)),
      cAge = (Age - 75) / 50,
      # PC1 = PC1_b0 + PC1_bd15 * (ti <= 15 / 365.25) + PC1_bd30 * (ti <= 30 / 365.25) + PC1_ba1 * cAge + PC1_ba2 * cAge ^2,
      PC1 = PC1_b0 + PC1_bd15 * (ti <= 15 / 365.25) + PC1_bd30 * (ti <= 30 / 365.25),
      PC1 = 1 / (1 + exp(-PC1)),
      Prop_0 = PZ,
      Prop_1 = PC1 * (1 - PZ),
      Prop_2 = 1 - Prop_0 - Prop_1,
      ql = Q_0 * Prop_0 + Q_1 * Prop_1 + Q_2 * Prop_2,
      ql = (1 - ql) * p_hz
      # ql_pn = pmin(Q_0, norm) * Prop_0 + pmin(Q_1, norm) * Prop_1 + pmin(Q_2, norm) * Prop_2,
      # ql_pn = (norm - ql_pn) * p_hz
    ) %>% 
    group_by(Key, Age, AgeT) %>% 
    #select(Key, Age, AgeT, ti, ql_ph, ql_pn, rate, p_health, p_hz)
    summarise(
      across(c(p_hz, p_health, ql), sum)
    )
  
  return(df)
  
}


sim_ql <- function(df) {
  df %>% 
    group_by(Key) %>% 
    mutate(
      surv = cumprod(1 - 0),
      surv = (surv + c(1, surv[-n()])) / 2,
      d15 = (1 + 0.015) ^ -(AgeT - Age),
      d35 = (1 + 0.035) ^ -(AgeT - Age)
    ) %>% 
    group_by(Key, Age) %>% 
    summarise(
      #LE = sum(surv),
      # QALE1 = sum((norm - ql_pn) * surv),
      # QALE0 = sum(norm * surv),
      # QL00 = sum(ql_pn * surv),
      # QL15 = sum(ql_pn * surv * d15),
      # QL35 = sum(ql_pn * surv * d35),
      QL00 = sum(ql * surv),
      QL15 = sum(ql * surv * d15),
      QL35 = sum(ql * surv * d35)
    ) %>% 
    ungroup()
}


simulate_shortfall <- function(pars, vset = "uk", age0 = 50, age1 = 99, year = 2024, mod = c("k", "f", "b")) {
  # pars <- tar_read(pars_shortfall, 1)
  # data_norm <- tar_read(data_norm)
  
  # pars_demo <- data_norm %>% 
  #   filter(Year == year) %>% 
  #   ungroup() %>% 
  #   select(Age, mr = mortality, norm = norm_leoss) 
  
  mod <- match.arg(mod)
  sim <- bind_rows(lapply(age0:age1, \(a0) {
    if (mod == "k") {
      df <- sim_pre_ql_k(age0 = a0, pars = pars)
    } else if (mod == "f") {
      df <- sim_pre_ql_f(age0 = a0, pars = pars, dt = 0.05)
    } else {
      df <- sim_pre_ql_b(age0 = a0, pars = pars, dt = 0.05)
    }
    
    df %>% sim_ql()
  }))
  
  # write_csv(sim, here::here("posteriors", paste0("QALY_loss_sims_", vset, ".csv")))
  
  return(sim)
}


summarise_shortfall <- function(sim, vset) {
  stats <- sim %>% 
    select(-Key) %>% 
    group_by(Age) %>% 
    summarise_all(amlu) %>% 
    pivot_longer(-Age, names_to = c("Index", "name"), names_sep = "_") %>% 
    pivot_wider()
  

  # write_csv(stats, here::here("docs", "tabs", paste0("QALY_loss_stats_", vset, ".csv")))
  
  return(stats)
}

