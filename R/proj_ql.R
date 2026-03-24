
sim_ql_s1 <- function(pars, age0, age_until = 5, dt = 0.1) {
  n_sim <- max(pars$Key)
  
  crossing(Key = 1:n_sim, Age = age0, ti = seq(0, age_until, dt)) %>% 
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
      PZ = PZ_b0 + PZ_bd15 * (ti <= 15 / 365.25) + PZ_bd30 * (ti <= 30 / 365.25),
      PZ = 1 / (1 + exp(-PZ)),
      cAge = (Age - 75) / 50,
      PC1 = PC1_b0 + PC1_bd15 * (ti <= 15 / 365.25) + PC1_bd30 * (ti <= 30 / 365.25),
      PC1 = 1 / (1 + exp(-PC1)),
      Prop_0 = PZ,
      Prop_1 = PC1 * (1 - PZ),
      Prop_2 = 1 - Prop_0 - Prop_1,
      d15 = exp(- 0.015 * ti),
      d35 = exp(- 0.035 * ti)
    )
}


simulate_ql <- function(pars, age0 = 50, age1 = 99, age_until = 5, dt = 0.1) {
  lapply(age0:age1, \(a0) {
    sim_ql_s1(pars, a0, age_until, dt) %>%
      mutate(
        ql = Q_0 * Prop_0 + Q_1 * Prop_1 + Q_2 * Prop_2,
        ql = (1 - ql) * p_hz,
        ql15 = ql * d15,
        ql35 = ql * d35,
      ) %>% 
      group_by(Key, Age, AgeT) %>% 
      summarise(across(c(p_hz, p_health, ql, ql15, ql35), sum)) %>% 
      group_by(Key) %>% 
      summarise(
        QL00 = sum(ql),
        QL15 = sum(ql15),
        QL35 = sum(ql35)
      ) %>% 
      mutate(Age = a0)
  }) %>% 
    bind_rows() %>% 
    ungroup()
}


get_norm_wt <- function(data_raw_pn, data_norm) {
  wts <- data_raw_pn %>% 
    group_by(study = SID) %>% 
    summarise(wt = n())
  
  data_norm %>% 
    left_join(wts) %>% 
    group_by(Age) %>% 
    summarise(
      across(c(norm, norm_leoss), \(x) weighted.mean(x, w = wt))
    ) %>% 
    mutate(
      study = "pooled"
    )
}


simulate_ql_pn <- function(pars, data_norm, age0, age1 = 99, age_until = 5, dt = 0.1) {
  lapply(age0:age1, \(a0) {
    sim_ql_s1(pars, a0, age_until, dt) %>% 
      left_join(data_norm, by = "Age", relationship = "many-to-many") %>% 
      mutate(
        Q_1 = 1 - pmax(norm_leoss - Q_1, 0),
        Q_2 = 1 - pmax(norm_leoss - Q_2, 0),
        ql = Q_0 * Prop_0 + Q_1 * Prop_1 + Q_2 * Prop_2,
        ql = (1 - ql) * p_hz,
        ql15 = ql * d15,
        ql35 = ql * d35,
      ) %>% 
      group_by(study, Key, Age, AgeT) %>% 
      summarise(across(c(p_hz, p_health, ql, ql15, ql35), sum)) %>% 
      group_by(study, Key)  %>% 
      summarise(
        QL00 = sum(ql),
        QL15 = sum(ql15),
        QL35 = sum(ql35)
      ) %>%  
      mutate(Age = a0)
  }) %>% 
    bind_rows() %>% 
    ungroup()
}


simulate_ql_0 <- function(pars, pars_baseline, age0, age1 = 99, age_until = 5, dt = 0.1) {
  lapply(age0:age1, \(a0) {
    sim_ql_s1(pars, a0, age_until, dt) %>%
      left_join(pars_baseline$p_base0 %>% rename(Age = age), by = "Age") %>% 
      left_join(pars_baseline$q, by = "Key") %>% 
      mutate(
        cap = ifelse(runif(n()) < p_ph, 1, EQ5D0),
        Q_1 = 1 - pmax(cap - Q_1, 0),
        Q_2 = 1 - pmax(cap - Q_2, 0),
        ql = Q_0 * Prop_0 + Q_1 * Prop_1 + Q_2 * Prop_2,
        ql = (1 - ql) * p_hz,
        ql15 = ql * d15,
        ql35 = ql * d35,
      ) %>% 
      group_by(Key, Age, AgeT) %>% 
      summarise(across(c(p_hz, p_health, ql, ql15, ql35), sum)) %>% 
      group_by(Key)  %>% 
      summarise(
        QL00 = sum(ql),
        QL15 = sum(ql15),
        QL35 = sum(ql35)
      ) %>% 
      mutate(Age = a0)
  }) %>% 
    bind_rows() %>% 
    ungroup()
}


summarise_ql <- function(pars_qloss, sim_qloss) {
  ages <- sim_qloss %>% pull(Age) %>% unique()
  
  pars_qloss %>% 
    crossing(Age = ages) %>% 
    mutate(
      rate = r0 * exp(Age * ba1),
      dur = 1 / rate
    ) %>% 
    group_by(Age) %>% 
    summarise(
      across(c(rate, dur), amlu)
    ) %>% left_join(
      sim_qloss %>% 
        group_by(Age) %>% 
        summarise(
          across(starts_with("QL"), amlu)
        )
    )
}



summarise_ql_pn <- function(sim_qloss_pn) {
  sim_qloss_pn %>% 
    group_by(study, Age) %>% 
    summarise(
      across(starts_with("QL"), amlu)
    )
}



bind_ql_tabs <- function(tab_qloss_ph, tab_qloss_pn, tab_qloss_0, ages = seq(50, 90, 10)) {
  bind_rows(
    tab_qloss_ph %>% mutate(Cap = "Perfect health"),
    tab_qloss_pn %>% mutate(Cap = "Population norm"),
    tab_qloss_0 %>% mutate(Cap = "Pre-HZ baseline")
  ) %>% 
    filter(Age %in% seq(50, 90, 10)) %>% 
    mutate(
      across(starts_with("dur"), \(x) round(x * 365.25)),
      Dur = sprintf("%s (%s - %s)", dur_M, dur_L, dur_U),
      QL00 = sprintf("%.3f (%.3f - %.3f)", QL00_M, QL00_L, QL00_U),
      QL15 = sprintf("%.3f (%.3f - %.3f)", QL15_M, QL15_L, QL15_U),
      QL35 = sprintf("%.3f (%.3f - %.3f)", QL35_M, QL35_L, QL35_U)
    ) %>% 
    select(Cap, Age, Dur, QL00, QL15, QL35)
  
}


simplify_ql_pn_sub_tabs <- function(tab_qloss_pn_sub, ages = seq(50, 90, 10)) {
  tab_qloss_pn_sub %>% 
    filter(Age %in% seq(50, 90, 10)) %>% 
    mutate(
      QL00 = sprintf("%.3f (%.3f - %.3f)", QL00_M, QL00_L, QL00_U),
      QL15 = sprintf("%.3f (%.3f - %.3f)", QL15_M, QL15_L, QL15_U),
      QL35 = sprintf("%.3f (%.3f - %.3f)", QL35_M, QL35_L, QL35_U)
    ) %>% 
    select(study, Age, QL00, QL15, QL35)
  
}


