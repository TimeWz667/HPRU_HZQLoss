
calc_timing <- \(df) {
  df %>% 
    group_by(SID, PID) %>% 
    select(SID, PID, Age, Agp, ti, EQ5D) %>% 
    mutate(
      loss = 1 - EQ5D,
      End = cumsum(loss) == sum(loss) & loss == loss[n()],
      EndZeros = End * (loss == 0)
    ) %>% 
    summarise(
      QL_act = mean(loss[!End]),
      QL_end = mean(loss[End]),
      T_last = max(ti[!End]),
      T_last = ifelse(is.infinite(T_last), 0, T_last),
      T_end = min(ti[End]),
      T_end = ifelse(T_last == 0, 0, T_end),
      T_obs = max(ti),
      T_evt = (T_last + T_end) / 2,
      Recovered = sum(EndZeros) > 0
    ) %>% 
    mutate(
      T_evt = ifelse(Recovered, T_evt, T_obs)
    ) %>% 
    ungroup()
}


calc_hazard <- \(df, t0 = 0, t1 = 2, dt= 1/12) {
  bind_rows(lapply(seq(t0, t1, dt), \(ti) {
    from <- ti
    to <- ti + dt
    
    df %>%
      filter(pmin(T_evt, T_obs) > from) %>% 
      summarise(
        ti = ti,
        AtRisk = n(),
        Evt = sum(T_evt <= to),
        Haz = Evt / AtRisk
      )
  }))
}


empirical_surv <- function(df, freq = 1 / 12) {
  haz_full <- df %>% calc_hazard(dt = 1/12)
  
  haz_season_full <- haz_full %>% 
    mutate(
      ti = floor(ti * 4) / 4
    ) %>% 
    group_by(ti) %>% 
    summarise(Haz = mean(Haz))
  
  haz <- with(haz_season_full %>% filter(!is.na(Haz)), {
    as.data.frame(approx(x = ti, y = Haz, xout = seq(0, 5, 1 / 12), rule = 2)) %>% 
      rename(ti = x, Haz = y)
  }) %>% 
    mutate(
      Haz = ifelse(ti == max(ti), 1, Haz),
      Surv = c(1, cumprod(1 - Haz)[-n()])
    )
  
  hf <- approxfun(x = haz$ti, y = haz$Surv, rule = 2)
  inv_hf <- approxfun(x = haz$Surv, y = haz$ti, rule = 2)
  
  sim_tte <- \(t0) inv_hf(hf(t0) * runif(1))
  
  list(
    fn = sim_tte,
    haz = haz
  )
  
  return(sim_tte)
}


format_tte <- function(df) {
  require(tidyverse)
  
  dat_tte <- df %>% calc_timing()
  
  sim_tte <- empirical_surv(dat_tte)
  
  ## Test simulation
  # tte <- sapply(1:1000, \(i) sim_tte(0))
  # plot(survival::Surv(tte))
  # points(haz$ti, haz$Surv, col = 2, pch = 20, cex = 2)
  # 
  # 
  # tte <- sapply(1:1000, \(i) sim_tte(0.5))
  # plot(survival::Surv(tte))
  # haz_t <- haz %>% filter(ti >= 0.5)
  # points(haz_t$ti, haz_t$Surv / hf(0.5), col = 2, pch = 20, cex = 2)
  # 
  # 
  
  ## Insert event points to those censored
  dat_filled <- dat_tte %>% 
    mutate(
      T_evt2 = ifelse(QL_end == 0, T_evt, sapply(T_obs, sim_tte))
    ) %>% 
    left_join(
      df %>% 
        select(PID, SID, Age, Agp) %>% 
        distinct()
    )
  
  
  #write_csv(dat_filled, here::here("data", "tte_uk.csv"))
  
  return(dat_filled)
}


format_qol <- function(df) {
  dat <- df %>% 
    arrange(PID, ti) %>% 
    group_by(PID) %>% 
    filter(!all(EQ5D == 1)) %>% 
    mutate(
      loss = 1 - EQ5D,
      End = cumsum(loss) == sum(loss) & loss == loss[n()],
      EndZeros = End * (loss == 0)
    ) %>% 
    filter(EndZeros == 0) %>% 
    ungroup() %>% 
    select(PID, SID, ti, Age, Agp, EQ5D, loss) %>% 
    mutate(
      Key = 1:n()
    )
  
  return(dat)
}


