describe_basic <- function(dat_raw, dat_tte, pars_qol, vset) {
  require(tidyverse)
  
  stats <- list()
  
  stats$Counts <- bind_rows(
    dat_raw %>% 
      group_by(SID, PID) %>% 
      summarise(N = n(), Age = mean(Age)) %>% 
      group_by(SID) %>% 
      summarise(
        N_record = sum(N),
        N_patient = n(),
        N_mean = mean(N),
        N_max = max(N),
        N_min = min(N)
      ),
    dat_raw %>% 
      group_by(PID) %>% 
      summarise(N = n(), Age = mean(Age)) %>% 
      ungroup() %>% 
      summarise(
        N_record = sum(N),
        N_patient = n(),
        N_mean = mean(N),
        N_max = max(N),
        N_min = min(N)
      )
  )
  
  stats$Agp <- dat_raw %>% 
    select(PID, Age, Agp) %>% 
    distinct() %>% 
    mutate(
      Gp = case_when(
        Age < 50 ~ "< 50",
        Age < 60 ~ "50-59",
        Age < 70 ~ "60-69",
        Age < 80 ~ "70-79",
        T ~ "80+"
      )
    ) %>% 
    group_by(Gp) %>% 
    summarise(N = n())
  
  stats$Perfect <- dat_raw %>% 
    summarise(p = sum(EQ5D == 1))
  
  ## End points
  
  stats$Endpoint <- dat_tte %>% 
    summarise(
      N = n(),
      N_rec = sum(Recovered),
      N_nr = N - N_rec,
      P_nr = N_nr / N
    )

  stats$Labels <- pars_qol$Labeled %>% 
    group_by(PID) %>% 
    summarise(
      has_severe = any(cluster == 2),
      has_severe90 = any(cluster == 2 & ti > 90 / 365.25),
      has_tw = any(cluster == 0)
    ) %>% 
    ungroup() %>% 
    select(-PID) %>% 
    summarise_all(list(n = sum, p = mean))
  

  jsonlite::write_json(stats, here::here("docs", "tabs", "stats_basic_" + glue::as_glue(vset) + ".json"))
  
  return(stats)
}
