
get_data_norm <- function(file_norm) {
  require(tidyverse)
  require(readxl)
  
  
  demo_proj <- read_xls(here::here("data", "external", "uk_pop_proj_2020.xls"), sheet = "Population") %>% 
    pivot_longer(-c(Sex, Age), names_to = "Year", values_to = "N") %>% 
    mutate(
      Age = ifelse(Age %in% as.character(0:100), Age, "100"),
      Age = as.numeric(Age),
      Year = as.numeric(Year),
      Sex = ifelse(Sex == 1, "m", "f")
    ) %>% 
    group_by(Year, Sex, Age) %>% 
    summarise(N = sum(N)) %>% 
    ungroup()
  
  
  demo_pre <- bind_rows(lapply(2012:2020, function(yr) {
    raw_pre <- read_xlsx(here::here("data", "external", "uk_pop_data.xlsx"), sheet = as.character(yr))
    raw_pre %>% 
      filter(geogcode %in% c("E92000001", "K02000001")) %>% 
      pivot_longer(starts_with(c("m_", "f_"))) %>% 
      separate(name, c("Sex", "Year", "Age")) %>% 
      filter(Age != "al") %>% 
      mutate(
        Age = as.numeric(Age),
        Year = 2000 + as.numeric(Year)
      ) %>% 
      select(Location = variable, Year, Age, Sex, N = value)
  })) %>% 
    mutate(
      Location = ifelse(Location == "ENGLAND", "England", "UK")
    )
  
  
  demo_pre_old <- demo_pre %>% 
    filter(Age == 90) %>% 
    select(- Age) %>% 
    full_join(crossing(Year = 2012:2020, Age = 90:100), relationship = "many-to-many") %>% 
    group_by(Location, Year, Sex) %>% 
    mutate(
      k = 100 - Age,
      k = k / sum(k),
      N = N * k
    ) %>% 
    ungroup() %>% 
    select(Location, Year, Sex, Age, N)
  
  
  demo_pre <- bind_rows(
    demo_pre %>% filter(Age != 90),
    demo_pre_old
  )
  
  
  rat <- demo_pre %>% 
    group_by(Year, Location, Sex) %>%
    summarise(N = sum(N)) %>% 
    pivot_wider(names_from = Location, values_from = N) %>% 
    mutate(rat = England / UK) %>% 
    filter(Year > 2016) %>% 
    group_by(Sex) %>% 
    summarise(rat = mean(rat))
  
  
  mor_proj <- bind_rows(
    read_xlsx(here::here("data", "external", "uk_mor.xlsx"), sheet = "males period qx", skip = 4) %>% 
      pivot_longer(-age, names_to = "Year", values_to = "mortality") %>% 
      rename(Age = age) %>% 
      mutate(Sex = "m"),
    read_xlsx(here::here("data", "external", "uk_mor.xlsx"), sheet = "females period qx", skip = 4) %>% 
      pivot_longer(-age, names_to = "Year", values_to = "mortality") %>% 
      rename(Age = age) %>% 
      mutate(Sex = "f")
  ) %>% 
    mutate(
      mortality = mortality * 1e-5,
      Year = as.numeric(Year)
    ) 
  
  
  demo_ons <- mor_proj %>% 
    inner_join(
      bind_rows(
        demo_pre %>% filter(Location == "England") %>% select(-Location),
        demo_proj %>% filter(Year > 2020) %>% left_join(rat) %>% 
          mutate(N = round(N * rat)) %>% 
          select(-rat)
      )
    )
  
  
  ##
  ## Source: https://github.com/bitowaqr/shortfall
  ##
  
  pn <- read_csv(file_norm) %>% 
    mutate(Sex = ifelse(sex == "female", "f", "m")) %>% 
    select(Sex, cw, co, age5_str, age5_start) %>% 
    distinct() %>% 
    group_by(Sex) %>% 
    mutate(
      age5_end = c(age5_start[-1] - 1, 100),
      age5_mid = (age5_start + age5_end) / 2
    ) %>% 
    ungroup()
  
  
  ## Combine them and aggregate across sexes
  
  sup_demo_s <- demo_ons  %>% 
    mutate(
      age5_str = floor(Age * 0.2) /0.2,
      age5_str = paste(age5_str, age5_str + 4, sep = "-"),
      age5_str = case_when(
        Age < 18 ~ "16-17",
        Age < 20 ~ "18-19",
        Age >= 90 ~ "90+",
        T ~ age5_str
      )
    ) %>% 
    left_join(pn) %>% 
    select(Age, Year, mortality, Sex, N, age5_str, norm = co)
  
  
  sup_demo <- sup_demo_s %>% 
    group_by(Year, age5_str) %>% 
    mutate(
      norm = weighted.mean(norm, N)
    ) %>% 
    group_by(Year, Age) %>% 
    summarise(
      mortality = weighted.mean(mortality, wt = N),
      norm = weighted.mean(norm, wt = N),
      N = sum(N)
    ) %>% 
    group_by(Year) %>% 
    mutate(
      norm_leoss = predict(loess(norm ~ Age))
    )
  
  
  save(sup_demo_s, sup_demo, file = here::here("data", "processed", "sup_demo.rdata"))
  
  return(sup_demo)
}


get_data_qol <- function(file) {
  ## QoL with uk value sets
  raw <- read_csv(file)
  
  max_qol <- max(raw$EQ5D)
  min_qol <- min(raw$EQ5D)
  
  dat_qol <- raw %>% 
    mutate(ti = time_points / 365.25) %>% 
    select(SID, PID, Age, Agp, ti = ti, EQ5D) %>% 
    mutate(
      Q_rescaled = (EQ5D - min_qol) / (max_qol - min_qol)
    )
  
  return(dat_qol)
}


get_data_pn <- function(data_raw, data_pn) {
  # data_raw <- tar_read(data_raw, 1)
  # data_pn <- tar_read(data_norm, 1)
  
  dat_qol_pn <- data_raw %>% 
    filter(!is.na(Age)) %>% 
    left_join(data_pn %>% select(SID = study, Age, norm_leoss), by = c("SID", "Age")) %>% 
    rename(EQ5D1 = EQ5D) %>% 
    mutate(
      EQ5D = pmin(norm_leoss, EQ5D1) + (1 - norm_leoss),
      loss = 1 - EQ5D,
      offset = EQ5D - EQ5D1
    )
  
  return(dat_qol_pn)
  
}


get_data_qol_shift <- function(data_raw, data_baseline) {
  ## QoL with uk value sets
  # data_baseline <- tar_read(data_baseline, 1)
  dat_qol_shift <- data_raw %>% 
    inner_join(data_baseline %>% select(PID = Patient.ID, EQ5D0 = EQ5D)) %>%
    rename(EQ5D1 = EQ5D) %>% 
    mutate(
      EQ5D = 1 - pmax(EQ5D0 - EQ5D1, 0),
      EQ5D2 = pmin(EQ5D0, EQ5D1) + (1 - EQ5D0),
      loss = 1 - EQ5D,
      offset = EQ5D - EQ5D1
    )
  
  return(dat_qol_shift)
}



