library(tidyverse)


list_studies <- read_csv(here::here("data", "list_studies.csv"))



root <- here::here("data", "inputs_for_gen")

dir.create(root, showWarnings = F)



## TTE sub ----
ds <- dir(here::here("data", "inputs_r1_post"))
ds <- ds[startsWith(ds, "stats_tte_sub")]

ds

pars_tte_r1 <- bind_rows(lapply(ds, \(f) {
  if (f != "RamCRI") {
    tab <- read_csv(here::here("data", "inputs_r1_post", f)) %>% 
      filter(Var %in% c("r0", "ba1")) %>% 
      mutate(SID = f)
  } else {
    tab <- read_csv(here::here("data", "inputs_r1_post", "stats_ph_tte.csv")) %>% 
      filter(Var %in% c("r0", "ba1")) %>% 
      mutate(SID = f)
  }
  
tab
})) %>% 
  tidyr::extract(SID, "SID", "stats_tte_sub_(\\w+).csv")


src_tte_gsk <- read_csv(here::here("data", "inputs_gsk_post", "Summary_tte.csv"))

pars_tte_gsk <- src_tte_gsk %>% 
  filter(Var %in% c("r0", "ba1")) %>% 
  filter(SID %in% c("DEU", "ESP", "JPN", "ITA"))


pars_tte <- bind_rows(pars_tte_r1, pars_tte_gsk) %>% 
  select(SID, Var, mean, sd)



write_csv(pars_tte, here::here(root, "pars_tte_ind.csv"))



pars_tte_r1_ph <- read.csv(here::here("data", "inputs_r1_post", "stats_ph_tte.csv"))


pars_tte <- bind_rows(
    pars_tte_r1 %>% 
      select(SID, Var) %>% 
      left_join(pars_tte_r1_ph), 
    pars_tte_gsk %>% 
      select(SID, Var) %>% 
      left_join(src_tte_gsk %>% filter(SID == "GSK") %>% select( - SID))
  ) %>% 
  select(SID, Var, mean, sd)

write_csv(pars_tte, here::here(root, "pars_tte_2gp.csv"))




## TTE slopes
src_slopes_r1 <- read_csv(here::here("data", "inputs_r1_post", "stats_tte_slopes.csv")) %>% 
  select(SID, Var, m) %>% 
  pivot_wider(names_from = Var, values_from = m)

src_slopes_gsk <- read_csv(here::here("data", "inputs_gsk_post", "Summary_slopes.csv")) %>% 
  tidyr::extract(file, "SID", "Post_tte_(\\w+).csv")


pars_slopes <- bind_rows(
  src_slopes_r1 %>% 
    filter(SID != "RamCRI") %>% 
    mutate(SID = ifelse(SID == "Combined", "RamCRI", SID)),
  src_slopes_gsk %>% 
    filter(SID %in% c("DEU", "ESP", "JPN", "ITA"))
)


write_csv(pars_slopes, here::here(root, "pars_tte_slopes_ind.csv"))



pars_slopes <- bind_rows(
  src_slopes_r1 %>% select(SID) %>% 
    filter(SID != "Combined") %>% 
    bind_cols(src_slopes_r1 %>% filter(SID == "Combined") %>% select(-SID)),
  src_slopes_gsk %>% select(SID) %>% 
    filter(SID %in% c("DEU", "ESP", "JPN", "ITA")) %>% 
    bind_cols(src_slopes_gsk %>% filter(SID == "gsk") %>% select(-SID))
)

write_csv(pars_slopes, here::here(root, "pars_tte_slopes_2gp.csv"))





### QOL ----
pars_qol_r1 <- list_studies %>% 
  filter(StudyGroup == "r1") %>% 
  pull(SID) %>% 
  lapply(\(sid) {
    bind_rows(
      read_csv(here::here("data", "inputs_r1_post", paste0("stats_qol_b_orig_sub_", sid, ".csv"))) %>% 
        filter(Var != "lp__") %>% 
        mutate(Var = paste0(Model, "_", Var)),
      read_csv(here::here("data", "inputs_r1_post", paste0("stats_qol_b_bg_orig_sub_", sid, ".csv"))) %>% 
        filter(Var != "lp__") %>% 
        mutate(Var = gsub("\\[(\\d+)\\]", "", paste0(Model, "_", Var)))
    ) %>% 
      select(Model, Var, mean, sd) %>% 
      mutate(SID = sid)
  }) %>% 
  bind_rows() %>% 
  arrange(SID, Model, Var)


pars_qol_r1


pars_qol_gsk0 <- read_csv(here::here("data", "inputs_gsk_post", "Summary_qols.csv")) %>% 
  distinct() %>% 
  filter(Var != "lp__") %>% 
  mutate(Var = gsub("\\[1\\]", "", Var)) %>% 
  filter(SID %in% c("DEU", "ESP", "JPN", "ITA")) %>% 
  mutate(Var = paste0(Model, "_", Var)) %>% 
  select(SID, Model, Var, mean, sd) %>% 
  mutate(
    mean = case_when(
      Var == "PZ_bg" ~ 0,
      Var == "PC1_bg" ~ 0,
      T ~ mean
    ),
    sd = case_when(
      Var == "PZ_bg" ~ 0,
      Var == "PC1_bg" ~ 0,
      T ~ sd
    )
  )


pars_qol <- bind_rows(pars_qol_r1, pars_qol_gsk0)


write_csv(pars_qol, here::here(root, "pars_qol.csv"))


pars_qol_gskk <- read_csv(here::here("data", "inputs_gsk_post", "tab_qols.csv"))


pars_qol <- bind_rows(
  pars_qol_r1,
  pars_qol_gskk %>% 
    summarise(
      PZ = Prop[Cluster == 0],
      PZ_bg = 0,
      PZ_bd15 = 0,
      PZ_bd30 = 0,
      PZ_b0 = log(PZ / (1 - PZ)), 
      PC1 = Prop[Cluster == 1] / (Prop[Cluster == 1] + Prop[Cluster == 2]),
      PC1_b0 = log(PC1 / (1 - PC1)), 
      PC1_bg = 0,
      PC1_bd15 = 0,
      PC1_bd30 = 0,
      .by = "Country"
    ) %>% 
    pivot_longer(- Country, values_to = "mean", names_to = "Var") %>% 
    mutate(sd = 0) %>% 
    tidyr::extract(Var, "Model", "(\\w+)_", remove = F) %>% 
    filter(!is.na(Model)) %>% 
    rename(SID = Country),
  pars_qol_gskk %>% 
    filter(Cluster != 0) %>%
    mutate(Model = paste0("C", Cluster), bg = 0) %>% 
    rename(b0 = mu, sigma = std) %>% 
    select(Model, SID = Country, b0, sigma, bg) %>% 
    pivot_longer(- c(SID, Model), values_to = "mean", names_to = "Var") %>% 
    mutate(sd = 0, Var = paste0(Model, "_", Var)),
  
)

write_csv(pars_qol, here::here(root, "pars_qol_km.csv"))




## Age distribution

prof_a_r1 <- read_csv(here::here("data", "inputs_studies", "tab_agp_by_studies.csv"))
# prof_a_r1 <- prof_a_r1 %>% rename(study = SID) %>% left_join(list_studies) %>% select(-c(study, SIK)) %>% 
#   mutate(SID = factor(SID, list_studies$SID))

prof_a_gsk <- read_csv(here::here("data", "inputs_gsk_post", "tab_cnt.csv"))


prof_a <- bind_rows(
  prof_a_r1 %>% 
    filter(!is.na(Agp)) %>% 
    select(study = SID, Agp, n) %>% 
    tidyr::extract(Agp, c("A0", "A1"), "\\[(\\d+),(\\d+)\\)", convert = T, remove = F) %>% 
    left_join(list_studies %>% select(study, SID)) %>% 
    select(SID, Agp, A0, A1, n),
  prof_a_gsk %>% 
    mutate(
      SID = case_when(
        value_set == "Germany" ~ "DEU", 
        value_set == "Japan"~ "JPN",
        value_set == "Italy" ~ "ITA",
        value_set == "Spain" ~ "ESP",
        T ~ NA
      )
    ) %>% 
    select(SID, Agp, A0 = age_min, A1 = age_max, n = n_sub)
)



write_csv(prof_a, here::here(root, "prof_a.csv"))
