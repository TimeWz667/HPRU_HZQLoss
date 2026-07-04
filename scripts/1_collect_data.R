library(tidyverse)


list_studies <- read_csv(here::here("data", "list_studies.csv"))

root <- here::here("data", "inputs_gen")
dir.create(root, showWarnings = F)


## TTE ----
src_tte_ph_gsk <- read_csv(here::here("data", "inputs_gsk_post", "Summary_tte.csv"))
src_tte_pn_gsk <- read_csv(here::here("data", "inputs_gsk_post", "Summary_tte_pn.csv"))

src_tte_ph_r1 <- read.csv(here::here("data", "inputs_r1_post", "stats_ph_tte.csv"))
src_tte_pn_r1 <- read.csv(here::here("data", "inputs_r1_post", "stats_pn_tte.csv"))
src_tte_base_r1 <- read.csv(here::here("data", "inputs_r1_post", "stats_base_tte.csv"))


src_slopes_r1 <- read_csv(here::here("data", "inputs_r1_post", "stats_tte_slopes.csv")) %>% 
  select(SID, Var, m) %>% 
  pivot_wider(names_from = Var, values_from = m)

src_slopes_ph_gsk <- read_csv(here::here("data", "inputs_gsk_post", "Summary_slopes.csv")) %>% 
  tidyr::extract(file, "SID", "Post_tte_(\\w+).csv")

src_slopes_pn_gsk <- read_csv(here::here("data", "inputs_gsk_post", "Summary_slopes_pn.csv")) %>% 
  tidyr::extract(file, "SID", "Post_tte_(\\w+).csv")



pars_tte_slopes <- bind_rows(
  src_slopes_r1 %>% 
    mutate(
      Set = "SetA",
      Base = case_when(
        SID == "Combined_pn" ~ "PN",
        SID == "Combined_ph" ~ "PH",
        T ~ "Baseline"
      )
    ) %>% 
    filter(startsWith(SID, "Combined")),
  src_slopes_ph_gsk %>% mutate(Set = "SetB", Base = "PH") %>% filter(SID == "gsk"),
  src_slopes_pn_gsk %>% mutate(Set = "SetB", Base = "PN") %>% filter(SID == "pn")
) %>% 
  pivot_longer(c(b0, b_lr0, rsd), names_to = "Var", values_to = "mean") %>% 
  select(- SID) %>% mutate(
    sd = 0
  )



pars_tte <- bind_rows(
  src_tte_ph_gsk %>% 
    mutate(Set = "SetB", Base = "PH") %>% filter(SID == "GSK"),
  src_tte_pn_gsk%>% 
    mutate(Set = "SetB", Base = "PN") %>% filter(SID == "GSK"),
  src_tte_ph_r1 %>% 
    mutate(Set = "SetA", Base = "PH"),
  src_tte_pn_r1 %>% 
    mutate(Set = "SetA", Base = "PN"),
  src_tte_base_r1 %>% 
    mutate(Set = "SetA", Base = "Baseline")
)  %>% 
  filter(Var %in% c("r0", "ba1")) %>% 
  select(Set, Base, Var, mean, sd) %>% 
  bind_rows(pars_tte_slopes) %>% 
  arrange(Set, Base, Var)



write_csv(pars_tte, here::here(root, "pars_tte.csv"))



## Quality of life ----

pars_qol_r1 <- list_studies %>% 
  filter(StudyGroup == "r1") %>% 
  pull(SID) %>% 
  lapply(\(sid) {
    bind_rows(
      read_csv(here::here("data", "inputs_r1_post", paste0("stats_qol_b_orig_sub_", sid, ".csv"))) %>% 
        filter(Var != "lp__") %>% 
        mutate(Var = paste0(Model, "_", Var), VSet = "Orig"),
      read_csv(here::here("data", "inputs_r1_post", paste0("stats_qol_b_bg_orig_sub_", sid, ".csv"))) %>% 
        filter(Var != "lp__") %>% 
        mutate(Var = gsub("\\[(\\d+)\\]", "", paste0(Model, "_", Var)), VSet = "Orig"),
      read_csv(here::here("data", "inputs_r1_post", paste0("stats_qol_b_uk_sub_", sid, ".csv"))) %>% 
        filter(Var != "lp__") %>% 
        mutate(Var = paste0(Model, "_", Var), VSet = "UK"),
      read_csv(here::here("data", "inputs_r1_post", paste0("stats_qol_b_bg_uk_sub_", sid, ".csv"))) %>% 
        filter(Var != "lp__") %>% 
        mutate(Var = gsub("\\[(\\d+)\\]", "", paste0(Model, "_", Var)), VSet = "UK")
    ) %>% 
      select(Model, Var, VSet, mean, sd) %>% 
      mutate(SID = sid)
  }) %>% 
  arrange(VSet, SID, Model, Var)


pars_qol_r1


pars_qol_gsk <- read_csv(here::here("data", "inputs_gsk_post", "Summary_qols_b.csv")) %>% 
  distinct() %>% 
  mutate(VSet = ifelse(value_set == "orig", "Orig", "UK")) %>% 
  filter(Var != "lp__") %>% 
  mutate(Var = gsub("\\[1\\]", "", Var)) %>% 
  filter(SID %in% c("DEU", "ESP", "JPN", "ITA")) %>% 
  mutate(Var = paste0(Model, "_", Var)) %>% 
  select(VSet, SID, Model, Var, mean, sd) %>% 
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
  ) %>% 
  arrange(VSet, SID, Model, Var)


pars_qol <- bind_rows(pars_qol_r1, pars_qol_gsk)


write_csv(pars_qol, here::here(root, "pars_qol.csv"))


# pars_qol_gskk <- read_csv(here::here("data", "inputs_gsk_post", "tab_qols.csv"))
# 
# 
# pars_qol <- bind_rows(
#   pars_qol_r1,
#   pars_qol_gskk %>% 
#     summarise(
#       PZ = Prop[Cluster == 0],
#       PZ_bg = 0,
#       PZ_bd15 = 0,
#       PZ_bd30 = 0,
#       PZ_b0 = log(PZ / (1 - PZ)), 
#       PC1 = Prop[Cluster == 1] / (Prop[Cluster == 1] + Prop[Cluster == 2]),
#       PC1_b0 = log(PC1 / (1 - PC1)), 
#       PC1_bg = 0,
#       PC1_bd15 = 0,
#       PC1_bd30 = 0,
#       .by = "Country"
#     ) %>% 
#     pivot_longer(- Country, values_to = "mean", names_to = "Var") %>% 
#     mutate(sd = 0) %>% 
#     tidyr::extract(Var, "Model", "(\\w+)_", remove = F) %>% 
#     filter(!is.na(Model)) %>% 
#     rename(SID = Country),
#   pars_qol_gskk %>% 
#     filter(Cluster != 0) %>%
#     mutate(Model = paste0("C", Cluster), bg = 0) %>% 
#     rename(b0 = mu, sigma = std) %>% 
#     select(Model, SID = Country, b0, sigma, bg) %>% 
#     pivot_longer(- c(SID, Model), values_to = "mean", names_to = "Var") %>% 
#     mutate(sd = 0, Var = paste0(Model, "_", Var)),
#   
# )
# 
# write_csv(pars_qol, here::here(root, "pars_qol_km.csv"))




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





