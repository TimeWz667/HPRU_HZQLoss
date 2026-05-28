library(tidyverse)



source(here::here("R", "gen_analyses.R"))



root_out <- here::here("docs", "tabs", "gen")
dir.create(root_out, showWarnings = F)



sets <- c("2gp", "2gp_km", "ind", "ind_km")


res_a <- list()
res_t <- list()


for (set in sets) {
  root_src <- here::here("out", paste0("gen_", set))
  last_v <- read_file(here::here(root_src, "last_version.txt")) %>% glue::as_glue()
  
  cat("Folder: ", root_src, "\n")
  cat("Version: ", last_v, "\n")
  
  
  batches <- dir(here::here(root_src), pattern = paste0("post_qol_s(\\d+)_", last_v, ".csv")) %>% 
    str_extract("post_qol_(s\\d+)_", 1)
  
  
  summ_qls <- lapply(batches, \(batch) {
    pars_qol <- read_csv(here::here(root_src, paste0("post_qol_", batch, "_", last_v, ".csv")))
    pars_tte <- read_csv(here::here(root_src, paste0("post_tte_exp_", batch, "_", last_v, ".csv"))) %>% 
      mutate(Key = 1:n())
    
    summ <- post_pars <- pars_qol %>% 
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
        ages = seq(50, 90, 5),
        times = c(0, 14, 30, 60, 90)
      )
    
    summ <- lapply(summ, \(x) x %>% mutate(Batch = batch, Set = set))
    
    return(summ)
  })
  
  
  summ_qls <- list(
    tab_time = bind_rows(lapply(summ_qls, \(summ) summ$tab_time)),
    tab_age = bind_rows(lapply(summ_qls, \(summ) summ$tab_age))
  )
  
  write_csv(summ_qls$tab_time, here::here(root_out, paste0("stats_qloss_time_", set, ".csv")))
  write_csv(summ_qls$tab_age, here::here(root_out, paste0("stats_qloss_age_", set, ".csv")))
  
  res_t[[length(res_t) + 1]] <- summ_qls$tab_time
  res_a[[length(res_a) + 1]] <- summ_qls$tab_age
}



res_t %>% bind_rows() %>% 
  select(- Batch) %>% 
  summarise(
    across(everything(), mean),
    .by = c("Time", "Set")
  ) %>% 
  write_csv(here::here(root_out, "stats_qloss_time.csv"))


res_a %>% bind_rows() %>% 
  select(- Batch) %>% 
  summarise(
    across(everything(), mean),
    .by = c("Age", "Set")
  ) %>% 
  write_csv(here::here(root_out, "stats_qloss_age.csv"))



summ_a_all <- res_a %>% bind_rows() %>% 
  select(- Batch) %>% 
  summarise(
    across(everything(), mean),
    .by = c("Age", "Set")
  ) 





root_src <- here::here("data", "inputs_gsk_post")
pars_qol <- read_csv(here::here(root_src, "Summary_qols.csv")) %>% 
  mutate(Var = paste0(Model, "_", Var)) %>% 
  filter(SID == "GSK") %>% 
  select(Var, mean, sd) %>% 
  filter(Var != "lp__")

pars_tte <- read_csv(here::here(root_src, "Summary_tte.csv")) %>% 
  filter(SID == "GSK") %>% 
  select(Var, mean, sd) %>% 
  filter(Var != "lp__")

pars_tte_slopes <- read_csv(here::here(root_src, "Summary_slopes.csv")) %>% 
  filter(file == "Post_tte_gsk.csv") %>% 
  select(- file)



boot_pars <- revive_pars(n_mc = 200, pars_qol, pars_tte, pars_tte_slopes)
summ_gsk_orig <- boot_pars %>% summary_proj(
  ages = seq(50, 90, 5),
  times = c(0, 14, 30, 60, 90)
)



pars_qol_gskk <- read_csv(here::here(root_src, "tab_qols.csv"))

pars_qol <- bind_rows(
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
    mutate(sd = 0, Var = paste0(Model, "_", Var))
)  %>% 
  summarise(
    across(c(mean, sd), mean), .by = Var
  )


boot_pars <- revive_pars(n_mc = 200, pars_qol, pars_tte, pars_tte_slopes)
summ_gsk_uk <- boot_pars %>% summary_proj(
  ages = seq(50, 90, 5),
  times = c(0, 14, 30, 60, 90)
)


theme_set(theme_bw())

summ_a_gsk <- bind_rows(
  summ_gsk_uk$tab_age %>% mutate(Value_Set = "Value set (GSK): UK"),
  summ_gsk_orig$tab_age %>% mutate(Value_Set = "Value set (GSK): Original")
) 



g_comparing_qloss <- summ_a_all %>% 
  mutate(
    Value_Set = ifelse(endsWith(Set, "km"), "Value set (GSK): UK", "Value set (GSK): Original"),
    Disuti_Group = ifelse(startsWith(Set, "2gp"), "Disutility: Two groups", "Disutility: Individual")
  ) %>% 
  ggplot(aes(x = Age)) +
  geom_ribbon(aes(ymin = QL35_L, ymax = QL35_U), alpha = 0.2) +
  geom_line(aes(y = QL35_M)) +
  #geom_ribbon(data = summ_r1$tab_age, aes(ymin = QL35_L, ymax = QL35_U, fill = "Set R1"), alpha = 0.2) +
  geom_line(data = summ_r1$tab_age, aes(y = QL35_M, colour = "Set R1"), linewidth = 2) +
  geom_line(data = summ_r1$tab_age, aes(y = QL35_L, colour = "Set R1"), linetype = 2, linewidth = 1.5) +
  geom_line(data = summ_r1$tab_age, aes(y = QL35_U, colour = "Set R1"), linetype = 2, linewidth = 1.5) +
  #geom_ribbon(data = summ_gsk$tab_age, aes(ymin = QL35_L, ymax = QL35_U, fill = "Set GSK"), alpha = 0.2) +
  geom_line(data = summ_a_gsk, aes(y = QL35_M, colour = "Set GSK"), linewidth = 2) +
  geom_line(data = summ_a_gsk, aes(y = QL35_L, colour = "Set GSK"), linetype = 2, linewidth = 1.5) +
  geom_line(data = summ_a_gsk, aes(y = QL35_U, colour = "Set GSK"), linetype = 2, linewidth = 1.5) +
  facet_grid(Value_Set~Disuti_Group) +
  scale_y_continuous("QALY loss") +
  expand_limits(y = 0)


ggsave(g_comparing_qloss, filename = here::here("docs", "figs", "g_comparing_gen.pdf"), width = 12, height = 10)



g_comparing_disuti <- summ_a_all %>% 
  mutate(
    Value_Set = ifelse(endsWith(Set, "km"), "Value set (GSK): UK", "Value set (GSK): Original"),
    Disuti_Group = ifelse(startsWith(Set, "2gp"), "Disutility: Two groups", "Disutility: Individual")
  ) %>% 
  filter(Value_Set == "Value set (GSK): Original") %>% 
  ggplot(aes(x = Age)) +
  geom_ribbon(aes(ymin = disuti_L, ymax = disuti_U), alpha = 0.2) +
  geom_line(aes(y = disuti_M)) +
  #geom_ribbon(data = summ_r1$tab_age, aes(ymin = QL35_L, ymax = QL35_U, fill = "Set R1"), alpha = 0.2) +
  geom_line(data = summ_r1$tab_age, aes(y = disuti_M, colour = "Set R1"), linewidth = 2) +
  geom_line(data = summ_r1$tab_age, aes(y = disuti_L, colour = "Set R1"), linetype = 2, linewidth = 1.5) +
  geom_line(data = summ_r1$tab_age, aes(y = disuti_U, colour = "Set R1"), linetype = 2, linewidth = 1.5) +
  #geom_ribbon(data = summ_gsk$tab_age, aes(ymin = QL35_L, ymax = QL35_U, fill = "Set GSK"), alpha = 0.2) +
  geom_line(data = summ_a_gsk %>% filter(Value_Set == "Value set (GSK): Original"), aes(y = disuti_M, colour = "Set GSK"), linewidth = 2) +
  geom_line(data = summ_a_gsk %>% filter(Value_Set == "Value set (GSK): Original"), aes(y = disuti_L, colour = "Set GSK"), linetype = 2, linewidth = 1.5) +
  geom_line(data = summ_a_gsk %>% filter(Value_Set == "Value set (GSK): Original"), aes(y = disuti_U, colour = "Set GSK"), linetype = 2, linewidth = 1.5) +
  facet_grid(Value_Set~Disuti_Group) +
  scale_y_continuous("Disutility period, month", labels = scales::number_format(scale = 12)) +
  expand_limits(y = 0)


ggsave(g_comparing_disuti, filename = here::here("docs", "figs", "g_comparing_gen_disuti.pdf"), width = 12, height = 6)





