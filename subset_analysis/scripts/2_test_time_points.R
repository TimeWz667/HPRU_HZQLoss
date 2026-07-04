library(tidyverse)
library(rstan)

theme_set(theme_bw())

options(mc.cores = 8)
rstan_options(auto_write = TRUE)


for (file in dir("R")) source(here::here("R", file))

model_tte <- stan_model(here::here("models", "time2zero_surv_age.stan"))
model_exp <- stan_model(here::here("models", "time2zero_age.stan"))


## Background population

list_studies <- read_csv(here::here("data", "list_studies.csv")) %>% 
  select(study, SID, SIK)
list_studies


prof_a <- read_csv(here::here("data", "inputs_studies", "tab_agp_by_studies.csv"))
prof_a <- prof_a %>% rename(study = SID) %>% left_join(list_studies) %>% select(-c(study, SIK)) %>% 
  mutate(SID = factor(SID, list_studies$SID))

prof_t <- read_csv(here::here("data", "inputs_studies", "tab_t_by_studies.csv"))
prof_t <- prof_t %>% rename(study = SID) %>% left_join(list_studies) %>% select(-c(study, SIK)) %>% 
  mutate(SID = factor(SID, list_studies$SID))


slopes <- read_csv(here::here("data", "inputs_r1_post", "stats_tte_slopes.csv")) %>% 
  filter(SID == "Combined")

stats_tte <- read_csv(here::here("data", "inputs_r1_post", "stats_ph_tte.csv"))

## Functions

repli <- function(ext) {
  ext %>% 
    crossing(Age = 50:95)  %>% 
    mutate(
      rate = r0 * exp(Age * ba1),
      disuti = rexp(n(), rate = rate)
    ) %>% 
    group_by(Age) %>% 
    summarise(
      dur_M = median(disuti),
      dur_L = quantile(disuti, 0.025),
      dur_U = quantile(disuti, 0.975),
      durr_M = median(1 / rate),
      durr_L = quantile(1 / rate, 0.025),
      durr_U = quantile(1 / rate, 0.975)
    )
}

## Simulate age distribution and tte

sims_a <- prof_a %>% 
  filter(!is.na(Agp)) %>% 
  select(SID, Agp, n) %>% 
  uncount(n) %>% 
  tidyr::extract(Agp, c("A0", "A1"), "\\[(\\d+),(\\d+)\\)", convert = T, remove = F)  %>% 
  mutate(
    PID = paste0(SID, "_", 1:n()),
    A0 = pmax(A0, 18),
    Age = floor(runif(n(), A0, A1))
  ) %>% 
  select(-c("A0", "A1")) %>% 
  crossing(
    stats_tte %>% 
      select(Var, mean, sd) %>% 
      filter(Var != "lp__")
  ) %>% 
  mutate(value = rnorm(n(), mean, sd)) %>% 
  select(-c(mean, sd)) %>% 
  pivot_wider(names_from = Var, values_from = value) %>% 
  bind_cols(
    slopes %>%
      select(Var, m) %>% 
      pivot_wider(names_from = Var, values_from = m)
  ) %>% 
  mutate(
    r0 = pmax(r0, 0.0001),
    ba1 = b0 + log(r0) * b_lr0 + rnorm(n(), 0, rsd),
    rate = r0 * exp(Age * ba1)
  )



### Distributed with data 



summary_orig <- lapply(paste0("Var", 1:20), function(v) {
  sims_at <- sims_a %>% 
    mutate(disuti = rexp(n(), rate = rate)) %>% 
    left_join(prof_t %>% select(SID, visit, t_l, t_u, Prop), relationship = "many-to-many") %>% 
    mutate(
      time_points = floor(runif(n(), t_l, t_u + 1)),
      ti = time_points / 365.25,
      mis = 1 - rbinom(n(), 1, Prop),
      EQ5D = case_when(
        ti > disuti ~ 1, 
        T ~ 0.5
      )
    ) %>% 
    ungroup() %>% 
    filter(!mis) %>% 
    select(
      SID, PID, Agp, Age, visit, ti, time_points, EQ5D, rate, disuti
    )
  
  
  
  data_tte <- format_tte(sims_at)
  pars_tte <- fit_tte(model_tte, data_tte, n_iter = 5000)
  
  pars_tte$Ext %>% repli %>% mutate(Var = v)
}) %>% 
  bind_rows()

write_csv(summary_orig, here::here("docs", "test", "sims_t_orig.csv"))



summary_nomiss <- lapply(paste0("Var", 1:20), function(v) {
  sims_at <- sims_a %>% 
    mutate(disuti = rexp(n(), rate = rate)) %>% 
    left_join(prof_t %>% select(SID, visit, t_l, t_u, Prop), relationship = "many-to-many") %>% 
    mutate(
      time_points = floor(runif(n(), t_l, t_u + 1)),
      ti = time_points / 365.25,
      mis = 1 - rbinom(n(), 1, Prop),
      EQ5D = case_when(
        ti > disuti ~ 1, 
        T ~ 0.5
      )
    ) %>% 
    ungroup() %>% 
    select(
      SID, PID, Agp, Age, visit, ti, time_points, EQ5D, rate, disuti
    )
  
  
  
  data_tte <- format_tte(sims_at)
  pars_tte <- fit_tte(model_tte, data_tte, n_iter = 5000)
  
  pars_tte$Ext %>% repli %>% mutate(Var = v)
}) %>% 
  bind_rows()

write_csv(summary_nomiss, here::here("docs", "test", "sims_t_nomiss.csv"))



### Distributed with data 


summary_freq <- lapply(paste0("Var", 1:20), function(v) {
  sims_at <- sims_a %>%
    mutate(disuti = rexp(n(), rate = rate)) %>% 
    crossing(
      tibble(visit = 1:54, time_points = 5 * (1:54))
    ) %>% 
    mutate(
      ti = time_points / 365.25,
      EQ5D = case_when(
        ti > disuti ~ 1, 
        T ~ 0.5
      )
    ) %>% 
    ungroup() %>% 
    select(
      SID, PID, Agp, Age, visit, ti, time_points, EQ5D, rate, disuti
    )
  
  
  data_tte <- format_tte(sims_at)
  pars_tte <- fit_tte(model_tte, data_tte, n_iter = 5000)
  
  pars_tte$Ext %>% repli %>% mutate(Var = v)
}) %>% 
  bind_rows()

write_csv(summary_freq, here::here("docs", "test", "sims_t_freq.csv"))


### Distributed with data, but shorter 


summary_pure <- lapply(paste0("Var", 1:20), function(v) {
  ds <- sims_a %>%
    mutate(disuti = rexp(n(), rate = rate))
  
  ds <- list(
    N = nrow(ds),
    Ts = ds$disuti,
    As = ds$Age
  )
  
  
  pars_tte <- sampling(model_exp, data = ds, iter = 5000, warmup = 4500)
  pars_tte <- restructure_stan(pars_tte)
  
  pars_tte$Ext %>% repli %>% mutate(Var = v)
}) %>% 
  bind_rows()


write_csv(summary_pure, here::here("docs", "test", "sims_t_pure.csv"))


### Results


summary_orig <- read_csv(here::here("docs", "test", "sims_t_orig.csv"))
summary_nomiss <- read_csv(here::here("docs", "test", "sims_t_nomiss.csv"))
summary_freq <- read_csv(here::here("docs", "test", "sims_t_freq.csv"))
summary_pure <- read_csv(here::here("docs", "test", "sims_t_pure.csv"))


sims0 <- read_csv(here::here("data", "inputs_r1_post", "summary_qloss_ph_orig.csv")) %>% 
  filter(Age <= 95)




g <- ggpubr::ggarrange(
  summary_orig %>% 
    ggplot() +
    geom_ribbon(aes(x = Age, ymin = durr_L, ymax = durr_U, group = Var), alpha = 0.03, fill = 2) +
    geom_line(aes(x = Age, y = durr_M, group = Var), colour = 2) +
    geom_ribbon(data = sims0, aes(x = Age, ymin = dur_L, ymax = dur_U), alpha = 0.3) +
    geom_line(data = sims0, aes(x = Age, y = dur_M)) +
    scale_y_continuous("Disutility period, year") +
    expand_limits(y = c(0, 1)) + 
    labs(subtitle = "Original time-point distribution"),
  summary_freq %>% 
    ggplot() +
    geom_ribbon(aes(x = Age, ymin = durr_L, ymax = durr_U, group = Var), alpha = 0.03, fill = 2) +
    geom_line(aes(x = Age, y = durr_M, group = Var), colour = 2) +
    geom_ribbon(data = sims0, aes(x = Age, ymin = dur_L, ymax = dur_U), alpha = 0.3) +
    geom_line(data = sims0, aes(x = Age, y = dur_M)) +
    scale_y_continuous("Disutility period, year") +
    expand_limits(y = c(0, 1)) + 
    labs(subtitle = "Visits every five days"),
  summary_nomiss %>% 
    ggplot() +
    geom_ribbon(aes(x = Age, ymin = durr_L, ymax = durr_U, group = Var), alpha = 0.03, fill = 2) +
    geom_line(aes(x = Age, y = durr_M, group = Var), colour = 2) +
    geom_ribbon(data = sims0, aes(x = Age, ymin = dur_L, ymax = dur_U), alpha = 0.3) +
    geom_line(data = sims0, aes(x = Age, y = dur_M)) +
    scale_y_continuous("Disutility period, year") +
    expand_limits(y = c(0, 1)) + 
    labs(subtitle = "Original time-point distribution, no missing"),
  summary_pure %>% 
    ggplot() +
    geom_ribbon(aes(x = Age, ymin = durr_L, ymax = durr_U, group = Var), alpha = 0.03, fill = 2) +
    geom_line(aes(x = Age, y = durr_M, group = Var), colour = 2) +
    geom_ribbon(data = sims0, aes(x = Age, ymin = dur_L, ymax = dur_U), alpha = 0.3) +
    geom_line(data = sims0, aes(x = Age, y = dur_M)) +
    scale_y_continuous("Disutility period, year") +
    expand_limits(y = c(0, 1)) + 
    labs(subtitle = "Fully observable disutility period")
)


ggsave(g, filename = here::here("docs", "test", "test_timepoints.pdf"), width = 12, heigh = 10)


