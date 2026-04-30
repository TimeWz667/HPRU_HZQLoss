



res_prev <- read_csv(here::here("docs", "tabs", "summary_qloss_ph.csv"))

res_2gp <- read_csv(here::here("docs", "full", "sim_qloss_2gp.csv"))

res_ind <- read_csv(here::here("docs", "full", "sim_qloss_ind.csv"))



g_ql <- res_ind %>% 
  ggplot(aes(x = Age)) + 
  geom_ribbon(aes(ymin = QL35_L, ymax = QL35_U, group = Var), alpha = 0.01) +
  geom_line(aes(y = QL35_M * exp(-0.002 * Age), group = Var)) +
  geom_ribbon(data = res_prev, aes(ymin = QL35_L, ymax = QL35_U, fill = "Old ver."), alpha = 0.2) +
  geom_line(data = res_prev, aes(y = QL35_M, colour = "Old ver."))  + expand_limits(y = 0)

g_ql

g_ql <- res_2gp %>% 
  ggplot(aes(x = Age)) + 
  geom_ribbon(aes(ymin = QL35_L, ymax = QL35_U, group = Var), alpha = 0.01) +
  geom_line(aes(y = QL35_M, group = Var)) +
  geom_ribbon(data = res_prev, aes(ymin = QL35_L, ymax = QL35_U, fill = "Old ver."), alpha = 0.2) +
  geom_line(data = res_prev, aes(y = QL35_M, colour = "Old ver.")) + expand_limits(y = 0)



g_ql


pars_qol <- read_csv(here::here("data", "inputs_for_gen", "pars_qol.csv"))

## Background population
list_studies <- read_csv(here::here("data", "list_studies.csv")) %>% 
  select(study, SID, SIK)
list_studies


root_src <- here::here("data", "inputs_for_gen")

## Shared inputs

prof_a <- read_csv(here::here(root_src, "prof_a.csv"))


pars_tte <- read_csv(here::here("data", "inputs_for_gen", "pars_tte_ind.csv"))
pars_tte_slopes <- read_csv(here::here("data", "inputs_for_gen", "pars_tte_slopes_ind.csv"))

du_ind <- prof_a %>% 
  uncount(n) %>% 
  mutate(
    PID = paste0(SID, "_", 1:n()),
    A0 = pmax(A0, 18),
    Age = floor(runif(n(), A0, A1))
  ) %>% 
  select(-c("Agp", "A0", "A1")) %>% 
  full_join(bind_rows(pars_tte, pars_qol) %>% select(- Model), by = "SID", relationship = "many-to-many")  %>% 
  mutate(value = rnorm(n(), mean, sd)) %>% 
  select(-c(mean, sd)) %>% 
  pivot_wider(names_from = Var, values_from = value) %>% 
  left_join(pars_tte_slopes) %>% 
  mutate(
    r0 = pmax(r0, 0.0001),
    ba1_sim = b0 + log(r0) * b_lr0, # + rnorm(n(), 0, rsd),
    rate = r0 * exp(Age * ba1_sim),
    disuti = 1 / rate
  ) %>% 
  group_by(Age) %>% 
  summarise(
    DU_M = median(disuti),
    DU_L = quantile(disuti, 0.025),
    DU_U = quantile(disuti, 0.975)
  ) %>% 
  mutate(
    Type = "Ind"
  )



pars_tte <- read_csv(here::here("data", "inputs_for_gen", "pars_tte_2gp.csv"))
pars_tte_slopes <- read_csv(here::here("data", "inputs_for_gen", "pars_tte_slopes_2gp.csv"))

du_2gp <- prof_a %>% 
  uncount(n) %>% 
  mutate(
    PID = paste0(SID, "_", 1:n()),
    A0 = pmax(A0, 18),
    Age = floor(runif(n(), A0, A1))
  ) %>% 
  select(-c("Agp", "A0", "A1")) %>% 
  full_join(bind_rows(pars_tte, pars_qol) %>% select(- Model), by = "SID", relationship = "many-to-many")  %>% 
  mutate(value = rnorm(n(), mean, sd)) %>% 
  select(-c(mean, sd)) %>% 
  pivot_wider(names_from = Var, values_from = value) %>% 
  left_join(pars_tte_slopes) %>% 
  mutate(
    r0 = pmax(r0, 0.0001),
    ba1_sim = b0 + log(r0) * b_lr0, # + rnorm(n(), 0, rsd),
    rate = r0 * exp(Age * ba1_sim),
    disuti = 1 / rate
  ) %>% 
  group_by(Age) %>% 
  summarise(
    DU_M = median(disuti),
    DU_L = quantile(disuti, 0.025),
    DU_U = quantile(disuti, 0.975)
  ) %>% 
  mutate(
    Type = "2GP"
  )



bind_rows(du_ind, du_2gp) %>% 
  ggplot() +
  geom_ribbon(aes(x = Age, ymin = DU_L, ymax = DU_U, fill = Type), alpha = 0.1) +
  geom_line(aes(x = Age, y = DU_M, colour = Type)) + 
  geom_ribbon(data = res_prev, aes(x = Age, ymin = dur_L, ymax = dur_U, fill = "Prev"), alpha = 0.1) +
  geom_line(data = res_prev, aes(x = Age, y = dur_M, colour = "Prev")) 








src_slopes_r1 <- read_csv(here::here("data", "inputs_r1_post", "stats_tte_slopes.csv")) %>% 
  filter(SID == "Combined") %>% 
  select(SID, Var, m) %>% 
  pivot_wider(names_from = Var, values_from = m)


pars_tte_r1_ph <- read.csv(here::here("data", "inputs_r1_post", "stats_ph_tte.csv"))


pars_prev <- pars_tte_r1_ph %>% 
  filter(Var == "r0") %>% 
  crossing(Key = 1:1000) %>% 
  mutate(
    r0 = rnorm(n(), mean, sd)
  ) %>% 
  select(Key, r0) %>% 
  bind_cols(src_slopes_r1) %>% 
  mutate(
    ba1_sim = b0 + log(r0) * b_lr0 + rnorm(n(), 0, rsd)
  )


src_slopes_gsk <- read_csv(here::here("data", "inputs_gsk_post", "Summary_slopes.csv")) %>% 
  tidyr::extract(file, "SID", "Post_tte_(\\w+).csv")

src_tte_gsk <- read_csv(here::here("data", "inputs_gsk_post", "Summary_tte.csv"))



pars_full <- src_tte_gsk %>% filter(SID == "FULL") %>% filter(Var == "r0") %>% 
  crossing(Key = 1:1000) %>% 
  mutate(
    r0 = rnorm(n(), mean, sd)
  ) %>% 
  select(Key, SID, r0) %>% 
  bind_cols(src_slopes_gsk %>% filter(SID == "full") %>% select(- SID)) %>% 
  mutate(
    ba1_sim = b0 + log(r0) * b_lr0 + rnorm(n(), 0, rsd)
  )



pars_gsk <- src_tte_gsk %>% filter(SID == "GSK") %>% filter(Var == "r0") %>% 
  crossing(Key = 1:1000) %>% 
  mutate(
    r0 = rnorm(n(), mean, sd)
  ) %>% 
  select(Key, SID, r0) %>% 
  bind_cols(src_slopes_gsk %>% filter(SID == "gsk") %>% select(- SID)) %>% 
  mutate(
    ba1_sim = b0 + log(r0) * b_lr0 + rnorm(n(), 0, rsd)
  )


plot(pars_full$r0, pars_full$ba1_sim)
plot(pars_gsk$r0, pars_gsk$ba1_sim)



sims <- bind_rows(
  pars_full %>% 
    crossing(Age = 50:95) %>% 
    mutate(
      rate = r0 * exp(Age * ba1_sim),
      dur = 1 / rate
    ) %>% 
    summarise(across(c(rate, dur), list(
      M = \(x) median(x),
      L = \(x) quantile(x, 0.025),
      U = \(x) quantile(x, 0.975)
    )), .by = Age) %>% 
    mutate(
      Group = "Full"
    ),
  pars_gsk %>% 
    crossing(Age = 50:95) %>% 
    mutate(
      rate = r0 * exp(Age * ba1_sim),
      dur = 1 / rate
    ) %>% 
    summarise(across(c(rate, dur), list(
      M = \(x) median(x),
      L = \(x) quantile(x, 0.025),
      U = \(x) quantile(x, 0.975)
    )), .by = Age) %>% 
    mutate(
      Group = "GSK"
    ),
  pars_prev %>% 
    crossing(Age = 50:95) %>% 
    mutate(
      rate = r0 * exp(Age * ba1_sim),
      dur = 1 / rate
    ) %>% 
    summarise(across(c(rate, dur), list(
      M = \(x) median(x),
      L = \(x) quantile(x, 0.025),
      U = \(x) quantile(x, 0.975)
    )), .by = Age) %>% 
    mutate(
      Group = "Prev"
    )
)


sims %>% 
  ggplot() +
  geom_ribbon(aes(x = Age, ymin = dur_L, ymax = dur_U, fill = Group), alpha = 0.1) +
  geom_line(aes(x = Age, y = dur_M, colour = Group))


res_prev <- read_csv(here::here("docs", "tabs", "summary_qloss_ph.csv"))


res <- read_csv(here::here("docs", "full", "sim_qloss_2gp.csv"))
res_prev %>% 
  ggplot() +
  geom_ribbon(aes(x = Age, ymin = dur_L, ymax = dur_U), alpha = 0.1) +
  geom_line(aes(x = Age, y = dur_M), alpha = 0.4)  +
  geom_ribbon(data = sims, aes(x = Age, ymin = dur_L, ymax = dur_U, fill = Group), alpha = 0.1) +
  geom_line(data = sims, aes(x = Age, y = dur_M, colour = Group))+
  expand_limits(y = 0)



gen <- read_csv(here::here("out", "gen_ind", "generated_s1.csv"))

gen %>% group_by(SID) %>% 
  filter(ti < disuti) %>% 
  summarise(
    mean(rate),
    mean(EQ5D)
  )











