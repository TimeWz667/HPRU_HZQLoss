library(tidyverse)


theme_set(theme_bw())




n_mc <- 2000



## Round 1

post_tte_r1 <- read_csv(here::here("posteriors", "post_ph_tte.csv"))

post_tte_r1 %>% 
  crossing(Age = seq(50, 80, 10)) %>% 
  mutate(
    rate = r0 * exp(Age * ba1)
  ) %>% 
  summarise(
    M = median(rate),
    L = quantile(rate, 0.025),
    U = quantile(rate, 0.975),
    .by = "Age"
  )



## GSK

src_slopes_gsk <- read_csv(here::here("data", "inputs_gsk_post", "Summary_slopes.csv")) %>% 
  tidyr::extract(file, "SID", "Post_tte_(\\w+).csv") %>% 
  mutate(SID = str_to_upper(SID))

src_tte_gsk <- read_csv(here::here("data", "inputs_gsk_post", "Summary_tte.csv"))

pars_tte_gsk <- src_tte_gsk %>% 
  filter(SID == "GSK") %>% 
  crossing(PID = 1:n_mc) %>% 
  mutate(
    value = rnorm(n(), mean, sd)
  ) %>% 
  select(name = Var, value, PID) %>% 
  pivot_wider() %>% 
  bind_cols(src_slopes_gsk %>% filter(SID == "GSK"))  %>% 
  mutate(
    ba1 = log(r0) * b_lr0 + b0 + rnorm(n(), 0, rsd)
  ) %>% 
  select(r0, ba1)

  


pars_tte_full <- src_tte_gsk %>% 
  filter(SID == "FULL") %>% 
  crossing(PID = 1:n_mc) %>% 
  mutate(
    value = rnorm(n(), mean, sd)
  ) %>% 
  select(name = Var, value, PID) %>% 
  pivot_wider() %>% 
  bind_cols(src_slopes_gsk %>% filter(SID == "FULL"))  %>% 
  mutate(
    ba1 = log(r0) * b_lr0 + b0 + rnorm(n(), 0, rsd)
  ) %>% 
  select(r0, ba1)

  
  
  
  
sims <- bind_rows(
  post_tte_r1 %>% mutate(Case = "R1"), 
  pars_tte_gsk %>% mutate(Case = "GSK"),
  pars_tte_full %>% mutate(Case = "Full")
) %>% 
  crossing(Age = seq(50, 90, 1)) %>% 
  mutate(
    rate = r0 * exp(Age * ba1)
  ) %>% 
  summarise(
    M = median(1 / rate),
    L = quantile(1 / rate, 0.025),
    U = quantile(1 / rate, 0.975),
    .by = c("Age", "Case")
  )


sims %>% 
  ggplot() + 
  geom_ribbon(aes(x = Age, ymin = L, ymax = U, fill = Case), alpha = 0.2) +
  geom_line(aes(x = Age, y = M, colour = Case)) +
  facet_grid(.~Case) + 
  scale_y_continuous("Disutility period, months", labels = scales::number_format(scale = 12)) +
  expand_limits(y = 0)


  

root_src <- here::here("out", "gen_2gp")


post <- read_csv(here::here(root_src, "post_tte_exp_s1.csv"))


sims_gen <- bind_rows(
  read_csv(here::here(root_src, "post_tte_exp_s1.csv")) %>% mutate(Ver = "S1"), 
  read_csv(here::here(root_src, "post_tte_exp_s2.csv")) %>% mutate(Ver = "S2"), 
  read_csv(here::here(root_src, "post_tte_exp_s3.csv")) %>% mutate(Ver = "S3"), 
  read_csv(here::here(root_src, "post_tte_exp_s4.csv")) %>% mutate(Ver = "S4"), 
  read_csv(here::here(root_src, "post_tte_exp_s5.csv")) %>% mutate(Ver = "S5"), 
) %>% 
  crossing(Age = seq(50, 90, 1)) %>% 
  mutate(
    rate = r0 * exp(Age * ba1)
  ) %>% 
  summarise(
    M = median(1 / rate),
    L = quantile(1 / rate, 0.025),
    U = quantile(1 / rate, 0.975),
    .by = c("Age", "Ver")
  )

sims_gen %>% 
  ggplot() + 
  geom_ribbon(data = sims %>% filter(Case == "Full"), aes(x = Age, ymin = L, ymax = U), alpha = 0.2) +
  geom_line(data = sims %>% filter(Case == "Full"), aes(x = Age, y = M)) +
  geom_ribbon(aes(x = Age, ymin = L, ymax = U, fill = Ver), alpha = 0.2) +
  geom_line(aes(x = Age, y = M, colour = Ver)) +
  
  scale_y_continuous("Disutility period, months", labels = scales::number_format(scale = 12)) +
  expand_limits(y = 0)


