library(tidyverse)




pars <- list(
  qol0 = 0.9,
  qol_diff = 0,
  qol_drop = 0.3,
  rec = 0.10,
  err = 0.05
)


sims <- with(pars, {
  qol_final = qol0 - qol_diff
  qol_dz = qol0 - qol_drop
  red = qol_final - qol_dz
  
  crossing(
    Key = 1:100,
    ti = -1:10
  ) %>% 
  mutate(
    qol = qol_final - red * exp(- rec * (ti - 0)),
    qol = ifelse(ti < 0, qol0, qol),
    qol = rnorm(n(), qol, err),
    .by = Key
  )
  
}) 


sims %>% 
  ggplot() +
  geom_line(aes(x = ti, y = qol, group = Key))



library(rstan)




model <- rstan::stan_model(here::here("models", "response_shift.stan"))


dat <- list(
  N = nrow(sims),
  Ts = sims$ti,
  Qs = sims$qol
)


post <- rstan::sampling(model, data = dat)


su <- summary(post)$summary


fitted <- bind_cols(
  Var = rownames(su),
  su
) %>% 
  filter(Var %in% names(pars)) %>% 
  mutate(TV = sapply(Var, \(x) pars[[x]]))


fitted %>% 
  ggplot() +
  geom_pointrange(aes(x = 0, y = `50%`, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_linerange(aes(x = 0, ymin = `25%`, ymax = `75%`), linewidth = 1.3) +
  geom_point(aes(x = 0, y = TV), size = 10, alpha = 0.2) +
  facet_wrap(~Var, scales = "free_y") +
  theme(axis.text.x = element_blank())







