library(tidyverse)
library(rstan)

theme_set(theme_bw())



## Loading stan models and settings
options(mc.cores = 8)
rstan_options(auto_write = TRUE)

model_t2z <- rstan::stan_model(here::here("models", "time2zero_age.stan"))



## Setup
n_iter <- 1000
n_warmup <- 900
n_batch <- 10

pss <- tibble(r0 = c(3, 5), ba1 = c(-0.025, -0.05))

ns_sample <- seq(100, 5000, 200)



## Run
results <- lapply(1:nrow(pss), \(k) {
  pars <- as.list(pss[k, ])
  print(pars)
  
  lapply(ns_sample, \(n_sample) {
    lapply(1:n_batch, \(run) {
      sims0 <- tibble(
        As = floor(runif(n_sample, 40, 80))
      ) %>% 
        mutate(
          rate = pars$r0 * exp(As * pars$ba1),
          Ts = rexp(n(), rate = rate)    
        )
      
      dat <- list(
        N = nrow(sims0),
        Ts = sims0$Ts,
        As = sims0$As
      )
      
      post <- rstan::sampling(model_t2z, data = dat, iter = n_iter, warmup = n_warmup)
      
      data.frame(as.table(summary(post)$summary)) %>% 
        rename(parameter = Var1, stats = Var2) %>% 
        filter(parameter %in% names(pars)) %>% 
        pivot_wider(names_from = stats, values_from = Freq) %>% 
        left_join(pss[k, ] %>% pivot_longer(everything(), names_to = "parameter", values_to = "true_value"), by = "parameter") %>% 
        mutate(Scenario = k, Batch = run, N_Sample = n_sample)
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows()


write_csv(results, here::here("docs", "experiments", "res_sample_size.csv"))


results %>% 
  group_by(parameter, Scenario, N_Sample) %>%
  summarise(
    m = mean(mean - true_value),
    l = mean(`25%` - true_value),
    u = mean(`75%` - true_value)
  ) %>% 
  ggplot() + 
  geom_pointrange(aes(x = N_Sample, y = m, ymin = l, ymax = u), position = position_dodge2(30)) +
  geom_hline(yintercept = 0) + 
  scale_y_continuous("Errors") +
  facet_wrap(parameter~Scenario, scale = "free_y")

g <- results %>% 
  ggplot() + 
  geom_pointrange(aes(x = N_Sample - true_value, y = mean - true_value, ymin = `25%` - true_value, ymax = `75%` - true_value), position = position_dodge2(30)) +
  geom_hline(yintercept = 0) + 
  scale_y_continuous("Errors") +
  facet_wrap(parameter~Scenario, scale = "free_y")


g



