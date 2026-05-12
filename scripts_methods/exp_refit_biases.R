library(tidyverse)
library(rstan)

theme_set(theme_bw())



model_t2z <- rstan::stan_model(here::here("models", "time2zero_age.stan"))
model_t2z_cen <- rstan::stan_model(here::here("models", "time2zero_age.stan"))





n_sample <- 2000
n_iter <- 2000
n_warmup <- 1800


pss <- crossing(r0 = c(1, 3, 5), ba1 = c(0, -0.025, -0.05))


results <- lapply(1:nrow(pss), \(k) {
  pars <- as.list(pss[k, ])
  print(pars)
  
  
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
  
  data.frame(as.table(summary(post)$c_summary)) %>% 
    filter(parameter %in% names(pars)) %>% 
    pivot_wider(names_from = stats, values_from = Freq) %>% 
    left_join(pss[k, ] %>% pivot_longer(everything(), names_to = "parameter", values_to = "true_value")) %>% 
    mutate(Scenario = k)
  
}) %>% 
  bind_rows()


results %>% 
  ggplot() + 
  geom_pointrange(aes(x = Scenario, y = mean, ymin = `25%`, ymax = `75%`, colour = chains), position = position_dodge2(0.3)) +
  geom_linerange(aes(xmin = Scenario - 0.3, xmax = Scenario + 0.3, y = true_value)) +
  facet_wrap(.~parameter, scales = "free")















