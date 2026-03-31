library(tidyverse)
library(rstan)


theme_set(theme_bw())


options(mc.cores = 8)
rstan_options(auto_write = TRUE)


for (file in dir("R")) source(here::here("R", file))



model_tte <- stan_model(here::here("models", "time2zero_surv_age.stan"))

n_mc <- 1000


list_studies <- read_csv(here::here("data", "list_studies.csv")) %>% 
  select(study, SID, SIK)
list_studies




for (sid in c("RamKOR", "RamMEX", "RamTWN", "RamTHA", "ScottUK")) {
  collector <- list()
  
  for (v in c("s1", "s2", "s3", "s4",  "s5")) {
    cat(sid, v, "\n")
    d0 <- get_data_qol(here::here("data", "generated", "data_generated_" + glue::as_glue(v) + ".csv"))
    
    data_raw <- d0 %>% filter(SID == sid)
    
    data_tte <- format_tte(data_raw)
    data_qol <- format_qol(data_raw)
    
    pars_tte <- fit_tte(model_tte, data_tte, n_iter = 5000)
    
    pars_qol <- fit_qol_bayes(data_qol, ind = T)
    
    write_csv(pars_tte$Ext, here::here("out", "post_tte.csv"))
    write_csv(pars_qol$Ext, here::here("out", "post_qol.csv"))

    sims <- boot_pars_bayes(here::here("out", "post_tte.csv"), 
                            here::here("out", "post_qol.csv"), n_sim = n_mc) %>% 
      simulate_ql(age0 = 50, age1 = 95, age_until = 5, dt = 0.01) %>% 
      summarise(
        M = median(QL35),
        L = quantile(QL35, 0.025),
        U = quantile(QL35, 0.975), .by = Age
      ) %>% 
      mutate(Var = v, SID = sid)
    
    collector[[length(collector) + 1]] <- sims
  }
  
  
  res <- bind_rows(collector)
  
  write_csv(res, here::here("docs", "test", "sim_qloss_" + glue::as_glue(sid) + ".csv"))
  
  g <- res %>% 
    ggplot() +
    geom_ribbon(aes(x = Age, ymin = L, ymax = U, fill = Var), alpha = 0.2) +
    geom_line(aes(x = Age, y = M, colour = Var)) +
    scale_y_continuous("QALY loss") +
    expand_limits(y = 0)
  
  ggsave(g, file = here::here("docs", "test", "sim_qloss_" + glue::as_glue(sid) + ".png"), width = 6, height = 4)
  
  
}


sid <- "all"
collector <- list()

for (v in c("s1", "s2", "s3", "s4",  "s5")) {
  cat(sid, v, "\n")
  d0 <- get_data_qol(here::here("data", "generated", "data_generated_" + glue::as_glue(v) + ".csv"))
  
  data_raw <- d0
  
  data_tte <- format_tte(data_raw)
  data_qol <- format_qol(data_raw)
  
  pars_tte <- fit_tte(model_tte, data_tte, n_iter = 5000)
  
  pars_qol <- fit_qol_bayes(data_qol, ind = T)
  
  write_csv(pars_tte$Ext, here::here("out", "post_tte.csv"))
  write_csv(pars_qol$Ext, here::here("out", "post_qol.csv"))
  
  sims <- boot_pars_bayes(here::here("out", "post_tte.csv"), 
                          here::here("out", "post_qol.csv"), n_sim = n_mc) %>% 
    simulate_ql(age0 = 50, age1 = 95, age_until = 5, dt = 0.01) %>% 
    summarise(
      M = median(QL35),
      L = quantile(QL35, 0.025),
      U = quantile(QL35, 0.975), .by = Age
    ) %>% 
    mutate(Var = v, SID = sid)
  
  collector[[length(collector) + 1]] <- sims
}

res <- bind_rows(collector)

write_csv(res, here::here("docs", "test", "sim_qloss_" + glue::as_glue(sid) + ".csv"))



sims0 <- boot_pars_bayes(here::here("posteriors", "post_ph_tte.csv"), 
                         here::here("posteriors", "post_qol_b_orig.csv"), n_sim = n_mc) %>% 
  simulate_ql(age0 = 50, age1 = 95, age_until = 5, dt = 0.01) %>% 
  summarise(
    M = median(QL35),
    L = quantile(QL35, 0.025),
    U = quantile(QL35, 0.975), .by = Age
  )



sims_sub <- lapply(list_studies$SID, function(sid) {

  sims0 <- boot_pars_bayes(here::here("posteriors", "post_ph_tte.csv"), 
                           here::here("posteriors", "post_qol_b_orig_sub_" + glue::as_glue(sid) + ".csv"), n_sim = n_mc) %>% 
    simulate_ql(age0 = 50, age1 = 95, age_until = 5, dt = 0.01) %>% 
    summarise(
      M = median(QL35),
      L = quantile(QL35, 0.025),
      U = quantile(QL35, 0.975), .by = Age
    ) %>% 
  mutate(
    SID = sid
  )
  
  sims0
  
})
sims_sub <- sims_sub %>% bind_rows()


for (sid in list_studies$SID) {
  if (sid != "RamCRI") {
    res <- read_csv(here::here("docs", "test", "sim_qloss_" + glue::as_glue(sid) + ".csv"))
    
    ss <- sims_sub %>% filter(SID == sid)
    g <- res %>% 
      ggplot() +
      geom_ribbon(aes(x = Age, ymin = L, ymax = U, fill = Var), alpha = 0.2) +
      geom_line(aes(x = Age, y = M, colour = Var)) +
      geom_ribbon(data = ss, aes(x = Age, ymin = L, ymax = U, fill = "Sub", colour = "Sub"), alpha = 0.7, linewidth = 1.5, linetype = 2) +
      geom_line(data = ss, aes(x = Age, y = M, colour = "Sub"), linewidth = 1.5) +
      geom_ribbon(data = sims0, aes(x = Age, ymin = L, ymax = U, fill = "Full", colour = "Full"), alpha = 0.7, linewidth = 1.5, linetype = 2) +
      geom_line(data = sims0, aes(x = Age, y = M, colour = "Full"), linewidth = 1.5) +
      scale_y_continuous("QALY loss") +
      expand_limits(y = 0) + 
      guides(fill = guide_none()) +
      labs(subtitle = sid)
    
    g
    ggsave(g, file = here::here("docs", "test", "sim_qloss_" + glue::as_glue(sid) + ".png"), width = 6, height = 4)
    
    
  }

}



d0 <- read_csv(here::here("data", "generated", "data_generated_" + glue::as_glue(v) + ".csv"))


d0 %>% 
  filter(visit == 1) %>% 
  ggplot() +
  geom_point(aes(x = Age, y = 1 / rate, colour = SID)) +
  geom_point(aes(x = Age, y = disuti)) +
  facet_wrap(.~SID)


d0 %>% 
  filter(visit == 1) %>% 
  summarise(tte = mean(disuti), .by = c(SID, Age)) %>% 
  ggplot() +
  geom_smooth(aes(x = Age, y = tte, colour = SID), method = "lm") +
  geom_point(data = d0 %>% select(Age, rate, SID) %>% distinct(), aes(x = Age, y = 1 / rate), alpha = 0.2) +
  facet_wrap(.~SID)


d0 %>% 
  filter(visit == 1) %>% 
  ggplot() +
  geom_smooth(aes(x = Age, y = disuti, colour = SID), method = "lm") +
  geom_point(data = d0 %>% select(Age, rate, SID) %>% distinct(), aes(x = Age, y = 1 / rate), alpha = 0.2) +
  facet_wrap(.~SID)


