library(tidyverse)
library(rstan)

theme_set(theme_bw())

options(mc.cores = 8)
rstan_options(auto_write = TRUE)


for (file in dir("R")) source(here::here("R", file))

model_tte <- stan_model(here::here("models", "time2zero_surv_age.stan"))

n_mc <- 1000


## Scenario: individual tte
root_src <- here::here("out", "gen_ind")


for(v in paste0("s", 1:100)) {
  cat(v, "\n")
  data_raw <- get_data_qol(here::here(root_src, "generated_" + glue::as_glue(v) + ".csv"))
  
  tryCatch({
    #data_tte <- format_tte(data_raw)
    data_qol <- format_qol(data_raw) %>% sample_n(3000)
    
    #pars_tte <- fit_tte(model_tte, data_tte, n_iter = 2000)
    pars_qol <- fit_qol_bayes(data_qol, ind = F, n_iter = 2000)
    
    #write_csv(pars_tte$Ext, here::here(root_src, "post_tte_" + glue::as_glue(v) + ".csv"))
    write_csv(pars_qol$Ext, here::here(root_src, "post_qol_" + glue::as_glue(v) + ".csv"))

      #when it throws an error, the following block catches the error
  }, error = function(msg){
    return(NULL)
  })
  
}


## Scenario: individual tte
root_src <- here::here("out", "gen_2gp")

for(v in paste0("s", 1:100)) {
  cat(v, "\n")
  data_raw <- get_data_qol(here::here(root_src, "generated_" + glue::as_glue(v) + ".csv"))
  
  tryCatch({
    #data_tte <- format_tte(data_raw)
    data_qol <- format_qol(data_raw) %>% sample_n(3000)
    
    #pars_tte <- fit_tte(model_tte, data_tte, n_iter = 2000)
    pars_qol <- fit_qol_bayes(data_qol, ind = F, n_iter = 2000)
    
    #write_csv(pars_tte$Ext, here::here(root_src, "post_tte_" + glue::as_glue(v) + ".csv"))
    write_csv(pars_qol$Ext, here::here(root_src, "post_qol_" + glue::as_glue(v) + ".csv"))
    
    #when it throws an error, the following block catches the error
  }, error = function(msg){
    return(NULL)
  })
  
}




## Simulate QALY loss

for(gp in c("ind", "2gp")) {
  gp <- glue::as_glue(gp)
  
  root_src <- here::here("out", "gen_" + gp)
  vs <- dir(root_src)
  vs <- vs[startsWith(vs, "post_tte")]
  vs <- gsub("post_tte_", "", vs)
  vs <- gsub(".csv", "", vs)
  
  
  res <- bind_rows(lapply(vs, \(v) {
    sims <- boot_pars_bayes(here::here(root_src, "post_tte_" + glue::as_glue(v) + ".csv"), 
                            here::here(root_src, "post_qol_" + glue::as_glue(v) + ".csv"), n_sim = n_mc) %>% 
      simulate_ql(age0 = 50, age1 = 95, age_until = 5, dt = 0.01) %>% 
      summarise(across(starts_with("QL"), list(
        M = \(x) median(x),
        L = \(x) quantile(x, 0.025),
        U = \(x) quantile(x, 0.975)
      )), .by = Age) %>% 
      mutate(Var = v)
    
    sims
  }))
  
  write_csv(res, here::here("docs", "full", "sim_qloss_" + gp + ".csv"))

  
}



g_ql <- res %>% 
  ggplot(aes(x = Age)) + 
  geom_ribbon(aes(ymin = QL35_L, ymax = QL35_U, group = Var), alpha = 0.2) +
  geom_line(aes(y = QL35_M, group = Var))



