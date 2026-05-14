library(tidyverse)
library(rstan)

theme_set(theme_bw())

## Loading stan models and settings
options(mc.cores = 8)
rstan_options(auto_write = TRUE)


for (file in dir("R")) source(here::here("R", file))

model_t2z <- rstan::stan_model(here::here("models", "time2zero_age.stan"))

## Setup
n_iter <- 2000
n_warmup <- 1500


root_src <- here::here("out", "gen_2gp")


for(v in paste0("s", 1:10)) {
  cat(v, "\n")
  
  data_tte <- read_csv(here::here(root_src, "generated_" + glue::as_glue(v) + ".csv")) %>% 
    filter(visit == 1) %>% 
    select(Age, TTe = disuti)
  
  
  dat <- list(
    N = nrow(data_tte),
    As = data_tte$Age,
    Ts = data_tte$TTe
  )

  
  tryCatch({
    post <- rstan::sampling(model_t2z, data = dat, iter = n_iter, warmup = n_warmup)
    
    ext <- rstan::extract(post, pars = c("ba1", "r0")) %>% data.frame()
    
    write_csv(ext, here::here(root_src, "post_tte_exp_" + glue::as_glue(v) + ".csv"))

    #when it throws an error, the following block catches the error
  }, error = function(msg){
    return(NULL)
  })
  
}




