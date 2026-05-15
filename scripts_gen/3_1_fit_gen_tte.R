library(tidyverse)
library(rstan)

theme_set(theme_bw())

## Loading stan models and settings
options(mc.cores = 8)
rstan_options(auto_write = TRUE)


for (file in dir("R")) source(here::here("R", file))

model_t2z <- rstan::stan_model(here::here("models", "time2zero_age.stan"))

 
list_studies <- read_csv(here::here("data", "list_studies.csv"))
list_studies

## Functions

fit_tte <- function(data_tte, n_iter, n_warmup) {
  dat <- list(
    N = nrow(data_tte),
    As = data_tte$Age,
    Ts = data_tte$TTE
  )
  
  
  tryCatch({
    post <- rstan::sampling(model_t2z, data = dat, iter = n_iter, warmup = n_warmup)
    ext <- rstan::extract(post, pars = c("ba1", "r0")) %>% data.frame()
    return(ext)
    #when it throws an error, the following block catches the error
  }, error = function(msg){
    return(NULL)
  })
}




## Setup
n_iter <- 1000
n_warmup <- 900

check_existence <- TRUE


for (type in c("ind", "2gp")) {
  type <- glue::as_glue(type)
  print(type)

  root_src <- here::here("out", "gen_" + type)
  
  last_version <- read_file(here::here(root_src, "last_version.txt")) %>% glue::as_glue()
  
  gens <- dir(root_src)
  gens <- gens[endsWith(gens, last_version + ".csv")]
  gens <- gens[startsWith(gens, "generated")]
  
  for (gen in gens) {
    v <- gsub("generated_", "", gen)
    v <- gsub("_" + last_version + ".csv", "", v)
    cat(type, v, "\n")
    out_file <- here::here(root_src, "post_tte_exp_" + glue::as_glue(v) + "_" + last_version + ".csv")
    
    if (!(check_existence & file.exists(out_file))) {
      ext <- fit_tte(data_tte, n_iter = n_iter, n_warmup = n_warmup)
      write_csv(ext, here::here(root_src, "post_tte_exp_" + glue::as_glue(v) + "_" + last_version + ".csv"))
    }
  }
}


## Check all done
for (type in c("ind", "2gp")) {
  type <- glue::as_glue(type)
  root_src <- here::here("out", "gen_" + type)
  
  last_version <- read_file(here::here(root_src, "last_version.txt")) %>% glue::as_glue()
  
  gens <- dir(root_src)
  gens <- gens[endsWith(gens, last_version + ".csv")]
  gens <- gens[startsWith(gens, "generated")]
  
  done <- sapply(gens, \(gen) {
    v <- gsub("generated_", "", gen)
    v <- gsub("_" + last_version + ".csv", "", v)
    
    out_file <- here::here(root_src, "post_tte_exp_" + glue::as_glue(v) + "_" + last_version + ".csv")
    
    return(file.exists(out_file))
  })
  
  cat("Type: ", type, "-> ", ifelse(all(done), "Completed", "Incompleted"), "\n")
}






