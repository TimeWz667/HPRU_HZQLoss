library(tidyverse)
library(rstan)

theme_set(theme_bw())

options(mc.cores = 8)
rstan_options(auto_write = TRUE)


for (file in dir("R")) source(here::here("R", file))



# Setup
n_iter <- 2000
n_warmup <- 1500

check_existence <- TRUE


fit_qol <- function(data_raw, n_iter, n_warmup) {
  
  max_qol <- max(data_raw$EQ5D)
  min_qol <- min(data_raw$EQ5D)
  
  data_qol <- data_raw %>% 
    filter(ti < disuti) %>% 
    select(SID, PID, Age, Agp, ti, EQ5D) %>% 
    mutate(
      Q_rescaled = (EQ5D - min_qol) / (max_qol - min_qol)
    )
  
  
  tryCatch({
    pars_qol <- fit_qol_bayes(data_qol, ind = F, n_iter = n_iter)
    return(pars_qol$Ext)
  }, error = function(msg){
    return(NULL)
  })
}


### Fitting
for (type in c("ind", "2gp", "ind_km", "2gp_km")) {
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
    
    data_raw <- read_csv(here::here(root_src, gen))
    
    out_file <- here::here(root_src, "post_qol_" + glue::as_glue(v) + "_" + last_version + ".csv")
    
    if (!(check_existence & file.exists(out_file))) {
      cat(type, v, "\n")
      ext <- fit_qol(data_raw, n_iter = n_iter, n_warmup = n_warmup)
      if (!is.null(ext)) write_csv(ext, out_file)
    }
  }
}


## Check all done
for (type in c("ind", "2gp", "ind_km", "2gp_km")) {
  type <- glue::as_glue(type)
  root_src <- here::here("out", "gen_" + type)
  
  last_version <- read_file(here::here(root_src, "last_version.txt")) %>% glue::as_glue()
  
  gens <- dir(root_src)
  gens <- gens[endsWith(gens, last_version + ".csv")]
  gens <- gens[startsWith(gens, "generated")]
  
  done <- sapply(gens, \(gen) {
    v <- gsub("generated_", "", gen)
    v <- gsub("_" + last_version + ".csv", "", v)
    
    out_file <- here::here(root_src, "post_qol_" + glue::as_glue(v) + "_" + last_version + ".csv")
    
    return(file.exists(out_file))
  })
  
  cat("Type: ", type, "-> ", ifelse(all(done), "Completed", "Incompleted"), "\n")
}

