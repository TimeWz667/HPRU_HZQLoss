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
    
    out_file <- here::here(root_src, "post_qol_" + glue::as_glue(v) + "_" + last_version + ".csv")
    
    return(file.exists(out_file))
  })
  
  cat("Type: ", type, "-> ", ifelse(all(done), "Completed", "Incompleted"), "\n")
}




# 
# ## Simulate QALY loss
# 
# for(gp in c("ind", "2gp")) {
#   gp <- glue::as_glue(gp)
#   
#   root_src <- here::here("out", "gen_" + gp)
#   vs <- dir(root_src)
#   vs <- vs[startsWith(vs, "post_tte")]
#   vs <- gsub("post_tte_", "", vs)
#   vs <- gsub(".csv", "", vs)
#   
#   
#   res <- bind_rows(lapply(vs, \(v) {
#     sims <- boot_pars_bayes(here::here(root_src, "post_tte_" + glue::as_glue(v) + ".csv"), 
#                             here::here(root_src, "post_qol_" + glue::as_glue(v) + ".csv"), n_sim = n_mc) %>% 
#       simulate_ql(age0 = 50, age1 = 95, age_until = 5, dt = 0.01) %>% 
#       summarise(across(starts_with("QL"), list(
#         M = \(x) median(x),
#         L = \(x) quantile(x, 0.025),
#         U = \(x) quantile(x, 0.975)
#       )), .by = Age) %>% 
#       mutate(Var = v)
#     
#     sims
#   }))
#   
#   write_csv(res, here::here("docs", "full", "sim_qloss_" + gp + ".csv"))
# 
#   
# }
# 
# 
# 
# g_ql <- res %>% 
#   ggplot(aes(x = Age)) + 
#   geom_ribbon(aes(ymin = QL35_L, ymax = QL35_U, group = Var), alpha = 0.2) +
#   geom_line(aes(y = QL35_M, group = Var))
# 


