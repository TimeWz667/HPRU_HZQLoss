format_qol_gen <- function(data_raw) {
  max_qol <- max(data_raw$EQ5D)
  min_qol <- min(data_raw$EQ5D)
  
  dat_qol <- data_raw %>% 
    filter(ti < disuti) %>% 
    select(SID, PID, Age, Agp, ti, EQ5D) %>% 
    mutate(
      Q_rescaled = (EQ5D - min_qol) / (max_qol - min_qol)
    )
  return(dat_qol)
}


fit_qol <- function(data_raw, n_iter, n_warmup, n_attempts = 10) {
  
  data_qol <- format_qol_gen(data_raw)
  i <- 0
  while (i < n_attempts) {
    ext <- tryCatch({
      pars_qol <- fit_qol_bayes(data_qol, ind = F, n_iter = n_iter, n_collect = n_iter - n_warmup)
      return(pars_qol$Ext)
    }, error = function(msg){
      print(msg)
      cat("re-try\n")
      return(NULL)
    })
    
    if (!is.null(ext)) {
      break
    }
    i <- i + 1
  }
  cat("failed\n")
  return(ext)
  
}



fit_tte_exact <- function(data_raw, n_iter, n_warmup, n_attempts = 10) {
  model_t2z <- rstan::stan_model(here::here("models", "time2zero_age.stan"))
  
  data_tte <- data_raw %>% 
    filter(visit == 1) %>% 
    select(Age, TTE = disuti)
  
  dat <- list(
    N = nrow(data_tte),
    As = data_tte$Age,
    Ts = data_tte$TTE
  )
  
  i <- 0
  while (i < n_attempts) {
    ext <- tryCatch({
      post <- rstan::sampling(model_t2z, data = dat, iter = n_iter, warmup = n_warmup)
      ext <- rstan::extract(post, pars = c("ba1", "r0")) %>% data.frame()
      return(ext)
    }, error = function(msg){
      print(msg)
      cat(", re-try\n")
      return(NULL)
    })
    
    if (!is.null(ext)) {
      break
    }
    i <- i + 1
  }
  cat("failed\n")
  return(ext)
}



