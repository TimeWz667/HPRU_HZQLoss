fit_tte <- function(model, dat_tte, n_iter = 2e4, n_collect = 500, n_chains = 4) {
  n_warmup <- floor(n_iter - n_collect)
  

  ds <- local({
    #dat_tte <- tar_read(data_tte) #%>% head(200)
    dat_tte <- dat_tte %>% filter(T_evt > 0) %>% filter(T_last < T_end)
    
    dat_int <- dat_tte %>% filter(Recovered)
    dat_cen <- dat_tte %>% filter(!Recovered) %>% filter(!is.na(Age))
    
    list(
      N = nrow(dat_int),
      Ts0 = dat_int$T_last,
      Ts1 = dat_int$T_end,
      As = dat_int$Age,
      N_Cen = nrow(dat_cen),
      Ts_Cen = dat_cen$T_end,
      As_Cen = dat_cen$Age
    )
  })
  
  post <- sampling(model, data = ds, pars = c("r0", "ba1"), init = \() { list(r0 = exp(-2.35), ba1 = 0.01)},
                   chains = n_chains, iter = n_iter, warmup = n_warmup)
  

  return(restructure_stan(post))
  
}

