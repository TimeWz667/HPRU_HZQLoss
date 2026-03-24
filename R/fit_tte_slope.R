fit_tte_slope <- function(pars_tte, suf = "") {
  pp <- pars_tte[[1]]
  
  fit <- lm(ba1 ~ log(r0), data = pp)
  
  tab <- bind_cols(
    M = coef(fit),
    confint(fit)
  ) %>% 
    mutate(Pars = names(coef(fit))) %>% 
    relocate(Pars)
  
  f <-  here::here("docs", "tabs", paste0("slope_sub_tte_", suf, ".csv"))
  write_csv(tab, f)
  f
}