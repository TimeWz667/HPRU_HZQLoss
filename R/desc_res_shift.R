
desc_disuti <- function(pars_tte, pars_tte0) {
  
  pars0 <- pars_tte0
  pars1 <- pars_tte 
  
  res <- pars0$Ext %>%
    mutate(Key = 1:n()) %>% 
    crossing(Age = c(50, 60, 70, 80, 90)) %>% 
    mutate(
      rate0 = r0 * exp(Age * ba1),
      dur0 = 365.25 / rate0
    ) %>% 
    select(Key, Age, rate0, dur0) %>% 
    arrange(Key) %>% 
    left_join(pars1$Ext %>% mutate(Key = 1:n())) %>% 
    mutate(
      rate1 = r0 * exp(Age * ba1),
      dur1 = 365.25 / rate1,
      diff = 1 - dur0 / dur1
    ) %>% 
    group_by(Age) %>% 
    summarise(across(c(dur0, dur1, diff), amlu))
  
}
