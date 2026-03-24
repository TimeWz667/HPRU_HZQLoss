
get_pars_baseline <- function(data_baseline, n_sim = 1e3) {
  ss <- glm((EQ5D == 1) ~ poly(age, 2), data = data_baseline %>% filter(!is.na(age)), family = binomial)

  p_base0 <- tibble(
    age = 18:100,
    p_ph = predict(ss, newdata = data.frame(age = 18:100), type = "response")
  )

  dens <- density(data_baseline %>% filter(!is.na(age)) %>% filter(EQ5D < 1) %>% pull(EQ5D))
  
  fn <- approxfun(cumsum(dens$y) / sum(dens$y), dens$x)

  
  q <- tibble(EQ5D0 = pmin(fn(runif(n_sim)), 1)) %>% 
    mutate(Key = 1:n())
  
  return(list(
    p_base0 = p_base0,
    q = q
  ))
}





