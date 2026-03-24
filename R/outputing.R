
output_pars_shortfall <- function(pars_sf, data_norm, vset = "uk") {
  require(jsonlite)
  
  pars_demo <- data_norm %>% 
    filter(Year == 2023) %>% 
    ungroup() %>% 
    select(Age, r_death = mortality, norm = norm_leoss, N = N) 
  
  js <- list(
    ValueSet = vset,
    pars_demo = pars_demo,
    pars_shortfall = pars_sf
  )
  
  f <- here::here("pars", paste0("pars_ql_", vset, ".csv"))
  write_csv(pars_sf, f)
  
  f <- here::here("pars", paste0("pars_qol_", vset, ".json"))
  jsonlite::write_json(js, f, auto_unbox = T)

  return(f)
}
