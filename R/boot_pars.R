
boot_pars <- function(file_pars_tte, file_pars_qol, n_sim = 1000) {
  boot_qol <- read_csv(file_pars_qol) %>% 
    mutate(
      Parameter = recode(Parameter, 
                         `(Intercept)` = "b0", 
                         `poly(Age, 2, raw = T)1` = "ba1", 
                         `poly(Age, 2, raw = T)2` = "ba2",
                         d15TRUE = "bd15",
                         d30TRUE = "bd30"),
      Model = recode(Model, `Pr(C0)` = "PZ", `Pr(C1|~C0)` = "PC1")
    ) %>% 
    crossing(Key = 1:n_sim) %>% 
    mutate(
      value = rnorm(n(), Estimate, Std),
      value = case_when(
        Parameter %in% c("ba1", "ba2") ~ Estimate,
        T ~ value
      ),
      name = paste0(Model, "_", Parameter)
    ) %>% 
    select(Key, value, name) %>% 
    pivot_wider()
   
  # boot_qol <- read_csv(file_pars_qol) %>%
  #   filter(Agp == "All") %>% 
  #   crossing(Key = 1:n_sim) %>% 
  #   mutate(
  #     std = ifelse(is.na(std), 0, std),
  #     Q = rnorm(n(), mu, std),
  #     Q = pmin(Q, 1)
  #   ) %>% 
  #   relocate(Key) %>% 
  #   select(Key, Q, Prop, Cluster) %>% 
  #   pivot_wider(names_from = Cluster, values_from = c(Q, Prop))
  
  boot_tte <- read_csv(file_pars_tte) %>% select(r0, ba1) %>% mutate(ID = 1:n()) 
  
  boot_tte <- tibble(Key = 1:n_sim, ID = sample.int(nrow(boot_tte), n_sim, replace = n_sim > nrow(boot_tte))) %>% 
    left_join(boot_tte, by = "ID") %>% 
    select(- ID)
  
  pars <- merge(boot_qol, boot_tte)
  return(pars)
}


boot_pars_bayes <- function(file_pars_tte, file_pars_qol, n_sim = 1000) {
  boot_qol <- read_csv(file_pars_qol) %>% 
    select(- lp__)
  
  boot_qol <- boot_qol  %>% 
    group_by(Model) %>% 
    sample_n(n()) %>% 
    mutate(Key = 1:n()) %>% 
    pivot_longer(-c(Key, Model)) %>% 
    mutate(name = paste0(Model, "_", name)) %>% 
    filter(!is.na(value)) %>% 
    ungroup() %>% 
    select(-Model) %>% 
    pivot_wider()
  
  boot_tte <- read_csv(file_pars_tte) %>% select(r0, ba1) %>% mutate(ID = 1:n()) 
  
  boot_tte <- tibble(Key = 1:n_sim, ID = sample.int(nrow(boot_tte), n_sim, replace = n_sim > nrow(boot_tte))) %>% 
    left_join(boot_tte, by = "ID") %>% 
    select(- ID)
  
  pars <- boot_tte %>% left_join(boot_qol)
  return(pars)
}


