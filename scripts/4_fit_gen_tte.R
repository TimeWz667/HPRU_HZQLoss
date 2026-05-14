library(tidyverse)
library(rstan)

theme_set(theme_bw())

## Loading stan models and settings
options(mc.cores = 8)
rstan_options(auto_write = TRUE)


for (file in dir("R")) source(here::here("R", file))

model_t2z <- rstan::stan_model(here::here("models", "time2zero_age.stan"))


## Functions

fit_tte <- function(data_tte, n_iter, n_warmup) {
  dat <- list(
    N = nrow(data_tte),
    As = data_tte$Age,
    Ts = data_tte$TTe
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


root_src <- here::here("out", "gen_2gp")


for(v in paste0("s", 1:10)) {
  cat(v, "\n")
  
  data_tte <- read_csv(here::here(root_src, "generated_" + glue::as_glue(v) + ".csv")) %>% 
    filter(visit == 1) %>% 
    select(Age, TTe = disuti)
  
  
  dat <- list(
    N = nrow(data_tte),
    As = data_tte$Age,
    Ts = data_tte$TTe
  )

  
  tryCatch({
    post <- rstan::sampling(model_t2z, data = dat, iter = n_iter, warmup = n_warmup)
    
    ext <- rstan::extract(post, pars = c("ba1", "r0")) %>% data.frame()
    
    write_csv(ext, here::here(root_src, "post_tte_exp_" + glue::as_glue(v) + ".csv"))

    #when it throws an error, the following block catches the error
  }, error = function(msg){
    return(NULL)
  })
  
}





list_studies <- read_csv(here::here("data", "list_studies.csv"))
list_studies

v <- "s1"
root_ind <- here::here("out", "gen_ind")
root_2gp <- here::here("out", "gen_2gp")


data_tte_ind <- read_csv(here::here(root_ind, "generated_" + glue::as_glue(v) + ".csv")) %>% 
  filter(visit == 1) %>% 
  select(SID, Age, TTe = disuti)

data_tte_2gp <- read_csv(here::here(root_2gp, "generated_" + glue::as_glue(v) + ".csv")) %>% 
  filter(visit == 1) %>% 
  select(SID, Age, TTe = disuti)


ext_ind <- fit_tte(data_tte_ind, n_iter = n_iter, n_warmup = n_warmup) %>% 
  mutate(Scenario = "All from independent pars")

ext_2gp <- fit_tte(data_tte_2gp, n_iter = n_iter, n_warmup = n_warmup) %>% 
  mutate(Scenario = "All")



ext_minus_bind <- lapply(list_studies$SID, \(sid) {
  print(sid)
  
  data_tte_bind <- bind_rows(
    data_tte_2gp %>% filter(SID != sid),
    data_tte_ind %>% filter(SID == sid)
  )
  
  ext_bind <- fit_tte(data_tte_bind, n_iter = n_iter, n_warmup = n_warmup) %>% 
    mutate(Scenario = paste0("All - ", sid))
}) %>% bind_rows()


ext_plus_bind <- lapply(list_studies$SID, \(sid) {
  print(sid)
  
  data_tte_bind <- bind_rows(
    data_tte_2gp %>% filter(SID == sid),
    data_tte_ind %>% filter(SID != sid)
  )
  
  ext_bind <- fit_tte(data_tte_bind, n_iter = n_iter, n_warmup = n_warmup) %>% 
    mutate(Scenario = paste0("Ind + ", sid))
}) %>% bind_rows()





ext <- bind_rows(
  ext_ind, ext_2gp,
  ext_minus_bind,
  ext_plus_bind
)


root_src <- here::here("out", "gen_bind")


write_csv(ext, here::here(root_src, "post_tte_exp_" + glue::as_glue(v) + ".csv"))



disuti <- ext %>% 
  crossing(Age = seq(50, 80, 5)) %>% 
  mutate(
    rate = r0 * exp(Age * ba1)
  ) %>% 
  summarise(
    M = median(1 / rate),
    L = quantile(1 / rate, 0.025),
    U = quantile(1 / rate, 0.975),
    .by = c("Age", "Scenario")
  )


disuti$Scenario %>% unique()

disuti_minus <- disuti %>% 
  filter(startsWith(Scenario, "All - ")) %>% 
  tidyr::extract(Scenario, "SID", "\\- (\\w+)") %>% 
  left_join(
    disuti %>% filter(Scenario == "All") %>% 
      select(Age, S1 = M)
  ) %>% 
  left_join(
    disuti %>% filter(Scenario == "All from independent pars") %>% 
      select(Age, S0 = M)
  ) 


disuti_plus <- disuti %>% 
  filter(startsWith(Scenario, "Ind")) %>% 
  tidyr::extract(Scenario, "SID", "\\+ (\\w+)") %>% 
  left_join(
    disuti %>% filter(Scenario == "All") %>% 
      select(Age, S1 = M)
  ) %>% 
  left_join(
    disuti %>% filter(Scenario == "All from independent pars") %>% 
      select(Age, S0 = M)
  ) 


g <- disuti_plus %>% 
  ggplot() +
  geom_segment(aes(x = S0, xend = M, y = SID, colour = "plus")) + 
  geom_pointrange(aes(x = M, xmin = L, xmax = U, y = SID, colour = "plus"), linewidth = 1.4) + 
  geom_segment(data = disuti_minus, aes(x = S1, xend = M, y = SID, colour = "minus")) + 
  geom_pointrange(data = disuti_minus, aes(x = M, xmin = L, xmax = U, y = SID, colour = "minus"), linewidth = 1.4) + 
  geom_vline(aes(xintercept = S0)) +
  geom_vline(aes(xintercept = S1)) +
  facet_grid(Age~.) +
  scale_x_continuous("Disutility period, month", labels = scales::number_format(scale = 12)) +
  scale_color_discrete("Baseline", labels = c(minus = "Two groups minus", plus = "Ind plus")) + 
  expand_limits(x = 0)



ggsave(g, filename = here::here("docs", "experiments", "g_+-.png"), width = 12, height = 20) 



sid_ord <- disuti_plus %>% 
  filter(Age == 60) %>% 
  mutate(Diff = abs(S0 - M)) %>% 
  arrange(Diff) %>% 
  pull(SID)


sid_ord



ext_stepwise_bind <- lapply(1:length(sid_ord[1:4]), \(k) {
  sid_sel <- sid_ord[1:k]

  print(sid_sel)
  
  data_tte_bind <- bind_rows(
    data_tte_2gp %>% filter(SID %in% sid_sel),
    data_tte_ind %>% filter(!SID %in% sid_sel)
  )
  
  ext_bind <- fit_tte(data_tte_bind, n_iter = n_iter, n_warmup = n_warmup) %>% 
    mutate(Scenario = paste0(" + ", sid_ord[k]))
}) %>% bind_rows()





disuti_step <- ext_stepwise_bind %>% 
  crossing(Age = seq(50, 80, 5)) %>% 
  mutate(
    rate = r0 * exp(Age * ba1)
  ) %>% 
  summarise(
    M = median(1 / rate),
    L = quantile(1 / rate, 0.025),
    U = quantile(1 / rate, 0.975),
    .by = c("Age", "Scenario")
  ) %>% 
  left_join(
    disuti %>% filter(Scenario == "All") %>% 
      select(Age, S1 = M)
  ) %>% 
  left_join(
    disuti %>% filter(Scenario == "All from independent pars") %>% 
      select(Age, S0 = M)
  ) 


disuti_step %>% 
  ggplot() +
  geom_segment(aes(x = S0, xend = M, y = Scenario, colour = "plus")) + 
  geom_point(aes(x = M, y = Scenario, colour = "plus"), linewidth = 1.4) + 
  geom_vline(aes(xintercept = S0)) +
  geom_vline(aes(xintercept = S1)) +
  facet_grid(Age~.) +
  scale_x_continuous("Disutility period, month", labels = scales::number_format(scale = 12)) +
  expand_limits(x = 0)





data_tte_ind %>% 
  group_by(SID, Age) %>% 
  summarise(M = median(TTe)) %>% 
  ggplot() +
  geom_smooth(aes(x = Age, y = M, colour = SID))



pars_tte %>%
  filter(SID != "RamCRI") %>% 
  select(- sd) %>% 
  pivot_wider(names_from = Var, values_from = mean) %>% 
  ggplot() +
    geom_point(aes(x = r0, y = ba1, colour = SID)) +
  geom_text(aes(x = r0, y = ba1, label = SID))


