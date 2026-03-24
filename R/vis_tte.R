visualise_tte <- function(dat_tte, pars_tte) {
  require(tidyverse)
  require(tidybayes)
  
  theme_set(theme_bw())
  
  gs <- list()
  
  dat_tte <- dat_tte %>% filter(!is.na(Age))
  
  rates <- pars_tte$Ext %>% 
    crossing(Age = 10:100) %>% 
    mutate(rate = r0 * exp(Age * ba1)) %>% 
    group_by(Age) %>% 
    summarise(rate = mean(rate), duration = mean(1/ rate)) %>% 
    ungroup() 

  imp <- dat_tte %>% 
    left_join(rates) %>% 
    mutate(
      rd = runif(n()),
      rd = ifelse(Recovered,  rd * (1 - exp(- rate * (T_end - T_last))), rd),
      rd = - log(1 - rd) / rate,
      T_imp = ifelse(Recovered, rd + T_last, rd + T_obs)
    ) %>% 
    filter(T_evt > 0) %>% 
    arrange(- T_imp) %>% 
    mutate(
      Key = 1:n() / n()
    ) %>% 
    arrange(- T_imp)
  
  imp_sel <- imp %>% 
    filter(1:nrow(imp) %in% round(seq(1, nrow(imp), length = 50)))
  
  gs$g_imp <- imp %>% 
    ggplot() +
    geom_line(aes(x = T_imp, y = Key)) +
    geom_linerange(data = imp_sel %>% filter(Recovered), aes(xmin = T_last, xmax = T_end, y = Key, colour = "Interval")) +
    geom_point(data = imp_sel %>% filter(!Recovered), aes(x = T_obs, y = Key, colour = "Last point")) +
    scale_color_discrete("Type of data") +
    scale_y_continuous("Probability of unrecovered, %", labels = scales::percent) +
    scale_x_continuous("Time since the onset of rash", breaks = seq(0, 2.5, 0.25), 
                       labels = scales::number_format(scale = 12, suffix = " mo."),
                       limits = c(0, 1.5))
  
  
  pars <- pars_tte$Ext %>% mutate(Key = 1:n())
  
  sims <- dat_tte %>% 
    mutate(
      Key = sample(pars$Key, n(), rep = T)
    ) %>% 
    left_join(pars) %>% 
    mutate(
      rate = (r0 * exp(Age * ba1)),
      T_Sim = rexp(n(), rate)
    )
  
  exp_sims <- pars %>%
    mutate(
      rd = runif(n()),
      rd = -log(1 - runif(n()))
    ) %>% 
    crossing(Age = min(sims$Age):100) %>% 
    mutate(
      rate = (r0 * exp(Age * ba1)),
      T_Sim = rd / rate,
      T_Exp = 1 / rate
    )
  
  exp_cri <- exp_sims %>%
    group_by(Age) %>%
    summarise(
      m = mean(T_Exp),
      l = quantile(T_Exp, 0.025),
      u = quantile(T_Exp, 0.975)
    )

  
  exp_ci <- exp_cri %>% 
    left_join(sims %>% group_by(Age) %>% summarise(n = n())) %>% 
    filter(n >= 10) %>% 
    mutate(
      l = 1 / (1 / m * (1 - 1.96 / sqrt(n))),
      u = 1 / (1 / m * (1 + 1.96 / sqrt(n)))
    ) 
  
  exp_range <- bind_rows(
    exp_cri %>% mutate(Source = "Simulated"), 
    exp_ci %>% mutate(Source = "Data")
  )

  gs$g_tte_data <- sims %>% 
    group_by(Age) %>% 
    summarise(
      Data = mean(T_evt2),
      Simulated = mean(T_Sim)
    ) %>% 
    pivot_longer(-Age, names_to = "Source") %>% 
    ggplot() +
    geom_point(aes(x = Age, y = value, colour = Source)) +
    geom_ribbon(data = exp_ci, aes(x = Age, ymin = l, ymax = u), alpha = 0.2) +
    geom_line(data = exp_ci, aes(x = Age, y = m)) + 
    scale_y_continuous("Time to recovery, months", breaks = seq(0, 2.5, 0.5), 
                       labels = scales::number_format(scale = 12, suffix = " mo."),
                       limits = c(0, 1.5)) + 
    scale_x_continuous("Age", limits = c(0, 100))+ 
    scale_colour_brewer(palette = "Accent")
  
  
  gs$g_tte_pred <- exp_cri %>% 
    ggplot() +
    geom_ribbon(aes(x = Age, ymin = l, ymax = u), alpha = 0.2) +
    geom_line(aes(x = Age, y = m)) +
    scale_y_continuous("Time to recovery, months", breaks = seq(0, 2.5, 0.5), 
                       labels = scales::number_format(scale = 12, suffix = " mo."),
                       limits = c(0, 1.5)) +
    scale_x_continuous("Age", limits = c(0, 100))
  
  
  gs$g_tte_bind1 <- ggpubr::ggarrange(
    gs$g_imp + labs(subtitle = "(A)") + theme(legend.position = c(1, 1), legend.justification = c(1.1, 1.1)),
    gs$g_tte_data + labs(subtitle = "(B)") + theme(legend.position = c(1, 1), legend.justification = c(1.1, 1.1)),
    gs$g_tte_pred + labs(subtitle = "(C)"),
    ncol = 1
  )
  
  gs$g_tte_bind2 <- ggpubr::ggarrange(
    gs$g_tte_data + labs(subtitle = "(A)") + theme(legend.position = c(1, 1), legend.justification = c(1.1, 1.1)),
    gs$g_tte_pred + labs(subtitle = "(B)"),
    ncol = 1
  )
  

  return(gs)
}


output_vis_tte <- function(gs, folder, modifier, ext = ".png") {
  modifier <- glue::as_glue(modifier)
  ext <- glue::as_glue(ext)
  
  ggsave(gs$g_imp, file = here::here(folder, "g_tte_imputed_" + modifier + ext), width = 6, height = 4)
  ggsave(gs$g_tte_data, file = here::here(folder, "g_tte_data_" + modifier + ext), width = 6, height = 4)
  ggsave(gs$g_tte_pred, filename = here::here(folder, "g_tte_cri_" + modifier + ext), width = 7, height = 4)
  ggsave(gs$g_tte_bind1, file = here::here(folder, "g_tte_bind1_" + modifier + ext), width = 9, height = 12)
  ggsave(gs$g_tte_bind2, file = here::here(folder, "g_tte_bind2_" + modifier + ext), width = 9, height = 8)
  return(here::here(folder, "g_tte_bind" + modifier + ext))
}


vis_tte_comparing <- function(tabs_tte_ph, tabs_tte_pn, tabs_tte_0, folder, ext = ".png", scale = 1) {
  sim <- bind_rows(
    tibble(Age = seq(50, 90, 5)) %>% cross_join(tabs_tte_ph$Ext) %>% mutate(Cap = "Perfect health"),
    tibble(Age = seq(50, 90, 5)) %>% cross_join(tabs_tte_pn$Ext) %>% mutate(Cap = "Population norm"),
    tibble(Age = seq(50, 90, 5)) %>% cross_join(tabs_tte_0$Ext) %>% mutate(Cap = "Pre-HZ baseline")
  ) %>% 
    mutate(
      rate = (r0 * exp(Age * ba1)),
      T_Exp = 1 / rate
    ) %>% 
    group_by(Age, Cap) %>%
    summarise(
      m = mean(T_Exp) * scale,
      l = quantile(T_Exp, 0.025) * scale,
      u = quantile(T_Exp, 0.975) * scale
    )
  
  g <- sim %>% 
    ggplot() +
    geom_ribbon(aes(x = Age, ymin = l, ymax = u, fill = Cap), alpha = 0.2) +
    geom_line(aes(x = Age, y = m, colour = Cap)) +
    scale_y_continuous("Time to recovery, months", breaks = seq(0, 2.5, 0.25), 
                       labels = scales::number_format(scale = 12, suffix = " mo."),
                       limits = c(0, 1)) +
    scale_x_continuous("Age", limits = c(50, 90))
  
  
  ggsave(g, filename = here::here(folder, "g_tte_compare.png"), height = 5, width = 7)
  
  return(sim)
}



