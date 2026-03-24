
vis_qol_t <- function(pars, pars_demo, vset, age = 50) {
  theme_set(theme_bw())
  
  norm0 <- pars_demo %>% filter(Year == 2023 & Age == age) %>% pull(norm)
  
  sim_qol <- pars %>% 
    filter(Key < 100) %>% 
    crossing(ti = seq(0, 1, 0.01)) %>% 
    mutate(
      Age = age, 
      rate = r0 * exp(Age * ba1),
      pr_hz = 1 - pexp(ti, rate),
      Cluster = runif(n()),
      Cluster = case_when(
        Cluster < Prop_0 ~ "0",
        Cluster < Prop_0 + Prop_1 ~ "1",
        T ~ "2",
      ),
      Q = case_when(
        Cluster == "0" ~ Q_0,
        Cluster == "1" ~ Q_1,
        T ~ Q_2
      ),
      Q_norm = pmin(Q, norm0),
      qol = (Q_0 * Prop_0 + Q_1 * Prop_1 + Q_2 * Prop_2),
      qol_norm = (pmin(Q_0, norm0) * Prop_0 + pmin(Q_1, norm0) * Prop_1 + pmin(Q_2, norm0) * Prop_2),
      Cluster = factor(Cluster)
    )
  
  gs <- list()
  
  gs$g_qol <- sim_qol %>% 
    ggplot() +
    geom_point(aes(x = ti, y = Q, colour = Cluster, alpha = pr_hz), pch = 19) +
    geom_hline(yintercept = norm0, linetype = 2) +
    geom_text(x = 1, y = norm0, vjust = -.5, hjust = 1.1, label = "Population norm") +
    scale_alpha("Probability of \ndisutility", range = c(0, 0.5), label = scales::percent) +
    scale_y_continuous("Health-related quality of life", breaks = seq(-0.5, 1, 0.5)) +
    scale_x_continuous("Months since rash onset", label = scales::number_format(scale = 12)) +
    expand_limits(y = c(-0.5, 1)) +
    labs(caption = paste0("Age: ", age))
  
  
  gs$g_qol_a <- sim_qol %>% 
    mutate(qol_t = qol * pr_hz + (1 - pr_hz)) %>% 
    ggplot() +
    geom_point(aes(x = ti, y = qol_t, colour = Cluster), alpha = 0.1, pch = 19) +
    geom_hline(yintercept = norm0, linetype = 2) +
    geom_text(x = 1, y = norm0, vjust = -.5, hjust = 1.1, label = "Population norm") +
    scale_y_continuous("Health-related quality of life", breaks = seq(-0.5, 1, 0.5)) +
    scale_x_continuous("Months since rash onset", label = scales::number_format(scale = 12)) +
    expand_limits(y = c(-0.5, 1)) +
    labs(caption = paste0("Age: ", age))
  
  
  gs$g_qol_t <- sim_qol %>% 
    mutate(qol_t = qol * pr_hz + (1 - pr_hz)) %>% 
    ggplot() +
    stat_lineribbon(aes(x = ti, y = qol_t)) +
    geom_hline(yintercept = norm0, linetype = 2) +
    geom_text(x = 1, y = norm0, vjust = -.5, hjust = 1.1, label = "Population norm") +
    #geom_point(aes(x = ti, y = qol_t), pch = 19) +
    scale_y_continuous("Health-related quality of life", breaks = seq(-0.5, 1, 0.5)) +
    scale_fill_brewer() +
    scale_x_continuous("Months since rash onset", label = scales::number_format(scale = 12)) +
    expand_limits(y = c(-0.5, 1)) +
    labs(caption = paste0("Age: ", age))
  
  
  gs$g_qol_norm <- sim_qol %>% 
    mutate(pr_hz = ifelse(Q <= norm0, pr_hz, 0)) %>% 
    ggplot() +
    geom_point(aes(x = ti, y = Q_norm, colour = Cluster, alpha = pr_hz), pch = 19) +
    geom_hline(yintercept = norm0, linetype = 2) +
    geom_text(x = 1, y = norm0, vjust = -.5, hjust = 1.1, label = "Population norm") +
    scale_alpha("Probability of \ndisutility", range = c(0, 0.5), label = scales::percent) +
    scale_y_continuous("Health-related quality of life", breaks = seq(-0.5, 1, 0.5)) +
    scale_x_continuous("Months since rash onset", label = scales::number_format(scale = 12)) +
    expand_limits(y = c(-0.5, 1)) +
    labs(caption = paste0("Age: ", age))
  
  
  gs$g_qol_a_norm <- sim_qol %>% 
    mutate(qol_t = qol_norm * pr_hz + norm0 * (1 - pr_hz)) %>% 
    ggplot() +
    geom_point(aes(x = ti, y = qol_t, colour = Cluster), alpha = 0.1, pch = 19) +
    geom_hline(yintercept = norm0, linetype = 2) +
    geom_text(x = 1, y = norm0, vjust = -.5, hjust = 1.1, label = "Population norm") +
    scale_y_continuous("Health-related quality of life", breaks = seq(-0.5, 1, 0.5)) +
    scale_x_continuous("Months since rash onset", label = scales::number_format(scale = 12)) +
    expand_limits(y = c(-0.5, 1)) +
    labs(caption = paste0("Age: ", age))
  
  
  gs$g_qol_t_norm <- sim_qol %>% 
    mutate(qol_t = qol_norm * pr_hz + norm0 * (1 - pr_hz)) %>% 
    ggplot() +
    stat_lineribbon(aes(x = ti, y = qol_t)) +
    geom_hline(yintercept = norm0, linetype = 2) +
    geom_text(x = 1, y = norm0, vjust = -.5, hjust = 1.1, label = "Population norm") +
    #geom_point(aes(x = ti, y = qol_t), pch = 19) +
    scale_y_continuous("Health-related quality of life", breaks = seq(-0.5, 1, 0.5)) +
    scale_fill_brewer() +
    scale_x_continuous("Months since rash onset", label = scales::number_format(scale = 12)) +
    expand_limits(y = c(-0.5, 1)) +
    labs(caption = paste0("Age: ", age))
  
  ggsave(gs$g_qol, filename = here::here("docs", "figs", paste0("g_qol_", age, "_", vset, ".png")), width = 8, height = 5)
  ggsave(gs$g_qol_a, filename = here::here("docs", "figs", paste0("g_qol_avg_", age, "_", vset, ".png")), width = 8, height = 5)
  ggsave(gs$g_qol_t, filename = here::here("docs", "figs", paste0("g_qol(", age, ")_", vset, ".png")), width = 8, height = 5)
  ggsave(gs$g_qol_norm, filename = here::here("docs", "figs", paste0("g_qol_", age, "_norm_", vset, ".png")), width = 8, height = 5)
  ggsave(gs$g_qol_a_norm, filename = here::here("docs", "figs", paste0("g_qol_avg_", age, "_norm_", vset, ".png")), width = 8, height = 5)
  ggsave(gs$g_qol_t_norm, filename = here::here("docs", "figs", paste0("g_qol(", age, ")_norm_", vset, ".png")), width = 8, height = 5)
  
  return(gs)
}
