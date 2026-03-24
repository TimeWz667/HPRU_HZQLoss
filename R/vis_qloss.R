

vis_shortfall <- function(sim, stats, vset, age0 = 50, age1 = 98) {
  gs <- list()
  gs$g_ql <- stats %>% 
    filter(startsWith(Index, "QL")) %>% 
    filter(endsWith(Index, "15")) %>% 
    filter(Age <= age1) %>% 
    filter(Age >= age0) %>% 
    mutate(Index = factor(Index, c("QLH15", "QL15"))) %>% 
    ggplot(aes(x = Age)) +
    geom_ribbon(aes(ymin = L, ymax = U), alpha = 0.2) +
    geom_line(aes(y = M)) +
    scale_y_continuous("QALY loss") +
    scale_x_continuous("Age at onset of HZ rash") +
    facet_wrap(.~Index, labeller = labeller(Index = c(QL15 = "From population norm", QLH15 = "From perfect health"))) +
    expand_limits(y = 0) +
    labs(caption = "discounting: 1.5%") +
    theme_bw()
  
  
  gs$g_ql_grad <- sim %>% 
    select(Age, QL15, QLH15) %>% 
    filter(Age <= age1) %>%  
    filter(Age >= age0) %>%  
    pivot_longer(-Age, names_to = "Index") %>% 
    mutate(Index = factor(Index, c("QLH15", "QL15"))) %>% 
    ggplot(aes(x = Age, y = value)) +
    stat_lineribbon() +
    scale_fill_brewer() +
    scale_y_continuous("QALY loss") +
    scale_x_continuous("Age at onset of HZ rash") +
    facet_wrap(.~Index, labeller = labeller(Index = c(QL15 = "From population norm", QLH15 = "From perfect health"))) +
    expand_limits(y = 0) +
    labs(caption = "discounting: 1.5%") +
    theme_bw()
  
  
  ggsave(gs$g_ql, filename = here::here("docs", "figs", paste0("g_ql_", vset, ".png")), width = 7, height = 4)
  ggsave(gs$g_ql_grad, filename = here::here("docs", "figs", paste0("g_ql_grad_", vset, ".png")), width = 7.3, height = 4)
  
  return(gs)
  
}
