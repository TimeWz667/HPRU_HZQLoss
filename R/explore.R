


explore_qol <- function(data_qol) {
  
  require(tidyverse)
  require(tidybayes)
  
  theme_set(theme_bw())
  
  
  data_qol <- tar_read(data_qol, 2)
  
  
  res <- list()
  
  res$stats_range <- data_qol %>% 
    group_by(SID) %>% 
    summarise(
      Age0 = min(Age),
      Age1 = max(Age),
      EQ5D0 = min(EQ5D),
      EQ5D1 = max(EQ5D)
    )

  
  res$g_hist_all <- data_qol %>% 
    ggplot(aes(x = EQ5D)) +
    #geom_histogram() +
    geom_density()
  
  
  res$g_hist_sid <- data_qol %>% 
    ggplot(aes(x = EQ5D)) +
    #geom_histogram() +
    geom_density() +
    facet_grid(.~SID)
  
  ## labelling
  dat_cluster <- data_qol %>% 
    left_join(
      data_qol %>% 
        filter(EQ5D != 1) %>% 
        mutate(
          Cluster = kmeans(EQ5D, 2)$cluster
        )
    ) %>% 
    mutate(
      Cluster = ifelse(is.na(Cluster), 0, Cluster),
      Cluster = as.character(Cluster)
    )

  
  res$stats_cluster_all <- dat_cluster %>% 
    group_by(Cluster) %>% 
    summarise(
      mu = mean(EQ5D),
      Q1 = quantile(EQ5D, 0.25),
      Q3 = quantile(EQ5D, 0.75),
      std = sd(EQ5D),
      N = n()
    ) %>% 
    ungroup() %>% 
    mutate(Prop = N / sum(N))
    
  
  res$stats_cluster_sid <- dat_cluster %>% 
    group_by(Cluster, SID) %>% 
    summarise(
      mu = mean(EQ5D),
      Q1 = quantile(EQ5D, 0.25),
      Q3 = quantile(EQ5D, 0.75),
      std = sd(EQ5D),
      N = n()
    ) %>% 
    group_by(SID) %>% 
    mutate(Prop = N / sum(N)) %>% 
    arrange(SID)
  
  
  res$stats_cluster_agp <- dat_cluster %>% 
    group_by(Cluster, Agp) %>% 
    summarise(
      mu = mean(EQ5D),
      Q1 = quantile(EQ5D, 0.25),
      Q3 = quantile(EQ5D, 0.75),
      std = sd(EQ5D),
      N = n()
    ) %>% 
    group_by(Agp) %>% 
    mutate(Prop = N / sum(N)) %>% 
    arrange(Agp)
  
  
  res$stats_cluster_tip <- dat_cluster %>% 
    mutate(
      tip = case_when(
        ti <= 0.25 ~ "0-3 mo",
        ti <= 0.5 ~ "4-6 mo",
        ti <= 1 ~ "1 year",
        T ~ ">= 1 year"
      ),
      tip = factor(tip, c("0-3 mo", "4-6 mo", "1 year", ">= 1 year"))
    )%>% 
    group_by(Cluster, tip) %>% 
    summarise(
      mu = mean(EQ5D),
      Q1 = quantile(EQ5D, 0.25),
      Q3 = quantile(EQ5D, 0.75),
      std = sd(EQ5D),
      N = n()
    ) %>% 
    group_by(tip) %>% 
    mutate(Prop = N / sum(N))

  
  res$stats_cluster_sid %>% 
    left_join(res$stats_cluster_all %>% select(Cluster, Prop0 = Prop)) %>% 
    ungroup() %>% 
    mutate(
      w = sqrt(N) / (Prop * (1 - Prop)),
      q = w * (Prop - Prop0) ^ 2
    ) %>% 
    summarise(q = sum(q), q = (q - n() + 1) / q) %>% pull(q)
    
  
  res$g_range_sid <- res$stats_cluster_sid %>% 
    #filter(Cluster != "0") %>% 
    ggplot() +
    geom_pointinterval(aes(x = mu, xmin = Q1, xmax = Q3, y = SID, colour = Cluster))
  
  
  res$g_range_sid <- res$stats_cluster_sid %>% 
    #filter(Cluster != "0") %>% 
    ggplot() +
    geom_pointinterval(aes(x = mu, xmin = Q1, xmax = Q3, y = SID, colour = Cluster))
  
  res$stats_cluster_sid %>% 
    #filter(Cluster != "0") %>% 
    ggplot() +
    geom_bar(aes(x = Prop, y = SID, fill = Cluster), stat = "identity")
  
  res$stats_cluster_agp %>% 
    #filter(Cluster != "0") %>% 
    ggplot() +
    geom_bar(aes(x = Prop, y = Agp, fill = Cluster), stat = "identity")
  
  res$stats_cluster_tip %>% 
    #filter(Cluster != "0") %>% 
    ggplot() +
    geom_bar(aes(x = Prop, y = tip, fill = Cluster), stat = "identity")
  
  
  res$g_labelled_q_sid <- dat_cluster %>% 
    ggplot(aes(x = ti, y = EQ5D)) + 
    geom_point(aes(colour = Cluster)) +
    geom_hline(data = res$stats_cluster_sid, aes(yintercept = mu)) +
    geom_hline(data = res$stats_cluster_sid, aes(yintercept = Q1), linetype = 2) +
    geom_hline(data = res$stats_cluster_sid, aes(yintercept = Q3), linetype = 2) +
    facet_wrap(.~SID)
  
  res$g_labelled_q_agp <- dat_cluster %>% 
    ggplot(aes(x = ti, y = EQ5D)) + 
    geom_point(aes(colour = Cluster)) +
    geom_hline(data = res$stats_cluster_agp, aes(yintercept = mu)) +
    geom_hline(data = res$stats_cluster_agp, aes(yintercept = Q1), linetype = 2) +
    geom_hline(data = res$stats_cluster_agp, aes(yintercept = Q3), linetype = 2) +
    facet_wrap(Agp~.)

  
  return(res)
}



