
restructure_stan <- function(post) {
  su <- summary(post)$summary
  su <- data.frame(su) %>% 
    mutate(Var = row.names(su)) %>% 
    relocate(Var)
  rownames(su) <- NULL
  
  ext <- data.frame(rstan::extract(post))
  
  list(
    Ext = ext,
    Summary = su
  )
}


output_posterior <- function(res, tag) {
  tag <- glue::as_glue(tag)
  write_csv(res$Summary, file = here::here("docs", "tabs", "stats_" + tag + ".csv"))
  write_csv(res$Ext, file = here::here("posteriors", "post_" + tag + ".csv"))
  return(here::here("posteriors", "post_" + tag + ".csv"))
}

