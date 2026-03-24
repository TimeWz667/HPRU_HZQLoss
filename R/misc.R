amlu <- list(
  A = mean,
  M = median,
  L = \(x) quantile(x, 0.025),
  U = \(x) quantile(x, 0.975)
)


save_tabs <- function(ss, prefix = "", folder = NA) {
  require(tidyverse)
  
  if (!is.na(folder)) {
    root <- here::here("docs", "tabs", folder)
    dir.create(root, showWarnings = F)
  } else {
    root <- here::here("docs", "tabs")
  }
  
  if (prefix != "") {
    prefix <- paste0(prefix, "_")
  }
  
  for (key in names(ss)) {
    write_csv(ss[[key]], here::here(root, paste0(prefix, key, ".csv")))
  }
}


save_tab <- function(ss, key, prefix = "", folder = NA) {
  require(tidyverse)
  
  if (!is.na(folder)) {
    root <- here::here("docs", "tabs", folder)
    dir.create(root, showWarnings = F)
  } else {
    root <- here::here("docs", "tabs")
  }
  
  if (prefix != "") {
    prefix <- paste0(prefix, "_")
  }
  
  write_csv(ss, here::here(root, paste0(prefix, key, ".csv")))
}

