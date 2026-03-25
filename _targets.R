library(targets)
library(tarchetypes)
library(rstan)
library(tidyverse)


tar_option_set(
  packages = c("tibble", "tidyverse", "tidybayes", "rstan", "readxl", "glue", "lme4")
)
options(dplyr.summarise.inform = FALSE)
options(mc.cores = 8)
rstan_options(auto_write = TRUE)
theme_set(theme_bw())


tar_source()




fig_path <- here::here("docs", "figs")
fig_ext <- ".png"
tab_path <- here::here("docs", "tabs")

dir.create("posteriors/", showWarnings = F)

dir.create("docs/", showWarnings = F)
dir.create(fig_path, showWarnings = F)
dir.create(tab_path, showWarnings = F)


n_mc <- 2000
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_target(file_qol, here::here("data", "generated", "data_generated_11667.csv"), format = "file"),
  tar_target(file_norm, here::here("data", "processed", "pn_mapped.csv"), format = "file"),
  
  ## extraction
  tar_target(data_norm, read_csv(file_norm)),
  
  ### Perfect health baseline
  tar_target(data_raw, get_data_qol(file_qol)),
  tar_target(data_tte, format_tte(data_raw)),
  tar_target(data_qol, format_qol(data_raw)),

  ### Population norm baseline
  tar_target(data_raw_pn, get_data_pn(data_raw, data_norm)),
  tar_target(data_tte_pn, format_tte(data_raw_pn)),
  tar_target(data_qol_pn, format_qol(data_raw_pn)),

  
  # modelling, TTE
  tar_target(file_model_tte, here::here("models", "time2zero_surv_age.stan"), format = "file"),
  tar_target(model_tte, stan_model(file_model_tte)),
  
  tar_target(pars_tte, fit_tte(model_tte, data_tte, n_iter = n_mc)),
  tar_target(gs_tte, visualise_tte(data_tte, pars_tte)),
  tar_target(file_posterior_tte, output_posterior(pars_tte, "ph_tte"), format = "file"),
  tar_target(file_vis_tte, output_vis_tte(gs_tte, folder = fig_path, "ph", ext = fig_ext)),
  
  tar_target(pars_tte_pn, fit_tte(model_tte, data_tte_pn, n_iter = n_mc)),
  tar_target(gs_tte_pn, visualise_tte(data_tte_pn, pars_tte_pn)),
  tar_target(file_posterior_tte_pn, output_posterior(pars_tte_pn, "pn_tte"), format = "file"),
  tar_target(file_vis_tte_pn, output_vis_tte(gs_tte_pn, folder = fig_path, "pn", ext = fig_ext)),
  
  # modelling, QoL
  tar_target(pars_qol, fit_qol_bayes(data_qol)),
  tar_target(file_pars_qol, summarise_qol_bayes(pars_qol, "full"), format = "file"),


  # simulate QALY loss
  tar_target(pars_qloss, boot_pars_bayes(file_posterior_tte, file_pars_qol, n_sim = n_mc)),
  tar_target(sim_qloss, simulate_ql(pars_qloss, age0 = 50, age1 = 99, age_until = 5, dt = 0.01)),
  tar_target(tab_qloss, summarise_ql(pars_qloss, sim_qloss)),
  tar_target(out_gloss, save_tab(tab_qloss, key = "summary_qloss_ph"))

  # tar_target(pars_qloss_pn, boot_pars_bayes(file_posterior_tte_pn, file_pars_qol, n_sim = n_mc)),
  # tar_target(sim_qloss_pn, simulate_ql_pn(pars_qloss_pn, data_norm, age0 = 50, age1 = 99, age_until = 5, dt = 0.01)),
  # # tar_target(tab_qloss_pn_sub, summarise_ql_pn(sim_qloss_pn)),
  # # tar_target(out_gloss_pn_sub, save_tab(tab_qloss_pn_sub, key = "summary_qloss_pn_sub")),
  # tar_target(tab_qloss_pn_sub_short, simplify_ql_pn_sub_tabs(tab_qloss_pn_sub)),
  # tar_target(out_gloss_pn_sub_short, save_tab(tab_qloss_pn_sub_short, key = "summary_qloss_pn_sub")),
  # 
  # tar_target(data_norm_pooled, get_norm_wt(data_raw_pn, data_norm)),
  # tar_target(sim_qloss_pn_pooled, simulate_ql_pn(pars_qloss_pn, data_norm_pooled, age0 = 50, age1 = 99, age_until = 5, dt = 0.01)),
  # tar_target(tab_qloss_pn, summarise_ql(pars_qloss_pn, sim_qloss_pn_pooled)),
  # tar_target(out_gloss_pn, save_tab(tab_qloss_pn, key = "summary_qloss_pn")),
  # 
  # 
  # tar_target(tab_qloss_bind, bind_ql_tabs(tab_qloss_ph = tab_qloss, tab_qloss_pn, tab_qloss_pn, ages = seq(50, 90, 10))),
  # tar_target(out_gloss_bind, save_tab(tab_qloss_bind, key = "summary_qloss_bind"))
)
