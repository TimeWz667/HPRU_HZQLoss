library(targets)
library(tarchetypes) # Load other packages as needed.







# Set target options:
tar_option_set(
  packages = c("tidyverse", "rstan")
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_target(
    name = data,
    command = tibble(x = rnorm(100), y = rnorm(100))
    # format = "qs" # Efficient storage for general data objects.
  ),
  tar_target(
    name = model,
    command = coefficients(lm(y ~ x, data = data))
  )
)
