## script to set up parallel processing with {future} and {furrr}

library(future)
library(furrr)

plan(multisession, workers = 6)


print("Parallel ready to go")
