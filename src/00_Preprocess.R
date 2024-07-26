library(here)
library(readr)

# Add the missing "End;" line at the end of the beast tree files ---------------
files <- list.files(here("data/real"), full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "\\.trees$")) |> 
  write_file("End;", append = TRUE)


