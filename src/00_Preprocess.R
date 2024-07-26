library(here)
library(readr)

# Add the missing "End;" line at the end of the beast tree files ---------------

write_file("End;",
           here("data/real/bantu_ctmc-strict-bd/bantu_ctmc-strict-bd.trees"),
           append = TRUE
)
write_file("End;",
           here("data/real/bantu_ctmc-strict-bd-subsample/bantu_ctmc-strict-bd-subsample.trees"),
           append = TRUE
)
write_file("End;",
           here("data/real/bantu_ctmc-strict-bd-subsample2/bantu_ctmc-strict-bd-subsample2.trees"),
           append = TRUE
)
write_file("End;",
           here("data/real/iecor_ctmc-strict-M1/iecor_ctmc-strict-M1.trees"),
           append = TRUE
)
write_file("End;",
           here("data/real/st_ctmc-strict-fbd/st_ctmc-strict-fbd.trees"),
           append = TRUE
)
write_file("End;",
           here("data/real/st_ctmc-strict-fbd-by-sens/st_ctmc-strict-fbd.trees"),
           append = TRUE
)
write_file("End;", here("data/real/tea_ctmc-strict-fbd-constrained/tea_ctmc-strict-fbd-constrained.trees"),
           append = TRUE
)

