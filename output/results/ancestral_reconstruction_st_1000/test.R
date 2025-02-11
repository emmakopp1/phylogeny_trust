library(parallel)

# Nombre de cœurs disponibles
ncl = 30
cl = makeCluster(ncl, type="PSOCK")
clusterSetRNGStream(cl)
       
sink("/home/users/kopp/work/ancestral_reconstruction/test.out", append = T)
t1 <- Sys.time()
res_list <- parLapply(cl, 1:20, function(k) {k^2})
t2 <- Sys.time()
cat("For ", ncl, "nodes and type PSOCK, the time of execution is", difftime(t2, t1, units = "secs"), "\n")
stopCluster(cl)
sink()






