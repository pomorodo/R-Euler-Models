#Euler Model
allele_freq<-function(k, Ne, tempo_sim, deltat){
  N<-round(tempo_sim/deltat)
  z<-matrix(NA, nrow=round(N)+1, ncol = k) 
  z[1,]<-rep(1/k,k)           #E(z(0))=1/k=z(0)
  ztilde<-z[, 1:(k-1)]        #z restriction
  #Ã generated as a time-dep hypermatrix,
  Atilde <- array(NA, dim = c(N, k-1, k-1))
  for (t in 1:N) {
    for (i in 1:(k-1)) {
      for (j in 1:(k-1)) {
        delta_ij <- ifelse(i == j, 1, 0)
        Atilde[t, i, j] <- (1 / (2 * Ne)) * ztilde[t, i] * (delta_ij - ztilde[t, j])
      }
    }
    #Writing the Euler model of the sde
    autov<-eigen(Atilde[t, , ], symmetric = TRUE)
    V_sqrt<-sqrt(pmax(autov$values, 0))  #Decomposizione spettrale invece di cholesky
    Q<-autov$vectors                     #ancora dubbi su valori negativi, come chol()
    
    L <- Q %*% diag(V_sqrt, nrow = k-1)
    
    dW <- rnorm(k - 1, mean = 0, sd = sqrt(deltat))
    ztilde[t+1, ] <- ztilde[t, ] + L %*% dW
    #if(sum(ztilde[t+1,])>1){break}
  }
  z[,1:(k-1)]<-ztilde
  for (h in 1:(N+1)) {
    z[h,k]<-1-sum(ztilde[h,])
  }
  
  return(ts(z, start = 0, deltat = deltat))
}
estinzioni<-function(num_sim=1000, k, Ne, tempo_sim, deltat){
  ext_times<-matrix(NA, nrow = num_sim, ncol = k) #Here we collect the extinction times
  
  for (s in 1:num_sim) {
    s_sim<-as.matrix(allele_freq(k, Ne, tempo_sim, deltat))
    #Looking for the first times in which the variants touch zero (here: < 10^-9)
    for (v in 1:K) {
      x <- which(as.matrix(s_sim)[, v] < 1e-9)[1] #x = first index | variants under threshold
      if (!is.na(x)) {
        ext_times[s, v] <- x * deltat 
      }
    }
    if (s %% 100 == 0) cat(sprintf("Simulazione %d / %d\n", s, num_sim))
  }
  colnames(ext_times) <- paste("Variant #", 1:k)
  return(ext_times)
}

#-------------------Simulazione--------------------------------
set.seed(111)
deltat <- 0.01
K <- 4
NE<-10
TSIM<-100
NSIM<-1000
#evoluzione <- allele_freq(k=K, Ne=NE, tempo_sim=tsim, deltat=deltat)
#colnames(evoluzione)<-paste("Variant #", 1:K)
extinctions<-estinzioni(num_sim = NSIM,k=K, Ne=NE, tempo_sim=TSIM, deltat=deltat)

#Plots---------------------------------------------------------
par(mfrow=c(2,2))

for (i in 1:K){
  ext_times<-extinctions[,i]
  ext_times<-ext_times[!is.na(ext_times)]
  if(length(ext_times)>0){
    hist(ext_times, breaks = 30, col = i, border = "white",
         main = paste("Variant #", i),
         xlab = "time of extinction", ylab = "frequence")
  }
  else{plot.new()
       title(paste("No extinction of variant #",i))}
  
}
par(mfrow=c(1,1))
#Plot of extinctions per variant
n_estinzioni <- colSums(!is.na(extinctions))
barplot(n_estinzioni,
        col    = 1:K,
        width  = 0.8,          
        space  = 0.3,
        names.arg = paste("V#", 1:K),
        main   = "Overall extinctions per variant",
        xlab   = "Variant", ylab   = "Number of extinctions",
        border = "white")

#Textual responses------------------------------------------------
# Extinctions Ratios
for (v in 1:K) {
  n_ext <- sum(!is.na(extinctions[, v]))         
  perc   <- round(100 * n_ext / NSIM, 1)           
  cat(sprintf("Variant #%d: extinct in %d/%d simulations  (%.1f%%)\n",
              v, n_ext, NSIM, perc))
}


# Amount of sim that had at least one extinction
sim_con_ext <- sum(rowSums(!is.na(extinctions)) > 0)
cat(sprintf("\n %d/%d (%.1f%%) simulations contain at least one extinction\n",
            sim_con_ext, NSIM, 100 * sim_con_ext / NSIM))


























