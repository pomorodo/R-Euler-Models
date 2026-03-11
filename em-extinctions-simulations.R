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
  }
  colnames(ext_times) <- paste("Variant #", 1:k)
  return(ext_times)
}


set.seed(111)
deltat <- 0.01
K <- 4
NE<-10
tsim<-10
NSIM<-10
#evoluzione <- allele_freq(k=K, Ne=NE, tempo_sim=tsim, deltat=deltat)
#colnames(evoluzione)<-paste("Variant #", 1:K)
extinctions<-estinzioni(num_sim = NSIM,k=K, Ne=NE, tempo_sim=tsim, deltat=deltat)

#Plot
hist(extinctions,
     breaks = 40, col = "firebrick",
     main = "",
     xlab = "time", ylab = "number of simulations")
par(mfrow=c(2,2))







