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
    if(sum(ztilde[t+1,])>1){break}
  }
  z[,1:(k-1)]<-ztilde
  for (h in 1:(N+1)) {
    z[h,k]<-1-sum(ztilde[h,])
  }
  
  return(ts(z, start = 0, deltat = deltat))
}
set.seed(111)
K<-4
evoluzione<-allele_freq(k=K, Ne=1000, tempo_sim = 100, deltat = 0.01)
colnames(evoluzione)<-paste("Variant #", 1:K)
#Plot 1 (each variant progression)
plot(evoluzione, col = 1:K, lwd=2,
     main = "k-allele variants progress",
     ylab = "Frequence", xlab = "Time")
#Plot 2 (multiprogression with legenda)
plot(evoluzione, plot.type = "single", col = 1:K, lwd=2,
     main = "k-allele variants progress",
     ylab = "Frequence", xlab = "Time")
legend("bottomleft",
       legend = colnames(evoluzione),
       col = 1:K,
       lwd = 2)
