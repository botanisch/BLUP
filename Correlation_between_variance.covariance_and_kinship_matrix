
   bv      <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/BLUP/master/raw_data/K-BLUP_estimates.csv"),header=T)      
   g.kin   <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/BLUP/master/raw_data/genomic_kinship_matrix.CSV"),row.names=1)
   var_cov <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/BLUP/master/raw_data/variance-covariance.csv"),row.names=1)

   g.kin <- g.kin[match(bv[,1],rownames(g.kin)),match(bv[,1],rownames(g.kin))] 
   colnames(p.kin) <- colnames(g.kin)
   rownames(p.kin) <- rownames(g.kin)


   cnt   <- 0
   vec_1 <- NULL
   vec_2 <- NULL
   for (i in 1:107){
     for (j in i:107){
        cnt        <- cnt + 1
        vec_1[cnt] <- g.kin[i,j]
        vec_2[cnt] <- p.kin[i,j]
   }}
   plot(vec_2,vec_1,xlab="Var(u) calculated by Equation 7",ylab="Var(u) calculated by λK")
