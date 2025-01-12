                                          
   install.packages("RCurl")                  
   require(RCurl) 

# Figure 1(A) Correlation between Naive-BLUP and K-BLUP   

   dat1 <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/BLUP/master/raw_data/Naive-BLUP_estimates.csv"),header=T)
   dat2 <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/BLUP/master/raw_data/K-BLUP_estimates.csv"),header=T)   
   plot(dat1[,2],dat2[,2],xlab="Naive-BLUP",ylab="K-BLUP")
   
   
# Figure 1(B) Correlation between Naive-BLUP and Average_vector   
   
   dat1 <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/BLUP/master/raw_data/Naive-BLUP_estimates.csv"),header=T)
   dat2 <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/BLUP/master/raw_data/pheno.csv"),header=T)   
  
   tmp <- NULL
   for (i in 1:dim(dat1)[1]){
      tmp[i] <- mean(as.numeric(dat2[which(dat2[,1] %in% dat1[i,1]),3])) 
   }
   plot(tmp,dat1[,2],xlab="Average-vector",ylab="Naive-BLUP")
 
 
# Figure 1(C) Correlation between K-BLUP and Average_vector   

   dat1 <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/BLUP/master/raw_data/K-BLUP_estimates.csv"),header=T)
   dat2 <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/BLUP/master/raw_data/pheno.csv"),header=T)   
   tmp <- NULL
   for (i in 1:dim(dat1)[1]){
      tmp[i] <- mean(as.numeric(dat2[which(dat2[,1] %in% dat1[i,1]),3])) 
   }
   plot(tmp,dat1[,2],xlab="Average-vector",ylab="K-BLUP")
