
# This computation used GWASpro code. 
# Please cite the GWASpro paper to use this code. Below is the citation.
# Kim, Bongsong, Xinbin Dai, Wenchao Zhang, Zhaohong Zhuang, Darlene L. Sanchez, Thomas Lübberstedt, Yun Kang, Michael K. Udvardi, William D. Beavis, Shizhong Xu*, Patrick X. Zhao* "GWASpro: A High-Performance Genome-Wide Association Analysis Server" Bioinformatics (2018): doi.org/10.1093/bioinformatics/bty989.
# GWASpro URL: https://bioinfo.noble.org/GWASPRO


 install.packages("RCurl")                  
 require(RCurl) 

 dat <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/BLUP/master/raw_data/pheno.csv"),header=T)           
 x   <- model.matrix(~1 + as.factor(dat$loc)) 
 y <- dat[,3]
 z   <- model.matrix(~dat[,1]-1)
              
 kk  <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/BLUP/master/raw_data/genomic_kinship_matrix.CSV"),row.names=1)
 kk  <- kk[match(sort(unique(dat[,1])),rownames(kk)),match(sort(unique(dat[,1])),rownames(kk))] 

 id_n <- rownames(kk)

  library(MASS)
  print(Sys.time()) 
  c11 <- t(x)%*%x
  c12 <- t(x)%*%z
  c21 <- t(c12)
  c22 <- t(z)%*%z
  w1  <- t(x)%*%y                           
  w2  <- t(z)%*%y                                    
  n   <- nrow(x)
  m   <- ncol(z)
  q   <- ncol(x)
 
  kk <- as.matrix(kk) + diag(c(rep(1,m))) * 1e-6                                     
  b  <<- 0
  u  <<- 0 
  ve <<- 1
  vg <<- 1     
  parm0<<-rep(1,q+2)
  ki  <- solve(kk)
  c22_i <- solve(c22)
  max.iter = 500
  min.err  = 1e-8  
     qq1 <- eigen(c22_i %*% ki,symmetric=F)
     u1  <- qq1$vectors
     d1  <- qq1$values
  
     m1  <- c12 %*% u1
     m2  <- t(u1) %*% c22_i %*% c21
     m3  <- t(u1) %*% c22_i %*% w2
     m4  <- t(u1) %*% c22_i 

 
  iter <<- 0                           
  err  <<- 1e8
  vcom <<- NULL
  while(iter < max.iter & err > min.err){ 
     h1   <- diag(vg/(vg+ve*d1))
     b    <- solve(c11 - m1 %*% h1 %*% m2) %*% (w1 - m1 %*% h1 %*% m3)    
     u    <- u1 %*% h1 %*% m4 %*% (w2 - c21 %*% b)                             
     r    <- y - x%*%b-z%*%u
     d22   <- solve((ve/vg)*diag(rep(1,dim(kk)[1])) + kk %*% (c22 - c21 %*% solve(c11) %*% c12))  
     vg   <<- (drop(t(u)%*%ki%*%u) + ve*sum(diag(d22))) /m   
     ve   <<- sum(r*y)/(n-q) 
     parm <- c(b,vg,ve)
     err  <-  sum((parm-parm0)^2)
     parm0<- parm
     iter <<- iter + 1
     cat("Iteration"," #",iter,":\tVariance component for fixed effect=",ve,"\tVariance component for random effect=",vg)
     cat("\n")        
     vcom  <<- rbind(vcom,paste(iter,vg,ve,Sys.time(),sep=","))
  } 
  head <- paste("Iteratioin","Vg","Ve","time",sep=",")
  vcom <- rbind(head,vcom)
  tail <- paste("Heritability:",vg / (vg + ve),"(H^2 = Vg /(Vg + Ve))",sep=",")
  vcom <- rbind(vcom,tail)

  u  <- cbind(id_n,u)
  colnames(u) <- c("ID","Breeding values")
  u  <- data.frame(u)

write.table(u, file = "K-BLUP.csv",quote=F,sep=",",row.names=F,col.names=T)
write.table(vcom, file = "variance_component.csv",quote=F,row.names=F,col.names=F,sep=",")
