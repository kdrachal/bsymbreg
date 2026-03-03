
lbsymbreg <- function(y,x,alpha=0.99,lambda=0.99,V=1,W=1,kappa=NULL,mods.incl=0,av=0,
                      N.pop=NULL,p.mut=0.05,p.cross=0.90,verbose=FALSE) 
  {
    if (is.null(N.pop)) { N.pop <- round(0.25*2^ncol(x)) }
    if (!class(mods.incl)=="matrix")
      {
        mods.incl.type <- mods.incl
        if (mods.incl.type==0)
          {
            mods.incl <- expand.grid(rep.int(list(0:1),ncol(x)))
            mods.incl <- as.matrix(cbind(rep.int(1,nrow(mods.incl)),mods.incl))
            mods.incl <- mods.incl[sample(x=1:nrow(mods.incl),size=N.pop),]
          }
        if (mods.incl.type==1) 
          {
            mods.incl <- sample(c(0,1),ncol(x)*N.pop,replace=TRUE)
            mods.incl <- matrix(mods.incl,nrow=N.pop,byrow=TRUE)
            mods.incl <- as.matrix(cbind(rep.int(1,nrow(mods.incl)),mods.incl))
            mods.incl <- unique(mods.incl)
          }
        if (mods.incl.type==2)
          {
            mods.incl <- as.matrix(rbind(rep.int(0,ncol(x)),diag(1,nrow=ncol(x),ncol=ncol(x))))
            mods.incl <- as.matrix(cbind(rep.int(1,ncol(x)+1),mods.incl))
          }
      }
    yhat <- vector()
    N.mods <- vector()
    thetas.exp <- as.vector(rep.int(0,ncol(mods.incl)))
    rvi <- as.vector(rep.int(0,ncol(mods.incl)))

    for (T in 1:nrow(x))
      {
        if (verbose==TRUE) { print(paste('round:',T)) }
        
        T.dma <- f.dma(y=y[1:T,,drop=FALSE],x=x[1:T,,drop=FALSE],alpha=alpha,lambda=lambda,V=V,W=W,kappa=kappa,
                       mods.incl=mods.incl,av=av)
        
        yhat[T] <- T.dma$y.hat[T]
        mods.old <- mods.incl
        N.mods[T] <- nrow(mods.incl)
        thetas.exp <- rbind(thetas.exp,T.dma$coef[nrow(T.dma$coef),,drop=FALSE])  
        rvi <- rbind(rvi,T.dma$rvi[nrow(T.dma$rvi),,drop=FALSE])
        
        p.mutate <- rvi[nrow(rvi),-1,drop=FALSE]
        
        if (nrow(mods.incl) > N.pop)
          {
            ind.red <- sort(as.vector(T.dma$w[nrow(T.dma$w),,drop=FALSE]),decreasing=TRUE,index.return=TRUE)$ix
            ind.red <- ind.red[1:N.pop]
            mods.incl <- mods.incl[ind.red,,drop=FALSE]          
          }
        
        ind.mutate <- rbinom(n=nrow(mods.incl),size=1,prob=p.mut)
        if (!all(ind.mutate==0))
          {
            ind.mutate <- which(ind.mutate==1)
            for (i.mut in 1:length(ind.mutate))
              {
                mut.pos <- sample(x=1:ncol(x),size=1,prob=p.mutate)
                mods.incl[ind.mutate[i.mut],1+mut.pos] <- 1
              }
            mods.incl <- unique(mods.incl)
          }
        
        ind.cross <- rbinom(n=nrow(mods.incl),size=1,prob=p.cross)
        if (!all(ind.cross==0))
          {
            ind.cross <- which(ind.cross==1)
            if (length(ind.cross)>1)
              {
                if (!(length(ind.cross) %% 2 == 0))
                  {
                    ind.cross <- ind.cross[-sample(x=1:length(ind.cross),size=1)]
                  }
                ind.split <- split(sample(1:length(ind.cross)),rep(1:(length(ind.cross)/2),each=2))
                for (i.cross in 1:length(ind.split))
                  {
                    cross.pos <- sample(x=1:(ncol(x)-1),size=1)
                    c1 <- cbind(mods.incl[ind.split[[i.cross]][1],1:(1+cross.pos),drop=FALSE],
                                mods.incl[ind.split[[i.cross]][2],(1+cross.pos+1):ncol(mods.incl),drop=FALSE])
                    c2 <- cbind(mods.incl[ind.split[[i.cross]][2],1:(1+cross.pos),drop=FALSE],
                                mods.incl[ind.split[[i.cross]][1],(1+cross.pos+1):ncol(mods.incl),drop=FALSE])
                    mods.incl <- rbind(mods.incl,c1,c2)
                  }
                mods.incl <- unique(mods.incl)
              }
          }
      }
        
    thetas.exp <- thetas.exp[-1,]
    rvi <- rvi[-1,]
    out <- list(yhat,T.dma$y.hat.all,mods.old,N.mods,as.vector(T.dma$w[nrow(T.dma$w),]),thetas.exp,rvi)
    names(out) <- c("y.hat","y.hat.all","mods","N.mods","w","coef","rvi")
    
    return(out)
  }
