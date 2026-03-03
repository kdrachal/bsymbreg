
f.dma <- function(y,x,alpha=0.99,lambda=0.99,V=1,W=1,kappa=NULL,mods.incl=NULL,av=0)
  {
    if (is.null(mods.incl))
      {
        mods.incl <- expand.grid(rep.int(list(0:1),ncol(x)))
        mods.incl <- as.matrix(cbind(rep.int(1,nrow(mods.incl)),mods.incl))
      }
    c <- 0.001 / (2^nrow(mods.incl))
    pi1 <- as.vector(rep.int(1/(nrow(mods.incl)),nrow(mods.incl)))
    w <- as.vector(rep.int(0,nrow(mods.incl)))
    yhat.all <- matrix(0,ncol=nrow(mods.incl),nrow=1)
    ydma <- vector()
    thetas.exp <- as.vector(rep.int(0,ncol(mods.incl)))

    f.param <- function(i)
      {
        if (length(which(mods.incl[i,-1,drop=FALSE]==1))>0)
          {
            xx <- x[,which(mods.incl[i,-1,drop=FALSE]==1),drop=FALSE]
            if (mods.incl[i,1]==1)
              {
                c.incl <- TRUE
              }
            else
              {
                c.incl <- FALSE
              }
          }
        else
          {
            xx <- matrix(,ncol=0,nrow=length(y))
            c.incl <- TRUE
          }
        
        out <- f.tvp(y=y,x=xx,V=V,W=W,lambda=lambda,kappa=kappa,c=c.incl)
            
        return(out)
      }

    est.models <- lapply(seq(nrow(mods.incl)),f.param)

    f.pdens <- function(i)
      {
        return(.subset2(.subset2(est.models,i),3)[t])
      }

    f.mse <- function(i)
      {
        out.mse <- mean((y[1:t] - .subset2(.subset2(est.models,i),1)[1:t])^2)
        if (out.mse==0) { out.mse <- c }
        return(1/out.mse)
      }

    f.yhat <- function(i)
      {
        return(.subset2(.subset2(est.models,i),1)[t])
      }

    f.thetas <- function(i)
      {
        theta.i.tmp <- as.vector(mods.incl[i,])
        theta.i.tmp[theta.i.tmp==1] <- .subset2(.subset2(est.models,i),2)[t,]
        return(theta.i.tmp)
      }

    for (t in 1:nrow(x))
      {
        pi2 <- (pi1^alpha + c) / (sum((pi1)^alpha + c))
        w <- rbind(w,pi2) 

        yhat <- unlist(lapply(seq(nrow(mods.incl)),f.yhat))
        yhat.all <- rbind(yhat.all,yhat) 
        ydma[t] <- crossprod(pi2,yhat)
        thetas <- t(sapply(seq(nrow(mods.incl)),f.thetas))
        thetas.exp <- rbind(thetas.exp,pi2 %*% thetas)

        if (av==0) { pdens <- unlist(lapply(seq(nrow(mods.incl)),f.pdens)) }
        if (av==1) { pdens <- unlist(lapply(seq(nrow(mods.incl)),f.mse)) } 
        
        pi1 <- (pi2 * pdens) / as.numeric(crossprod(pi2,pdens))
     }

    thetas.exp <- thetas.exp[-1,,drop=FALSE]
    w <- as.matrix(w[-1,,drop=FALSE])
    yhat.all <- yhat.all[-1,,drop=FALSE]
    rvi <- w %*% mods.incl

    out <- list(ydma,yhat.all,w,thetas.exp,rvi)
    names(out) <- c("y.hat","y.hat.all","w","coef","rvi")

    return(out)
  }
  