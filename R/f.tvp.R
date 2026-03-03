
f.tvp <- function(y,x,V=1,W=1,lambda=0.99,kappa=NULL,c=TRUE)
  {
    x <- as.matrix(x)
    if (ncol(x)==0) { c <- TRUE }

    xx <- cbind(rep.int(1,nrow(x)),x)

    if (c==TRUE)
      {
        theta <- matrix(0,ncol=1,nrow=ncol(xx))
      }
    else
      {
        theta <- matrix(0,ncol=1,nrow=ncol(xx)-1)
      }

    E <- diag(W,ncol=1+ncol(x),nrow=1+ncol(x))

    if (c==FALSE)
      {
        xx <- xx[,-1,drop=FALSE]
        E <- E[-1,-1,drop=FALSE]
      }

    tvpcppout <- kalman(xx,y,theta,V,E,lambda,kappa)

    thetas <- tvpcppout[[1]]
    y.tvp <- tvpcppout[[2]]
    pdensi <- tvpcppout[[3]]

    thetas <- t(thetas[,-ncol(thetas)])
    if ((c==FALSE && ncol(x)==1) || ncol(x)==0) { thetas <- t(thetas) }

    out <- list(as.vector(y.tvp),thetas,pdensi)
    names(out) <- c("y.hat","coef","pdens")
    
    return(out)
  }
