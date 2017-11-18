# ---- newVB.funcs ----

new_varbvsnormupdate <- function (X, sigma, sa, xy, d, alpha0, mu0, Xr0) {

  # Get the number of samples (n) and variables (p).
  n <- nrow(X)
  p <- ncol(X)

  L = nrow(alpha0) # alpha0 and mu0 must be L by p
  pi = rep(1,p)

  # Check input X.
  if (!is.double(X) || !is.matrix(X))
    stop("Input X must be a double-precision matrix")

  # Check inputs sigma and sa.
  if (length(sigma) != 1 | length(sa) != 1)
    stop("Inputs sigma and sa must be scalars")


  # Check input Xr0.
  if (length(Xr0) != n)
    stop("length(Xr0) must be equal to nrow(X)")


  # Initialize storage for the results.
  alpha <- alpha0
  mu    <- mu0
  Xr    <- Xr0

  # Repeat for each effect to update
  for (l in 1:L) {
    # remove lth effect
    Xr = Xr - X %*% (alpha[l,]*mu[l,])

    # Compute the variational estimate of the posterior variance.
    s <- sa*sigma/(sa*d + 1)

    # Update the variational estimate of the posterior mean.
    mu[l,] <- s/sigma * (xy - t(X) %*% Xr)

    # Update the variational estimate of the posterior inclusion
    # probability. This is basically prior (pi) times BF.
    # The BF here comes from the normal approx - could be interesting
    # to replace it with a t version that integrates over sigma?
    alpha[l,] <- pi*exp((log(s/(sa*sigma)) + mu[l,]^2/s)/2)

    alpha[l,] <- alpha[l,]/sum(alpha[l,])

    # Update Xr by adding back in the $l$th effect
    Xr <- Xr + X %*% (alpha[l,]*mu[l,])
  }

  return(list(alpha = alpha,mu = mu,Xr = Xr,s=s))
}

# Just repeated applies those updates
# X is an n by p matrix of genotypes
# Y a n vector of phenotypes
# sa the variance of the prior on effect sizes (actually $\beta \sim N(0,sa sigma)$ where the residual variance sigma here is fixed to the variance of Y based on a small effect assumption.)
fit = function(X,Y,sa=1,sigma=NULL,niter=100,L=5,calc_elbo=FALSE){
  if(is.null(sigma)){
    sigma=var(Y)
  }
  p =ncol(X)
  xy = t(X) %*% Y
  d = colSums(X * X)
  alpha0= mu0 = matrix(0,nrow=L,ncol=p)
  Xr0 = X %*% colSums(alpha0*mu0)
  elbo = rep(NA,niter)
  for(i in 1:niter){
    res = new_varbvsnormupdate(X, sigma, sa, xy, d, alpha0, mu0, Xr0)
    alpha0 = res$alpha
    mu0 = res$mu
    Xr0 = res$Xr
    if(calc_elbo){
      elbo[i] = elbo(X,Y,sigma,sa,mu0,alpha0)
    }
  }
  return(c(res,list(elbo=elbo)))
}

#this is for scaled prior in which effect prior variance is sa * sigma
elbo = function(X,Y,sigma,sa,mu,alpha){
  L = nrow(alpha)
  n = nrow(X)
  p = ncol(X)

  Xr = (alpha*mu) %*% t(X)
  Xrsum = colSums(Xr)

  d = colSums(X*X)
  s <- sa*sigma/(sa*d + 1)

  postb2 = alpha * t(t(mu^2) + s)

  Eloglik =  -(n/2) * log(2*pi* sigma) -
    (1/(2*sigma)) * sum(Y^2) +
    (1/sigma) * sum(Y * Xrsum) -
    (1/(2*sigma)) * sum(Xrsum^2) +
    (1/(2*sigma)) * sum((Xr^2)) -
    (1/(2*sigma)) * sum(d*t(postb2))

  KL1 = sum(alpha * log(alpha/(1/p)))
  KL2 = - 0.5* sum(t(alpha) * (1 + log(s)-log(sigma*sa)))
  + 0.5 * sum(postb2)/(sigma*sa)

  return(Eloglik - KL1 - KL2)
}


# This computes the average lfsr across SNPs for each l, weighted by the
# posterior inclusion probability alpha
lfsr_fromfit = function(res){
  pos_prob = pnorm(0,mean=t(res$mu),sd=sqrt(res$s))
  neg_prob = 1-pos_prob
  1-rowSums(res$alpha*t(pmax(pos_prob,neg_prob)))
}

#find how many variables in the 95% CI
# x is a probability vector
n_in_CI_x = function(x){
  sum(cumsum(sort(x,decreasing = TRUE))<0.95)+1
}

# return binary vector indicating if each point is in CI
# x is a probability vector
in_CI_x = function(x){
  n = n_in_CI_x(x)
  o = order(x,decreasing=TRUE)
  result = rep(0,length(x))
  result[o[1:n]] = 1
  return(result)
}

# Return binary matrix indicating which variables are in CI of each
# of effect
in_CI = function(res){
  t(apply(res$alpha,1,in_CI_x))
}

n_in_CI = function(res){
  apply(res$alpha,1,n_in_CI_x)
}

# computes z score for association between each
# column of X and y
calc_z = function(X,y){
  z = rep(0,ncol(X))
  for(i in 1:ncol(X)){
    z[i] = summary(lm(y ~ X[,i]))$coeff[2,3]
  }
  return(z)
}


# plot p values of data and color in the 95% CIs
# for simulated data, specify b = true effects (highlights in red)
pplot = function(X,y,res,pos=NULL,b=NULL,CImax = 400,...){
  z = calc_z(X,y)
  zneg = -abs(z)
  logp = log10(pnorm(zneg))
  if(is.null(b)){b = rep(0,ncol(X))}
  if(is.null(pos)){pos = 1:ncol(X)}
  plot(pos,-logp,col="grey",xlab="",ylab="-log10(p)",...)
  points(pos[b!=0],-logp[b!=0],col=2,pch=16)
  for(i in 1:nrow(res$alpha)){
    if(n_in_CI(res)[i]<CImax)
      points(pos[which(in_CI(res)[i,]>0)],-logp[which(in_CI(res)[i,]>0)],col=i+2)
  }

}
