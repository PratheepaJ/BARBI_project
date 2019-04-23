set.seed(789)
N <- 5
probs <- c(.5, .3, .2 )
alpha0 <- .5
alpha <- alpha0*probs
alpha <- matrix( alpha, nrow=N, ncol=length(alpha), byrow=TRUE)
x <- sirt::dirichlet.simul(alpha)
dirichlet.mle(x)