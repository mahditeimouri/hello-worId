dst0 <- function(x, mu, sigma, lambda, nu){z <- (x - mu)/sigma; 2/sigma*dt(z, df = nu)*pt(lambda*z*sqrt( (1 + nu)/(z^2 + nu) ), df = nu + 1 )}
mlest <- function(x){
    n <- length(x)
    N <- 1000
    cri <- 10e-5;
    j <- 2
    eps <- 1
    s2 <- s3 <- Mu <- M <- K <- A <- d <- PDF <- tau <- rep(NA, n)
    out <- matrix(NA, ncol = 4, nrow = N)
    mu <- median(x)
    m3 <- mean( (x - mu)^3 )
    sigma <- sqrt( var(x) )[1]
    lambda <- sign( m3/sigma^3 )[1]
    nu <- 2
    out[1, ] <- c(mu, sigma, lambda, nu)
    Del <- lambda/sqrt( 1 + lambda^2 )*sigma
    Gam <- (1 - lambda^2/(1 + lambda^2) )*sigma^2
    while( eps > 0.5 & j < N ){
    d <- ( ( x - mu )/sigma )^2
    A <- lambda*( x - mu )/sigma
    Mu <- Del/(Gam + Del^2)*( x - mu )
    M <- sqrt(Gam/(Gam + Del^2))
    PDF <- dst0(x, mu, sigma, lambda, nu)
        K <- 4*nu^(nu/2)*gamma( (nu + 3)/2 )*(d + nu)^( -(nu + 3)/2 )*pt( sqrt( (nu + 3)/(d + nu) )*A, nu + 3)/(PDF*gamma(nu/2)*sqrt(pi)*sigma)
        tau <- 2*nu^(nu/2)*gamma( (nu + 2)/2 )*( d + nu + A^2 )^( -(nu + 2)/2 )/( PDF*gamma(nu/2)*pi*sigma )
        s2 <- K*Mu + M*tau
        s3 <- K*Mu^2 + M^2 + Mu*M*tau
      mu <- sum(K*x - Del*s2)/sum(K)
      Del <- sum((x - mu)*s2)/sum(s3)
      Gam <- mean( K*(x - mu)^2 - 2*(x - mu)*Del*s2 + Del^2*s3 )
     sigma <- sqrt(Del^2 + Gam)
      lambda <- Del/sqrt(Gam)
      obj <- function(par){ sum( log( dst0(x, mu, sigma, lambda, par[1]) ) )}
      nu <- suppressWarnings( optimize(obj, c(1, 100), tol = 0.000001, maximum = TRUE)$maximum )
      out[j, ] <- c(mu, sigma, lambda, nu)
      if (sum( abs( out[j - 1, ] - out[j, ] ) ) < cri || j >= (N - 1) ){
        eps <- 0
      }else{ j <- j + 1 }
    }
    return( list(mu = out[j - 1, 1], sigma = out[j - 1, 2], sigma = out[j - 1, 2], lambda = out[j - 1, 3], nu = out[j - 1, 4]) )
  }
