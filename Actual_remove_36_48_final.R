# Define the parameters
n <- 42  # Number of trials (total mutations)
p <- 0.5  # Probability of success (proportion of shared mutations)

# Generate the probability values
probabilities <- dbinom(0:42, size = 42, prob = p)

# Update the first 7 rows with the specified probability value
probabilities[1:6] <- 0.001
probabilities[7] <- 0.137
# Convert the vector to a column vector
probabilities
# Create the data frame
probability <- data.frame(Probability = probabilities)


##############
# Define the number of rows and columns
num_rows <- 43
num_cols <- 47 

# Create a vector of zeros
zero_vector <- rep(0, num_rows * num_cols)

# Reshape the vector into a matrix with the desired dimensions
zero_matrix <- matrix(zero_vector, nrow = num_rows, ncol = num_cols)

# Create the data frame from the matrix
mutation <- as.data.frame(zero_matrix)


##############


# Concatenate probability and mutation data frames
details <- cbind(probability, mutation)

# Rename the data frame
details <- setNames(details, c("probability", paste0("Tumour Pair", 1:47)))#49

# View the resulting data frame



details[, 2] <- c(rep(0, 7), rep(2, 16), rep(0, nrow(details) - 23))
details[, 3] <- c(rep(0, 7), rep(2, 6), rep(0, nrow(details) - 13))
details[, 4] <- c(rep(0, 7), rep(2, 9), rep(0, nrow(details) - 16))
details[, 5] <- c(rep(0, 7), rep(2, 18), rep(0, nrow(details) - 25))
details[, 6] <- c(rep(0, 7), rep(2, 11), rep(0, nrow(details) - 18))
details[, 7] <- c(rep(1, 1), rep(0, 6),rep(2, 6), rep(0, nrow(details) - 13))
details[, 8] <- c(rep(0, 7), rep(2, 5), rep(0, nrow(details) - 12))
details[, 9] <- c(rep(0, 7), rep(2, 17), rep(0, nrow(details) - 24))
details[, 10] <- c(rep(0, 7), rep(2, 6), rep(0, nrow(details) - 13))
details[, 11] <- c(rep(0, 7), rep(2, 13), rep(0, nrow(details) - 20))
details[, 12] <- c(rep(0, 7), rep(2, 18), rep(0, nrow(details) - 25))
details[, 13] <- c(rep(0, 7), rep(2, 14), rep(0, nrow(details) - 21))
details[, 14] <- c(rep(0, 7), rep(2, 10), rep(0, nrow(details) - 17))
details[, 15] <- c(rep(0, 7), rep(2, 7), rep(0, nrow(details) - 14))
details[, 16] <- c(rep(0, 7), rep(2, 14), rep(0, nrow(details) - 21))
details[, 17] <- c(rep(0, 7), rep(2, 7), rep(0, nrow(details) - 14))
details[, 18] <- c(rep(0, 7), rep(2, 10), rep(0, nrow(details) - 17))
details[, 19] <- c(rep(0, 7), rep(2, 11), rep(0, nrow(details) - 18))
details[, 20] <- c(rep(0, 7), rep(2, 9), rep(0, nrow(details) - 16))
details[, 21] <- c(rep(0, 7), rep(2, 4), rep(0, nrow(details) - 11))
details[, 22] <- c(rep(0, 7), rep(2, 11), rep(0, nrow(details) - 18))
details[, 23] <- c(rep(0, 7), rep(2, 11), rep(0, nrow(details) - 18))
details[, 24] <- c(rep(0, 7), rep(2, 9), rep(0, nrow(details) - 16))
details[, 25] <- c(rep(0, 7), rep(2, 10), rep(0, nrow(details) - 17))
details[, 26] <- c(rep(0, 7), rep(2, 9), rep(0, nrow(details) - 16))
details[, 28] <- c(rep(0, 7), rep(2, 10), rep(0, nrow(details) - 17))
details[, 29] <- c(rep(0, 7), rep(2, 11), rep(0, nrow(details) - 18))
details[, 30] <- c(rep(0, 7), rep(2, 9), rep(0, nrow(details) - 16))
details[, 31] <- c(rep(0, 7), rep(2, 8), rep(0, nrow(details) - 15))
details[, 32] <- c(rep(0, 7), rep(2, 30), rep(0, nrow(details) - 37))
details[, 33] <- c(rep(0, 7), rep(2, 7), rep(0, nrow(details) - 14))
details[, 35] <- c(rep(0, 7), rep(2, 12), rep(0, nrow(details) - 19))
details[, 36] <- c(rep(0, 7), rep(2, 7), rep(0, nrow(details) - 14))
details[, 37] <- c(rep(0, 7), rep(2, 7), rep(0, nrow(details) - 14))
details[, 38] <- c(rep(0, 7), rep(2, 5), rep(0, nrow(details) - 12))
details[, 39] <- c(rep(0, 7), rep(2, 8), rep(0, nrow(details) - 15))
details[, 40] <- c(rep(0, 6),rep(1, 1),rep(2, 10), rep(0, nrow(details) - 17))
details[, 41] <- c(rep(0, 7), rep(2, 9), rep(0, nrow(details) - 16))
details[, 42] <- c(rep(0, 7), rep(2, 36))
details[, 43] <- c(rep(0, 6),rep(1, 1),rep(2, 3), rep(0, nrow(details) - 10))
details[, 44] <- c(rep(0, 7), rep(2, 7), rep(0, nrow(details) - 14))
details[, 45] <- c(rep(0, 7), rep(2, 4), rep(0, nrow(details) - 11))
details[, 46] <- c(rep(0, 7), rep(2, 4), rep(0, nrow(details) - 11))
details[, 47] <- c(rep(0, 7), rep(2, 3), rep(0, nrow(details) - 10))
details[, 48] <- c(rep(0, 6),rep(1, 1),rep(2, 5), rep(0, nrow(details) - 12))
details[, 49] <- c(rep(0, 7), rep(2, 12), rep(0, nrow(details) - 19))
details[, 50] <- c(rep(0, 7), rep(2, 4), rep(0, nrow(details) - 11))


# Convert details data frame to a matrix
details_matrix <- as.matrix(details)


#########################
# calculates the density of a log-normal distribution at specified points in the vector xigrid
xidens <- function(pmu, psig, xigrid) { sapply(xigrid, function(x) dlnorm(-log(1-x), pmu, psig)/(1-x)) }

tolerance=0.000001
pmu1=-2;psig1=1.5
pmu2=-1;psig2=1
pmu3=0.7;psig3=0.3
# Define xigrid
xigrid <- seq(0.000001, 0.999, length.out = 100)  # Example: 100 grid points from 0 to 1


# Compute the density
density1 <- xidens(pmu1, psig1, xigrid)
density2 <- xidens(pmu2, psig2, xigrid)
density3 <- xidens(pmu3, psig3, xigrid)


# Plot the density
plot(xigrid, density1,col="blue" ,type = "l", xlab = "x", ylab = "Density", main = "Log-normal Density Plot")
lines(xigrid, density2, col = "red", lwd = 2)
lines(xigrid, density3, col = "purple", lwd = 2)


############################################################################

# grid.lik=calculate the likelihood of a set of mutations 
grid.lik <- function(xigrid, mutns, probamut) {
  # find the index of shared and private mutations for each case
  mutidx <- apply(mutns, 2, function(muts) {
    ishared <- which(muts==1)
    iprivate <- which(muts==2)
    list(ishared, iprivate)
  })
  sapply(xigrid, function(xi, probamut, mutidx) {
    # p is the conditional prob of shared mutation
    p <- {xi + {1-xi}*probamut}/{xi + {1-xi}*{2-probamut}}
    # helper vectors
    logp <- log(p)
    log1p <- log(1-p)                    # ii represents an element of mutidx, ii is automatically taken from mutidx
    # unlist =convert the resulting list into a vector.
    # lapply = apply a given function to each element of a list or vector.
    unlist(lapply(mutidx, function(ii, logp, log1p) {  # sub function computes the likelihood for each case based on the indices of shared and private mutations, and the corresponding logarithms of conditional probabilities.   
      exp(sum(logp[ii[[1]]]) + sum(log1p[ii[[2]]]))  # sum(logp[ii[[1]]]=sum of logarithms of the conditional probabilities for shared mutations.
    }, logp, log1p))                                 # sum(log1p[ii[[2]]])=sum of logarithms of the conditional probabilities for private mutations.
  }, probamut, mutidx)
}

mutmat=details_matrix
probamut <- mutmat[, 1]
mutns <- mutmat[, -1]



#############################################################################

#____________________       Auxiliary functions


clonEM <- function(mutmat, init.para, xigrid, conv.crit, niter) {
  
  probamut <- mutmat[, 1]
  mutns <- mutmat[, -1]
  # number of cases
  ncase <- ncol(mutns)
  # number of grid points
  ngrid <- length(xigrid)
  # initialize parameters (mu, sigma, pi)
  para <- init.para
  # probability of matches for given clonality signal
  pYj0 <- grid.lik(0, mutns, probamut)
  pYj <- grid.lik(xigrid, mutns, probamut)
  # clonality signal density for xigrid
  g <- xidens(para[1], para[2], xigrid)
  gj <- rep(0, ngrid)
  vxij <- exij <- vpij <- rep(NA_real_, ncase)
  zz <- log(-log(1-xigrid))
  
  # E-step
  for(j in 1:ncase) {
    # product density
    gj <- pYj[j,]*g
    # integral
    igj <- sum(gj)/ngrid
    # conditional estimate of case being clonal (## clonality probabilities (vpij))
    vpij[j] <- para[3]*igj/((1-para[3])*pYj0[j] + para[3]*igj)
    # expected value of zz 
    exij[j] <- sum(zz*gj)/sum(gj)
    # variance of zz 
    vxij[j] <- sum((zz^2)*gj)/sum(gj) - exij[j]^2
  }
  #M-step
  # estimate the parameters from maximizing the conditional expectation
  para[1] <- sum(vpij*exij)/sum(vpij)
  para[2] <- max(0.1, sqrt(sum(vpij*(exij - para[1])^2)/sum(vpij)))
  para[3] <- sum(vpij)/ncase
  # iterate through the parameters
  ii <- 1
  conv <- 0
  while(conv == 0 & ii <= niter) {
    prev <- para
    g <- xidens(para[1], para[2], xigrid)
    for(j in 1:ncase) { 
      gj <- pYj[j,]*g
      igj <- sum(gj)/ngrid
      vpij[j] <- para[3]*igj/((1-para[3])*pYj0[j] + para[3]*igj)
      exij[j] <- sum(zz*gj)/sum(gj)
    }
    para[1] <- sum(vpij*exij)/sum(vpij)
    # variance is expectation of variance + variance of expectation
    vxi <- sum(vpij*vxij)/sum(vpij) + sum(vpij*(exij - para[1])^2)/sum(vpij)
    # lower bound sigma by 0.1
    para[2] <- max(0.1, sqrt(vxi))
    para[3] <- sum(vpij)/ncase
    #convergence reached
    #if( max(abs(para - prev)) < conv.crit ) conv <- 1
    if (!any(is.na(para)) && !any(is.na(prev)) && max(abs(para - prev), na.rm = TRUE) < conv.crit) conv <- 1
    
    ii <- ii + 1
  }
  return( list(para = para, n.iter = (ii-1), convergence = conv) )
}

# Prepare initial parameters (mu, sigma, pi)
init.para <- c(0.7, 0.3, 0.1)  # You may need to adjust these initial values based on your data and prior knowledge.

# Define xigrid
xigrid <- seq(0.00001, 0.999, length.out = 500)

# Set convergence criterion and maximum number of iterations
conv.crit <- 0.001
niter <- 1000

# Call the clonEM function
result_1 <- clonEM(mutmat, init.para, xigrid, conv.crit, niter)
result_1

cat("estimated proportion of clonality after removing 36 and 48 cases =",result_1$para[3])


