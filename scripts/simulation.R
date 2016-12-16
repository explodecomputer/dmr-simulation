#' Simulate an outcome based on exposures and effects of those exposures on the outcome
#'
#' Can take a single or multiple exposures (vector or matrix).
#' Each exposure must have an effect. The exposure will be scaled to have unit variance and 
#' zero mean therefore the effect size is equivalent to the correlation coefficient.
#' 
#'
#' @param x Vector or matrix of exposure variables
#' @param r Vector of correlation coefficients, must be same length as ncol(x), or length of 1 if x is a vector
#'
#' @example
#'
#' # e.g If you want to simulate a variable y that is influenced by x1 and x2, and they each 
#' explain 10% of the variance, and x1 increases y, while x2 decreases y, then
#' x <- matrix(rnorm(1000 * 2), 1000, 2)
#' y <- sim_y_from_x(x, c(sqrt(0.1), -sqrt(0.1)))
#'
#' @export
#' @return Vector of y
sim_y_from_x <- function(x, r)
{
	x <- scale(x)
	stopifnot(ncol(x) == length(r))
	stopifnot(sum(r^2) <= 1)
	e <- scale(rnorm(nrow(x))) * sqrt(1-sum(r^2))
	g <- x %*% r
	y <- g + e
	return(y)
}



# Model 1

# Each CpG influenced by one SNP
# Each CpG has an independent effect on the outcome

nid <- 2000
ncpg <- nsnp <- 10

r_cpg_y <- sample(c(0.06, -0.06), ncpg, replace=TRUE)
r_geno_cpg <- rep(0.3,ncpg)

geno <- matrix(rbinom(nid * nsnp, 2, 0.5), nid, ncpg)

cpg <- matrix(0, nid, ncpg)
for(i in 1:ncpg)
{
	cpg[,i] <- sim_y_from_x(geno[,i], r_geno_cpg[i])
}

y <- sim_y_from_x(cpg, r_cpg_y)

summary(lm(y ~ cpg))



# Model 2

# Same as model 1 but the the SNP-CpG are more complex

nid <- 2000
ncpg <- 10
nsnp <- 20

r_geno_cpg <- matrix(rnorm(nsnp * ncpg, sd = 0.1), ncpg, nsnp)
r_geno_cpg[sample(1:length(r_geno_cpg), length(r_geno_cpg) * 0.5, replace=FALSE)] <- 0
r_cpg_y <- sample(c(0.06, -0.06), ncpg, replace=TRUE)

geno <- matrix(rbinom(nid * nsnp, 2, 0.5), nid, nsnp)

cpg <- matrix(0, nid, ncpg)
for(i in 1:ncpg)
{
	cpg[,i] <- sim_y_from_x(geno, r_geno_cpg[i,])
}

y <- sim_y_from_x(cpg, r_cpg_y)

summary(lm(y ~ cpg))
summary(lm(cpg[,5] ~ geno))


# Model 3

# G causes a latent variable (gene expression level) and this causes CpG and y

nid <- 2000
ncpg <- 10
nsnp <- 10

r_geno_l <- rep(0.1, nsnp)
geno <- matrix(rbinom(nid * nsnp, 2, 0.5), nid, nsnp)

r_l_cpg <- rnorm(ncpg, sd=0.3)

l <- sim_y_from_x(geno, r_geno_l)
summary(lm(l ~ geno))

cpg <- matrix(0, nid, ncpg)
for(i in 1:ncpg)
{
	cpg[,i] <- sim_y_from_x(l, r_l_cpg[i])
}

y <- sim_y_from_x(l, 0.2)
summary(lm(y ~ cpg))



# Model 4

# Mixture of models 3 and 4

nid <- 2000
ncpg_in <- 5
ncpg_out <- 5
nsnp <- 20

r_geno_cpg <- matrix(rnorm(nsnp * ncpg_in, sd = 0.1), ncpg_in, nsnp)
r_geno_cpg[sample(1:length(r_geno_cpg), length(r_geno_cpg) * 0.5, replace=FALSE)] <- 0
r_cpg_l <- sample(c(0.06, -0.06), ncpg_in, replace=TRUE)

geno <- matrix(rbinom(nid * nsnp, 2, 0.5), nid, nsnp)

cpg_in <- matrix(0, nid, ncpg_in)
for(i in 1:ncpg_in)
{
	cpg_in[,i] <- sim_y_from_x(geno, r_geno_cpg[i,])
}

l <- sim_y_from_x(cpg_in, r_cpg_l)

r_l_cpg <- rnorm(ncpg_out, sd=0.3)

cpg_out <- matrix(0, nid, ncpg_out)
for(i in 1:ncpg_out)
{
	cpg_out[,i] <- sim_y_from_x(l, r_l_cpg[i])
}

y <- sim_y_from_x(l, 0.2)

