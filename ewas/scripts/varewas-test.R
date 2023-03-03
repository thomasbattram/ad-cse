library(devtools)
install_github("jrs95/jlst")
library(jlst)


# Data  
x <- rbinom(1000, 1, 0.5)
y <- 0.5 + 0.025 * x + rnorm(1000, 0, sqrt(0.005 * x)) + rnorm(1000, 0, 0.1)

vartest(y, x = as.factor(x), type = 1) # Breusch-Pagan test
vartest(y, x = as.factor(x), type = 2) # Brown-Forsythe test --- THIS ONE

# Joint location-and-scale test using Fisher's method
jlsp(y, x = as.factor(x), var.type = 1) # Breusch-Pagan variance test
jlsp(y, x = as.factor(x), var.type = 2) # Brown-Forsythe variance test

## 

all_dat <- data.frame(id = 1:1000, 
					  ad = x, 
					  y = y, 
					  covarA = rnorm(1000), 
					  covarB = rnorm(1000))

vartest(all_dat$y, x = as.factor(all_dat$x), type = 2)

res <- vartest(y = all_dat$y, 
			   x = as.factor(all_dat$ad), 
			   covar = all_dat[, -c(1,2,3)],
			   covar.var = TRUE,
			   type = 2)

tidy_var_res <- function(var_res)
{
	CO <- var_res$coef
	CO <- CO[which(rownames(CO) == "x1"),]
	CO[names(CO) == "Estimate"]
	out <- data.frame(beta = CO[names(CO) == "Estimate"], 
					  se = CO[names(CO) == "Std. Error"], 
					  p = CO[names(CO) == "Pr(>|t|)"]
					  )
	rownames(out) <- NULL
	return(out)
}

tidy_var_res(res)