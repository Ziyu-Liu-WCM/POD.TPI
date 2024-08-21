
library(dfcrm)
library(BOIN)
library(truncdist)
library(MCMCpack)
library(poisbinom)
library(zoo)

Mode = function(x){
  ux = unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


pava = function(x, wt = rep(1, length(x))){
  n = length(x)
  if (n <= 1) 
    return(x)
  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets = (1:n)
  repeat {
    viol = (as.vector(diff(x)) < 0)
    if (!(any(viol))) 
      break
    i = min((1:(n - 1))[viol])
    lvl1 = lvlsets[i]
    lvl2 = lvlsets[i + 1]
    ilvl = (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] = sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] = lvl1
  }
  x
}

betavar = function(a, b){
  var = a * b / ((a+b)^2*(a+b+1))
  return(var)
}


GetMTD = function(p, p_T, epsilon){
  MTD = which((p + 1e-10 >= p_T - epsilon[1]) & (p - 1e-10 <= p_T + epsilon[2]))
  
  if(length(MTD) == 0){
    MTD = which(p <= p_T)
    if(length(MTD) > 0){
      MTD = max(MTD)
    } else{
      MTD = c(-2, -1, 0)
    }
  }
  
  return(MTD)
}


GetMTD_unique = function(p, p_T, epsilon){
  
  
  MTD = which((p + 1e-10 >= p_T - epsilon[1]) & (p - 1e-10 <= p_T + epsilon[2]))
  
  if(length(MTD) > 1){
    
    D = length(p)
    p = p + (1:D)*1E-10
    
    MTD = which.min(abs(p - p_T))
    if(length(MTD) > 1){
      p_diff = abs(p[MTD] - p_T)
      if(sum(p[MTD] == p_T - p_diff) >= 1){
        MTD = max(which(p == p_T - p_diff))
      } else{
        MTD = min(MTD)
      }
    }
  } else if(length(MTD) == 0){
    MTD = which(p <= p_T)
    if(length(MTD) > 0){
      MTD = max(MTD)
    } else{
      MTD = 0
    }
  }
  
  return(MTD)
}

beta_fun = function(t, h){
# return [beta(t, 1), ..., beta(t, K)]
# h is [h_1, ..., h_K] but does not include h_0, which is always 0
# h[K] = W

  K = length(h)
  beta_all = rep(0, K)
  h2 = c(0, h[1:(K-1)])
  
  beta_all[t > h] = 1
  beta_all[t > h2 & t <= h] = (t - h2[t > h2 & t <= h]) / (h[t > h2 & t <= h] - h2[t > h2 & t <= h])
  
  return(beta_all)
  
}


gen_weibull_params = function(p_d, W, W_star, q){
# generates the shape and scale parameters for a Weibull distribution such that
# T ~ Weibull(shape, scale), Pr(T <= W) = p_d, Pr(T <= W_star) = (1-q)*p_d (for W_star < W)
  
  shape = log(log(1 - p_d) / log(1 - p_d + q * p_d)) / (log(W / W_star))
  scale = W/((-log(1 - p_d))^(1 / shape))

  return(c(shape, scale))

}



gen_DLT_ACC_table = function(n, p, W, W_star, q, xi){

# return: a n * (D+1) dimensional matrix

# each row: the first D entries record the DLT times of a patient should the patient
#           be treated at dose d, and the last entry is the time of this patient 
#           becoming available for enrollment

  D = length(p)
  DLT_ACC_table = matrix(NA, n, D+1)
  # weibull_params_all: 2 * D matrix
  weibull_params_all = sapply(p, gen_weibull_params, W = W, W_star = W_star, q = q)
  # apply: returns n * D matrix
  DLT_ACC_table[ , 1:D] = apply(weibull_params_all, 2, function(x){rweibull(n = n, shape = x[1], scale = x[2])})
  DLT_ACC_table[ , D+1] = cumsum(c(0, rexp(n = n-1, rate = xi)))
  
  return(DLT_ACC_table)
}
