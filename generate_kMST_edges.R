library(MASS)
library(ade4)


generate_kMST_edges <- function(tau_, n, l, p,
                                rho1, beta1, epsilon1, nu11, nu12,
                                rho2, beta2, epsilon2, nu21, nu22,
                                sigma, k = 3, seed = 16) {
  
  set.seed(seed)
  
  generate_mat <- function(beta, epsilon, nu1, nu2, rho, l, p, sigma) {
    within_subject <- diag(1 - rho, l) + matrix(rho, nrow = l, ncol = l)
    covariance_matrix <- kronecker(within_subject, diag(p)) * sigma^2
    
    ak <- mvrnorm(n = 1, mu = beta, Sigma = diag(epsilon^2, p))
    mean_vector <- rep(ak, l)
    theta <- mvrnorm(n = 1, mu = mean_vector, Sigma = covariance_matrix)
    omega_u <- runif(1, nu1, nu2)
    Z <- mvrnorm(n = 1, mu = theta, Sigma = (omega_u^2) * diag(l*p))
    Z_mat <- matrix(Z, nrow = l, ncol = p, byrow = TRUE)
    return(Z_mat)
  }
  
  generate_n_mat <- function(n, beta, epsilon, nu1, nu2, rho, l, p, sigma){
    result_mat <- matrix(nrow = 0, ncol = p)
    for (i in 1:n) {
      mat <- generate_mat(beta, epsilon, nu1, nu2, rho, l, p, sigma)
      result_mat <- rbind(result_mat, mat)
    }
    return(result_mat)
  }
  
  find_index_row <- function(k, ncol = l) {
    row <- ((k - 1) %/% ncol) + 1
    return(row)
  }
  
  # generate data
  group1_matrix <- generate_n_mat(tau_, beta1, epsilon1, nu11, nu12, rho1, l, p, sigma)
  group2_matrix <- generate_n_mat(n - tau_, beta2, epsilon2, nu21, nu22, rho2, l, p, sigma)
  group_matrix <- rbind(group1_matrix, group2_matrix)
  
  distance_matrix <- as.matrix(dist(group_matrix))
  kmst_matrix <- ade4::mstree(as.dist(distance_matrix), ngmax=k)
  edges <- array(sapply(kmst_matrix, find_index_row), dim = c(nrow(kmst_matrix), 2))
  edges <- edges[order(edges[,1], edges[,2]), ]
  
  return(edges)
}