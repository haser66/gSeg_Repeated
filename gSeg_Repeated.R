
library(dplyr)


gseg1_repeated = function(n, l, edges, n0=0.05*n, n1=0.95*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100, alpha=1, kappa=1){
  
  r1 = list()
  n0 = ceiling(n0)
  n1 = floor(n1)
  Ebynode = vector("list", n)
  
  for(i in 1:n) Ebynode[[i]]=rep(0,0)
  for(i in 1:nrow(edges)){
    if(edges[i,1]==edges[i,2]){
      Ebynode[[edges[i,1]]] = c(Ebynode[[edges[i,1]]],edges[i,2])
    }else{
      Ebynode[[edges[i,1]]] = c(Ebynode[[edges[i,1]]],edges[i,2])
      Ebynode[[edges[i,2]]] = c(Ebynode[[edges[i,2]]],edges[i,1])
    }
  }
  
  n0_us = n0
  n1_us = n1
  
  if(n0<2){
    n0=2
  }
  if(n1>(n-2)){
    n1=n-2
  }
  
  r1$scanZ = gcp1bynode(n,l,Ebynode,n0,n1,alpha=alpha,kappa=kappa)
  
  if (pval.appr==TRUE){
    mypval1 = pval1(n,l,edges,Ebynode,r1$scanZ,skew.corr,n0,n1)
    r1$pval.appr = mypval1
  }
  if (pval.perm==TRUE){
    mypval2 = permpval1(n,l,Ebynode,r1$scanZ,B,n0,n1,alpha=alpha,kappa=kappa)
    r1$pval.perm = mypval2
  }
  
  
  
  cat("Repeated edge-count statistic: \n")
  if(n0_us<=1){
    cat("  Note: Starting index has been set to n0 = 2 as the repeated edge-count test statistic is not well-defined for t<2. \n")
  }
  if(n1_us>=n-1){
    cat("  Note: Ending index has been set to n1 =", n-2, " as the repeated edge-count test statistic is not well-defined for t>",n-2,". \n")
  }
    
  cat("  Estimated change-point location:", r1$scanZ$rmax$tauhat, "\n")
  cat("  Test statistic (M):", r1$scanZ$rmax$Zmax, "\n")
  if (pval.appr==TRUE){
    cat("  Final Approximated p-value:", r1$pval.appr$pval, "\n")
  }
  if (pval.perm==TRUE){
    cat("  Final p-value from", B, "permutations:", r1$pval.perm$pval, "\n")
  }
  
  
  return(r1)
}

## the Nu function
Nu = function(x){
  y = x/2
  (1/y)*(pnorm(y)-0.5)/(y*pnorm(y) + dnorm(y))
}

## gcp1bynode - single change-point
gcp1bynode = function(n, l, Ebynode, n0=ceiling(0.05*n), n1=floor(0.95*n), alpha=1, kappa=1){
  # "n" is the total number of nodes.
  # "l" is the number of repeated measures.
  # Ebynode[[i]] is the list of nodes that are connect to i by an edge.
  # The nodes are numbered by their order in the sequence.
  
  Dmat <- matrix(0, nrow = n, ncol = n)
  
  for (i in seq_len(n)) {
    neighbors <- Ebynode[[i]]
    for (j in neighbors) {
      Dmat[i, j] <- Dmat[i, j] + 1
    }
  }
  
  Duu <- diag(Dmat)
  Du <- colSums(Dmat) - diag(Dmat)
  Duv <- Dmat * (1 - diag(nrow(Dmat)))
  nodedeg <- colSums(Dmat) + diag(Dmat)
  
  sum_Du2 <- sum(Du^2)
  sum_Du_Duu <- sum(Du*Duu)
  sum_Duv2 <- sum(Duv^2)
  sum_Duu2 <- sum(Duu^2)
  
  G_in <- sum(Duu)
  G_out <- sum(Duv)/2
  
  sumEisq = sum(nodedeg^2)
  nE = sum(nodedeg)/2
  
  
  g = rep(1,n)
  R = rep(0,n)
  R1 = rep(0,n)
  R2 = rep(0,n)
  Ro1 = rep(0,n)
  Ro2 = rep(0,n)
  Ri1 = rep(0,n)
  Ri2 = rep(0,n)
  
  for(i in 1:n){
    g[i] = 0  # update g=group
    links = Ebynode[[i]]
    
    if(i==1){
      if(length(links)>0){
        R[i] = sum(rep(g[i],length(links)) != g[links])
        Ri1[i] = sum(links == i)
      } else {
        R[i] = 0
      }
      R1[i] = Ri1[i]
      Ro1[i] = R1[i] - Ri1[i]
      R2[i] = nE-length(links)
      Ri2[i] = G_in - Ri1[i]
      Ro2[i] = R2[i] - Ri2[i]
    } else {
      if(length(links)>0){
        add = sum(rep(g[i],length(links)) != g[links])
        subin = sum(links == i)
        subtract = length(links)-add-subin
        R[i] = R[i-1]+add-subtract
        Ri1[i] = Ri1[i-1]+subin
        R1[i] = R1[i-1]+subtract+subin
      } else {
        R[i] = R[i-1]
        Ri1[i] = R[i-1]
        R1[i] = R1[i-1]
      }
      Ro1[i] = R1[i]-Ri1[i]
      R2[i] = nE-R[i]-R1[i]
      Ri2[i] = G_in - Ri1[i]
      Ro2[i] = R2[i] - Ri2[i]
    }
    
  }
  
  tt = 1:n
  temp=n0:n1
  
  scanZ = list()
  
  Rw = ((n-tt-1)*R1+(tt-1)*R2)/(n-2)
  mu.Rw = nE*((n-tt-1)*tt*(tt-1)+(tt-1)*(n-tt)*(n-tt-1))/(n*(n-1)*(n-2))
    
  mu.R1 = nE*tt*(tt-1)/(n*(n-1))
  mu.R2 = nE*(n-tt)*(n-tt-1)/(n*(n-1))
    
  v11 = mu.R1*(1-mu.R1) + 2*(0.5*sumEisq-nE)*(tt*(tt-1)*(tt-2))/(n*(n-1)*(n-2)) + (nE*(nE-1)-2*(0.5*sumEisq-nE))*(tt*(tt-1)*(tt-2)*(tt-3))/(n*(n-1)*(n-2)*(n-3))
  v22 = mu.R2*(1-mu.R2) + 2*(0.5*sumEisq-nE)*((n-tt)*(n-tt-1)*(n-tt-2))/(n*(n-1)*(n-2)) + (nE*(nE-1)-2*(0.5*sumEisq-nE))*((n-tt)*(n-tt-1)*(n-tt-2)*(n-tt-3))/(n*(n-1)*(n-2)*(n-3))
    
  v12 = (nE*(nE-1)-2*(0.5*sumEisq-nE))*tt*(n-tt)*(tt-1)*(n-tt-1)/(n*(n-1)*(n-2)*(n-3)) - mu.R1*mu.R2
  var.Rw=((n-tt-1)/(n-2))^2*v11 + 2*((n-tt-1)/(n-2))*((tt-1)/(n-2))*v12+((tt-1)/(n-2))^2*v22
  Zw = -(mu.Rw-Rw)/sqrt(apply(cbind(var.Rw,rep(0,n)),1,max))
    
    
  mu.Ro1 = G_out*(tt*(tt-1))/(n*(n-1))
  mu.Ro2 = G_out*(n-tt)*(n-tt-1)/(n*(n-1))
  mu.Ri1 = G_in*tt/n
    
  V1 = (tt*(tt-1)*(n-tt)*(n-tt-1)/(n*(n-1)*(n-2)*(n-3))) * ((1/2)*sum_Duv2 + ((tt-2)/(n-tt-1))*(sum_Du2-4*((G_out)^2)/n) - 2*((G_out)^2)/(n*(n-1)))
  V2 = (tt*(tt-1)*(n-tt)*(n-tt-1)/(n*(n-1)*(n-2)*(n-3))) * ((1/2)*sum_Duv2 + ((n-tt-2)/(tt-1))*(sum_Du2-4*((G_out)^2)/n) - 2*((G_out)^2)/(n*(n-1)))
  V3 = (tt*(n-tt)/(n*(n-1))) * (sum_Duu2 - ((G_in)^2)/n)
    
  C12 = (tt*(tt-1)*(n-tt)*(n-tt-1)/(n*(n-1)*(n-2)*(n-3))) * ((1/2)*sum_Duv2 - (sum_Du2-4*((G_out)^2)/n) - 2*((G_out)^2)/(n*(n-1)))
  C13 = (tt*(n-tt)*(tt-1)/(n*(n-1)*(n-2))) * (sum_Du_Duu-2*G_in*G_out/n)
  C23 = -(tt*(n-tt)*(n-tt-1)/(n*(n-1)*(n-2))) * (sum_Du_Duu-2*G_in*G_out/n)
    
  var.Row = ((n-tt-1)^2)*V1 + ((tt-1)^2)*V2 + 2*(n-tt-1)*(tt-1)*C12
    
  Zow = ((n-tt-1)*Ro1 + (tt-1)*Ro2 - (n-tt-1)*mu.Ro1 - (tt-1)*mu.Ro2)/sqrt(apply(cbind(var.Row,rep(0,n)),1,max))
  Zod = abs((Ro1-Ro2 - (mu.Ro1-mu.Ro2))/sqrt(apply(cbind(V1 + V2 - 2*C12,rep(0,n)),1,max)))
  Zin = abs((Ri1-mu.Ri1)/sqrt(apply(cbind(V3,rep(0,n)),1,max)))
    
  Mout = apply(cbind(Zod,kappa*Zow),1,max)
  M = apply(cbind(Zin,alpha*Mout),1,max)
  tauhat = temp[which.max(M[n0:n1])]
  rmax = list(tauhat=tauhat, Zmax=M[tauhat], Moutmax=Mout[tauhat], Zowmax=Zow[tauhat], Zodmax=Zod[tauhat], Zinmax=Zin[tauhat],
              Zow=Zow, Zod=Zod, Zin=Zin, M=M, Mout=Mout, alpha=alpha, kappa=kappa)
  scanZ$rmax = rmax
  
  return(scanZ)
}



## rho_one = n h_G
rho_one = function(n, s, sumE, sumEisq){
  f1 = 4*(n-1)*(2*s*(n-s)-n)
  f2 = ((n+1)*(n-2*s)^2-2*n*(n-1))
  f3 = 4*((n-2*s)^2-n)
  f4 = 4*n*(s-1)*(n-1)*(n-s-1)
  f5 = n*(n-1)*((n-2*s)^2-(n-2))
  f6 = 4*((n-2)*(n-2*s)^2-2*s*(n-s)+n)
  n*(n-1)*(f1*sumE + f2*sumEisq - f3*sumE^2)/(2*s*(n-s)*(f4*sumE + f5*sumEisq - f6*sumE^2))
}
rho_one_Rw = function(n, t){
  -((2*t^2 - 2*n*t + n)*(n^2 - 3*n + 2)^4)/(2*t*(n - 1)^3*(n - 2)^4*(t - 1)*(n^2 - 2*n*t - n + t^2 + t))
}


## p-value approximation for single change-point, sub functions
pval1_sub_1 = function(n,b,r,x,lower,upper){
  if (b<1){        #b<0
    return(1)
  }
  theta_b = rep(0,n-1)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  for(i in 1:length(theta_b[pos])){
    if (is.na(theta_b[pos][i])==TRUE){
      theta_b[pos][i]=0
    }
  }
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = x*Nu(sqrt(2*b^2*x)) * ratio
  
  nn.l = ceiling(n/2)-length(which(1+2*r[1:ceiling(n/2)]*b>0))
  nn.r = ceiling(n/2)-length(which(1+2*r[ceiling(n/2):(n-1)]*b>0))
  if (nn.l>0.35*n || nn.r>0.35*n){
    return(0)
  }
  if (nn.l>=lower){
    neg = which(1+2*r[1:ceiling(n/2)]*b<=0)
    dif = c(diff(neg),n/2-nn.l)
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
  }
  if (nn.r>=(n-upper)){
    neg = which(1+2*r[ceiling(n/2):(n-1)]*b<=0 )
    id1 = min(neg+ceiling(n/2)-1,ceiling(n/2)-1)
    id2 = id1 - ceiling(0.03*n)
    id3 = id2 - ceiling(0.09*n)
    inc = (ratio[id3]-ratio[id2])/(id3-id2)
    ratio[id2:(n-1)] = ratio[id2-1]+inc*((id2:(n-1))-id2)
    ratio[ratio<0]=0
    a[(n/2):(n-1)] = (x*Nu(sqrt(2*b^2*x)) * ratio)[(n/2):(n-1)] # update a after extrapolation 
  }
  neg2 = which(a<0)
  a[neg2] = 0
  integrand = function(s){
    a[s]
  }
  result = try(2*dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
  return(result)
  
}
pval1_sub_2 = function(n,b,r,x,lower,upper){
  if (b<1){        #b<0
    return(1)
  }
  theta_b = rep(0,n-1)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = x*Nu(sqrt(2*b^2*x)) * ratio
  a_na = which(is.na(a)==TRUE )
  a[a_na] = 0
  nn = n-1-length(pos)
  if (nn>0.75*n){
    return(0)
  }
  if (nn>=(lower-1)+(n-upper)){
    neg = which(1+2*r*b<=0)
    dif = neg[2:nn]-neg[1:(nn-1)]
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
    a[(n/2+1):n] = a[(n/2):1]
    neg2 = which(a<0 | is.na(a)==TRUE)
    a[neg2] = 0
  }
  integrand = function(s){
    a[s]
  }
  result = try(dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
  return(result)
}


## make edges_cnt from Ebynode
edges_count <- function(Ebynode){
  edge_list <- list()
  
  for (i in seq_along(Ebynode)) {
    for (j in Ebynode[[i]]) {
      if (i != j) {
        a <- min(i, j)
        b <- max(i, j)
        edge_list <- append(edge_list, list(c(a, b)))
      }
    }
  }
  
  edge_mat <- do.call(rbind, edge_list)
  edge_df <- as.data.frame(edge_mat)
  colnames(edge_df) <- c("V1", "V2")
  
  edges_cnt <- edge_df %>%
    group_by(V1, V2) %>%
    summarise(count = n()/2, .groups = "drop") %>%
    as.matrix()
  
  colnames(edges_cnt) <- NULL
  return(edges_cnt)
}


## p value approximation for single change-point
pval1 = function(n, l, edges, Ebynode, scanZ, skew.corr=TRUE, lower=ceiling(0.05*n), upper=floor(0.95*n)){
  
  output = list()
  
  Dmat <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    neighbors <- Ebynode[[i]]
    for (j in neighbors) {
      Dmat[i, j] <- Dmat[i, j] + 1
    }
  }
  
  Duu <- diag(Dmat)
  Du <- colSums(Dmat) - diag(Dmat)
  Duv <- Dmat * (1 - diag(nrow(Dmat)))
  deg <- colSums(Dmat) + diag(Dmat)
  
  sum_Du2 <- sum(Du^2)
  sum_Du_Duu <- sum(Du*Duu)
  sum_Duv2 <- sum(Duv^2)
  sum_Duu2 <- sum(Duu^2)
  
  G_in <- sum(Duu)
  G_out <- sum(Duv)/2
  
  sumEisq = sum(deg^2)
  sumE = sum(deg)/2
  
  
  
  if (skew.corr==FALSE){

    b1 = scanZ$rmax$Zowmax
    if (b1>1){                      
      integrandR1 = function(t){
        x = rho_one_Rw(n,t)
        x*Nu(sqrt(2*b1^2*x))
      }
      pval.repeated1 = dnorm(b1)*b1*integrate(integrandR1, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    }else{
      pval.repeated1 = 1
    }
      
    b2 = scanZ$rmax$Zodmax
    if (b2>1){                      
      integrandR2 = function(t){
        x = n/(2*t*(n - t))
        x*Nu(sqrt(2*b2^2*x))
      }
      pval.repeated2 = 2*dnorm(b2)*b2*integrate(integrandR2, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    }else{
      pval.repeated2 = 1
    }
      
    b = scanZ$rmax$Moutmax
    kappa = scanZ$rmax$kappa
    alpha = scanZ$rmax$alpha
    if (b>1){
      integrand1 = function(t){
        x1 = n/(2*t*(n - t))
        x1*Nu(sqrt(2*b^2*x1))
      }
      integrand2 = function(t){
        x2 = rho_one_Rw(n,t)
        x2*Nu(sqrt(2*b^2*x2))
      }
      pval_u1 = 2*dnorm(b/alpha)*(b/alpha)*integrate(integrand1, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval_u2 = dnorm(b/(kappa*alpha))*(b/(kappa*alpha))*integrate(integrand2, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval.rmout = as.numeric(1-(1-min(pval_u1,1))*(1-min(pval_u2,1)))
    }else{
      pval.rmout = 1
    }
      
      
    b3 = scanZ$rmax$Zinmax
    if (b3>1){                      
      integrandR3 = function(t){
        x = n/(2*t*(n - t))
        x*Nu(sqrt(2*b3^2*x))
      }
      pval.repeated3 = 2*dnorm(b3)*b3*integrate(integrandR3, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    }else{
      pval.repeated3 = 1
    }
      
    pval1 = min(pval.repeated1,1)
    pval2 = min(pval.repeated2,1)
    pval3 = min(pval.repeated3,1)
    pvalMout = min(pval.rmout,1)
      
    pvals = c(pvalMout, pval3)
      
    w <- rep(1/2, 2)
    CCT_T <- sum(w * tan((0.5 - pmin(pvals, 0.9)) * pi))
    pCCT <- 0.5 - (1/pi) * atan(CCT_T)
      
    out.repeated = list(pval1 = pval1, pval2 = pval2, pvalMout=pvalMout, pval3 = pval3, pval = pCCT)
    output = out.repeated
    
    return(output)
  }
  
  
  
  #skewness
  edges_cnt <- edges_count(Ebynode)
  
  xo1 <- sum(Duv^2*(Du-Duv))
  
  xo2=0
  for (i in 1:n){
    selected_rows <- edges_cnt[edges_cnt[,1] == i | edges_cnt[,2] == i, , drop=FALSE]
    values <- selected_rows[, 3]
    if (length(values) >= 3) {
      xo2 = xo2 + 6*sum(apply(combn(values, 3), 2, prod))
    }
  }
  
  xo4 <- sum(Duv*(Du-Duv)*(G_out-Du))
  
  xo3 = 0
  xo5 = 0
  for (i in 1:nrow(edges_cnt)){
    j = edges_cnt[i,1]
    k = edges_cnt[i,2]
    s1 = sum(!(Ebynode[[j]] %in% c(j, k)))
    s2 = sum(!(Ebynode[[k]] %in% c(j, k)))
    xo3 = xo3 + 2*s1*s2*edges_cnt[i,3]
    
    ls <- setdiff(intersect(Ebynode[[j]], Ebynode[[k]]), c(j, k)) # j, k에 대하여, j, k 제외 (j,l), (k,l)이 있는 모든 l값
    row1 = which((edges_cnt[,1] == j & edges_cnt[,2] == k) | (edges_cnt[,1] == k & edges_cnt[,2] == j))
    row2 <- which((edges_cnt[,1] %in% ls & edges_cnt[,2] == k) | (edges_cnt[,1] == k & edges_cnt[,2] %in% ls))
    row3 <- which((edges_cnt[,1] %in% ls & edges_cnt[,2] == j) | (edges_cnt[,1] == j & edges_cnt[,2] %in% ls))
    xo5 = xo5 + 2*sum(edges_cnt[row1, 3]*edges_cnt[row2, 3]*edges_cnt[row3, 3])
  }
  
  xo6 <- sum(Duv^2*(G_out-Du))
  
  xof <- (G_out^3  - ((sum(Duv^3)/2) + 3*xo1 + (3/2)*(xo6-xo1) + xo2 + (3*xo3-3*xo5) + xo5 + (3*xo4+3*xo5-6*xo3)) )
  
  

  t = 1:(n-1)
  Ao1 = (sum(Duv^3)/2)*t*(t-1)/(n*(n-1)) + 3*xo1*t*(t-1)*(t-2)/(n*(n-1)*(n-2)) + (3/2)*(xo6-xo1)*t*(t-1)*(t-2)*(t-3)/(n*(n-1)*(n-2)*(n-3))  + xo2*t*(t-1)*(t-2)*(t-3)/(n*(n-1)*(n-2)*(n-3)) + (3*xo3-3*xo5)*(t*(t-1)*(t-2)*(t-3))/(n*(n-1)*(n-2)*(n-3)) + xo5*(t*(t-1)*(t-2))/(n*(n-1)*(n-2)) + (3*xo4+3*xo5-6*xo3)*t*(t-1)*(t-2)*(t-3)*(t-4)/(n*(n-1)*(n-2)*(n-3)*(n-4)) + xof*t*(t-1)*(t-2)*(t-3)*(t-4)*(t-5)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
    
  Bo1 = (1/2)*(xo6-xo1)*(t*(t-1)*(n-t)*(n-t-1))/(n*(n-1)*(n-2)*(n-3)) + (xo4+xo5-2*xo3)*(t*(t-1)*(t-2)*(n-t)*(n-t-1))/(n*(n-1)*(n-2)*(n-3)*(n-4)) + xof*t*(t-1)*(t-2)*(t-3)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
    
  Co1 = (1/2)*(xo6-xo1)*(n-t)*(n-t-1)*t*(t-1)/(n*(n-1)*(n-2)*(n-3)) + (xo4+xo5-2*xo3)*(n-t)*(n-t-1)*(n-t-2)*t*(t-1)/(n*(n-1)*(n-2)*(n-3)*(n-4)) + xof*t*(t-1)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
    
  Do1 = (sum(Duv^3)/2)*(n-t)*(n-t-1)/(n*(n-1)) + 3*xo1*(n-t)*(n-t-1)*(n-t-2)/(n*(n-1)*(n-2)) + (3/2)*(xo6-xo1)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3))  + xo2*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3)) + (3*xo3 - 3*xo5) *((n-t)*(n-t-1)*(n-t-2)*(n-t-3))/(n*(n-1)*(n-2)*(n-3)) + xo5*((n-t)*(n-t-1)*(n-t-2))/(n*(n-1)*(n-2)) + (3*xo4+3*xo5-6*xo3)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)*(n-t-4)/(n*(n-1)*(n-2)*(n-3)*(n-4)) + xof*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)*(n-t-4)*(n-t-5)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
    
  xi1 = sum(Duu)*sum(Duu^2) - sum(Duu^3)
  Ai1 = sum(Duu^3)*t/n + 3*xi1*t*(t-1)/(n*(n-1)) + (G_in^3 - 3*xi1 - sum(Duu^3))*t*(t-1)*(t-2)/(n*(n-1)*(n-2))
    
    
  E1 = G_out*t*(t-1)/(n*(n-1))
  E2 = G_out*(n-t)*(n-t-1)/(n*(n-1))
  E3 = G_in*t/n
    
  V1 = (t*(t-1)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3))) * ((1/2)*sum_Duv2 + ((t-2)/(n-t-1))*(sum_Du2-4*((G_out)^2)/n) - 2*((G_out)^2)/(n*(n-1)))
  V2 = (t*(t-1)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3))) * ((1/2)*sum_Duv2 + ((n-t-2)/(t-1))*(sum_Du2-4*((G_out)^2)/n) - 2*((G_out)^2)/(n*(n-1)))
  V3 = (t*(n-t)/(n*(n-1))) * (sum_Duu2 - ((G_in)^2)/n)
    
  C12 = (t*(t-1)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3))) * ((1/2)*sum_Duv2 - (sum_Du2-4*((G_out)^2)/n) - 2*((G_out)^2)/(n*(n-1)))
  C13 = (t*(n-t)*(t-1)/(n*(n-1)*(n-2))) * (sum_Du_Duu-2*G_in*G_out/n)
  C23 = -(t*(n-t)*(n-t-1)/(n*(n-1)*(n-2))) * (sum_Du_Duu-2*G_in*G_out/n)
  
  x1 = rho_one_Rw(n,t)
    
  q1=(n-t-1)/(n-2)
  p1=(t-1)/(n-2)
    
  mu1 = q1*E1 + p1*E2
  sig11 = q1^2*V1 + p1^2*V2 + 2*q1*p1*C12
  sig1 = sqrt(sig11)
  ER31 = q1^3*Ao1 + 3*q1^2*p1*Bo1 + 3*q1*p1^2*Co1 + p1^3*Do1
  r1 =  (ER31- 3*mu1*sig1^2 - mu1^3)/sig1^3
  br1 = scanZ$rmax$Zowmax
  result.ur_ow = pval1_sub_2(n,br1,r1,x1,lower,upper)
    
  r.Rw = r1
  x.Rw = x1
    
    
  x2 = n/(2*t*(n - t))
    
  q2=1
  p2=-1
    
  mu2 = q2*E1 + p2*E2
  sig12 = q2^2*V1 + p2^2*V2 + 2*q2*p2*C12
  sig2 = sqrt(sig12)
  ER32 = q2^3*Ao1 + 3*q2^2*p2*Bo1 + 3*q2*p2^2*Co1 + p2^3*Do1
  r2 =  (ER32- 3*mu2*sig2^2 - mu2^3)/sig2^3
  br2 = scanZ$rmax$Zodmax
  result.ur_od = pval1_sub_1(n,br2,r2,x2,lower,upper)
    
  r.Rd = r2
  x.Rd = x2
    
    
  x3 = n/(2*t*(n - t))
  r3 = (Ai1 - 3*E3*V3 - E3^3)/V3^(3/2)
  br3 = scanZ$rmax$Zinmax
  result.ur_in = pval1_sub_1(n,br3,r3,x3,lower,upper)
    
  r.Ri = r3
  x.Ri = x3
    
    
  br = scanZ$rmax$Moutmax
  kappa = scanZ$rmax$kappa
  alpha = scanZ$rmax$alpha
    
  result.ur1 = pval1_sub_1(n,(br/alpha),r2,x2,lower,upper)        #Zod
  result.ur2 = pval1_sub_2(n,(br/(kappa*alpha)),r1,x1,lower,upper)  #Zow
    
    
    
  if (is.numeric(result.ur_ow) && result.ur_ow > 0){
    pval.repeated1 = min(result.ur_ow,1)
  }else{
    if (result.ur_ow ==0){
      cat("Extrapolation for skewness-corrected p-value approximation (out, weighted) could not be performed. \n")
    }
    cat("Repeated edge-count statistic: p-value approximation (out, weighted) without skewness correction is reported.\n")
    b = scanZ$rmax$Zowmax
    if (b>1){
      integrandW = function(t){
        x = rho_one_Rw(n,t)
        x*Nu(sqrt(2*b^2*x))
      }
      pval.repeated1 = dnorm(b)*b*integrate(integrandW, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    }else{
      pval.repeated1 = 1
    }
  }
      
  if (is.numeric(result.ur_od) && result.ur_od > 0){
    pval.repeated2 = min(result.ur_od,1)
  }else{
    if (result.ur_od ==0){
      cat("Extrapolation for skewness-corrected p-value approximation (out, diff) could not be performed. \n")
    }
    cat("Repeated edge-count statistic: p-value approximation (out, diff) without skewness correction is reported.\n")
    b = scanZ$rmax$Zodmax
    if (b>1){
      integrandW = function(t){
        x = n/(2*t*(n - t))
        x*Nu(sqrt(2*b^2*x))
      }
      pval.repeated2 = 2*dnorm(b)*b*integrate(integrandW, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          
    }else{
      pval.repeated2 = 1
    }
  }
      
  if (is.numeric(result.ur1) && is.numeric(result.ur2) && result.ur1 > 0 && result.ur2 > 0){
    pval.rmout = 1-(1-min(result.ur1,1))*(1-min(result.ur2,1))
  }else{
    if(result.ur1 ==0 || result.ur2 == 0){
      cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
    }
    cat("Max-type (out) edge-count statistic: p-value approximation without skewness correction is reported.\n")
    br = scanZ$rmax$Moutmax
    if (br>1){
      integrand1 = function(t){
        x1 = n/(2*t*(n - t))
        x1*Nu(sqrt(2*br^2*x1))
      }
      integrand2 = function(t){
        x2 = rho_one_Rw(n,t)
        x2*Nu(sqrt(2*br^2*x2))
      }
      pval_ur1 = 2*dnorm(br/alpha)*(br/alpha)*integrate(integrand1, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval_ur2 = dnorm(br/(kappa*alpha))*(br/(kappa*alpha))*integrate(integrand2, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval.rmout = as.numeric(1-(1-min(pval_ur1,1))*(1-min(pval_ur2,1)))
    }else{
      pval.rmout = 1
    }
  }
      
  if (is.numeric(result.ur_in) && result.ur_in > 0){
    pval.repeated3 = min(result.ur_in,1)
  }else{
    if (result.ur_in ==0){
      cat("Extrapolation for skewness-corrected p-value approximation (in) could not be performed. \n")
    }
    cat("Repeated edge-count statistic: p-value approximation (in) without skewness correction is reported.\n")
    b = scanZ$rmax$Zinmax
    if (b>1){
      integrandW = function(t){
        x = n/(2*t*(n - t))
        x*Nu(sqrt(2*b^2*x))
      }
      pval.repeated3 = 2*dnorm(b)*b*integrate(integrandW, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    }else{
      pval.repeated3 = 1
    }
  }
      
      
  pval1 = min(pval.repeated1,1)
  pval2 = min(pval.repeated2,1)
  pval3 = min(pval.repeated3,1)
  pvalMout = min(pval.rmout,1)
      
  pvals = c(pvalMout, pval3)
  w <- rep(1/2, 2)
  CCT_T <- sum(w * tan((0.5 - pmin(pvals, 0.9)) * pi))
  pCCT <- 0.5 - (1/pi) * atan(CCT_T)
      
  output = list(pval1 = pval1, pval2 = pval2, pvalMout=pvalMout, pval3 = pval3, pval = pCCT)
      
  
  return(output)
}

## p value from permutation for single change point
permpval1 = function(n, l, Ebynode, scanZ, B, n0=ceiling(0.05*n), n1=floor(0.95*n), alpha=alpha, kappa=kappa){
  
  Z1.repeated = Z2.repeated = Z3.repeated = Zmout.repeated = Zmax.repeated = matrix(0,B,n)
  for(b in 1:B){
    if(b%%1000 ==0) {
      cat(b, "permutations completed.\n")
    }
    perm = sample(n)
    permmatch = rep(0,n)
    for(i in 1:n) permmatch[perm[i]] = i
    Ebnstar =  vector("list", n)
    for(i in 1:n){
      oldlinks = Ebynode[[permmatch[i]]]
      Ebnstar[[i]] = perm[oldlinks]
    }
    
    gcpstar=gcp1bynode(n,l,Ebnstar,n0,n1,alpha=alpha,kappa=kappa)
    
    Z1.repeated[b,] = gcpstar$rmax$Zow
    Z2.repeated[b,] = gcpstar$rmax$Zod
    Z3.repeated[b,] = gcpstar$rmax$Zin
    Zmout.repeated[b,] = gcpstar$rmax$Mout
    Zmax.repeated[b,] = gcpstar$rmax$M
    
  }
  
  output = list()
  p=1-(0:(B-1))/B
  
  
  maxZ1 = apply(Z1.repeated[,n0:n1],1,max)
  maxZ1s = sort(maxZ1)
  maxZ2 = apply(Z2.repeated[,n0:n1],1,max)
  maxZ2s = sort(maxZ2)
  maxZ3 = apply(Z3.repeated[,n0:n1],1,max)
  maxZ3s = sort(maxZ3)
    
  maxZout = apply(Zmout.repeated[,n0:n1],1,max)
  maxZouts = sort(maxZout)
    
    
  pval1 = length(which(maxZ1s>=scanZ$rmax$Zowmax))/B
  pval2 = length(which(maxZ2s>=scanZ$rmax$Zodmax))/B
  pval3 = length(which(maxZ3s>=scanZ$rmax$Zinmax))/B
  pvalMout = length(which(maxZouts>=(scanZ$rmax$Moutmax/scanZ$rmax$alpha)))/B
    
  pvals = c(pvalMout, pval3)
  w <- rep(1/2, 2)
  CCT_T <- sum(w * tan((0.5 - pmin(pvals, 0.9)) * pi))
  pCCT <- 0.5 - (1/pi) * atan(CCT_T)
    
  output = list(pval1=pval1, pval2=pval2, pvalMout=pvalMout, pval3=pval3, pval=pCCT,
                maxZ1s=maxZ1s, maxZ2s=maxZ2s, maxZ3s=maxZ3s, maxZouts=maxZouts,
                Z1=Z1.repeated, Z2=Z2.repeated, Z3=Z3.repeated, Zmout=Zmout.repeated)
  
  
  return(output)
}




