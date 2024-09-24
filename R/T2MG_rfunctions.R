B_Omega = function(A){
  
  n <- nrow(A)
  B <- matrix(0, n, n)
  Omega <- matrix(0, n, n)
  
  # Find directed edges and set corresponding entries in B
  directed_edges <- which(A != 0 & A != t(A), arr.ind = TRUE)
  B[directed_edges] <- 1
  
  # Find bi-directed edges and set corresponding entries in Omega
  bidirected_edges <- which(A != 0 & A == t(A), arr.ind = TRUE)
  Omega[bidirected_edges] <- 1
  
  # Set diagonal elements to 1 in Omega
  diag(Omega) <- 1
  
  # Return the separated matrices
  return(list(B = B, Omega = Omega))
}



SEM.fit = function(p,ntotal,Data,pop,A,graph.type="MG"){
  if (graph.type == "DAG") {
    if (p == ncol(A)) {
      p.model = p
      # p: the # of regression! not the number of nodes in the gene expression data
      Qest = pmatrix = matrix(0,p.model,p.model) # estimated Q's
      Epsi = matrix(0,ntotal,p.model) # estimated Epsilon matrices
      Rdiag = numeric()
      for (i in 1:p.model)
      {
        non0i = which(A[i,]!=0) # index for the xi's that are included in the i-th regression:
        if (length(non0i)>0) # if there exist directed edges from xj's to this xi
        {
          formulai=formula(paste0(paste(paste0(pop,i), paste(paste0(pop,non0i),collapse="+"), sep="~"), " -1"))
          fiti = lm(formulai,data = Data)
          Qest[i,non0i] = fiti$coefficients
          pmatrix[i,non0i] = summary(fiti)$coef[,'Pr(>|t|)']
          Epsi[,i] = Data[,i] - as.matrix(Data[,non0i])%*%t(t(fiti$coefficients))
          #### assume diagonal:
        }
        if (length(non0i)==0) # if there is no directed edge from any xj to this xi
        {
          Epsi[,i] = Data[,i]
        }
        Rdiag[i] = sum((Epsi[,i]-mean(Epsi[,i]))^2)/(ntotal - 2 -length(non0i))#var(Epsi[,i]) # estimated R's
      }
      #Rinvest = diag(1/Rdiag) # estimated Rinverse's
      return(list(Qest=Qest,Rest=Rdiag))
    }
  }
  if (graph.type == "MG" || missing(graph.type)) {
    B = B_Omega(A)[[1]]
    Omega = B_Omega(A)[[2]]
    Y_new = t(data.matrix(Data))
    out = ricf(B, Omega, Y_new)
    return(list(Qest=out$BHat,Rest=out$OmegaHat,Sest=out$SigmaHat))
  }
}


Pathway_testing = function(X, Y, A, graph_type = "MG") {
  
  ################### Process with the adjacency matrix A ######################
  A <- data.matrix(A)
  diag(A) = 0
  abssum = function(x) sum(abs(x))
  nonzero.A = apply(A,1,abssum)
  d = max(nonzero.A)
  p = nrow(A)
  p0 = sum(nonzero.A > 0)
  sparsity = sum(abs(A) > 0)/p^2
  
  if (graph_type == "DAG") {
    if (identical(A, t(A))) {
      A.dag <- lower.tri(A) * A
    } else {
      a = adj.remove.cycles(A,maxlength=p)
      A.dag = a$adjmat.acyclic
    }
  }
  # n.circles = sum(a$adjmat.removed)
  
  ###################### Process with the data X and Y #########################
  n = list(X=nrow(X),Y=nrow(Y))
  nx = n[['X']]; ny = n[['Y']]
  N = n[['X']]+n[['Y']] ### in BS paper 'n' = N-2
  nmin = min(nx,ny) # minimum sample size
  n.confounder = 0
  
  Z = rbind(scale(X, center = T, scale = F),
            scale(Y, center = T, scale = F))
  Z = as.data.frame(Z)
  colnames(Z) = c(paste0("Z", 1:(p+n.confounder)))#,paste0('M',1:n.confounder))
  meandiff = colMeans(X)-colMeans(Y)
  
  pval=rep(NA,2) # calculated p values
  rej=rep(NA,2) # 1: H0 rejected, 0: H0 not rejected
  test_name_clx = paste("T2", graph_type, "CLX", sep="_")
  test_name_z = paste("T2", graph_type, "Z", sep="_")
  names(pval) = names(rej) = c(test_name_clx, test_name_z)
  

  if (graph_type == "MG" || missing(graph_type)) {
    
    if (missing(graph_type)) warning("graph type not specified; using mixed by default")
    
    ###################### T2MG: ######################
    SEM.results.MG = SEM.fit(p,N,Z,'Z', A[1:(p+n.confounder),1:(p+n.confounder)],"MG")
    Sigma.hat.MG = SEM.results.MG$Sest
    Omega.hat.MG = solve(Sigma.hat.MG)
    T.graph.MG = (n[["X"]]*n[["Y"]])/(N)*(t(colMeans(X)-colMeans(Y))%*% solve(Sigma.hat.MG[1:p,1:p]) %*%t(t(colMeans(X)-colMeans(Y))))
    
    pval[test_name_clx] = pchisq(q=T.graph.MG,df=p,lower.tail=FALSE)
    # rej[test_name_clx] = ifelse(pval[test_name_clx]<alpha,1,0)
    
    T.graph2.MG = (T.graph.MG - p)/sqrt(2*p)
    pval[test_name_z] =  2*pnorm(q=abs(T.graph2.MG),mean=0,sd=1,lower.tail=FALSE)
    # rej[test_name_z] = ifelse(pval[test_name_z]<alpha,1,0)
    
    Sigma.hat = Sigma.hat.MG # For returning result
    Omega.hat = Omega.hat.MG # For returning result

  } else {
    
    ###################### T2DAG: ######################
    SEM.results = SEM.fit(p,N,Z,'Z', A.dag[1:(p+n.confounder),1:(p+n.confounder)],"DAG")
    Qest = SEM.results$Qest
    Rest = SEM.results$Rest
    Sigma.hat = solve(diag(1,p+n.confounder)-Qest)%*%diag(Rest)%*%solve(t(diag(1,p+n.confounder)-Qest))
    Omega.hat = t(diag(1,p+n.confounder)-Qest) %*% diag(1/Rest) %*% (diag(1,p+n.confounder)-Qest)
    T.graph = (n[["X"]]*n[["Y"]])/(N)*(t(colMeans(X)-colMeans(Y))%*% solve(Sigma.hat[1:p,1:p]) %*%t(t(colMeans(X)-colMeans(Y))))
    
    pval[test_name_clx] = pchisq(q=T.graph,df=p,lower.tail=FALSE)
    # rej[test_name_clx] = ifelse(pval[test_name_clx]<alpha,1,0)
    
    T.graph2 = (T.graph - p)/sqrt(2*p)
    pval[test_name_z] =  2*pnorm(q=abs(T.graph2),mean=0,sd=1,lower.tail=FALSE)
    # rej[test_name_z] = ifelse(pval[test_name_z]<alpha,1,0)
  }
  basic.info = c(as.character(round(c(unlist(n),p,d,p0),0)),signif(sparsity,3))
  names(basic.info) = c('n_x', 'n_y','p','d','p0','sparisty')
  # print(paste0('Pathway ', pathway.index, ': ', pathwayID,' Completed'))
  return(list(basic.info=basic.info,pval=pval,Sigma.hat=Sigma.hat,Omega.hat=Omega.hat))
}