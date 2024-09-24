fisher_lda_posterior <- function(X1, X2, Z_obs, Omega.hat) {
  
  # # Prior probabilities
  # prior1 <- n1 / (n1 + n2)
  # prior2 <- n2 / (n1 + n2)
  prior1 <- prior2 <- 0.5
  
  # Means for each group
  mu1 <- colMeans(X1)
  mu2 <- colMeans(X2)
  
  # Discriminant function
  discriminant_score <- function(Z, mu, prior, Omega.hat) {
    t(mu) %*% Omega.hat %*% t(Z) - (1/2) * t(mu) %*% Omega.hat %*% mu + log(prior)
  }
  
  # Calculate discriminant scores for Z_obs
  # Z_obs should be p*1, not 1*p
  score1 <- discriminant_score(Z_obs, mu1, prior1, Omega.hat)
  score2 <- discriminant_score(Z_obs, mu2, prior2, Omega.hat)
  
  # Posterior probabilities
  max_score <- max(score1, score2)
  post_prob1 <- exp(score1 - max_score) / (exp(score1 - max_score) + exp(score2 - max_score))
  post_prob2 <- exp(score2 - max_score) / (exp(score1 - max_score) + exp(score2 - max_score))
  
  # Prediction
  class_prediction <- ifelse(post_prob1 > post_prob2, 1, 2)
  
  # Return a list containing the posterior probabilities and prediction
  return(list(probability_X1 = as.numeric(post_prob1),
              probability_X2 = as.numeric(post_prob2),
              class = as.numeric(class_prediction)))
}



AdaLDA <- function(xt, yt, ztest, lambda0 = 1) {
  
  n_xt <- nrow(xt)
  n_yt <- nrow(yt)
  n <- min(n_xt, n_yt)
  p <- ncol(xt)
  
  hatmux <- colMeans(xt)
  hatmuy <- colMeans(yt)
  mu <- cbind(hatmux, hatmuy)
  
  hatdelta <- hatmux - hatmuy
  hatSigma <- ((n_xt - 1) * cov(xt) + (n_yt - 1) * cov(yt)) / (n_xt + n_yt - 2)
  
  # Step 1 of the AdaLDA procedure
  d <- diag(hatSigma)
  a <- sqrt(log(p)/n) * sqrt(d)
  B <- lambda0 * sqrt(log(p)/n) * sqrt(d) %*% t(hatdelta)
  f <- rep(1, 2*p)
  
  CoeffM <- rbind(cbind(hatSigma-B, -(hatSigma-B)), cbind(-(hatSigma+B), hatSigma+B))
  # CoeffM <- matrix(c(hatSigma-B, -(hatSigma-B), -(hatSigma+B), hatSigma+B), byrow = T)
  Coeffb <- c(a + hatdelta, a - hatdelta)
  
  uv <- lp(direction = "min", f, CoeffM, "<=", Coeffb, all.bin = FALSE)$solution
  beta0 <- uv[1:p] - uv[(p+1):(2*p)]
  
  # Step 2 of the AdaLDA procedure
  lambda <- sapply(1:p, function(k) {
    sqrt(log(p)/n) * sqrt(lambda0 * hatSigma[k,k] * (abs(sum(beta0 * hatdelta)) + 1))
  })
  
  CoeffMnew <- rbind(cbind(hatSigma, - hatSigma), cbind(- hatSigma, hatSigma))
  Coeffbnew <- c(lambda + hatdelta, lambda - hatdelta)
  
  uv_new <- lp(direction = "min", f, CoeffMnew, "<=", Coeffbnew, all.bin = FALSE)$solution
  beta <- uv_new[1:p] - uv_new[(p+1):(2*p)]
  
  # Evaluation on the testing data
  # IDX <- as.numeric(rowSums(sweep(ztest, 2, colMeans(mu), "-") %*% matrix(beta)) <= 1e-06) + 1
  # IDX <- ifelse(((ztest - matrix(rep(rowMeans(mu), nrow(ztest)), ncol = p, byrow = TRUE)) %*% beta) <= 1e-06, 1, 2)
  # if(length(ztest) != p) {
  #   score <- ((ztest - matrix(1, nrow(ztest), 1) %*% t(as.matrix(colMeans(t(mu))))) %*% as.matrix(beta))
  #   IDX <- as.numeric(score <= 1e-06) + 1
  # } else {
  #   score <- ((t(as.matrix(ztest)) - matrix(1, nrow(t(as.matrix(ztest))), 1) %*% t(as.matrix(colMeans(t(mu)))))
  #             %*% as.matrix(beta))
  #   IDX <- as.numeric(score <= 1e-06) + 1
  # }
  score <- as.numeric(as.matrix(ztest - matrix(1, nrow(ztest), 1) %*% t(as.matrix(rowMeans(mu)))) %*% as.matrix(beta))
  class_prediction <- as.numeric(score <= 1e-06) + 1
  
  # error <- sum(abs(IDX - label_z)) / nrow(ztest)
  
  # Apply Sigmoid function to transform the score into a probability
  probability <- 1 / (1 + exp(score))
  
  # Modify the return statement to include probability
  return(list(class = class_prediction, score = score, probability = probability))
}



AdaGRN <- function(xt, yt, ztest, lambda0 = 1, Sigma.hat) {
  
  n_xt <- nrow(xt)
  n_yt <- nrow(yt)
  n <- min(n_xt, n_yt)
  p <- ncol(xt)
  
  hatmux <- colMeans(xt)
  hatmuy <- colMeans(yt)
  mu <- cbind(hatmux, hatmuy)
  
  hatdelta <- hatmux - hatmuy
  hatSigma <- Sigma.hat
  
  # Step 1 of the AdaLDA procedure
  d <- diag(hatSigma)
  a <- sqrt(log(p)/n) * sqrt(d)
  B <- lambda0 * sqrt(log(p)/n) * sqrt(d) %*% t(hatdelta)
  f <- rep(1, 2*p)
  
  CoeffM <- rbind(cbind(hatSigma-B, -(hatSigma-B)), cbind(-(hatSigma+B), hatSigma+B))
  Coeffb <- c(a + hatdelta, a - hatdelta)
  
  uv <- lp(direction = "min", f, CoeffM, "<=", Coeffb, all.bin = FALSE)$solution
  beta0 <- uv[1:p] - uv[(p+1):(2*p)]
  
  lambda <- sapply(1:p, function(k) {
    sqrt(log(p)/n) * sqrt(lambda0 * hatSigma[k,k] * (abs(sum(beta0 * hatdelta)) + 1))
  })
  
  CoeffMnew <- rbind(cbind(hatSigma, - hatSigma), cbind(- hatSigma, hatSigma))
  Coeffbnew <- c(lambda + hatdelta, lambda - hatdelta)
  
  uv_new <- lp(direction = "min", f, CoeffMnew, "<=", Coeffbnew, all.bin = FALSE)$solution
  beta <- uv_new[1:p] - uv_new[(p+1):(2*p)]
  
  score <- as.numeric(as.matrix(ztest - matrix(1, nrow(ztest), 1) %*% t(as.matrix(rowMeans(mu)))) %*% as.matrix(beta))
  class_prediction <- as.numeric(score <= 1e-06) + 1
  
  # Apply Sigmoid function to transform the score into a probability
  probability <- 1 / (1 + exp(score))
  
  return(list(class = class_prediction, score = score, probability = probability))
}



permute.rows <- function(x){
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = T)
}


balanced.folds <- function(y, nfolds = min(min(table(y)), 10)){
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)
  # makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)
  # nice way to get the ids in a list, split by class
  ###Make a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
  for(i in seq(totals)) {
    bigmat[seq(totals[i]), i] <- sample(yids[[i]])
  }
  smallmat <- matrix(bigmat, nrow = nfolds) # reshape the matrix
  ### Now do a clever sort to mix up the NAs
  smallmat <- permute.rows(t(smallmat)) ### Now a clever unlisting
  x <- apply(smallmat, 2, function(x) x[!is.na(x)])
  if(is.matrix(x)){
    xlist <- list()
    for(i in 1:ncol(x)){
      xlist[[i]] <- x[,i]
    }
    return(xlist)
  }
  return(x)
}



MakeYMat <- function(y){
  return(diag(length(unique(y)))[y,])
}


MakeMeanVecs <- function(x,y){
  Y <- MakeYMat(y)
  return(solve(msqrt(t(Y)%*%Y))%*%t(Y)%*%x)
}

msqrt <- function(mat){
  if(sum((mat-t(mat))^2)>1e-8) stop("Msqrt function only works if mat is symmetric....")
  redo <- TRUE
  while(redo){
    eigenmat <- eigen(mat)
    d <- eigenmat$values
    d[abs(d)<(1e-12)] <- 0
    a <- eigenmat$vectors%*%diag(sqrt(d))%*%t(eigenmat$vectors)
    if(sum(is.na(a))==0) redo <- FALSE
    if(redo) print('did one loop') 
  }
  return(a)
}

soft <- function(mat,lam){
  return(sign(mat)*pmax(abs(mat)-lam, 0))
}


Penalty <- function(v,lambda,type,chrom, lambda2){
  if(type=="standard") return(lambda*sum(abs(v)))
  if(type=="ordered"){
    tots <- lambda*sum(abs(v))
    for(chr in sort(unique(chrom))){
      tots <- tots+lambda2*sum(abs(diff(v[chrom==chr])))
    }
    return(tots)
  }
}


PenalizedPCACrit <- function(x, P, v, lambda, d, type, chrom, lambda2){
  return(t(v)%*%t(x)%*%P%*%x%*%v-d*Penalty(v, lambda, type, chrom, lambda2))
}



PenalizedPCA <- function(x, lambda, K, type="standard",  chrom=NULL, lambda2=NULL, maxiter=30, trace=FALSE){
  # Notice that this K is the number of components desired, NOT the number of classes in the classification problem.
  # Here, x is (# of classes) \times p
  
  # The criterion is maximize_b (b' Sigmabet b) - P(b) s.t. b' Sigmawit b = 1
  # Where Sigmawit=I and where P(b) = lambda||b||_1 or P(b) = lambda||b||_1 + lambda2 ||b_i - b_{i-1}||_1
  # We take a MINORIZATION approach to this problem.
  if(type=="ordered" && is.null(chrom)) chrom <- rep(1, ncol(x))
  if(is.null(lambda2)) lambda2 <- lambda
  crits <-  NULL
  betas <- matrix(0, nrow=ncol(x), ncol=K)
  critslist <- list()
  for(k in 1:K){
    if(trace) cat("Starting on component ", k, fill=TRUE)
    if(k>1){
      svda <- svd(x%*%betas)
      u <- svda$u[,svda$d>(1e-10)]
      P <- diag(nrow(x)) - u%*%t(u)
    }
    if(k==1) P <- diag(nrow(x))
    svdx <- svd(t(x)%*%P)
    d <- svdx$d[1]^2
    beta <- svdx$u[,1]
    crits <- c(crits, PenalizedPCACrit(x, P, beta, lambda, d, type, chrom=chrom, lambda2))
    for(iter in 1:maxiter){
      if((length(crits)<4 || abs(crits[length(crits)]-crits[length(crits)-1])/max(1e-3, crits[length(crits)]) > (1e-6)) && sum(abs(beta))>0){
        if(trace) cat(iter,fill=FALSE)
        tmp <- (t(x)%*%P)%*%(x%*%beta)
        if(type=="standard") beta <- soft(tmp, d*lambda/2)
        if(type=="ordered"){
          for(chr in sort(unique(chrom))){
            beta[chrom==chr] <- as.numeric(flsa(tmp[chrom==chr],  d*lambda/2,  d*lambda2/2))
          }
        }
        beta <- beta/l2n(beta)
        beta[is.na(beta)] <- 0
        crits <- c(crits, PenalizedPCACrit(x, P, beta, lambda, d, type, chrom=chrom, lambda2))
      }
    }
    if(trace) cat(fill=TRUE)
    betas[,k] <- beta#cbind(betas, beta)
    critslist[[k]] <- crits
    if(min(diff(crits))<(-1e-6)) stop("min diff crits is too small!!!")
    crits <- NULL
  }
  return(list(v=betas, crits=as.vector(critslist)))
}



l2n <- function(vec){
  return(sqrt(sum(vec^2)))
}



diag.disc <-function(x, centroids, prior) {
  dd <- t(x) %*% centroids
  dd0 <- (rep(1, nrow(centroids)) %*% (centroids^2))/2 - log(prior)
  scale(dd, as.numeric(dd0), FALSE) # this is -.5*||x_i - mu_k||^2+log(pi_k)
}


softmax <- function(x) {
  exp_x <- exp(x - apply(x, 1, max)) # Subtract max for numerical stability
  probabilities <- exp_x / rowSums(exp_x)
  return(probabilities)
}



Classify <- function(xtr,xte,ytr,equalpriors=FALSE){ # I introduced unequal priors on 02/22/2010
  prior <- rep(1/length(unique(ytr)), length(unique(ytr)))
  if(!equalpriors){             
    for(k in 1:length(unique(ytr))) prior[k] <- mean(ytr==k)
  }
  # classify test obs to nearest training centroid.
  if(is.matrix(xtr) && ncol(xtr)>1){
    mus <- matrix(0, nrow=ncol(xtr), ncol=length(unique(ytr)))
    for(k in 1:length(unique(ytr))){
      mus[,k] <- apply(xtr[ytr==k,], 2, mean)
    }
  } else {
    mus <- matrix(NA, nrow=1, ncol=length(unique(ytr)))
    for(k in 1:length(unique(ytr))) mus[1,k] <- mean(xtr[ytr==k])
  }
  negdists <- diag.disc(t(xte), mus, prior)
  negdists_attrless <- negdists
  attr(negdists_attrless, "scaled:center") <- NULL
  probabilities <- softmax(negdists_attrless)
  return(list(Classification = apply(negdists,1,which.max), Probs = probabilities, Probability = probabilities[,2]))
}  



wcsd <- function(vec, y){
  K <- length(unique(y))
  n <- length(vec)
  tots <- 0
  for(k in unique(y)){
    tots <- tots + sum((vec[y==k]-mean(vec[y==k]))^2)
  }
  return(sqrt(tots/n))
}


wcsd.matrix <- function(x,Y){
  n <- nrow(x)
  return(sqrt((1/n)*apply(((diag(n)-Y%*%diag(1/apply(Y,2,sum))%*%t(Y))%*%data.matrix(x))^2,2,sum)))
  # return(sqrt((1/n)*apply(((diag(n)-Y%*%diag(1/apply(Y,2,sum))%*%t(Y))%*%x)^2,2,sum)))
}



PenalizedLDA_Prob <- function(x, y, xte=NULL,  type="standard",  lambda,  K=2, 
                              chrom=NULL, lambda2=NULL, standardized=FALSE, wcsd.x=NULL, 
                              ymat=NULL, maxiter=20, trace=FALSE){
  if(sum(1:length(unique(y)) != sort(unique(y)))>0) stop("y must be a numeric vector, with values as follows: 1, 2, ....")
  if(sum(is.na(x))>0 || sum(is.na(y))>0 || (!is.null(xte) && sum(is.na(xte))>0)) stop("No missing values allowed!!!")
  if(K>=length(unique(y))) stop("Can have at most K-1 components of K unique classes")
  yclass <- y
  y <- ymat
  if(is.null(ymat)) y <- MakeYMat(yclass)
  if(type=="ordered" && is.null(lambda2)) stop("For type 'ordered', lambda2 must be specified.")
  xorig <- x
  
  if(!standardized){
    if(is.null(wcsd.x)){
      if(length(y)<=200){
        wcsd.x <- wcsd.matrix(x,y)#apply(x, 2, wcsd, yclass)
      } else {
        wcsd.x <- apply(x,2,wcsd,yclass)
      }
      if(min(wcsd.x)==0) stop("Some features have 0 within-class standard deviation.")
    }
    if(!is.null(xte)) xte <- scale(xte, center=apply(x,2,mean), scale=wcsd.x)
    x <- scale(x, T, scale=wcsd.x)
  }
  sqrt.sigma.bet <- t(scale(y, F, sqrt(apply(y, 2, sum))))%*%x/sqrt(nrow(x))
  while(sum(is.na(sqrt.sigma.bet))>0){
    sqrt.sigma.bet <- t(scale(y, F, sqrt(apply(y, 2, sum))))%*%x/sqrt(nrow(x))
    cat("retrying", fill=TRUE)
  }
  penpca <- PenalizedPCA(x=sqrt.sigma.bet, lambda=lambda,  K=K, type=type, chrom=chrom, lambda2=lambda2, maxiter=maxiter, trace=trace)
  Usparse <- penpca$v#matrix(penpca$v, ncol=K)
  if(K==1) Usparse <- matrix(Usparse,ncol=1)
  if(sum(is.na(Usparse))>0){
    Usparse[is.na(Usparse)] <- 0
    #Usparse <- matrix(Usparse,ncol=K)
    if(K==1) Usparse <- matrix(Usparse,ncol=1)
  }
  xtranssparse <- x%*%Usparse# matrix(x%*%Usparse, ncol=K)
  if(K==1) xtranssparse <- matrix(xtranssparse, ncol=1)
  if(!is.null(xte)){
    #xtetranssparse <- matrix(xte%*%Usparse, ncol=K)
    xtetranssparse <- xte%*%Usparse
    if(K==1) xtetranssparse <- matrix(xtetranssparse,ncol=1)
    ypredsparsemat <- yprobsparsemat <- matrix(NA, ncol=K, nrow=nrow(xte))
    for(k in 1:K){
      classify_output  <- Classify(matrix(xtranssparse[,1:k],ncol=k),matrix(xtetranssparse[,1:k],ncol=k),yclass)
      ypredsparsemat[,k] <- classify_output$Classification
      yprobsparsemat[,k] <- classify_output$Probability
    }
    obj <- list(
      ypred = ypredsparsemat,
      yprob = yprobsparsemat, # Add this to store probabilities
      discrim = Usparse,
      xproj = xtranssparse,
      teproj = xtetranssparse,
      K = K,
      crits = penpca$crits,
      type = type,
      lambda = lambda,
      lambda2 = lambda2,
      wcsd.x = wcsd.x,
      x = xorig,
      y = yclass
    )
    class(obj) <- "penlda"
    return(obj)
  } else {
    obj <- list(
      discrim = Usparse,
      xproj = xtranssparse,
      K = K,
      crits = penpca$crits,
      type = type,
      lambda = lambda,
      lambda2 = lambda2,
      wcsd.x = wcsd.x,
      x = xorig,
      y = yclass
    )
    class(obj) <- "penlda"
    return(obj)
  }
}



# Function to calculate ROC AUC for each lambda and select the best one
calculate_best_lambda <- function(result, lambda, true_labels) {
  # Determine the number of rows and lambda values
  n_lambda <- length(lambda)
  
  # Initialize a vector to store the AUC for each lambda
  auc_values <- numeric(n_lambda)
  
  # Iterate over each lambda value to calculate the AUC
  for (j in 1:n_lambda) {
    # Extract the pred_prob for the current lambda
    probs <- sapply(result, function(x) x$pred_prob[j])
    
    # Calculate the ROC AUC
    roc_obj <- roc(true_labels, probs, levels=c(1, 2), direction = "<")
    auc_values[j] <- auc(roc_obj)
  }
  
  # Identify the lambda with the highest AUC
  best_lambda_index <- which.max(auc_values)
  best_lambda <- lambda[best_lambda_index]
  best_auc <- auc_values[best_lambda_index]
  
  # Extract the pred_class and pred_prob for the best lambda
  best_class <- sapply(result, function(x) x$pred_class[best_lambda_index])
  best_prob <- sapply(result, function(x) x$pred_prob[best_lambda_index])
  
  # Create a dataframe with the best results
  result_df <- data.frame(
    pred_class = best_class,
    pred_prob = best_prob
  )
  
  # Return the dataframe, best lambda, and best AUC
  return(list(result_df = result_df, best_lambda = best_lambda, best_auc = best_auc))
}



Pathway_classification <- function(i, Sigma, Omega, X, labels, case_labels = 2, C_type = "LDA", lambda) {
  
  X_test = X[i, , drop = F]
  X_train = X[-i,]
  labels_test = labels[i]
  labels_train = labels[-i]
  
  X1 <- X_train[labels_train == 1,]
  X2 <- X_train[labels_train == 2,]
  
  if (C_type == "LDA") {
    result_updated = fisher_lda_posterior(X1, X2, X_test, Omega)
    pred_class_updated <- result_updated$class
    if (case_labels == 2) {
      prob_updated = result_updated$probability_X2
    } else {
      prob_updated = result_updated$probability_X1
    }
    
    lda_fit <- lda(X_train, grouping = labels_train)
    prediction <- predict(lda_fit, newdata = X_test)
    if (case_labels == 2) {
      prob_classic = as.numeric(prediction$posterior)[2]
    } else {
      prob_classic = as.numeric(prediction$posterior)[1]
    }
    pred_class_classic = prediction$class
    
    return(list(pred_class_updated = as.numeric(pred_class_updated), 
                pred_class_benchmk = pred_class_classic, 
                prob_updated = as.numeric(prob_updated), 
                prob_benchmk = prob_classic))
  } else if (C_type == "AdaLDA") {
    predicted_class = probabilities = vector()
    for (j in 1:length(lambda)) {
      ada_fit <- AdaLDA(X1, X2, X_test, lambda[j])
      predicted_class[j] <- ada_fit$class
      probabilities[j] <- ada_fit$probability
    }
    return(list(pred_class = predicted_class, pred_prob = probabilities))
  } else if (C_type == "AdaGRN") {
    predicted_class = probabilities = vector()
    for (j in 1:length(lambda)) {
      ada_fit <- AdaGRN(X1, X2, X_test, lambda[j], Sigma)
      predicted_class[j] <- ada_fit$class
      probabilities[j] <- ada_fit$probability
    } 
    return(list(pred_class = predicted_class, pred_prob = probabilities))
  } else if (C_type == "PenalizedLDA") {
    predicted_class = probabilities = vector()
    for (j in 1:length(lambda)){
      lda_result <- PenalizedLDA_Prob(x = X_train, y = labels_train, xte = X_test, 
                                      type = "standard", lambda = lambda[j], K = 1)
      predicted_class[j] <- lda_result$ypred
      probabilities[j] <- lda_result$yprob
    }
    return(list(pred_class = predicted_class, pred_prob = probabilities))
  } else if (C_type == "PenalizedGLM") {
    predicted_class = probabilities = vector()
    
    glm_fit <- glmnet(X_train, labels_train, family = "binomial", alpha = 1, lambda = lambda)
    
    probabilities <- predict(glm_fit, newx = X_test, type = "response")
    
    predicted_class <- ifelse(probabilities > 0.5, 2, 1)
    
    return(list(pred_class = predicted_class, pred_prob = probabilities))
  }
}