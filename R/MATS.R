#' Multi-Ancestry TwaS
#'
#' This function implements Multi-Ancestry TWAS (MATS), based on the framework outlined in Knutson & Pan, 2022
#' @param Y a vector of either continuous or binary trait values for n subjects.
#' @param xhat a numeric vector of imputed gene expression for n subjects.
#' @param groups a vector giving ancestral group labels for the n subjects. Must include at least 2 groups.
#' @param C an optional matrix of covariates for the n subjects (n rows). If any of the covariates are categorical, specify their columnn names in the categorical.vars parameter. Dp not include a column for the intercept.
#' @param P a matrix of genetic principal component on the n subjects (n rows). Must have at least nP PCs (at least nP columns). Columns should be orders PC1, PC2, ... and should not include any other columns aside from PCs (i.e., no subject ID column).
#' @param ev a numeric vector of ordered eigenvalues for the genetic principal component. Should be the same length as the number of columns in P.
#' @param K an n x n genetic covariance  matrix.
#' @param categorical.vars a character vector of categorical variable column names in C (i.e., c("variable1", "variable2"))
#' @param trait_type specificy either "binary" or "continuous" for the trait (Y) type.
#' @param nP the number of principal components to include in all models to adjust for population stratification, and the number used in the full PC-interaction model.
#' @param n The sample size used to estimate the disease GWAS summary data
#' @param trait_type Either "Continuous" or "Binary", for continous or binary disease traits
#' @param n_case The number of cases used in the disease GWAS sample. Required for trait_type = "Binary", default for Continuous trait is NULL
#' @param n_control The number of controls used in the disease GWAS sample. Required for trait_type = "Binary", default for Continuous trait is NULL
#' @import dplyr
#' @import lmtest
#' @import bigQF
#' @import fastDummies
#' @import Rcpp
#' @import RcppArmadillo
#' @export
#' @examples
#' MATS()

MATS <- function(Y, xhat, groups, C, P, ev, K, categorical.vars = NULL, trait_type = "binary", np = 5){

  if(sum(c(missing(Y), missing(xhat), missing(groups), missing(P), missing(ev), missing(K))) > 0){
    cat("At least one of the arguments Y, xhat, group, P, ev, or K is missing \n")
  }else if(sum(is.na(Y)) != 0){
    cat("Y cannot contain NA values.", "\n")
  }else if(sum(is.na(xhat)) != 0){
    cat("xhat cannot contain NA values.", "\n")
  }else if(sum(is.na(P)) != 0){
    cat("P cannot contain NA values.", "\n")
  }else if(!(trait_type %in% c("continuous", "binary"))){
    cat("Trait Type must be specified as continuous or binary \n")
  }else if(class(K) != "matrix"){
    cat("Correlation Matrix must be a matrix", "\n")
  }else if(is.vector(ev) & class(ev) != "numeric" & length(ev) < nP){
    cat("ev must be a numeric vector of length >= nP", "\n")
  }else if(!is.vector(Y)){
    cat("Y must be a vector", "\n")
  }else if(!is.vector(groups)){
    cat("groups must be a vector", "\n")
  }else if(!is.vector(xhat)  & class(xhat) != "numeric"){
    cat("xhat must be a numeric vector", "\n")
  }else if(length(as.factor(levels(groups))) > 1){
    cat("You must include more than one ancestral group!", "\n")
  }else if(ncol(P) < nP){
    cat("nP must be <= ncol(P)", "\n")
  }else if(var(c(length(Y), length(xhat), nrow(P), nrow(K), ncol(K), length(groups))) != 0){
    cat("length of Y, xhat, P, and the dimensions of K (rows and columns) must be equal.", "\n")
  }

  cat("Did you confirm the the subject order of Y, xhat, groups, K (cols and rows), P, and C are all the same? If not, results will be misleading!")

  names(P) <- paste0("P", 1:ncol(P))
  gmmat.data  <- data.frame(trait = Y, C, G = groups, P[,1:nP])
  dummy_groups <- (as.data.frame(fastDummies::dummy_cols(data.frame(G = groups), select_columns = 'G')))[,-c(1:2)]

  for(var in categorical.vars){
    x <- as.data.frame(fastDummies::dummy_cols(data.frame(as.factor(C[,var])))[,-c(1:2)]); names(x) <- paste0(var, ".", levels(C[,var])[-1])
    dummy_groups <- cbind(dummy_groups, x)
  }

  null.mod.C <- apply(list(intercept = rep(1, nrow(gmmat.data)), gmmat.data %>% dplyr::select(-G, -categorical.vars), dummy_groups) %>% do.call("cbind", .), 2, as.numeric)

  if(trait_type == "binary"){
    mod0 <- glm(trait ~ 1 + ., data = gmmat.data, family = binomial(link = "logit"))
  }else{
    mod0 <- glm(trait ~ 1 + ., data = gmmat.data, family = gaussian(link = "identity"))
  }

  n <- dim(C)[1] #number of samples
  I <- diag(1, n) #identity matrix of dimension n
  mu.mod0 <- mod0$fitted.value #mu tilde under H0: alpha0 = 0 and tau = 0
  if(trait_type == "binary"){
    mod0.var <- mu.mod0 * (1 - mu.mod0) #estimated variances under H0: alpha0 = 0 and tau = 0
  }else{
    mod0.var <- rep((1/n)*sum((Y - mu.mod0)^2), n) #estimated variances under H0: alpha0 = 0 and tau = 0
  }
  mod0.residuals <- Y - mu.mod0 #residuals under H0: alpha0 = 0 and tau = 0
  mod0.var <- diag(mod0.var) #estimated diagonal covariance matrix of Y under H0: alpha = 0 and tau = 0
  tCD0C.alpha <- t(null.mod.C) %*% mod0.var %*% null.mod.C
  inv.tCD0C.alpha <- solve_mm(tCD0C.alpha)$Imp
  D0C <- arma_mm(mod0.var, null.mod.C)
  Delta.alpha1 <- mod0.var - (arma_mm(arma_mm(D0C, inv.tCD0C.alpha), t(D0C)))

  E <- matrix(xhat, ncol = 1); colnames(E) = "E"
  M <- data.frame(gmmat.data, E = E[,1])
  Ediag <- diag(nrow(E)); diag(Ediag) <- E #set up matrix for GLMM version

  if(trait_type == "binary"){
    null.mod.pop0 <- glm(trait ~ 1 + ., data = gmmat.data %>% mutate(E = E), family = binomial(link = "logit"))
    null.mod.pop1 <- glm(trait ~ 1 + . + E*groups, data = gmmat.data %>% mutate(E = E), family = binomial(link = "logit"))
  }else{
    null.mod.pop0 <- glm(trait ~ 1 + ., data = gmmat.data %>% mutate(E = E), family = gaussian(link = "identity"))
    null.mod.pop1 <- glm(trait ~ 1 + . + E*groups, data = gmmat.data %>% mutate(E = E), family = gaussian(link = "identity"))
  }

  pop.index.total <- c(grep(":E", names(null.mod.pop1$coefficients)), which( names(null.mod.pop1$coefficients) == "E"))
  parameter.vector.total <- matrix(null.mod.pop1$coefficients[pop.index.total], ncol = 1)
  V.pop.total <- (vcov(null.mod.pop1))[pop.index.total,pop.index.total]

  wald.total <- t(parameter.vector.total) %*% solve(V.pop.total) %*% parameter.vector.total
  pop.p.total <- pchisq(wald.total, df = length(parameter.vector.total), lower.tail = FALSE)
  pop.het.total <-  data.frame(pop.het.wald.total = wald.total, pop.het.df.total = length(parameter.vector.total), pop.het.p.total = pop.p.total)

  pop.lik <- lmtest::lrtest(null.mod.pop0, null.mod.pop1)
  pop.het.LRT <-  data.frame(pop.het.chisq.LRT = pop.lik$Chisq, pop.het.df.LRT = pop.lik$Df, pop.het.p.LRT = pop.lik$`Pr(>Chisq)`)[2, ]

  pop.index <- grep("E:", names(null.mod.pop1$coefficients))
  parameter.vector <- matrix(null.mod.pop1$coefficients[pop.index], ncol = 1)
  V.pop <- (vcov(null.mod.pop1))[pop.index,pop.index]

  wald <- t(parameter.vector) %*% solve(V.pop) %*% parameter.vector
  pop.p <- pchisq(wald, df = length(parameter.vector), lower.tail = FALSE)
  pop.het <-  data.frame(pop.het.wald = wald, pop.het.df = length(parameter.vector), pop.het.p = pop.p)

  X.shared <- as.matrix(cbind(E, E*dummy_groups[,3:ncol(dummy_groups)]))
  Delta.alpha <- arma_mm(arma_mm(t(E), Delta.alpha1), E)
  S.alpha <- arma_mm(t(E ), matrix((Y - mu.mod0), ncol = 1))
  U.alpha <- t(S.alpha^2)/Delta.alpha
  p.U.alpha <- 1 - stats::pchisq(U.alpha, df = dim(E)[2])

  mu.0.tau <- null.mod.pop1$fitted.value #mu hat under H0: tau = 0
  d.0.tau <- mu.0.tau * (1 - mu.0.tau) # estimated variances under H0: tau = 0
  D.0.tau <- diag(d.0.tau) #estimated diagonal covariance matrix of Y under H0: tau = 0 , n x
  res.0.tau <- resid <- Y - mu.0.tau #residuals under H0: tau = 0
  U.tau <- (1/n)*arma_mm(arma_mm(matrix(E*resid, nrow = 1), K), matrix(E*resid, ncol = 1))
  WE <- diag(E*sqrt(diag(D.0.tau)))
  tmp <- sweep(K, MARGIN=1, E*sqrt(diag(D.0.tau)), `*`)

  Mat <- (1/n)*sweep(tmp, MARGIN=2, E*sqrt(diag(D.0.tau)), `*`)
  p.S.tau <- bigQF::pQF(U.tau, Mat, neig=100,convolution.method="integration", method = "satterthwaite")

  q.fisher <- -2*( log( p.S.tau ) + log( pop.p.total ) )
  MATS.p.overall <- 1-pchisq( q.fisher, df=4 )

  #full model with PC interactions
  full.mod.data <- cbind(gmmat.data %>% mutate(E = E), P[, 1:nP])

  if(trait_type == "binary"){
    reduced_PCadj.mod <- glm(trait ~ 1 + . + E*groups, data = data.frame(full.mod.data), family = binomial(link = "logit"))
  }else{
    reduced_PCadj.mod <- glm(trait ~ 1 + . + E*groups, data = data.frame(full.mod.data), family = gaussian(link = "identity"))
  }


  w.inters <- list()
  for(col in paste0("P", 1:nP)){
    w.inters[[length(w.inters) + 1]] <- full.mod.data[,"E"] * full.mod.data[,col] * sqrt(ev[as.numeric(gsub("P", "", col))])
  }
  w.inters <- w.inters %>% do.call("cbind", .); colnames(w.inters) <- paste0("w.inter", 1:nP)
  w.full.mod <- glm(trait ~ 1 + . + E*groups, data = data.frame(full.mod.data,w.inters), family = binomial(link = "logit"))

  w.lrt <- lmtest::lrtest(reduced_PCadj.mod, w.full.mod)

  if(trait_type == "binary"){
    twas <- stats::glm(formula = Y ~ ., data = data.frame(Y, C, E, P[,1:nP]), family = binomial(link = "logit"))
  }else{
    twas <- stats::glm(formula = Y ~ ., data = data.frame(Y, C, E, P[,1:nP]), family = gaussian(link = "identity"))
  }

  final.MATS.summary = data.frame(overall.expression.p = MATS.p.overall, indiv.het.S = U.tau, indiv.het.p = p.S.tau, pop.het.S = wald, pop.het.p = pop.p)

  if(MATS.p.overall > 0.05){
    recommended.model = "There is no expression effect."
  }else if(p.S.tau < 0.05 & pop.p < 0.05){
      recommended.model = "There is evidence of individual and population-level heterogeneity. Use the full model with PC-interactions."
  }else if(p.S.tau < 0.05 & pop.p > 0.05){
      recommended.model = "There is evidence of individual-level heterogeneity. Use the full model with PC-interactions to assess whether population-level heterogeneity is also present."
  }else if(p.S.tau > 0.05 & pop.p < 0.05){
      recommended.model = "There is evidence of population-level heterogeneity but not individual-level heterogeneity. Use the reduced model with population-specific expression effects."
  }else if(p.S.tau > 0.05 & pop.p > 0.05){
      recommended.model = "There is no evidence heterogeneity. Use standard TWAS."
  }

  cat("Done! \n")
  cat(recommended.model, "\n")

  return(list(test.summary = final.MATS.summary, recommended.model = recommended.model, full.model.pc.interactions = list(model = w.full.mod, lrt.pc.interaction.effect = w.lrt), reduced.model.pop.specific.effects = null.mod.pop1, twas = twas))

}

