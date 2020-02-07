#' Linear mixed model with heteroscedastisity
#'
#' \code{lmmHeterosced} Linear mixed model with the log error variance as linear function of data.
#'
#' @param formula model formula lmer style
#' @param data a data frame with the data
#' @param heterosced_formula an optional heteroscadisity formula
#'
#' @return \code{lmmHeterosced} a list with elements:
#' Fixef (fixed effects estaimtes),
#' Heterosced (parameters of the heteroscedasticity function, a linear function of the log residual variance),
#' Ranef_log_Var (random effect variances on log),
#' Ranef (random effect),
#' fit (nlminb output),
#' obj (MakeADFun ouptut),
#' optTime (Optimization time)
#'
#' @author Geir H. Bolstad
#'
#' @examples
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom lme4 lFormula lmerControl
#' @importFrom Matrix Matrix t
#' @useDynLib lmmHeterosced
#'
#' @export


lmmHeterosced <- function(formula, data, heterosced_formula = ~ 1){

  mf <- lme4::lFormula(formula, data)

  Y <- as.matrix(mf$fr[1])
  X <- Matrix::Matrix(mf$X, sparse = T)
  Z <- Matrix::t(Matrix::Matrix(mf$reTrms$Zt))
  n_random <- length(mf$reTrms$flist)
  n_u <- sapply(mf$reTrms$flist, function(x) length(unique(x)))
  u_start <- c(0, cumsum(n_u))
  u_end <- cumsum(n_u)

  # The model matrix log residual variance
    Q <- model.matrix(heterosced_formula, data = data)


  #----------------
  # Starting values
  #----------------
  param <- list(b = matrix(0, ncol = 1, nrow = ncol(X)),
                u = matrix(0, ncol = 1, nrow = ncol(Z)),
                b_ln_R = matrix(0, ncol = 1, nrow = ncol(Q)),
                ln_G = rep(0, n_random)
                )
  random <- c("b", "u") # the fixed effects are added to get REML estimates

  #------------
  # Data to TMB
  #------------
  dt <- list(Y = Y, X = X, Z = Z,
             n_random = n_random,
             u_start = u_start,
             u_end = u_end,
             Q = Q)


  #--------------
  # Model fitting
  #--------------
  obj <- TMB::MakeADFun(dt, param, DLL = "lmmHeterosced", random=random, silent = TRUE)
  optTime <- system.time(
    fit <- with(obj, nlminb(start=par, objective=fn, gradient=gr))
  )

  #-------
  # Return
  #-------

  results <- summary(TMB::sdreport(obj))

  Fixef <- matrix(results[which(rownames(results)=="b"),], ncol = 2,
                  dimnames=list(dimnames(X)[[2]], c("Estimate", "Std. Error")))
  Heterosced <- matrix(results[which(rownames(results)=="b_ln_R"),], ncol = 2,
                       dimnames=list(dimnames(Q)[[2]], c("Estimate", "Std. Error")))
  Ranef_log_Var <-  matrix(results[which(rownames(results)=="ln_G"),], ncol = 2,
                           dimnames=list(colnames(mf$reTrms$flist), c("Estimate", "Std. Error")))
  Ranef <- matrix(results[which(rownames(results)=="u"),], ncol = 2,
                  dimnames=list(dimnames(Z)[[2]], c("Estimate", "Std. Error")))

  Residuals <- matrix(results[which(rownames(results)=="residuals"),], ncol = 2,
                  dimnames=list(NULL, c("Estimate", "Std. Error")))

  list(Fixef=Fixef, Heterosced=Heterosced, Ranef_log_Var=Ranef_log_Var, Ranef=Ranef, Residuals = Residuals,
       fit = fit, obj = obj, optTime = optTime)

}

