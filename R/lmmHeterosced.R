#' Linear mixed model with heteroscedasticity
#'
#' \code{lmmHeterosced} Linear mixed model with the log error variance as linear function of data.
#'
#' @param formula model formula lmer style
#' @param heterosced_formula an optional heteroscedasticity formula
#' @param SE an optional vector of standard errors of each value in the response variable
#' @param data a data frame with the data
#' @param REML if restricted maximum likelihood should be used
#'
#' @return \code{lmmHeterosced} returns a list with elements:
#' Fixef (fixed effects estimates),  
#' Heterosced (parameters of the heteroscedasticity function, a linear function of the log residual variance),  
#' Ranef_log_Var (random effect variances on log),  
#' Ranef (random effect),
#' Residuals,
#' data (the original data),
#' error_var_matrix (variance matrix of the estimation error)
#' fit (nlminb output),
#' obj (MakeADFun output),
#' optTime (Optimization time)
#' 
#' @details 
#' When providing standard errors in the SE argument, measurement error will be takein into account when 
#' estimating the residual variance and the heteroscedasticity relationship. The measurement error is assumed to 
#' be normally distributed. The vector of standard error need to have the same length and order as the response.
#'
#' @author Geir H. Bolstad
#'
#' @examples
#' # See the vignette 'lmmHeterosced'.
#' @importFrom TMB MakeADFun sdreport 
#' @importFrom lme4 lFormula lmerControl
#' @importFrom Matrix Matrix t
#' @useDynLib lmmHeterosced
#'
#' @export


lmmHeterosced <- function(formula, heterosced_formula = ~ 1, SE = NULL, data, REML = TRUE){

  # Fixed and random effect design matrices
  mf <- try(lme4::lFormula(formula, data), silent = TRUE)
  if(mf[1]=="Error : No random effects terms specified in formula\n"){
    modtype <- "lm"
  }else{
    modtype <- "lmm"
  }
  
  if(modtype == "lm"){
    Y <- as.matrix(model.frame(formula, data = data)[,1])
    X <- Matrix::Matrix(model.matrix(formula, data = data), sparse = T)
    # dummy values to make the code run without too many ifelse:
      Z <- X
      n_random <- 1
      n_u <- ncol(Z)
  }else{ # modtype == "lmm"
    Y <- as.matrix(mf$fr[1])
    X <- Matrix::Matrix(mf$X, sparse = T)
    Z <- Matrix::t(Matrix::Matrix(mf$reTrms$Zt))
    n_random <- length(mf$reTrms$flist)
    n_u <- sapply(mf$reTrms$flist, function(x) length(unique(x)))
  }
  u_start <- c(0, cumsum(n_u))
  u_end <- cumsum(n_u)
  
  # The model matrix for the log residual variance
  Q <- model.matrix(heterosced_formula, data = data)

  # Measurement error
  if(is.null(SE)){
    SE   <- rep(0, nrow(Y))
  } 
    
  #----------------
  # Starting values
  #----------------
    param <- list(b = matrix(0, ncol = 1, nrow = ncol(X)),
                u = matrix(0, ncol = 1, nrow = ncol(Z)),
                b_ln_R = matrix(0, ncol = 1, nrow = ncol(Q)),
                ln_G = rep(0, n_random)#,
                #m = matrix(0, ncol = 1, nrow(Y))
                )
    if(REML){
      #random <- c("b", "u", "m")
      random <- c("b", "u")
    }else{
      #random <- c("u", "m")
      random <- c("u")
    }
    
  #----
  # map
  #----
    if(modtype == "lm"){
      map <- list(u = factor(rep(NA, length(param$u))), ln_G = factor(NA))
    } else{
      map <- list()
    }
    
    if(is.null(SE)){
      #map$m <- factor(rep(NA, length(param$m)))
      #SE   <- rep(1, nrow(Y))
      SE   <- rep(0, nrow(Y))
    } 

  #------------
  # Data to TMB
  #------------
  d <- list(Y = Y, X = X, Z = Z,
             n_random = n_random,
             u_start = u_start,
             u_end = u_end,
             Q = Q,
             SE = SE)


  #--------------
  # Model fitting
  #--------------
  obj <- TMB::MakeADFun(d, param, DLL = "lmmHeterosced", random=random, map = map, silent = TRUE)
  optTime <- system.time(
    fit <- with(obj, nlminb(start=par, objective=fn, gradient=gr))
  )

  #-------
  # Return
  #-------

  SD_report <- TMB::sdreport(obj, getJointPrecision = TRUE)
  results <- summary(SD_report)
  
  Fixef <- matrix(results[which(rownames(results)=="b"),], ncol = 2, 
                  dimnames=list(dimnames(X)[[2]], c("Estimate", "Std. Error")))
  Heterosced <- matrix(results[which(rownames(results)=="b_ln_R"),], ncol = 2, 
                       dimnames=list(dimnames(Q)[[2]], c("Estimate", "Std. Error")))
  if(modtype == "lmm"){
  Ranef_log_Var <-  matrix(results[which(rownames(results)=="ln_G"),], ncol = 2,
                           dimnames=list(colnames(mf$reTrms$flist), c("Estimate", "Std. Error")))
  Ranef <- matrix(results[which(rownames(results)=="u"),], ncol = 2,
                  dimnames=list(dimnames(Z)[[2]], c("Estimate", "Std. Error")))
  }
  Residuals <- obj$report()
  
  
  if(modtype == "lm" & REML == FALSE){
    Fixef_param <- which(colnames(SD_report$cov.fixed)=="b")
    Heterosced_param <- which(colnames(SD_report$cov.fixed)=="b_ln_R")
    error_var <- as.matrix(SD_report$cov.fixed[c(Fixef_param, Heterosced_param), c(Fixef_param, Heterosced_param)])
  }else{
    Fixef_param <- which(colnames(SD_report$jointPrecision)=="b")
    Heterosced_param <- which(colnames(SD_report$jointPrecision)=="b_ln_R")
    error_var <- as.matrix(solve(SD_report$jointPrecision)[c(Fixef_param, Heterosced_param), c(Fixef_param, Heterosced_param)])
  }
  colnames(error_var) <- rownames(error_var) <- c(rownames(Fixef), rownames(Heterosced))
  
  
  logLik <- -fit[["objective"]]
  n = nrow(Y)
  if(modtype =="lm"){
    k = nrow(Fixef)+nrow(Heterosced)
  }else{
    k = nrow(Fixef)+nrow(Heterosced)+nrow(Ranef_log_Var)
  }
  AICc <- 2*k - 2*logLik + 2*k*(k+1)/(n-k-1)

  
  if(modtype =="lm"){
    report <- list(Fixef=Fixef, Heterosced=Heterosced, Residuals = Residuals, AICc = AICc, 
                   logLik = logLik, data = data, response = names(model.frame(formula, data = data))[1],
                   error_var_matrix = error_var, fit = fit, obj = obj, optTime = optTime)
  }else{
    report <- list(Fixef=Fixef, Heterosced=Heterosced, Ranef_log_Var=Ranef_log_Var, 
                   Ranef=Ranef, Residuals = Residuals, AICc = AICc, logLik = logLik, data = data, response = names(mf$fr[1]),
                   error_var_matrix = error_var, fit = fit, obj = obj, optTime = optTime)
  }
  
  
  class(report) = "lmmHeterosced"
  return(report)
}

#' Plot of lmmHeterosced object
#'
#' \code{plot} method for class \code{'lmmHeterosced'}.
#'
#' @param x An object of class \code{'lmmHeterosced'}.
#' @param x_axis The name of the explanatory variable that should be plotted.
#' @param col_data The colour of the datapoints.
#' @param col_expectation The colour of the regression line showing the expectation
#' @param col_sd The colour of the regression lines showing +/- one standard deviation
#' @param xlab as in \code{\link{plot}}.
#' @param ylab as in \code{\link{plot}}.
#' @param col as in \code{\link{plot}}.
#' @param ... Additional arguments passed to \code{\link{plot}}.
#' @details Plots the regression fitted by the \code{\link{lmmHeterosced}}
#'   function. 
#' @return \code{plot} returns a plot of the lmmHeterosced regression.
#' @examples
#' # See the vignette 'lmmHeterosced'.
#' @author Geir H. Bolstad
#' @importFrom graphics plot lines
#' @export
plot.lmmHeterosced <- function(x, x_axis, col_data = "grey", col_expectation = "black", 
                               col_sd = "blue", xlab = NULL, ylab = NULL, ...) {
  y_axis <- x$response
  y <- x$data[,which(names(x$data)==y_axis)]
  xx <- x$data[,which(names(x$data)==x_axis)]
  X <- seq(min(xx), max(xx), length.out = 100)
  a <- x$Fixef[1,1]
  b <- x$Fixef[x_axis,1]
  aa <- x$Heterosced[1,1]
  bb <- x$Heterosced[x_axis,1]
  
  if(is.null(xlab)) xlab <- x_axis
  if(is.null(ylab)) ylab <- y_axis
  
  plot(xx, y, col = col_data, xlab = xlab, ylab = y_axis, ...)
  lines(X, a + b*X, col = col_expectation)
  lines(X, a + b*X + sqrt(exp(aa + bb*X)), col = col_sd)
  lines(X, a + b*X - sqrt(exp(aa + bb*X)), col = col_sd)
}

#' Plot of absolute y values
#'
#' \code{plot_abs_y} plots the expected change in the absolute values of y from an object of class \code{'lmmHeterosced'}.
#'
#' @param mod An object of class \code{'lmmHeterosced'}.
#' @param x The name of the explanatory variable that should be plotted.
#' @param xlab as in \code{\link{plot}}.
#' @param ylab as in \code{\link{plot}}.
#' @param col as in \code{\link{plot}}.
#' @param line_col The colour of the regression line.
#' @param cex.legend A character expansion factor relative to current par("cex")
#'   for the printed parameters.
#' @param ... Additional arguments passed to \code{\link{plot}}.
#' @return A plot.
#' @examples
#' # See the vignette 'Analyzing rates of evolution'.
#' @author Geir H. Bolstad
#' @importFrom graphics plot lines legend
#' @export
#' 
plot_abs_y <- function(mod, x, xlab = NULL, ylab = NULL, col = "grey", line_col = "black", cex.legend = 1, ...){
  fn_mu<-function(m,s){ # expectation of folded normal
    s*sqrt(2/pi)*exp((-1*m^2)/(2*s^2))+m*(1-2*pnorm(-1*m/s,0,1)) 
  }
  y <- mod$response
  yy <- mod$data[,which(names(mod$data)==y)]
  xx <- mod$data[,which(names(mod$data)==x)]
  X <- seq(min(xx), max(xx), length.out = 100)
  keep <- which(colnames(mod$error_var_matrix) %in% c("(Intercept)", x))
  errvarmat <- mod$error_var_matrix[keep, keep]
  V <- eigen(errvarmat)
  eigenv_samples <- sapply(V$values, function(q) rnorm(1000, 0, sqrt(q)))
  dist_param <- t(apply(eigenv_samples, 1, function(q){ # Distribution of parameter values
    q%*%t(V$vectors)+c(mod$Fixef[1,1], mod$Fixef[x,1],mod$Heterosced[1,1],mod$Heterosced[x,1])  
  }
  )
  )
  dist_absy <- t(apply(dist_param, 1, function(B) fn_mu(B[1] + B[2]*X, sqrt(exp(B[3] + B[4]*X)))))
  low_absy <- apply(dist_absy, 2, function(q) quantile(q, 0.025))
  high_absy <- apply(dist_absy, 2, function(q) quantile(q, 0.975))
  expectation_absy <- fn_mu(mod$Fixef[1,1] + mod$Fixef[x,1]*X, sqrt(exp(mod$Heterosced[1,1] + mod$Heterosced[x,1]*X)))
  Avg_change <- coef(lm(fn_mu(mod$Fixef[1,1] + mod$Fixef[x,1]*xx, sqrt(exp(mod$Heterosced[1,1] + mod$Heterosced[x,1]*xx)))~xx))[2]
  Avg_change_dist <- apply(dist_param, 1, function(B) coef(lm(fn_mu(B[1] + B[2]*xx, sqrt(exp(B[3] + B[4]*xx)))~xx))[2])
  result <- c(Avg_change, quantile(Avg_change_dist, 0.025), quantile(Avg_change_dist, 0.975))
  res <- format(signif(result, 2), trim=TRUE, scientific = FALSE)
  if(is.null(xlab)) xlab <- x
  if(is.null(ylab)) ylab <- y
  plot(xx,abs(yy), xlab = xlab, ylab = ylab, col=col, ...)
  lines(X, expectation_absy, col = line_col)
  lines(X, low_absy, lty = "dashed", col = line_col)
  lines(X, high_absy, lty = "dashed", col = line_col)
  legend("topleft", paste0("average change = ", res[1], " (", res[2], ", ", res[3], ")"),
         box.lty = 0, bg = "transparent", xjust = 0, cex = cex.legend)
}


