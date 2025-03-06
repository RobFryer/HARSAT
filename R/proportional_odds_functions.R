# 2_66 (OSPAR 2022)

# ctsm.VDS.cl - tidy up code dealing with Tritia (more needs to be done) and 
#  change factors to characters in output


#' ctsm.VDS.varlist
#' 
#' A list of names for functions and values necessary to export to
#' cluster prcesses for parallel computation
#' 
#' @export
ctsm.VDS.varlist <- c(
  "ctsm.VDS.p.calc", "ctsm.VDS.loglik.calc", "ordinal_theta_est", 
  "ordinal_theta_cl"
)

#' Detects the environment from the call
#' 
#' This is a utility function that detects the package environment.
#' If you have imported `harsat` as a package, it returns the package
#' environment. Otherwise, it returns the global environment. You can safely
#' export functions from the result of this, for example, when you are
#' setting up a cluster of child processes for parallel computation.
#' 
#' @export 
cstm.VDS.environment <- function() {
  environment(sys.function(sys.nframe()))
}

#' ctsm.VDS.p.calc
#' 
#' @param theta vector of values
#' @param cumulate a boolean, whether to use cumulative probabilities
#' @export
ctsm.VDS.p.calc <- function(theta, cumulate = FALSE) {
  if (cumulate) theta <- cumsum(theta)
  cumProb <- c(plogis(theta), 1)
  n <- length(cumProb)
  names(cumProb)[n] <- as.character(n-1)
  c(cumProb[1], diff(cumProb))
}

#' ctsm.VDS.loglik.calc
#' 
#' @param theta The maximum imposex stage
#' @param data Individual imposex data.
#' @param index.theta Allows for stage names that do not run from 0 to `theta` 
#'   (optional) 
#' @param minus.twice A logical specifying whehter to calculated the liklihood
#'   (FALSE) or the deviance (TRUE); default FALSE
#' @param cumulate A logical specifying whether to use cumulative probabilities; 
#'   default FALSE
#' @export
ctsm.VDS.loglik.calc <- function(
  theta, data, index.theta, minus.twice = FALSE, cumulate = FALSE) {

  vds <- with(data, table(indexID, VDS))
  
  if (missing(index.theta)) {
    index.theta <- rep(0, nrow(vds))
    names(index.theta) <- row.names(vds)
  }
  
  out <- sapply(row.names(vds), function(x) {
    theta[1] <- theta[1] + index.theta[x]
    dmultinom(c(vds[x,]), prob = ctsm.VDS.p.calc(theta, cumulate), log = TRUE)
  })
  if (all(is.finite(out))) 
    out <- sum(out)
  else 
    out <- -1e7
  
  if (minus.twice) out <- - 2 * out
  
  out
}

#' Estimates cut-points in ordinal model (imposex assessments)
#'
#' @param data individual ordinal data (e.g. for imposex assessments)
#' @param theta an (optional) set of intial parameter values for theta; defaults
#'   to NULL
#' @param ref_level an (optional) reference level for parameter estimation; 
#'   defaults to a index with intermediate levels of imposex 
#' @param calc_vcov logical specifying whether to calculate the covariance 
#'   matrix of the parameter estimates; defaults to FALSE
#'   
#' @export
ordinal_theta_est <- function(
    data, theta = NULL, ref_level = NULL, calc_vcov = FALSE) {

  # silence non-standard evaluation warnings
  .data <- NULL


  # ensure indexID is a character - remove when top level function is sorted
  # also ensure VDS is integer valued
  
  data$indexID <- as.character(data$indexID)
  data$VDS <- as.integer(data$VDS)
  
  # basic data checks
  
  counts <- table(data$VDS)
  n_theta <- max(data$VDS)
  
  if (!identical(names(counts), as.character(0:n_theta))) {
    stop("some categories have zero counts: contact harsat development team")
  }
  
  
  # initialise theta
  
  if (!is.null(theta)) {
    if (length(theta) != n_theta) {
      stop("length of theta is not compatible with the data")
    }
    if (any(theta[2:n_theta] <= 0)) {
      stop("theta cannot have negative values (apart from its first value)")
    }
  }
  
  if (is.null(theta)) {
    theta <- cumsum(counts) / sum(counts)
    theta <- theta[as.character(0:(n_theta - 1))]
    theta <- qlogis(theta)
    theta <- c(theta[1], diff(theta))
  }
  

  # deal with indices 

  index_all <- sort(unique(data$indexID))
  
  
  # identify indices where all values are zero - uninformative so exclude from 
  # the estimation problem

  index_est <- tapply(data$VDS, data$indexID, max)
  
  index_est <- names(index_est)[index_est > 0]

    
    
  
  if (length(index_est) == 0L) {
    stop("no informative data - cannot proceed")
  }
  
  data <- dplyr::filter(data, indexID %in% index_est)
  
  n_index <- length(index_est)

  

  # initialise ref_level
  
  if (!is.null(ref_level)) {
    if (!ref_level %in% names(index_all)) {
      stop("ref_level not found in data")
    }
    if (!ref_level %in% index_est) {
      warning("ref_level has no informative data - an alternative will be chosen")
    }
  }

  if (is.null(ref_level)) {
    ref_choose <- tapply(data$VDS, data$indexID, mean)
    ref_level <- names(sort(ref_choose))[n_index %/% 2]
  }
  
  
  # add in a single category 1 observation for each index that is all zeros
  # and a single category n - 1 observation for each index that is all maxed out
  
  wk.zero <- with(data, tapply(VDS, indexID, function(x) all(x == 0)))
  wk.id <- with(data, indexID %in% names(wk.zero)[wk.zero] & !duplicated(indexID))
  data[wk.id, "VDS"] <- 1
  
  wk.max <- with(data, tapply(VDS, indexID, function(x) all(x == n_theta)))
  wk.id <- with(data, indexID %in% names(wk.max)[wk.max] & !duplicated(indexID))
  data[wk.id, "VDS"] <- n_theta - 1

  
  # get starting values - need to exlude the reference level
  
  theta_index <- rep(0, n_index - 1)
  names(theta_index) <- setdiff(index_est, ref_level)

  in_par <- c(theta, theta_index)
  
  wk_optim <- function(par, data, n_theta) {
    
    cutID <- as.character(0:(n_theta - 1))
    
    theta <- par[cutID]
    
    par <- par[setdiff(names(par), cutID)]
    
    index.theta <- rep(0, dplyr::n_distinct(data$indexID))
    names(index.theta) <- sort(unique(data$indexID))
    index.theta[names(par)] <- par
    
    ctsm.VDS.loglik.calc(theta, data, index.theta, minus.twice = TRUE, cumulate = TRUE)
  }
  
  out <- optim(
    in_par, 
    wk_optim, 
    data = data, 
    n_theta = n_theta, 
    method = "L-BFGS-B", 
    lower = c(-Inf, rep(0, n_theta - 1), rep(-Inf, length(in_par) - n_theta)), 
    control = list(trace = 1, maxit = 500, REPORT = 10), 
    hessian = calc_vcov
  )
  
  
  out$ref_level <- ref_level
  
  if (calc_vcov) {

    out$vcov <- 2 * solve(out$hessian)
    
    out$summary <- data.frame(
      est = out$par, 
      se = sqrt(diag(out$vcov))
    )
    
    out$summary <- dplyr::mutate(
      out$summary, 
      t = .data$est / .data$se,
      p = round(2 * pnorm(abs(.data$t), lower.tail = FALSE), 4)
    )
    
    out$summary <- dplyr::relocate(out$summary, "p", .after = "t")
  }
  
  out$K <- n_theta

  out$index <- index_all
    
  out
}


#' Calculates confidence limits for imposex time series
#' 
#' @param fit The output from a call to `ordinal_theta_est` (sort of)
#' @param nsim The number of simulations on which each set of confidence limits
#' is based; default 1000
#'
#' @export
ordinal_theta_cl <- function(fit, nsim = 1000, annual_indices = NULL) {
    
  nCuts <- fit$K
  cutsID <- as.character(0:(nCuts-1))
  categories <- 0:nCuts
  
  indexID <- setdiff(names(fit$par), cutsID)

  # ad-hoc treatment for indices where all the individuals are zero
  
  extra <- setdiff(c(fit$index), c(indexID, fit$ref_level))
  
  if (length(extra) >= 1L){
    if (is.null(annual_indices)) {
      stop(
        "annual_indices must be supplied if there are zero indices\n", 
        "this is a temporary measure until the problem is fixed properly\n",
        "please contact the harsat development team", 
        call. = FALSE)
    }
    
    fit <- ordinal_update_fit(fit, indexID, extra, annual_indices)
    
    indexID <- c(indexID, extra)
  }
  

  set.seed(fit$seed)
  
  data <- MASS::mvrnorm(nsim, fit$par, fit$vcov)

  data.cuts <- data[, cutsID, drop = FALSE]
  if (nCuts > 1) data.cuts <- t(apply(data.cuts, 1, cumsum))

  data.index <- matrix(
    0, 
    nrow = nsim, 
    ncol = length(indexID) + 1, 
    dimnames = list(NULL, c(indexID, fit$ref_level))
  )
  
  data.index[, indexID] <- data[, indexID]
  data.index <- as.data.frame(data.index)

  cl <- sapply(data.index, FUN = function(i) {
    out <- data.cuts + i
    out <- sort(apply(out, 1, function(x) sum(ctsm.VDS.p.calc(x) * categories)))
    out <- out[round(nsim * c(0.05, 0.95))]
  }, simplify = FALSE)

  cl <- data.frame(do.call("rbind", cl))
  names(cl) <- c("lower", "upper")

  
  # ad-hoc treatment for indices where all the individuals are zero
  
  # extra <- setdiff(fit$index, row.names(cl))
  # 
  # if (length(extra > 0L)) {
  #   warning("ad-hoc estimation of upper cl for zero indices")
  #   extra <- data.frame(index_id = extra, lower = 0, upper = 0.5)
  #   extra <- tibble::column_to_rownames(extra, "index_id")
  #   cl <- dplyr::bind_rows(cl, extra)
  # }
  
  if (length(extra) >= 1L) {
    cl[extra, "lower"] <- 0
  }
  

  n_tail = 2L
  if (any(grepl("Tritia nitida (reticulata)", row.names(cl), fixed = TRUE))) {
    # warning("ad-hoc fix for Tritia nitida (reticulata)")
    n_tail = 3L
  }
  
  namesID <- strsplit(row.names(cl), " ", fixed = TRUE)
  cl$species <- sapply(namesID, function(x) paste(tail(x, n_tail), collapse = " "))
  cl$year <- as.numeric(sapply(namesID, function(x) x[length(x) - n_tail]))
  cl$station_code <- sapply(
    namesID, 
    function(x) paste(head(x, length(x) - n_tail - 1), collapse = " ")
  )

  cl
}


ordinal_update_fit <- function(fit, indexID, extraID, annual_indices) {
  
  if (!all(c(indexID, extraID) %in% annual_indices$indexID)) {
    stop("missing indices")
  }

  annual_indices <- tibble::column_to_rownames(annual_indices, "indexID")
    
  positive <- annual_indices[indexID, ]
  zero <- annual_indices[extraID, ]
  
  if (!all(zero$annual_index < 0.0001)) {
    stop("non-zero indices present where they shouldn't be")
  }

  positive$par <- fit$par[indexID]
  positive$var <- diag(fit$vcov[indexID, indexID])
  
  # predict par for zero index based on linear model on square root scale
  
  # lattice::xyplot(par ~ annual_index^0.5, data = positive, pch = 16)
  
  pred_model <- lm(par ~ sqrt(annual_index), data = positive)
  
  zero$par <- predict(pred_model, zero)


  # predict var for zero index based on upper quantile, since get relatively
  # less precise when based on many zeros

  # lattice::xyplot(I(sqrt(var * n)) ~ annual_index^0.5, data = positive, pch = 16)

  pred_sd <- quantile(sqrt(positive$var * positive$n), p = 0.9)

  zero$var <- pred_sd * pred_sd / zero$n
  
  
  # update fit$par and fit$vcov
  
  wk_n <- nrow(fit$vcov) + nrow(zero)
  wk_names <- c(names(fit$par), row.names(zero))
  
  par <- rep(0, wk_n)
  names(par) <- wk_names
  
  vcov <- matrix(
    0, nrow = wk_n, ncol = wk_n, dimnames = list(wk_names, wk_names)
  )
  
  par[names(fit$par)] <- fit$par
  par[row.names(zero)] <- zero$par
  
  vcov[names(fit$par), names(fit$par)] <- fit$vcov
  for (i in row.names(zero)) {
    vcov[i, i] <- zero[i, "var"]
  }

  fit$par <- par
  fit$vcov <- vcov

  return(fit)
}





# finds station year species combinations for which both individual and pooled 
# data have been submitted and checks for consistency

# ctsm.VDS.check <- function(ctsmOb) {
#   
#   data <- droplevels(subset(ctsmOb$data, determinand %in% determinands$Biota$imposex))
#   data <- data[c("seriesID", "station", "year", "species", "determinand", "concentration", "n_individual", 
#                  "%FEMALEPOP")]
#   
#   data[c("country", "region")] <- 
#     ctsmOb$stations[as.character(data$station), ][c("country", "region")]
#   data <- droplevels(data)
#   
#   data <- within(data, {
#     indexID <- factor(paste(station, determinand, species, year))
#     stopifnot(round(n_individual) == n_individual)
#     n_individual <- as.integer(n_individual)
#   })
#   
#   
#     # identify points with both individual and pooled data
#   
#   mixedData <- with(data, tapply(n_individual, indexID, function(x) any(x == 1) & any(x > 1)))
#   mixedData <- names(mixedData)[mixedData]
#   data <- droplevels(subset(data, indexID %in% mixedData))
#   
#   
#   # calculate VDSI based on individuals
#   
#   splitID <- factor(data$n_individual == 1, levels = c(FALSE, TRUE), labels = c("index", "stage"))
#   data <- split(data, splitID)
#   
#   stage <- with(data$stage, aggregate(concentration, by = list(indexID = indexID), function(x) 
#     c("SIndex" = mean(x), "SSum" = sum(x), "SN" = length(x))))
#   stage <- with(stage, cbind(indexID, as.data.frame(x)))
#   
#   data <- merge(data$index, stage, all = TRUE)
#   
#   names(data)[match("concentration", names(data))] <- "Index"
#   
#   data <- data[c("country", "station", "year", "species", "determinand", "n_individual", "%FEMALEPOP", 
#                  "Index", "SIndex", "SN", "SSum")]
#   names(data)[match("%FEMALEPOP", names(data))] <- "fprop"
#   
#   data <- within(data, {
#     okValue <- abs(Index - SIndex) <= 0.01  
#     okN <- round(n_individual * fprop / 100) == SN
#   })
#   
#   return(data)
# }



