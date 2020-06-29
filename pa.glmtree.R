require(partykit)

pa.glmfit <- function (y, x, start = NULL, weights = NULL, offset = NULL, 
                       cluster = NULL, ..., estfun = FALSE, object = FALSE, 
                       caseweights = TRUE, glm.weights=weights) 
{
  args <- list(...)
  ctrl <- list()
  for (n in c("epsilon", "maxit")) {
    if (n %in% names(args)) {
      ctrl[[n]] <- args[[n]]
      args[[n]] <- NULL
    }
  }
  args$control <- do.call("glm.control", ctrl)
  if (is.null(x)) 
    x <- matrix(1, nrow = NROW(y), ncol = 1L, dimnames = list(NULL, 
                                                              "(Intercept)"))
  args <- c(list(x = x, y = y, start = start, weights = glm.weights, 
                 offset = offset), args)
  z <- do.call("glm.fit", args)
  df <- z$rank
  if (z$family$family %in% c("gaussian", "Gamma", "inverse.gaussian")) 
    df <- df + 1
  if (substr(z$family$family, 1L, 5L) != "quasi") 
    objfun <- z$aic/2 - df
  else objfun <- z$deviance
  rval <- list(coefficients = z$coefficients, objfun = objfun, 
               estfun = NULL, object = NULL)
  if (estfun) {
    wres <- as.vector(z$residuals) * z$weights
    dispersion <- if (substr(z$family$family, 1L, 17L) %in% 
                      c("poisson", "binomial", "Negative Binomial")) {
      1
    }
    else {
      if (!is.null(glm.weights) && caseweights) {
        sum(wres^2/glm.weights, na.rm = TRUE)/sum(z$weights, 
                                              na.rm = TRUE)
      }
      else {
        sum(wres^2, na.rm = TRUE)/sum(z$weights, na.rm = TRUE)
      }
    }
    rval$estfun <- wres * x[, !is.na(z$coefficients), drop = FALSE]/dispersion
  }
  if (object) {
    class(z) <- c("glm", "lm")
    z$offset <- if (is.null(offset)) 
      0
    else offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- attr(x, "xlevels")
    cl <- as.call(expression(glm))
    cl$formula <- attr(x, "formula")
    if (!is.null(offset)) 
      cl$offset <- attr(x, "offset")
    z$call <- cl
    z$terms <- attr(x, "terms")
    if (!is.null(glm.weights) && caseweights) {
      z$df.null <- z$df.null - sum(glm.weights > 0) + sum(glm.weights)
      z$df.residual <- z$df.residual - sum(glm.weights > 0) + 
        sum(glm.weights)
    }
    rval$object <- z
  }
  return(rval)
}

pa.glmtree <- function (formula, data, subset, na.action, weights, offset, 
          cluster, family = gaussian, epsilon = 1e-08, maxit = 25, 
          glm.weights=weights,
          ...) 
{
  control <- mob_control(...)
  cl <- match.call(expand.dots = TRUE)
  f <- Formula::Formula(formula)
  if (length(f)[2L] == 1L) {
    attr(f, "rhs") <- c(list(1), attr(f, "rhs"))
    formula[[3L]] <- formula(f)[[3L]]
  }
  else {
    f <- NULL
  }
  if (inherits(family, "family")) {
    fam <- TRUE
  }
  else {
    fam <- FALSE
    if (is.character(family)) 
      family <- get(family)
    if (is.function(family)) 
      family <- family()
  }
  glmfit0 <- function(y, x, start = NULL, weights = NULL, offset = NULL, 
                      cluster = NULL, ..., estfun = FALSE, object = FALSE, 
                      caseweights = TRUE, glm.weights=weights) {
    pa.glmfit(y = y, x = x, start = start, weights = weights, 
           offset = offset, cluster = cluster, ..., estfun = estfun, 
           object = object, caseweights = control$caseweights,
           glm.weights=glm.weights)
  }
  m <- match.call(expand.dots = FALSE)
  if (!is.null(f)) 
    m$formula <- formula
  m$fit <- glmfit0
  m$control <- control
  m$epsilon <- epsilon
  m$maxit <- maxit
  if ("..." %in% names(m)) 
    m[["..."]] <- NULL
  if (!fam) 
    m$family <- family
  #m[[1L]] <- as.call(quote(partykit::mob))
  m[[1L]] <- pa.mob
  rval <- eval(m, parent.frame())
  rval$info$call <- cl
  rval$info$family <- family$family
  class(rval) <- c("glmtree", class(rval))
  return(rval)
}

pa.mob <- function (formula, data, subset, na.action, weights, offset, 
          cluster, fit, control = mob_control(), glm.weights=weights, ...) 
{
  fitargs <- names(formals(fit))
  if (!all(c("y", "x", "start", "weights", "offset") %in% fitargs)) {
    stop("no suitable fitting function specified")
  }
  if (!all(c("estfun", "object") %in% fitargs)) {
    afit <- function(y, x = NULL, start = NULL, weights = NULL, 
                     offset = NULL, cluster = NULL, ..., estfun = FALSE, 
                     object = FALSE) {
      obj <- if ("cluster" %in% fitargs) {
        fit(y = y, x = x, start = start, weights = weights, 
            offset = offset, cluster = cluster, ...)
      }
      else {
        fit(y = y, x = x, start = start, weights = weights, 
            offset = offset, ...)
      }
      list(coefficients = coef(obj), objfun = -as.numeric(logLik(obj)), 
           estfun = if (estfun) sandwich::estfun(obj) else NULL, 
           object = if (object) obj else NULL)
    }
  }
  else {
    if ("cluster" %in% fitargs) {
      afit <- fit
    }
    else {
      afit <- function(y, x = NULL, start = NULL, weights = NULL, 
                       offset = NULL, cluster = NULL, ..., estfun = FALSE, 
                       object = FALSE) {
        fit(y = y, x = x, start = start, weights = weights, 
            offset = offset, ..., estfun = estfun, object = object)
      }
    }
  }
  cl <- match.call()
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", 
               "offset", "cluster"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if (length(formula)[2L] < 2L) {
    formula <- Formula::Formula(formula(Formula::as.Formula(formula(formula), 
                                                            ~0), rhs = 2L:1L))
    xreg <- FALSE
  }
  else {
    if (length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1L:2L))
      warning("Formula must not have more than two RHS parts")
    }
    xreg <- TRUE
  }
  mf$formula <- formula
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- terms(formula, data = data)
  mtY <- terms(formula, data = data, rhs = if (xreg) 
    1L
    else 0L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- switch(control$ytype, vector = Formula::model.part(formula, 
                                                          mf, lhs = 1L)[[1L]], matrix = model.matrix(~0 + ., Formula::model.part(formula, 
                                                                                                                                 mf, lhs = 1L)), data.frame = Formula::model.part(formula, 
                                                                                                                                                                                  mf, lhs = 1L))
  X <- if (!xreg) 
    NULL
  else switch(control$xtype, matrix = model.matrix(mtY, mf), 
              data.frame = Formula::model.part(formula, mf, rhs = 1L))
  if (!is.null(X) && ncol(X) < 1L) {
    X <- NULL
    xreg <- FALSE
  }
  if (xreg) {
    attr(X, "formula") <- formula(formula, rhs = 1L)
    attr(X, "terms") <- mtY
    attr(X, "offset") <- cl$offset
  }
  Z <- Formula::model.part(formula, mf, rhs = 2L)
  n <- nrow(Z)
  nyx <- length(mf) - length(Z) - as.numeric("(weights)" %in% 
                                               names(mf)) - as.numeric("(offset)" %in% names(mf)) - 
    as.numeric("(cluster)" %in% names(mf))
  varindex <- match(names(Z), names(mf))
  weights <- model.weights(mf)
  if (is.null(weights)) 
    weights <- 1L
  if (length(weights) == 1L) 
    weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  offset <- if (xreg) 
    model.offset(mf)
  else NULL
  cluster <- mf[["(cluster)"]]
  if (!is.null(control$prune)) {
    if (is.character(control$prune)) {
      control$prune <- tolower(control$prune)
      control$prune <- match.arg(control$prune, c("aic", 
                                                  "bic", "none"))
      control$prune <- switch(control$prune, aic = {
        function(objfun, df, nobs) (2 * objfun[1L] + 
                                      2 * df[1L]) < (2 * objfun[2L] + 2 * df[2L])
      }, bic = {
        function(objfun, df, nobs) (2 * objfun[1L] + 
                                      log(n) * df[1L]) < (2 * objfun[2L] + log(n) * 
                                                            df[2L])
      }, none = {
        NULL
      })
    }
    if (!is.function(control$prune)) {
      warning("Unknown specification of 'prune'")
      control$prune <- NULL
    }
  }
  nodes <- pa.mob_partynode(Y = Y, X = X, Z = Z, weights = weights, 
                         offset = offset, cluster = cluster, fit = afit, control = control, 
                         varindex = varindex, ..., glm.weights = glm.weights)
  fitted <- fitted_node(nodes, data = mf)
  fitted <- data.frame(`(fitted)` = fitted, check.names = FALSE, 
                       row.names = rownames(mf))
  if (!identical(weights, rep.int(1L, n))) 
    fitted[["(weights)"]] <- weights
  if (!is.null(offset)) 
    fitted[["(offset)"]] <- offset
  if (!is.null(cluster)) 
    fitted[["(cluster)"]] <- cluster
  rval <- party(nodes, data = if (control$model) 
    mf
    else mf[0, ], fitted = fitted, terms = mt, info = list(call = cl, 
                                                           formula = oformula, Formula = formula, terms = list(response = mtY, 
                                                                                                               partitioning = mtZ), fit = afit, control = control, 
                                                           dots = list(...), nreg = max(0L, as.integer(xreg) * (nyx - 
                                                                                                                  NCOL(Y)))))
  class(rval) <- c("modelparty", class(rval))
  return(rval)
}

pa.mob_partynode <- function (Y, X, Z, weights = NULL, offset = NULL, cluster = NULL, 
                              fit, control = mob_control(), varindex = 1L:NCOL(Z), glm.weights=weights,...) 
{
  if (missing(X)) 
    X <- NULL
  xreg <- !is.null(X)
  n <- nrow(Z)
  if (is.null(weights)) 
    weights <- 1L
  if (length(weights) < n) 
    weights <- rep(weights, length.out = n)
  if(is.null(glm.weights))
    glm.weights <- weights
  minsize <- control$minsize
  if (!is.null(minsize) && !is.integer(minsize)) 
    minsize <- as.integer(minsize)
  verbose <- control$verbose
  rnam <- c("estfun", "object")
  terminal <- lapply(rnam, function(x) x %in% control$terminal)
  inner <- lapply(rnam, function(x) x %in% control$inner)
  names(terminal) <- names(inner) <- rnam
  w2n <- function(w) if (control$caseweights) 
    sum(w)
  else sum(w > 0)
  suby <- function(y, index) {
    if (control$ytype == "vector") 
      y[index]
    else y[index, , drop = FALSE]
  }
  subx <- if (xreg) {
    function(x, index) {
      sx <- x[index, , drop = FALSE]
      attr(sx, "contrasts") <- attr(x, "contrasts")
      attr(sx, "xlevels") <- attr(x, "xlevels")
      attr(sx, "formula") <- attr(x, "formula")
      attr(sx, "terms") <- attr(x, "terms")
      attr(sx, "offset") <- attr(x, "offset")
      sx
    }
  }
  else {
    function(x, index) NULL
  }
  subz <- function(z, index) z[index, , drop = FALSE]
  root.matrix <- function(X) {
    if ((ncol(X) == 1L) && (nrow(X) == 1L)) 
      return(sqrt(X))
    else {
      X.eigen <- eigen(X, symmetric = TRUE)
      if (any(X.eigen$values < 0)) 
        stop("Matrix is not positive semidefinite")
      sqomega <- sqrt(diag(X.eigen$values))
      V <- X.eigen$vectors
      return(V %*% sqomega %*% t(V))
    }
  }
  mob_grow_fluctests <- function(estfun, z, weights, obj = NULL, 
                                 cluster = NULL) {
    m <- NCOL(z)
    pval <- rep.int(NA_real_, m)
    stat <- rep.int(0, m)
    ifac <- rep.int(FALSE, m)
    mtest <- if (m <= control$mtry) 
      1L:m
    else sort(sample(1L:m, control$mtry))
    process <- as.matrix(estfun)
    ww0 <- (weights > 0)
    process <- process[ww0, , drop = FALSE]
    z <- z[ww0, , drop = FALSE]
    k <- NCOL(process)
    n <- NROW(process)
    nobs <- if (control$caseweights && any(weights != 1L)) 
      sum(weights)
    else n
    process <- process/sqrt(nobs)
    vcov <- control$vcov
    if (is.null(obj)) 
      vcov <- "opg"
    if (vcov != "opg") {
      bread <- vcov(obj) * nobs
    }
    if (vcov != "info") {
      meat <- if (is.null(cluster)) {
        crossprod(if (control$caseweights) 
          process/sqrt(weights)
          else process)
      }
      else {
        crossprod(as.matrix(apply(if (control$caseweights) 
          process/sqrt(weights)
          else process, 2L, tapply, as.numeric(cluster), 
          sum)))
      }
    }
    J12 <- root.matrix(switch(vcov, opg = chol2inv(chol(meat)), 
                              info = bread, sandwich = bread %*% meat %*% bread))
    process <- t(J12 %*% t(process))
    if (!is.null(control$parm)) 
      process <- process[, control$parm, drop = FALSE]
    k <- NCOL(process)
    from <- if (control$trim > 1) 
      control$trim
    else ceiling(nobs * control$trim)
    from <- max(from, minsize)
    to <- nobs - from
    lambda <- ((nobs - from) * to)/(from * (nobs - to))
    beta <- partykit:::mob_beta_suplm
    logp.supLM <- function(x, k, lambda) {
      if (k > 40L) {
        logp_estrella2003 <- function(x, k, lambda) -lgamma(k/2) + 
          k/2 * log(x/2) - x/2 + log(abs(log(lambda) * 
                                           (1 - k/x) + 2/x))
        p <- ifelse(x <= 1.5 * k, (x/(1.5 * k))^sqrt(k) * 
                      logp_estrella2003(1.5 * k, k, lambda), logp_estrella2003(x, 
                                                                               k, lambda))
      }
      else {
        nb <- ncol(beta) - 1L
        tau <- if (lambda < 1) 
          lambda
        else 1/(1 + sqrt(lambda))
        beta <- beta[(((k - 1) * 25 + 1):(k * 25)), ]
        dummy <- beta[, (1L:nb)] %*% x^(0:(nb - 1))
        dummy <- dummy * (dummy > 0)
        pp <- pchisq(dummy, beta[, (nb + 1)], lower.tail = FALSE, 
                     log.p = TRUE)
        if (tau == 0.5) {
          p <- pchisq(x, k, lower.tail = FALSE, log.p = TRUE)
        }
        else if (tau <= 0.01) {
          p <- pp[25L]
        }
        else if (tau >= 0.49) {
          p <- log((exp(log(0.5 - tau) + pp[1L]) + exp(log(tau - 
                                                             0.49) + pchisq(x, k, lower.tail = FALSE, 
                                                                            log.p = TRUE))) * 100)
          if (!is.finite(p)) 
            p <- mean(c(pp[1L], pchisq(x, k, lower.tail = FALSE, 
                                       log.p = TRUE)))
        }
        else {
          taua <- (0.51 - tau) * 50
          tau1 <- floor(taua)
          p <- log(exp(log(tau1 + 1 - taua) + pp[tau1]) + 
                     exp(log(taua - tau1) + pp[tau1 + 1L]))
          if (!is.finite(p)) 
            p <- mean(pp[tau1 + 0L:1L])
        }
      }
      return(as.vector(p))
    }
    for (i in mtest) {
      zi <- z[, i]
      if (length(unique(zi)) < 2L) 
        next
      if (is.factor(zi)) {
        oi <- order(zi)
        proci <- process[oi, , drop = FALSE]
        ifac[i] <- TRUE
        iord <- is.ordered(zi) & (control$ordinal != 
                                    "chisq")
        zi <- zi[oi]
        zi <- factor(zi, levels = unique(zi))
        segweights <- if (control$caseweights) 
          tapply(weights[oi], zi, sum)
        else table(zi)
        segweights <- as.vector(segweights)/nobs
        if (length(segweights) < 2L) {
          stat[i] <- 0
          pval[i] <- NA_real_
        }
        else if (iord) {
          proci <- apply(proci, 2L, cumsum)
          tt0 <- head(cumsum(table(zi)), -1L)
          tt <- head(cumsum(segweights), -1L)
          if (control$ordinal == "max") {
            stat[i] <- max(abs(proci[tt0, ]/sqrt(tt * 
                                                   (1 - tt))))
            pval[i] <- log(as.numeric(1 - mvtnorm::pmvnorm(lower = -stat[i], 
                                                           upper = stat[i], mean = rep(0, length(tt)), 
                                                           sigma = outer(tt, tt, function(x, y) sqrt(pmin(x, 
                                                                                                          y) * (1 - pmax(x, y))/((pmax(x, y) * 
                                                                                                                                    (1 - pmin(x, y)))))))^k))
          }
          else {
            proci <- rowSums(proci^2)
            stat[i] <- max(proci[tt0]/(tt * (1 - tt)))
            pval[i] <- log(strucchange::ordL2BB(segweights, 
                                                nproc = k, nrep = control$nrep)$computePval(stat[i], 
                                                                                            nproc = k))
          }
        }
        else {
          stat[i] <- sum(sapply(1L:k, function(j) (tapply(proci[, 
                                                                j], zi, sum)^2)/segweights))
          pval[i] <- pchisq(stat[i], k * (length(levels(zi)) - 
                                            1), log.p = TRUE, lower.tail = FALSE)
        }
      }
      else {
        oi <- if (control$breakties) {
          mm <- sort(unique(zi))
          mm <- ifelse(length(mm) > 1L, min(diff(mm))/10, 
                       1)
          order(zi + runif(length(zi), min = -mm, max = +mm))
        }
        else {
          order(zi)
        }
        proci <- process[oi, , drop = FALSE]
        proci <- apply(proci, 2L, cumsum)
        tt0 <- if (control$caseweights && any(weights != 
                                              1L)) 
          cumsum(weights[oi])
        else 1:n
        stat[i] <- if (from < to) {
          xx <- rowSums(proci^2)
          xx <- xx[tt0 >= from & tt0 <= to]
          tt <- tt0[tt0 >= from & tt0 <= to]/nobs
          max(xx/(tt * (1 - tt)))
        }
        else {
          0
        }
        pval[i] <- if (from < to) 
          logp.supLM(stat[i], k, lambda)
        else NA
      }
    }
    best <- which.min(pval)
    if (length(best) < 1L) 
      best <- NA
    rval <- list(pval = exp(pval), stat = stat, best = best)
    names(rval$pval) <- names(z)
    names(rval$stat) <- names(z)
    if (!all(is.na(rval$best))) 
      names(rval$best) <- names(z)[rval$best]
    return(rval)
  }
  mob_grow_findsplit <- function(y, x, zselect, weights, offset, 
                                 cluster,glm.weights, ...) {
    if (minsize > 0.5 & minsize < 1) 
      minsize <- 1 - minsize
    if (minsize < 0.5) 
      minsize <- ceiling(w2n(weights) * minsize)
    if (is.numeric(zselect)) {
      uz <- sort(unique(zselect))
      if (length(uz) == 0L) 
        stop("Cannot find admissible split point in partitioning variable")
      if (control$restart) {
        get_dev <- function(i) {
          zs <- zselect <= uz[i]
          if (w2n(weights[zs]) < minsize || w2n(weights[!zs]) < 
              minsize) {
            return(Inf)
          }
          else {
            fit_left <- fit(y = suby(y, zs), x = subx(x, 
                                                      zs), start = NULL, weights = weights[zs], 
                            offset = offset[zs], cluster = cluster[zs],
                            ..., glm.weights=glm.weights[zs])
            fit_right <- fit(y = suby(y, !zs), x = subx(x, 
                                                        !zs), start = NULL, weights = weights[!zs], 
                             offset = offset[!zs], cluster = cluster[!zs],
                             ..., glm.weights=glm.weights[!zs])
            return(fit_left$objfun + fit_right$objfun)
          }
        }
        dev <- unlist(control$applyfun(1L:length(uz), 
                                       get_dev))
      }
      else {
        dev <- vector(mode = "numeric", length = length(uz))
        start_left <- NULL
        start_right <- NULL
        for (i in 1L:length(uz)) {
          zs <- zselect <= uz[i]
          if (control$restart || !identical(names(start_left), 
                                            names(start_right)) || !identical(length(start_left), 
                                                                              length(start_right))) {
            start_left <- NULL
            start_right <- NULL
          }
          if (w2n(weights[zs]) < minsize || w2n(weights[!zs]) < 
              minsize) {
            dev[i] <- Inf
          }
          else {
            fit_left <- fit(y = suby(y, zs), x = subx(x, 
                                                      zs), start = start_left, weights = weights[zs], 
                            offset = offset[zs], cluster = cluster[zs], glm.weights=glm.weights[zs],
                            ...)
            fit_right <- fit(y = suby(y, !zs), x = subx(x, 
                                                        !zs), start = start_right, weights = weights[!zs], 
                             offset = offset[!zs], cluster = cluster[!zs], glm.weights=glm.weights[!zs],
                             ...)
            start_left <- fit_left$coefficients
            start_right <- fit_right$coefficients
            dev[i] <- fit_left$objfun + fit_right$objfun
          }
        }
      }
      if (all(!is.finite(dev))) {
        split <- list(breaks = NULL, index = NULL)
      }
      else {
        split <- list(breaks = if (control$numsplit == 
                                   "center") {
          as.double(mean(uz[which.min(dev) + 0L:1L]))
        } else {
          as.double(uz[which.min(dev)])
        }, index = NULL)
      }
    }
    else {
      if (!is.ordered(zselect) & control$catsplit == "multiway") {
        return(list(breaks = NULL, index = seq_along(levels(zselect))))
      }
      olevels <- levels(zselect)
      zselect <- factor(zselect)
      al <- partykit:::mob_grow_getlevels(zselect)
      get_dev <- function(i) {
        w <- al[i, ]
        zs <- zselect %in% levels(zselect)[w]
        if (w2n(weights[zs]) < minsize || w2n(weights[!zs]) < 
            minsize) {
          return(Inf)
        }
        else {
          if (nrow(al) == 1L) 
            1
          else {
            fit_left <- fit(y = suby(y, zs), x = subx(x, 
                                                      zs), start = NULL, weights = weights[zs], 
                            offset = offset[zs], cluster = cluster[zs],glm.weights=glm.weights[zs], 
                            ...)
            fit_right <- fit(y = suby(y, !zs), x = subx(x, 
                                                        !zs), start = NULL, weights = weights[!zs], 
                             offset = offset[!zs], cluster = cluster[zs], glm.weights=glm.weights[!zs],
                             ...)
            fit_left$objfun + fit_right$objfun
          }
        }
      }
      dev <- unlist(control$applyfun(1L:nrow(al), get_dev))
      if (all(!is.finite(dev))) {
        split <- list(breaks = NULL, index = NULL)
      }
      else {
        if (is.ordered(zselect)) {
          split <- list(breaks = match(levels(zselect)[which.min(dev)], 
                                       olevels), index = NULL)
        }
        else {
          ix <- structure(rep.int(NA_integer_, length(olevels)), 
                          .Names = olevels)
          ix[colnames(al)] <- !al[which.min(dev), ]
          ix <- as.integer(ix) + 1L
          split <- list(breaks = NULL, index = ix)
        }
      }
    }
    return(split)
  }
  mob_grow <- function(id = 1L, y, x, z, weights, offset, cluster, glm.weights,
                       ...) {
    if (verbose) {
      if (id == 1L) 
        cat("\n")
      cat(sprintf("-- Node %i %s\n", id, paste(rep("-", 
                                                   32 - floor(log10(id)) + 1L), collapse = "")))
      cat(sprintf("Number of observations: %s\n", w2n(weights)))
    }
    mod <- fit(y, x, weights = weights, offset = offset, 
               cluster = cluster, ..., estfun = TRUE, object = terminal$object | 
                 control$vcov == "info",glm.weights=glm.weights)
    mod$test <- NULL
    mod$nobs <- w2n(weights)
    mod$p.value <- NULL
    if (is.null(minsize)) 
      minsize <<- as.integer(ceiling(10L * length(mod$coefficients)/NCOL(y)))
    TERMINAL <- FALSE
    if (w2n(weights) < 2 * minsize) {
      if (verbose) 
        cat(sprintf("Too few observations, stop splitting (minsize = %s)\n\n", 
                    minsize))
      TERMINAL <- TRUE
    }
    if (depth >= control$maxdepth) {
      if (verbose) 
        cat(sprintf("Maximum depth reached, stop splitting (maxdepth = %s)\n\n", 
                    control$maxdepth))
      TERMINAL <- TRUE
    }
    if (TERMINAL) {
      return(partynode(id = id, info = mod))
    }
    test <- if (is.null(mod$estfun)) 
      NULL
    else try(mob_grow_fluctests(mod$estfun, z, weights, mod$object, 
                                cluster))
    if (!is.null(test) && !inherits(test, "try-error")) {
      if (control$bonferroni) {
        pval1 <- pmin(1, sum(!is.na(test$pval)) * test$pval)
        pval2 <- 1 - (1 - test$pval)^sum(!is.na(test$pval))
        test$pval <- ifelse(!is.na(test$pval) & (test$pval > 
                                                   0.001), pval2, pval1)
      }
      best <- test$best
      TERMINAL <- is.na(best) || test$pval[best] > control$alpha
      mod$p.value <- as.numeric(test$pval[best])
      if (verbose) {
        cat("\nParameter instability tests:\n")
        print(rbind(statistic = test$stat, p.value = test$pval))
        cat(sprintf("\nBest splitting variable: %s", 
                    names(test$stat)[best]))
        cat(sprintf("\nPerform split? %s", ifelse(TERMINAL, 
                                                  "no\n\n", "yes\n")))
      }
    }
    else {
      if (verbose && inherits(test, "try-error")) 
        cat("Parameter instability tests failed\n\n")
      TERMINAL <- TRUE
      test <- list(stat = NA, pval = NA)
    }
    mod$test <- rbind(statistic = test$stat, p.value = test$pval)
    if (TERMINAL) {
      return(partynode(id = id, info = mod))
    }
    else {
      zselect <- z[[best]]
      sp <- mob_grow_findsplit(y, x, zselect, weights, 
                               offset, cluster,glm.weights, ...)
      if (is.null(sp$breaks) & is.null(sp$index)) {
        if (verbose) 
          cat(sprintf("No admissable split found in %s\n\n", 
                      sQuote(names(test$stat)[best])))
        return(partynode(id = id, info = mod))
      }
      else {
        sp <- partysplit(as.integer(best), breaks = sp$breaks, 
                         index = sp$index)
        if (verbose) 
          cat(sprintf("Selected split: %s\n\n", paste(character_split(sp, 
                                                                      data = z)$levels, collapse = " | ")))
      }
    }
    kidids <- kidids_split(sp, data = z)
    depth <<- depth + 1L
    kids <- vector(mode = "list", length = max(kidids))
    for (kidid in 1L:max(kidids)) {
      nxt <- kidids == kidid
      if (kidid > 1L) {
        myid <- max(nodeids(kids[[kidid - 1L]]))
      }
      else {
        myid <- id
      }
      kids[[kidid]] <- mob_grow(id = myid + 1L, suby(y, 
                                                     nxt), subx(x, nxt), subz(z, nxt), weights[nxt], 
                                offset[nxt], cluster[nxt],glm.weights[nxt], ...)
    }
    depth <<- depth - 1L
    sp$varid <- as.integer(varindex[sp$varid])
    return(partynode(id = id, split = sp, kids = kids, info = mod))
  }
  mob_prune <- function(node) {
    nd <- as.list(node)
    if (is.null(control$prune)) 
      return(nd)
    id <- seq_along(nd)
    kids <- lapply(nd, "[[", "kids")
    tmnl <- sapply(kids, is.null)
    check <- sapply(id, function(i) !tmnl[i] && all(tmnl[kids[[i]]]))
    while (any(check)) {
      pnode <- which(check)
      objfun <- sapply(nd, function(x) x$info$objfun)
      pok <- sapply(pnode, function(i) control$prune(objfun = c(objfun[i], 
                                                                sum(objfun[kids[[i]]])), df = c(length(nd[[1]]$info$coefficients), 
                                                                                                length(kids[[i]]) * length(nd[[1]]$info$coefficients) + 
                                                                                                  as.integer(control$dfsplit)), nobs = c(nd[[i]]$info$nobs, 
                                                                                                                                         n)))
      pnode <- pnode[pok]
      if (length(pnode) < 1L) 
        break
      pkids <- sort(unlist(sapply(pnode, function(i) nd[[i]]$kids)))
      for (i in id) {
        nd[[i]]$kids <- if (nd[[i]]$id %in% pnode || 
                            is.null(kids[[i]])) {
          NULL
        }
        else {
          nd[[i]]$kids - sapply(kids[[i]], function(x) sum(pkids < 
                                                             x))
        }
      }
      nd[pkids] <- NULL
      id <- seq_along(nd)
      for (i in id) nd[[i]]$id <- i
      kids <- lapply(nd, "[[", "kids")
      tmnl <- sapply(kids, is.null)
      check <- sapply(id, function(i) !tmnl[i] && all(tmnl[kids[[i]]]))
    }
    return(nd)
  }
  depth <- 1L
  nodes <- mob_grow(id = 1L, Y, X, Z, weights, offset, cluster, glm.weights,
                    ...)
  if (verbose && !is.null(control$prune)) 
    cat("-- Post-pruning ---------------------------\n")
  nodes <- mob_prune(nodes)
  for (i in seq_along(nodes)) {
    if (is.null(nodes[[i]]$kids)) {
      nodes[[i]]$split <- NULL
      if (!terminal$estfun) 
        nodes[[i]]$info$estfun <- NULL
      if (!terminal$object) 
        nodes[[i]]$info$object <- NULL
    }
    else {
      if (!inner$estfun) 
        nodes[[i]]$info$estfun <- NULL
      if (!inner$object) 
        nodes[[i]]$info$object <- NULL
    }
  }
  as.partynode(nodes)
}
