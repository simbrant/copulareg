
##' predict
##' @name predict.copulareg
##' @aliases predict.copulareg
##' @description Computes predictions based on a fitted copulareg model.
##' @param object Model fit as returned by copulareg
##' @param new_x optional matrix of covariate values to compute the predicted
##' values of the outcome for. If not specified, the predicted values for the
##' training sample is returned.
##' @param eps Interval between each interpolation point when integrating to
##' evaluate the predicted value of y, in the case where y is continuous. If
##' y is discrete this parameter is ignored.
##' @param cont_method Specifies the method used to compute the expected
##' values. Can be specified as 'Localmedian' or 'Trapezoidalsurv'. The
##' first method divides the range of the observed values of y into
##' subintervals according to the argument 'eps', where the sub-integral
##' is approximated as the measure of the interval weighted by the local
##' median on the interval. The second method computes the integral by
##' integrating the survival function using the trapezoidal rule, by
##' transforming the outcome into a positive variable by adding a constant.
##' @param ... unused.
##' @rdname predict.copulareg
##' @export
predict.copulareg <- function(object, new_x = NULL, eps = 1E-2,
                            cont_method = "Localmedian", ...) {
  eval_at_u_y <- function(u_y, model, n_test, cond_x, y_type) {

    vinemat <- rvinecopulib::as_rvine_matrix(model$structure)
    h <- rep(0, n_test)
    u_ <- .weave_transformed(u_y,
                             get(.make_hfunc_key(vinemat[1, 1], c()),
                                 envir = cond_x))
    h <- .bicop_2_hbicop(u_, model$pair_copulas[[1]][[1]],
                         return_u_minus = (y_type == "d"))

    for (t in 2:(nrow(vinemat) - 1)) {
      u_ <- .weave_transformed(
        h, get(.make_hfunc_key(vinemat[t, 1], vinemat[1:(t - 1), 1]),
               envir = cond_x))

      h <- .bicop_2_hbicop(u_, model$pair_copulas[[t]][[1]],
                           return_u_minus = (y_type == "d"))
    }
    h
  }

  # Extract x and y types from model
  y_type <- object$model$var_types[object$model$structure$d]
  x_type <- object$model$var_types[1:(object$model$structure$d - 1)]

  # Transform covariates
  if (!is.null(new_x)) {

    n_obs <- nrow(new_x)
    u_x <- object$distr_x$transform(new_x, x_type)

    # Extract covariate model from model
    sub_mod <- object$model
    sub_mod$structure <- rvinecopulib::as_rvine_structure(
      rvinecopulib::as_rvine_matrix(
        sub_mod$structure
      )[1:(sub_mod$structure$d - 1),
        2:sub_mod$structure$d]
    )
    sub_mod$pair_copulas <- lapply(
      1:(sub_mod$structure$d - 1),
      function(l) sub_mod$pair_copulas[[l]][2:(sub_mod$structure$d - (l - 1))]
    )
    sub_mod$var_types <- x_type

    # Compute all the conditionals of the covariate model for the
    # data points where we want to compute the expectation
    cond_x <- .compute_conditionals(u_x, sub_mod)

  } else {
    cond_x <- object$conditionals
    n_obs <- object$model$nobs
  }

  if (y_type == "c") {

    ###################
    # Continuous case #
    ###################

    if (!(cont_method %in% c("Localmedian", "Trapezoidalsurv"))) {
      stop(paste0(c("Argument cont_method must be either 'Localmedian' or",
                    " 'Trapezoidalsurv'!"),
                  collapse = ""))
    }

    if (cont_method == "Localmedian") {
      u_y <- seq(eps, 1, eps)

      # The local median within each subinterval
      y_u <- stats::quantile(object$distr_y$margins[[1]], probs = u_y - eps / 2)

      # Compute the integral, compute the conditional CDF at the right endpoints
      # of the intervals corresponding to each evaluation point first (the last
      # has to be 1)
      uu <- cbind(
        sapply(u_y[-length(u_y)],
               function(u) eval_at_u_y(rep(u, n_obs), object$model, n_obs,
                                       cond_x, y_type)),
        rep(1, n_obs)
      )

      # Compute the
      udiff <- t(apply(uu, 1, diff))
      uu[, 2:length(u_y)] <- udiff

      uu %*% y_u

    } else {
      u_y <- seq(0, 1, eps)

      # Inverse CDF of Y at evaluation points
      y_u <- stats::quantile(object$distr_y$margins[[1]], probs = u_y)

      # Compute the conditional survival function at each evaluation point
      surv <- sapply(u_y, function(u) 1 - eval_at_u_y(rep(u, n_obs),
                                                      object$model, n_obs,
                                                      cond_x, y_type))

      (stats::quantile(y_u, 0) +
          sapply(1:(length(u_y) - 1),
                 function(j) (surv[, j] + surv[, j + 1]) / 2)
        %*% (y_u[-1] - y_u[-length(y_u)]))

    }
  } else {

    #################
    # Discrete case #
    #################

    if (is.null(object$y)) {
      stop("y-values used for fitting the model needed when y is discrete")
    }
    yval <- sort(unique(object$y))
    u_y <- object$distr_y$transform(yval, "d")
    uu <- vapply(1:(nrow(u_y)),
                 function(i) eval_at_u_y(matrix(rep(u_y[i, ], n_obs),
                                                ncol = 2, byrow = T),
                                         object$model, n_obs, cond_x,
                                         y_type = "d"),
                 FUN.VALUE = matrix(0, nrow = n_obs, ncol = 2))
    p_y <- sapply(seq(nrow(u_y)), function(k) apply(uu[, , k], 1,
                                                    function(x) diff(rev(x))))

    p_y %*% yval

  }
}
