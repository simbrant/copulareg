


.compute_expectation <- function(fit, new_x = NULL, eps = 1E-2) {

  eval_at_u_y <- function(u_y, model, n_test, cond_x, y_type) {

    vinemat <- rvinecopulib::as_rvine_matrix(model$structure)
    d <- model$structure$d
    h <- rep(0, n_test)
    u_ <- .weave_transformed(u_y,
                             get(.make_hfunc_key(vinemat[1, 1], c()),
                                 envir=cond_x))
    h <- .bicop_2_hbicop(u_, model$pair_copulas[[1]][[1]],
                         return_u_minus = (y_type == "d"))

    for (t in 2:(nrow(vinemat) - 1)){
      u_ <- .weave_transformed(
        h, get(.make_hfunc_key(vinemat[t, 1], vinemat[1:(t-1), 1]),
               envir=cond_x))

      h <- .bicop_2_hbicop(u_, model$pair_copulas[[t]][[1]],
                           return_u_minus = (y_type == "d"))
    }
    h
  }

  # Extract x and y types from model
  y_type <- fit$model$var_types[fit$model$structure$d]
  x_type <- fit$model$var_types[1:(fit$model$structure$d - 1)]

  # Transform covariates
  if (!is.null(new_x)){

    n_obs <- nrow(new_x)
    u_x <- fit$distr_x$transform(new_x, x_type)

    # Extract covariate model from model
    sub_mod <- fit$model
    sub_mod$structure <- rvinecopulib::as_rvine_structure(
      rvinecopulib::as_rvine_matrix(sub_mod$structure)[1:(sub_mod$structure$d - 1),
                                                       2:sub_mod$structure$d])
    sub_mod$pair_copulas <- lapply(1:(sub_mod$structure$d - 1),
                                   function(l) sub_mod$pair_copulas[[l]][
                                     2:(sub_mod$structure$d- (l - 1))])
    sub_mod$var_types <- x_type

    # Compute all the conditionals of the covariate model for the
    # data points where we want to compute the expectation
    cond_x <- .compute_conditionals(u_x, sub_mod)

  } else {
    cond_x <- fit$conditionals
    n_obs <- fit$model$nobs
  }

  if (y_type == "c"){

    ###################
    # Continuous case #
    ###################

    # Make a grid to integrate y on (on the distribution scale)
    u_y <- seq(eps, 1, eps)

    # Inverse CDF of Y at evaluation points (Try mean instead of median?)
    y_u <- quantile(fit$distr_y$margins[[1]], probs = u_y - eps/2)

    # Compute the integral, compute the conditional CDF at the right endpoints
    # of the intervals corresponding to each evaluation point first (the last
    # has to be 1)
    uu <- cbind(sapply(u_y[-length(u_y)],
                       function(u) eval_at_u_y(rep(u, n_obs), fit$model, n_obs,
                                               cond_x, y_type)),
                rep(1, n_obs))

    # Compute the
    udiff <- t(apply(uu, 1, diff))
    uu[, 2:length(u_y)] <- udiff

    uu %*%y_u

  } else {

    #################
    # Discrete case #
    #################

    if(is.null(y)){
      stop("Must supply y-values used for fitting the model, when y is discrete")
    }
    yval <- sort(unique(fit$y))
    u_y <- fit$distr_y$transform(yval, "d")
    uu <- vapply(1:(nrow(u_y)),
                 function(i) eval_at_u_y(matrix(rep(u_y[i, ], n_obs),
                                                ncol=2, byrow=T),
                                         fit$model, n_obs, cond_x, y_type = "d"),
                 FUN.VALUE = matrix(0, nrow=n_obs, ncol=2))
    p_y <- sapply(1:nrow(u_y), function(k) apply(uu[, , k], 1, function(x) diff(rev(x))))

    p_y %*% yval

  }
}
