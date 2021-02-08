##' predict
##' @name plot.copulareg
##' @aliases plot.copulareg
##' @description Plot the residuals against the fitted values for a copulareg
##' object, or predicted values against the prediction error
##' @param x Model fit as returned by copulareg
##' @param new_x optional matrix of covariate values to compute the predicted
##' values of the outcome for. If not specified, the fitted values for the
##' training sample are used.
##' @param new_y Optional vector if the plot should show predicted values and
##' prediction error.
##' @param ... additional parameters to plot.
##' @rdname plot.copulareg
##' @export
##' 
plot.copulareg <- function(x, new_x=NULL, new_y=NULL, ...){

  if (is.null(new_x) & is.null(new_y)){
    graphics::plot(copulareg::predict.copulareg(x),
                   x$y - copulareg::predict.copulareg(x), 
                   ylab = "Residuals", xlab = "Fitted values",
                   main = "Copula regression: Fitted values vs Residuals",
                   ...)
  } else if (!is.null(new_x) & !is.null(new_y)){
    graphics::plot(copulareg::predict.copulareg(x, new_x),
                   new_y - copulareg::predict.copulareg(x, new_x),
                   ylab = "Prediction error", xlab = "Predicted values",
                   main = "Copula regression: Predicted values vs Error",
                   ...)
  } else {
    stop("Must supply both new_x and new_y, not only one of the two.")
  }

}