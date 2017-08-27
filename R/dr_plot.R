#' Plot Delta R values from an object produced by dr_df()
#'
#' This function makes a plot of Delta R values from a list produced by the dr_df() function. It can produce two kinds of plots: either medians of computed Delta R values can be plotted along with their quantiles, or curves of Delta R densities can be plotted on the same plot.
#' @export
#' @param deltar_list a list produced by the dr_df() function.
#' @param name an identification code of a data set. Default is blank.
#' @param method a character vector. If “quantile” is passed to this argument, medians of computed Delta R values will be plotted along with their quantiles. The latter are defined by the “CI” argument during the dr_df() function call. If “density” is passed to this argument, curves of Delta R densities will be plotted on the same plot.
#' @param lim a numeric vector, which contains values used to expand x-axis limits. Default is zero.
#' @usage dr_plot(deltar_list, name = "", method = "quantile", lim = c(0,0))
#' @details This function creates a plot of Delta R values computed for a data set by the dr_df() function. If “quantile” is passed to the “method” argument,  Delta R values computed for each sample or pair are being ordered by their medians and then their medians are plotted along with intervals containing true Delta R value with a previously set probability. If “density” is passed to the “method” argument, curves of Delta R densities will be plotted on the same plot.
#' @return None. Invoked for a side effect (plot).
#' @seealso \code{\link{dr_shell}} \code{\link{dr_df}} \code{\link{dr_pair}}
#' @examples # Acquire "adak" data set
#' data(adak)
#' # Compute Delta R values and store them in the new object
#' adak_res <- dr_df(adak, make.plot = FALSE)
#' # Plot Delta R values with "quantile" method
#' dr_plot(adak_res, name = "Adak")
#' # Expand x-axis of the plot
#' dr_plot(adak_res, name = "Adak", lim = c(100, 200))
#' # Plot Delta R values with "density" method
#' dr_plot(adak_res, name = "Adak", method = "density")

dr_plot <- function(deltar_list, name = "", method = "quantile", lim = c(0,0))
{
  if (typeof(deltar_list) != "list") stop("deltar_list argument must be a list")
  else if (is.vector(lim) != TRUE || length(lim) != 2) stop("lim must be a vector of length 2")
  else
  {
    if (method == "quantile")
    {
      df <- deltar_list[[1]][order(deltar_list[[1]]$median),]
      n <- 1:nrow(df)
      labl <- row.names(df)
      plot(df$median, n, xlim = c(min(df$low.lim) - lim[1], max(df$up.lim) + lim[2]),
           main = name, xlab="Delta R, years", ylab = "", frame.plot = FALSE, yaxt = "n")
      arrows(df$low.lim, n, df$up.lim, n, code = 3, angle = 90, length = 0.1)
      text(df$up.lim, n, labels = labl, pos = 4)
    }else if (method == "density"){
      df <- deltar_list[[2]]
      y.lim <- numeric(length(names(df)))
      for(i in 1:length(names(df)))
      {
        y.lim[i] <- max(density(df[, i])$y)
      }
      plot(density(df[, 1]), type = "l", lty = 1, col = 1, main = name,
           xlab = "Delta R, years", xlim = c(min(unlist(deltar_list[[2]])) - lim[1],
                                           max(unlist(deltar_list[[2]])) + lim[2]), ylim = c(0, max(y.lim)))
      for(i in 2:length(names(df)))
      {
        lines(density(df[, i]), lty = i, col = i)
      }
      legend("topleft", legend = names(df), lty = 1:length(names(df)),
             col = 1:length(names(df)))
    } else stop("Specify method")
  }
}
