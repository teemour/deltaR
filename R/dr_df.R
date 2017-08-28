#' Compute Delta R values for a set of dates
#'
#' Computes Delta R values for a set of dates in a data frame either with dr_shell() or dr_pair() function. Multiple radiocarbon dates of marine samples with known collection dates or true ages identified in other ways can be downloaded as a data frame and Delta R for each column of the latter is computed
#' @export
#' @importFrom grDevices x11
#' @importFrom graphics arrows hist legend lines plot text
#' @importFrom stats density dnorm ks.test median quantile rnorm sd
#' @importFrom utils data
#' @md
#' @param df A data frame. Its first column indicates a type of dates constituting each sample or pair. If Delta R is computed for dates amalgamated in pairs, each pair is stored in its own column. The order of rows is essential: first, terrestrial (or uranium-thorium) date of a pair, then its standard deviation, then marine radiocarbon date, and its standard deviation (see "adak" and "coral" datasets as examples). If Delta R is computed for marine samples with known collection date, each sample is stored in its own column. The order of rows is essential: first, the collection date (in calendar years AD/BC), then radiocarbon date of a sample and its standard deviation (see "BSea" dataset as an example). Column names are further used as identification codes.
#' @param N Number of iterations. Must be a positive integer. Note that numbers too big can result in very long computation times. Default is 10000.
#' @param CI Probability with which true Delta R value is contained within computed interval. Must be between 0 and 1. Default is 0.95.
#' @param method A character vector. If Delta R is computed for dates amalgamated in pairs, "pair" should be supplied to this argument. If Delta R is computed for marine samples with known collection date, "shell" should be supplied to this argument.
#' @param calCurves A character vector to determine whether to calibrate the first date of a pair or not. If the first date is terrestrial radiocarbon date, it should be calibrated with "intcal13" calibration curve (or "shcal13" for Southern Hemisphere, see Reimer et al, 2013 and Hogg et al, 2013 for calibration curves details) and the same curve has to be supplied to the calCurves argument. If it is uranium-thorium date, calibration is not needed and "normal" curve, which treats ages as normally distributed with given standard deviation, should be passed to calCurves. Default is "intcal13".
#' @param make.plot A logical flag. Makes histograms of the Delta R values for each sample or pair if TRUE.
#' @usage dr_df(df, N = 10000, CI = 0.95, method = "pair",
#'                        calCurves = "intcal13", make.plot = FALSE)
#' @details The function computes Delta R values for each column of the data frame according to the chosen method. See details of dr_pair() and dr_shell() functions. See examples of data frame format in "adak" and "coral" data sets for method "pair" and in "BSea" data set for method "shell". The output of this function can be supplied to the dr_plot() function for plotting Delta R values. \cr \cr
#' If the make.plot argument is set to TRUE, the function produces histograms of Delta R values for each sample or a pair. The histograms have a total area of one and are supplemented with curves of probability densities of the corresponding normal distributions.
#' @return If the make.plot argument is set to TRUE, the function produces histograms of Delta R values obtained during each iteration. Histogram of each sample or a pair is produced in its own window. Regardless of the make.plot argument settings, function returns a list containing:
#'
#' * statistics — a data frame. Basic statistics computed for pairs or samples supplied to the function are the following: mean (mean Delta R value), median (median Delta R value), sd (standard deviation of Delta R value), low.lim and up.lim (the lower and upper limits of an interval containing true Delta R value with a probability set by the CI argument), p.value (p value of the Kolmogorov – Smirnov test of normality of computed Delta R values). These statistics are organized as columns, while pairs or samples supplied to the function are organized as rows.
#' * deltaR — a data frame. For each pair or a sample supplied to the function a numeric vector of length N containing a series of assessments of Delta R values obtained during each iteration is created. All of them are combined in this data frame with an identification code of a pair or a sample as a column name.
#'
#' @seealso \code{\link{dr_shell}} \code{\link{dr_pair}} \code{\link{dr_plot}}
#' @references Hogg AG, Hua Q, Blackwell PG, Niu M, Buck CE, Guilderson TP, Heaton TJ, Palmer JG,  Reimer PJ, Ron W, Turney CSM, Zimmerman SR. 2013. SHCal13 Southern Hemisphere calibration, 0–50,000 cal yr BP. Radiocarbon. 55 (4): 1889-1903. <\doi{10.2458/azu_js_rc.55.16783}> \cr \cr
#' Reimer PJ, Bard E, Bayliss A, Beck JW, Blackwell PG, Bronk Ramsey C, Buck CE, Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP, Haflidason H, Hajdas I, Hatté C, Heaton TJ, Hoffmann DL, Hogg AG, Hughen KA, Kaiser KF, Kromer B, Manning SW, Niu M, Reimer RW, Richards DA, Scott EM, Southon JR, Staff RA, Turney CSM, van der Plicht J. 2013. IntCal13 and Marine13 radiocarbon age calibration curves 0–50,000 years cal BP. Radiocarbon 55(4): 1869–87. <\doi{10.2458/azu_js_rc.55.16947}>

#' @examples
#' #' # Acquire "BSea" data set
#' data(BSea)
#' # Compute Delta R values and store them in the new object
#' # Note, that samples with known collection date need "shell" method
#' bsea_res <- dr_df(BSea, method = "shell")
#' # See basic statistics of the computed Delta R values
#' bsea_res$statistics
#' \dontrun{
#' # Acquire "adak" data set
#' data(adak)
#' # Compute Delta R values and store them in the new object
#' adak_res <- dr_df(adak)
#' # See basic statistics of the computed Delta R values
#' adak_res$statistics
#' # Save the results in the file "Adak.txt"
#' write.table(adak_res$statistics, "Adak.txt")
#' # Note, that it will be saved in the working directory
#' # See the path to the working directory
#' getwd()
#' # Acquire "coral" data set
#' data(coral)
#' # Compute Delta R values and store them in the new object
#' # Note, that 230Th dates do not need calibration
#' coral_res <- dr_df(coral, calCurves = "normal")
#' # See basic statistics of the computed Delta R values
#' coral_res$statistics
#'}
dr_df <- function (df, N = 10000, CI = 0.95, method = "pair",
                   calCurves = "intcal13", make.plot = FALSE)
{
  if (is.data.frame(df) != TRUE) stop("df argument must be a data frame")
  else if (typeof(N) != "double" || length(N) != 1 || N <= 0 || N%%1 != 0) stop("N must be a positive integer")
  else if (typeof(CI) != "double" || length(CI) != 1 || CI <= 0 || CI >= 1) stop ("CI must be a number between 0 and 1")
  else if (calCurves %in% c("intcal13", "shcal13", "normal") == FALSE) stop("calCurves must be either intcal13, shcal13 or normal")
  else
  {
    if (method == "pair") {
      ml <- lapply(df[, -1], dr_pair, N = N, CI = CI, calCurves = calCurves, make.plot = FALSE)
    } else if (method =="shell") {
      ml <- lapply(df[, -1], dr_shell, N = N, CI = CI, make.plot = FALSE)
    } else stop("Specify method")
    mean <- sapply(ml, "[[", 1)
    median <- sapply(ml, "[[", 2)
    sd <- sapply(ml, "[[", 3)
    low.lim <- sapply(ml, "[[", 4)[1, ]
    up.lim <- sapply(ml, "[[", 4)[2, ]
    p.value <- sapply(ml, "[[", 5)
    sum.stats <- data.frame(mean, median, sd, low.lim, up.lim, p.value)
    deltaR<-as.data.frame(sapply(ml, "[[", 6))
    out.list <- list(statistics = sum.stats, deltaR = deltaR)
    dfplot <- function(df)
    {
      for(i in 1:length(names(df)))
      {
        x11()
        hist(df[, i], freq = FALSE, ylim = c(0, max(c(density(df[, i])$y,
            dnorm(seq(min(df[, i]), max(df[, i]), length = 500), mean(df[, i]), sd(df[, i]))))),
            main = names(df)[i], xlab = "Delta R, years")
        lines(seq(min(df[, i]), max(df[, i]), length = 500),
            dnorm(seq(min(df[, i]), max(df[, i]), length = 500), mean(df[, i]), sd(df[, i])),
            col = "red")
      }
    }
    if (make.plot == FALSE) return(out.list)
    else
    {
      dfplot(out.list$deltaR)
      return(out.list)
    }
  }
}
