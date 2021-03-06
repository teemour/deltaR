#' Delta R from the radiocarbon date of the marine material with known collection date
#'
#' Computes Delta R value from the radiocarbon date of a marine sample with known collection date. In most cases such materials are represented by molluscan shells from museum collections.
#' @export
#' @md
#' @param dates A numeric vector, which supplies dates to the function. The order of dates in this vector is essential: first, a collection date (in calendar years AD/BC), then marine radiocarbon date and its standard deviation.
#' @param name Laboratory code of a radiocarbon date or any other sample ID. Default is blank.
#' @param N  Number of iterations. Must be a positive integer. Note that numbers too big can result in very long computation times. Default is 10000.
#' @param CI Probability with which true Delta R value is contained within computed interval. Must be between 0 and 1. Default is 0.95.
#' @param make.plot A logical flag. Makes a histogram of the Delta R values obtained during each iteration if TRUE.
#' @return If the "make.plot" argument is set to TRUE, the function produces a histogram of Delta R values obtained during each iteration. Regardless of the "make.plot" argument settings, the function returns a list containing:
#' * mean — mean Delta R value
#' * median — median Delta R value
#' * sd — standard deviation of Delta R value
#' * quantile — an interval containing true Delta R value with a probability set by the "CI" argument. Default is 0.95
#' * p.value — the p value of Kolmogorov – Smirnov test of normality of computed Delta R values
#' * delta — a numeric vector containing a series of assessments of Delta R value obtained during each iteration
#' @usage dr_shell(dates, name = "", N = 10000, CI = 0.95,  make.plot = FALSE)
#' @details Delta R or regional correction is defined as a difference between measured and modeled radiocarbon ages of a sample (Stuiver, Braziunas, 1993). Measured radiocarbon age of a sample is supplied to the function by the "dates" argument. To compute modeled radiocarbon age the function first transforms a collection date into the year BP and looks for a corresponding modeled radiocarbon age in the "marine13" data set (Reimer et al. 2013). For full details and references of the "marine13" calibration curve see http://www.radiocarbon.org/IntCal13.htm. For this and other functions of the 'deltar' package the "marine13" data set is made available via the 'Bchron' package (Haslett, Parnell, 2008; Parnell et al., 2008). \cr \cr
#' Then the function computes Delta R value in a series of iterations, number of which is defined by the "N" argument. During each iteration it chooses one year of measured radiocarbon age. The latter is supplied as its mean and standard deviation, thus forming normal distribution from which the year is chosen according to its probability. The function does the same with the modeled radiocarbon age and finds the difference between the two. This gives a series of assessments of Delta R values stored in a vector of length N. The function then computes basic statistics of these Delta R assessments and produces a histogram of Delta R probability densities. The histogram has a total area of one and is supplemented with a curve of probability densities of the corresponding normal distribution. \cr \cr
#' The results of the function call can be slightly different with the same arguments supplied to the function. This difference is due to the chosen numerical method of computations and  is insignificant with reasonable number of iterations.
#' @examples ## Calculation of Delta R for Pseudocardium sybillae shell
#' ## collected in Sakhalin in 1906 and radiocarbon dated
#' ## 826 ± 35 years BP, laboratory code TERRA-072305a15
#' ## (data from Yoneda et al., 2007, table 2)
#'
#' # Compute and store the results in the new object
#' Ps <- dr_shell(c(1906, 826, 35), name = "TERRA-072305a15")
#' # Mean value of Delta R
#' Ps$mean
#' # Median value of Delta R
#' Ps$median
#' # Standard deviation of Delta R
#' Ps$sd
#' # An interval containing true Delta R value with probability 0.95
#' Ps$quantile
#' # p value of Kolomogorov - Smirnov test
#' Ps$p.value
#' @seealso \code{\link{dr_pair}} \code{\link{dr_df}} \code{\link{dr_plot}}
#' @references Haslett J, Parnell AC. 2008. A simple monotone process with application to radiocarbondated depth chronologies. Journal of the Royal Statistical Society, Series C. 57: 399-418. <\doi{10.1111/j.1467-9876.2008.00623.x}>  \cr \cr
#' Parnell AC, Haslett J, Allen JRM, Buck CE, Huntley B. 2008. A flexible approach to assessing synchroneity of past events using Bayesian reconstructions of sedimentation history. Quaternary Science Reviews. 27(19-20): 1872-1885. <\doi{10.1111/j.1467-9876.2008.00623.x}> \cr \cr
#' Reimer PJ, Bard E, Bayliss A, Beck JW, Blackwell PG, Bronk Ramsey C, Buck CE, Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP, Haflidason H, Hajdas I, Hatté C, Heaton TJ, Hoffmann DL, Hogg AG, Hughen KA, Kaiser KF, Kromer B, Manning SW, Niu M, Reimer RW, Richards DA, Scott EM, Southon JR, Staff RA, Turney CSM, van der Plicht J. 2013. IntCal13 and Marine13 radiocarbon age calibration curves 0–50,000 years cal BP. Radiocarbon 55(4): 1869–87. <\doi{10.2458/azu_js_rc.55.16947}> \cr \cr
#' Stuiver M, Braziunas TF. 1993. Modeling atmospheric 14C influences and 14C ages of marine samples to 10,000 BC. Radiocarbon 35(1):137–89. <\doi{10.1017/S0033822200013874}> \cr \cr
#' Yoneda M, Uno H, Shibata Y, Suzuki R, Kumamoto Y, Yoshida K, Sasaki T, Suzuki A, Kawahata H. 2007. Radiocarbon marine reservoir ages in the western Pacific estimated by pre-bomb molluscan shells. Nuclear Instruments and Methods in Physics Research B 259(1): 432–7. <\doi{10.1016/j.nimb.2007.01.184}>

dr_shell <- function (dates, name = "", N = 10000, CI = 0.95, make.plot = FALSE)
{
  if (is.vector(dates) != TRUE || length(dates) != 3) stop("Check dates supplied to the function")
  else if (typeof(N) != "double" || length(N) != 1 || N <= 0 || N%%1 != 0) stop("N must be a positive integer")
  else if (typeof(CI) != "double" || length(CI) != 1 || CI <= 0 || CI >= 1) stop ("CI must be a number between 0 and 1")
  else
  {
    me <- new.env()
    data(marine13, package = "Bchron", envir = me)
    prm1 <- me$marine13[which.min(abs(me$marine13[, 1] - (1950 - dates[1]))), 2]
    prm2 <- me$marine13[which.min(abs(me$marine13[, 1] - (1950 - dates[1]))), 3]
    model <- rnorm(N, prm1, prm2)
    c14 <- rnorm(N, dates[2], dates[3])
    delta <- c14 - model
    out_list <- list(mean = mean(delta), median = median (delta), sd = sd(delta),
                     quantile = quantile(delta, c((1 - CI)/2, 1 - (1 - CI)/2)),
                     p.value = ks.test(delta,"pnorm", mean(delta), sd(delta))$p.value,
                     delta = delta)
    if (make.plot == FALSE) return(out_list)
    else
    {
      hist(delta, freq = FALSE, ylim = c(0, max(c(density(delta)$y,
                                                  dnorm(seq(min(delta), max(delta), length = 500), mean(delta), sd(delta))))),
           main = name, xlab = "Delta R, years")
      lines(seq(min(delta), max(delta), length = 500),
            dnorm(seq(min(delta), max(delta), length = 500), mean(delta), sd(delta)),
            col = "red")
      return(out_list)
    }
  }
}
