#' Calculate Delta R from a pair of dates
#'
#' Computes Delta R value for a radiocarbon dated marine sample, the true age of which was identified with one of the following methods: 1. by measuring other radioactive isotopes, mainly uranium-thorium (230Th/234U), ratio; 2. by radiocarbon dating of its terrestrial counterpart.
#' @export
#' @param dates A numeric vector, which supplies dates to the function. The order of dates in this vector is essential: first, terrestrial (or uranium-thorium) date of a pair, then its standard deviation, then marine radiocarbon date and its standard deviation.
#' @param name An identification code of a pair of dates. Default is blank.
#' @param N Number of iterations. Must be a positive integer. Note that numbers too big can result in very long computation times. Default is 10000
#' @param CI Probability with which true Delta R value is contained within computed interval. Must be between 0 and 1. Default is 0.95.
#' @param calCurves A character vector to determine whether to calibrate the first date or not. If the first date is terrestrial radiocarbon date, it should be calibrated with the “intcal13” calibration curve (or “shcal13” for Southern Hemisphere, see Reimer et al, 2013 and Hogg et al, 2013 for calibration curves details) and the same curve should be supplied to calCurves argument. If it is uranium-thorium date, calibration is not needed and the “normal” curve, which treats ages as normally distributed with given standard deviation, should be passed to this argument. Default is “intcal13”.
#' @param make.plot A logical flag. Makes a histogram of the Delta R values obtained during each iteration if TRUE.
#' @return If the make.plot argument is set to TRUE, the function produces a histogram of Delta R values obtained during each iteration. Regardless of the make.plot argument settings, the function returns a list containing:
#' \itemize{
#'   \item {mean} {- mean Delta R value}
#'   \item {median} {- median Delta R value}
#'   \item {sd}  {- standard deviation of the Delta R value}
#'   \item {quantile} {- an interval containing true Delta R value with a probability set by the “CI” argument. Default is 0.95}
#'   \item {p.value} {- the p value of Kolmogorov – Smirnov test of normality of computed Delta R values}
#'   \item {delta} {- a numeric vector of length N containing a series of assessments of Delta R value obtained during each iteration}
#' }
#' @details Delta R is defined as a difference between measured and modeled radiocarbon ages of a marine sample (Stuiver, Braziunas, 1993). Measured radiocarbon age of a sample is supplied to the function by the “dates” argument. To compute modeled radiocarbon age the function first looks at the “calCurves” argument. If it is supplied with “intcal13” or “shcal13” curve, the function calibrates terrestrial date of a pair with the BchronCalibrate() function from the “Bchron” package (Haslett, Parnell, 2008; Parnell et al., 2008) using “intcal13” or “shcal13” calibration curve (Reimer et al. 2013; Hogg et al, 2013). Calibration creates a grid of ages covering that part of the calibration curve, for which age probabilities are greater than 1e-05, and calculates the probability of each age from this range. \cr \cr
#' If the “calCurves” argument is supplied with the “normal” curve, BchronCalibrate() treats the first date of a pair as normally distributed with given standard deviation. For each age from the grid produced by BchronCalibrate() it is possible to look for a corresponding modeled radiocarbon age in the “marine13” data set (Reimer et al. 2013). For full details and references of the “marine13” calibration curve see http://www.radiocarbon.org/IntCal13.htm. For this and other functions of the “deltaR” package “marine13”, ”intcal13” and “shcal13” data sets are made available via the “Bchron” package (Haslett, Parnell, 2008; Parnell et al., 2008).\cr \cr
#' Then the function computes Delta R value in a series of iterations, number of which is defined by the “N” argument. During each iteration it chooses: 1. One year of measured radiocarbon age of the marine sample. The latter is supplied as its mean and standard deviation thus forming normal distribution from which the year is chosen according to its probability. 2. One year from the grid of ages produced by BchronCalibrate(). The function chooses a year according to its probability and then looks for a corresponding modeled radiocarbon age in the “marine13” data set. 3. Modeled radiocarbon age is also accompanied with its standard deviation, so the function again chooses one year from the corresponding normal distribution. Then it finds the difference between measured and modeled radiocarbon ages. This gives a series of assessments of Delta R value stored in a vector of length N. The function then computes basic statistics of these Delta R assessments and produces a histogram of Delta R probability densities. The histogram has a total area of one and is supplemented with a curve of probability densities of a corresponding normal distribution. \cr \cr
#' The results of the function call can be slightly different with the same arguments supplied to function. This difference is due to the chosen numerical method of computations and is insignificant with reasonable number of iterations.
#' @examples ## Calculation of Delta R for the last pair (pair 9)
#' ## from the “adak” data set (see Khasanov et al., 2015
#' ## for details). Terrestrial date of this pair is 835 ± 20
#' ## (laboratory code NUTA2-20937), marine date 1765 ± 20
#' ## (laboratory code NUTA2-20553).
#'
#' # Compute and store the results in the new object
#' p9 <- dr_pair(dates = c(835, 20, 1765, 20), name = "Adak, pair 9")
#' # Mean value of Delta R
#' p9$mean
#' # Median value of Delta R
#' p9$median
#' # Standard deviation of Delta R
#' p9$sd
#' # An interval containing true Delta R value with probability 0.95
#' p9$quantile
#' # p value of Kolomogorov - Smirnov test
#' p9$p.value
#'
#' ## Calculation of Delta R for the pair of terrestrial and
#' ## marine materials from the Lazaret Midden (Moreton Bay,
#' ## Australia). Charcoal sample was dated 500 ± 50
#' ## (laboratory code Wk-8009) and marine shell yielded 840 ± 50
#' ## (laboratory code Wk-8013). Data from Ulm et al, 2009, table 2.
#'
#' # Compute and store the results in the new object
#' # Note, that for 14C dates from the Southern Hemisphere
#' # "shcal13" curve is used
#' LM <- dr_pair(c(500, 50, 840, 50), name = "Lazaret Midden", calCurves = "shcal13")
#'
#' ## Calculation of Delta R for the coral M2-3 dated
#' ## 2170 ± 15 with 230Th and 2550 ± 30 with radiocarbon
#' ## (data from Yu et al, 2010; table 2)
#'
#' # Compute and store the results in the new object
#' # Note, that 230Th dates do not need calibration
#' M2_3 <- dr_pair(c(2170, 15, 2550, 30),  name = "M2-3", calCurves = "normal")
#' @usage dr_pair(dates, name = "", N = 10000, CI = 0.95,
#' calCurves = "intcal13", make.plot = FALSE)
#' @seealso \code{\link{dr_shell}} \code{\link{dr_df}} \code{\link{dr_plot}}
#' @references Haslett J, Parnell AC. 2008. A simple monotone process with application to radiocarbondated depth chronologies. Journal of the Royal Statistical Society, Series C. 57: 399-418.
#' Hogg AG, Hua Q, Blackwell PG, Niu M, Buck CE, Guilderson TP, Heaton TJ, Palmer JG,  Reimer PJ, Ron W, Turney CSM, Zimmerman SR. 2013. SHCal13 Southern Hemisphere calibration, 0–50,000 cal yr BP. Radiocarbon. 55 (4): 1889-1903. \cr \cr
#' Khasanov BF, Nakamura T, Okuno M, Gorlova EN, Krylovich OA, West DL, Hatfield V, Savinetsky AB. 2015. The Marine Radiocarbon Reservoir Effect on Adak Island (Central Aleutian Islands), Alaska. Radiocarbon. 57(5): 955-964 \cr \cr
#' Parnell AC, Haslett J, Allen JRM, Buck CE, Huntley B. 2008. A flexible approach to assessing synchroneity of past events using Bayesian reconstructions of sedimentation history. Quaternary Science Reviews. 27(19-20): 1872-1885. \cr \cr
#' Reimer PJ, Bard E, Bayliss A, Beck JW, Blackwell PG, Bronk Ramsey C, Buck CE, Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP, Haflidason H, Hajdas I, Hatté C, Heaton TJ, Hoffmann DL, Hogg AG, Hughen KA, Kaiser KF, Kromer B, Manning SW, Niu M, Reimer RW, Richards DA, Scott EM, Southon JR, Staff RA, Turney CSM, van der Plicht J. 2013. IntCal13 and Marine13 radiocarbon age calibration curves 0–50,000 years cal BP. Radiocarbon 55(4): 1869–87. \cr \cr
#' Stuiver M, Braziunas TF. 1993. Modeling atmospheric 14C influences and 14C ages of marine samples to 10,000 BC. Radiocarbon 35(1):137–89. \cr \cr
#' Ulm S, Petchey F, Ross A. 2009. Marine reservoir corrections for Moreton Bay, Australia. Archaeol. Oceania. 44: 160–168. \cr \cr
#' Yu K, Hua Q, Zhao J, Hodge E, Fink D, Barbetti M. 2010. Holocene marine 14C reservoir age variability: Evidence from 230Th‐dated corals in the South China Sea. Paleoceangraphy. 25: PA3205. \cr \cr


dr_pair <- function (dates, name = "", N = 10000, CI = 0.95,
                     calCurves = "intcal13", make.plot = FALSE)
{
  if (is.vector(dates) != TRUE || length(dates) != 4) stop("Check dates supplied to the function")
  else if (typeof(N) != "double" || length(N) != 1 || N <= 0 || N%%1 != 0) stop("N must be a positive integer")
  else if (typeof(CI) != "double" || length(CI) != 1 || CI <= 0 || CI >= 1) stop ("CI must be a number between 0 and 1")
  else if (calCurves %in% c("intcal13", "shcal13", "normal") == FALSE) stop("calCurves must be either intcal13, shcal13 or normal")
  else
  {
    me <- new.env()
    data(marine13, package = "Bchron", envir = me)
    cal.t <- Bchron::BchronCalibrate(dates[1], dates [2], calCurves)
    Year <- sample(cal.t$Date1$ageGrid, N, replace = TRUE,
                   prob = cal.t$Date1$densities)
    f <- function(Year)
    {
      prm1 <- me$marine13[which.min(abs(me$marine13[, 1] - Year)), 2]
      prm2 <- me$marine13[which.min(abs(me$marine13[, 1] - Year)), 3]
      model <- rnorm(1, prm1, prm2)
      c14 <- rnorm(1, dates[3], dates[4])
      delta <- c14 - model
    }
    deltaR <- sapply(Year, f)
    out_list <- list(mean = mean(deltaR), median = median (deltaR), sd = sd(deltaR),
                     quantile = quantile(deltaR, c((1 - CI)/2, 1-(1-CI)/2)),
                     p.value = ks.test(deltaR,"pnorm", mean(deltaR), sd(deltaR))$p.value,
                     delta = deltaR)
    if (make.plot == FALSE) return(out_list)
    else
    {
      hist(deltaR, freq = FALSE, ylim = c(0, max(c(density(deltaR)$y,
                                                   dnorm(seq(min(deltaR), max(deltaR), length = 500), mean(deltaR), sd(deltaR))))),
           main = name, xlab = "Delta R, years")
      lines(seq(min(deltaR), max(deltaR), length = 500),
            dnorm(seq(min(deltaR), max(deltaR), length = 500), mean(deltaR), sd(deltaR)),
            col = "red")
      return(out_list)
    }
  }
}
