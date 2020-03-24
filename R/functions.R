
#' diff_geo_means
#'
#' @param x .
#'
#' @return
#'
#' @examples
geo_mean <- function (x) {
  g = t(exp(apply(log(x), 2, mean)))
  return(g/sum(g))
}


#' diff_geo_means
#'
#' @param data .
#' @param i .
#'
#' @return
#'
#' @examples
diff_geo_means <- function(data, i){

  d <- data[i, ]

  d1 <- d[d$group == 1, 1:(ncol(data)-1)]
  g1 <- geo_mean(d1)

  d2 <- d[d$group == 2, 1:(ncol(data)-1)]
  g2 <- geo_mean(d2)

  return(as.numeric(log(g2 / g1)))
}


#' log_ratio_difference
#'
#' @description Returns the log-ratio difference for each component of a composition
#' between two groups. The confidence intervals of the difference are estimated
#' via bootstrapping.
#'
#' @param composition A composition (from the acomp function)
#' @param group A vector (single column) of the grouping variable.
#' This variable must only have two levels (e.g. 'Male' and 'Female').
#' @param n_rep Number of bootstrap iterations
#' @param alpha Alpha (i.e. p value) to use for confidence intervals
#' @param confint_type Type of confidence interval. One of ("norm","basic", "stud", "perc", "bca").
#' Default value is 'bca'
#'
#' @references Martín Fernández, Josep Antoni, Josep Daunis i Estadella, and Glòria Mateu i Figueras.
#' "On the interpretation of differences between groups for compositional data."
#' SORT: statistics and operations research transactions, 2015, vol. 39, núm. 2, p. 231-252 (2015).
#'
#' @return
#'
#' @export
#'
#' @examples
#'
#' @importFrom boot boot
#' @importFrom  boot boot.ci
log_ratio_difference <- function(composition, group, n_rep = 1000, alpha = 0.05, confint_type = 'bca') {

  group <- as.factor(group)

  if(nlevels(group) != 2)
    stop("The 'group' must have 2 levels")

  reference <- levels(group)[1]
  group <- as.integer(group)

  df <- as.data.frame(cbind(composition, group))

  bs <- boot(df, diff_geo_means, n_rep, strata = group)

  c_int <- data.frame(lower = numeric(), upper = numeric())

  for(i in 1:ncol(composition)) {
    c_int[nrow(c_int)+1,] <- boot.ci(bs, index = i, conf = 1 - alpha, type = confint_type)[[4]][4:5]
  }

  data.frame(name = factor(names(composition), levels = names(composition)), estimate = bs$t0, c_int, reference)
}


#' plot_log_ratio_difference
#'
#' @description Plots the log-ratio difference and confidence intervals that
#' were estimated from the \code{\link{log_ratio_difference}} function.
#'
#' @param data A dataframe output from \code{\link{log_ratio_difference}}
#'
#' @references Martín Fernández, Josep Antoni, Josep Daunis i Estadella, and Glòria Mateu i Figueras.
#' "On the interpretation of differences between groups for compositional data."
#' SORT: statistics and operations research transactions, 2015, vol. 39, núm. 2, p. 231-252 (2015).
#'
#' @return
#' @export
#'
#' @examples
#'
#' @import ggplot2
#'
plot_log_ratio_difference <- function(data) {

  ggplot(data, aes(x = name, y = estimate, ymin = lower, ymax = upper)) +
    geom_hline(yintercept = 0, lty=2) +
    geom_pointrange() +
    theme(legend.position = "none", axis.title = element_text()) +
    labs(x = "", y = "Log-ratio difference", title="", caption = paste0('Reference group: ', data$reference[1]))
}


#' plot_geo_means
#'
#' @description Creates a geometric mean barplot to visualise differences between
#' groups.
#'
#' @param composition A composition (from the acomp function)
#' @param group A vector (single column) of the grouping variable. This can be
#' any number of levels.
#'
#' @references Martín Fernández, Josep Antoni, Josep Daunis i Estadella, and Glòria Mateu i Figueras.
#' "On the interpretation of differences between groups for compositional data."
#' SORT: statistics and operations research transactions, 2015, vol. 39, núm. 2, p. 231-252 (2015).
#'
#' @return
#' @export
#'
#' @examples
#'
#' @import patchwork
#'
plot_geo_means <- function(composition, group) {

  overall_mean = geo_mean(composition)
  group = as.factor(group)

  result <- data.frame()

  for (i in levels(group)) {
    means <- log(geo_mean(composition[group==i,])/overall_mean)
    means <- data.frame(means, group = i)
    result <- rbind(result, means)
  }

  result <- result %>%
    pivot_longer(cols = -group)


  p1 <- ggplot(result, aes(x = group, y = value, fill = name)) +
    geom_col(position = position_dodge()) +
    labs(fill = 'Component', y = 'log(gk/G(X))', x = NULL)

  p2 <- ggplot(result, aes(x = name, y = value, fill = group)) +
    geom_col(position = position_dodge()) +
    labs(fill = 'Group', y = 'log(gk/G(X))', x = NULL)

  p1 / p2
}
