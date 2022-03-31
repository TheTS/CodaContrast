
#' geo_mean
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
#' @param percentage Return the difference as a percentage (default = `FALSE`)
#'
#' @return
#'
#' @examples
diff_geo_means <- function(data, i, percentage = FALSE){

  d <- data[i, ]

  d1 <- d[d$group == 1, 1:(ncol(data)-1)]
  g1 <- geo_mean(d1)

  d2 <- d[d$group == 2, 1:(ncol(data)-1)]
  g2 <- geo_mean(d2)

  lr <- as.numeric(log(g2 / g1))

  if (percentage) {
    return((exp(lr) - 1) * 100)
  }

  lr
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
#' @param reference_group The level of `group` that should be treated as the
#' reference group for the contrast
#' @param percentage Return the difference as a percentage (default = `FALSE`)
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
log_ratio_difference <- function(composition, group, n_rep = 1000, alpha = 0.05, confint_type = 'bca', reference_group = NULL, percentage = FALSE) {

  group <- as.factor(group)

  if(nlevels(group) != 2)
    stop("The 'group' must have 2 levels")

  if (is.null(reference_group)) {
    reference <- levels(group)[1]
  } else {
    if (!(reference_group %in% levels(group)))
      stop(paste("The 'reference_group' is not a level within 'group'. The available levels are:", paste(levels(group), collapse = ', ')))
    reference <- reference_group
  }

  group <- ifelse(group == reference, 1, 2)

  df <- as.data.frame(cbind(composition, group))

  bs <- boot(df, diff_geo_means, n_rep, strata = group, percentage = percentage)

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
#' @param h_line_color String. Color of the horozontal line crossing 0.
#' @param percentage Label y axis as percentage difference (default = `FALSE`)
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
plot_log_ratio_difference <- function(data, h_line_color = 'black', percentage = FALSE) {

  ggplot(data, aes(x = name, y = estimate, ymin = lower, ymax = upper)) +
    geom_hline(yintercept = 0, color = h_line_color) +
    geom_errorbar(width = 0.15) +
    geom_point(size = 2) +
    theme(legend.position = "none", axis.title = element_text()) -> p

  if (percentage) {

    p <- p + labs(x = "", y = "Difference (%)", title="", caption = paste0('Reference group: ', data$reference[1]))
  } else {
    p <- p +  labs(x = "", y = "Log-ratio difference", title="", caption = paste0('Reference group: ', data$reference[1]))
  }

  p

}


#' plot_geo_means
#'
#' @description Creates a geometric mean barplot to visualise differences between
#' groups. The logratio between the whole geometric mean and the geometric mean of
#' the group is calculated. These are presented on a log scale.
#'
#' @param composition A composition (from the acomp function)
#' @param group A vector (single column) of the grouping variable. This can be
#' any number of levels.
#' @param type The type of plot to return. Can be one of:
#' \itemize{
#'   \item "component" The colours of the bars represent the components of the composition
#'   \item "group" The colours of the bars represent the levels of `group`.
#'   \item "both" (default) Returns a 2-panel plot with both 'component' and 'group' plots.
#' }
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
plot_geo_means <- function(composition, group, type = 'both') {

  if(!require(patchwork)) {
    stop('The `plot_geo_means` function requires the `patchwork` R package. Please install it first and try again.')
  }

  overall_mean = geo_mean(composition)
  group = as.factor(group)

  result <- data.frame()

  for (i in levels(group)) {
    means <- log(geo_mean(composition[group==i,])/overall_mean)
    means <- data.frame(means, group = i)
    result <- rbind(result, means)
  }

  result <- result %>%
    pivot_longer(cols = -group) %>%
    mutate(name = factor(name, levels = names(composition)))


  p1 <- ggplot(result, aes(x = group, y = value, fill = name)) +
    geom_col(position = position_dodge()) +
    labs(fill = 'Component', y = 'log(gk/G(X))', x = NULL)

  p2 <- ggplot(result, aes(x = name, y = value, fill = group)) +
    geom_col(position = position_dodge()) +
    labs(fill = 'Group', y = 'log(gk/G(X))', x = NULL)

  if (type == 'both') {
    p1 / p2
  } else if (type == 'component') {
    p1
  } else if (type == 'group') {
    p2
  } else {
    stop("The 'type' must be one of: 'both', 'group', 'component'")
  }

}



#' pairwise_hotelling_test
#'
#' The Hotelling T-squared test is the multivariate equivalent of a t-test. This
#' function can be used as a post-hoc test to a MANOVA.
#'
#' @param comp A dataframe or composition object of multiple parts
#' @param groups A vector of the grouping factor. Must be same length as `comp`
#' @param adjust Method of p value adjustment (default ='holm'). See \code{\link[stats]{p.adjust.methods}}
#'
#' @return A dataframe of pairwise contrasts
#' @export
#'
#' @examples
#'
#' @importFrom Hotelling hotelling.test
#' @importFrom dplyr bind_rows
#' @importFrom dplyr case_when
#'
pairwise_hotelling_test <- function(comp, groups, adjust = c('holm')) {

  if(!require(Hotelling)) {
    stop('The `pairwise.hotelling.test` function requires the `Hotelling` R package. Please install it first and try again.')
  }

  groups <- as.factor(groups)

  if(length(groups) != nrow(comp)) {
    stop('The `comp` and `groups` variables must be the same length.')
  }

  lvls <- levels(groups)
  pair <- combn(1:length(lvls), 2)
  result <- list()

  for(i in 1:ncol(pair)) {

    t <- hotelling.test(comp ~ groups, pair = c(pair[,i][1], pair[,i][2]))

    result[[i]] <- data.frame(pair = paste(lvls[pair[,i]], collapse=" vs. "),
                              t_squared = round(t$stats$statistic * t$stats$m, 5),
                              df = paste0('(', t$stats$df[1], ', ', t$stats$df[2], ')'),
                              p.value = round(t$pval, 5))
  }

  result <- bind_rows(result)
  result$p.value.adj <- p.adjust(result$p.value, method = adjust)
  result$sig <- case_when(result$p.value.adj < 0.001 ~ '***',
                          result$p.value.adj < 0.01 ~ '**',
                          result$p.value.adj < 0.05 ~ '*',
                          result$p.value.adj < 0.1 ~ '.',
                          TRUE ~ '')
  result
}


