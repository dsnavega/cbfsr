#' Genetic Algorithm Feature Selection Parameters
#'
#' @author David Senhora Navega
#' @export
#' @import stats
#'
#' @param x a symmetric numeric matrix with inter-correlation (association) of
#' the features (input) under study. See Details.
#' @param y a vector with the correlation (association) of x (features) with the
#' class (output). See Details
#' @param control a list of parameters for genetic algorithm as created by
#' \code{\link{gafs.control}}
#'
#' @return a gafs object
#'
#' @details
#' In this implementation, and following Hall (1999), correlation is used its
#' more general sense. Is up to the user to compute the matrix x and vector y
#' using the correlation / association measure that best suite the data and
#' problem at hand.
#'
#' @references
#' Hall, M.A., 1999. Correlation-based feature selection for machine learning.
#'
gafs <- function(x, y, control = gafs.control()) {

  if (!is.gafs.control(control))
    stop("\n(-) control is not a gafs.control object")

  start <- Sys.time()

  same <- 0
  generation <- 0
  kept <- NULL
  fittest <- NULL
  generational <- NULL
  pressure_0 <- control$pressure
  crossover_0 <- control$crossover
  mutation_0 <- control$mutation

  ngeneration <- control$ngeneration
  warmup <-  control$warmup
  early <- control$early
  adaptive <- control$adaptive
  rate <- control$rate

  # Create Initial Population
  population <- create_binary_population(x = x, control = control)

  while (generation < ngeneration) {

    # Evaluate Fitness
    fitness <- sapply(1:nrow(population), function(ith) {

      chromossome <- as.logical(population[ith, ])
      x_rho <- x[chromossome, chromossome]
      y_rho <- y[chromossome]
      fitness.value <- correlation_merit(x_rho, y_rho)
      return(fitness.value)

    })

    # Keep Best
    best <- which.max(fitness)
    kept <- rbind(kept, population[best,, drop = F])

    # Evolve Population
    population <- evolve_binary_population(
      x = population, fitness = fitness, control = control
    )

    # Update Information
    fittest <- c(fittest, max(fitness, na.rm = T))
    generational <- c(generational, stats::median(fitness, na.rm = T))
    generation <- generation + 1


    if (adaptive & (generation > (warmup * ngeneration))) {
      correct <- warmup * ngeneration
      adapt <- exp(-rate * (generation - correct))
      control$pressure <- 1 - (1 - pressure_0) * adapt
      control$crossover <- 1 - (1 - crossover_0) * adapt
      control$mutation <- pmax(mutation_0 * adapt, 0.001)
    }

    # Early Stop Rule
    if (generation > ceiling(warmup * ngeneration)) {

      if (fittest[generation] == fittest[generation - 1]) {
        same <- same + 1
      } else {
        same <- 0
      }

    }

    if (same >= ceiling(early * ngeneration))
      break

  }

  solution <- gafs.test(
    fittest = kept, alpha = control$alpha, digits = control$digits
  )

  solution <- data.frame(feature = colnames(x), as.data.frame(solution))
  rownames(solution) <- NULL

  object <- structure(
    .Data = list(
      n = nrow(kept),
      solution = solution,
      fittest = kept,
      fitness = data.frame(
        Generation = seq_len(nrow(kept)),
        Population = generational,
        Fittest = fittest,
        Cardinality = rowSums(kept)
      ),
      control = control,
      runtime = Sys.time() - start
    ),
    class = c("gafs", "cbfsr")
  )

  return(object)

}

#' Genetic Algorithm Feature Selection Parameters
#'
#' @author David Senhora Navega
#' @export
#'
#' @param npopulation number of individuals in the initial population.
#' @param ngeneration number of generation
#' @param pressure selective pressure probability. Probability of selecting the
#' best individual in a binary tournament selection
#' @param crossover crossover probability.
#' @param mutation mutation probability.
#' @param elitism fraction of best individuals that are transferred directly to
#' the next generation.
#' @param warmup fraction of generations at which early stop and adaptive
#' behavior will trigger.
#' @param early fraction of consecutive generations at which the GA will be
#' terminated if no improvement is detected to the fitness function.
#' @param adaptive a logical stating if probabilities defined by pressure,
#' crossover and mutation should be adapted according to an exponential function
#' after the warm-up phase. Default is TRUE.
#' @param rate the rate of the exponential function to adapt the parameters of
#' the GA. Default is 0.01 (i.e. probability will increase or decrease by 1%)
#' @param alpha level of significance of binomial test performed to assess
#' if a feature should be accepted or rejected based the the best individual
#' chromosomes of all generations. Default is 0.01
#' @param digits an integer defining the precision of the round function.
#' Default is 3.
#'
gafs.control <- function(
  npopulation = 100, ngeneration = 200,
  pressure = 0.8, crossover = 0.8, mutation = 0.05, elitism = 0.05,
  warmup = 0.5, early = 0.05, adaptive = T, rate = 0.01, alpha = 0.01,
  digits = 3
) {

  control.list <- structure(
    .Data = list(
      npopulation = npopulation, ngeneration = ngeneration,
      crossover = crossover, mutation = mutation, pressure = pressure,
      elitism = elitism, warmup = warmup, early = early, adaptive = adaptive,
      rate = rate, alpha = alpha, digits = digits
    ),
    class = "gafs.control"
  )

  return(control.list)

}

#' @author David Senhora Navega
#' @noRd
is.gafs.control <- function(x) {
  return(inherits(x = x, what = "gafs.control"))
}

#' @author David Senhora Navega
#' @export
#'
#' @noRd
#'
print.gafs <- function(x, ...) {

  mfeature <- sum(x$solution$select)
  nfeature <- length(x$solution$select)

  cat("\n Genetic Algorithm for Feature Selection:")
  cat("\n Selected", mfeature, "of", nfeature, "features.")
  cat("\n Results: \n\n")

  solution <- x$solution
  print(solution[solution$select, ])

}

#' @author David Senhora Navega
#' @export
#' @noRd
#' @import tidyr
#'
plot.gafs <- function(x,...) {

  datum <- with(x,
    tidyr::gather(
      data = x$fitness, key = "Group", value = "Fitness", Fittest:Population
    )
  )

  p <- with(datum,
    ggplot2::ggplot(data = datum) +
      ggplot2::geom_point(
        mapping = aes(
          x = Generation,
          y = Fitness,
          shape = Group
        )
      ) +
      ggplot2::guides(
        shape = ggplot2::guide_legend(
          title.hjust = 0.5, title.position = "top",
          label.position = "right", label.hjust = 0,
          override.aes = list(size = 2)
        )
      ) +
      ggplot2::scale_x_continuous(
        limits = range(c(0,datum[[1]])),
        breaks = scales::extended_breaks(n = 10)
      ) +
      ggplot2::scale_y_continuous(
        limits = range(datum[[4]]),
        breaks = scales::extended_breaks(n = 10)
      ) +
      ggplot2::geom_vline(xintercept = 0.0, linetype = 2) +
      ggplot2::geom_hline(yintercept = max(datum[[4]]), linetype = 2) +
      ggplot2::theme_classic(base_line_size = 1) +
      ggplot2::theme(legend.position = "right") +
      ggplot2::ggtitle("Feature Selection Genetic Algorithm Fitness Evolution\n")
  )
  return(p)

}


#' @author David Senhora Navega
#' @noRd
#'
#' @note
#' Implementation is a combination of roulette wheel selection and binary
#' tournament selection.
#'
binary_selection <- function(x, fitness, probability = 0.75) {

  # Population Size
  n <- nrow(x)

  # Fitness Proportional Selection Probability (Wheel)
  wheel <- fitness / sum(fitness, na.rm = T)
  wheel[is.na(wheel)] <- 0

  select_individual <- function() {

    # Roullette Wheel Selection
    index <- sample(x = seq_len(n), size = 2, replace = F, prob = wheel)

    # Tournament Selection
    selected <- runif(n = 1, min = 0, max = 1) < probability
    if (selected) {
      index <- index[which.max(fitness[index])]
    } else {
      index <- index[which.min(fitness[index])]
    }
    return(index)

  }

  # Apply Operator
  index <- replicate(n = 2, expr = select_individual())
  x <-   x[index, ]

  return(x)

}

#' @author David Senhora Navega
#' @noRd
#'
binary_crossover <- function(x, probability = 0.8) {

  # Chromosome Length
  m <- ncol(x)

  # Duplicate
  z <- x

  # Create Uniform Crossover Mask
  mask <- sample(c(FALSE, TRUE), size = m, replace = T)

  # Recombination
  recombine <- runif(n = 1, min = 0, max = 1) < probability
  if (recombine) {
    x[1, mask] <- z[2, mask]
    x[2, mask] <- z[1, mask]
  } else {
    x <- z
  }

  return(x)

}

#' @author David Senhora Navega
#' @noRd
#'
binary_mutation <- function(x, probability = 0.05) {

  n <- nrow(x)
  p <- ncol(x)

  for (i in seq_len(n)) {
    for (j in seq_len(p)) {
      zero_one <- c(0, 1)
      if (runif(n = 1, min = 0, max = 1) < probability) {
        x[i, j] <- zero_one[x[i, j] != zero_one]
      }
    }
  }

  return(x)

}

#' @author David Senhora Navega
#' @noRd
#'
create_binary_population <- function(x, control) {

  npopulation <- control$npopulation
  nloci <- ncol(x)
  size <- nloci  * npopulation
  aleles <- sample(c(0, 1), size = size, replace = T, prob = c(0.5, 0.5))
  chromosomes <- matrix(aleles, nrow = npopulation, ncol = nloci)
  dimnames(chromosomes) <- list(seq_len(npopulation), colnames(x))

  return(chromosomes)

}

#' @author David Senhora Navega
#' @noRd
#'
evolve_binary_population <- function(x, fitness, control) {

  n <- nrow(x) / 2
  new_individuals <- function(x, fitness) {
    x <- binary_selection(x, fitness = fitness, probability = control$pressure)
    x <- binary_crossover(x, probability = control$crossover)
    x <- binary_mutation(x, probability = control$mutation)
    return(x)
  }

  evolved_population <- do.call(
    what = rbind,
    args = lapply(
      X = seq_len(length.out = n),
      FUN = function(i) {
        new_individuals(x = x, fitness = fitness)
      }
    )
  )

  if (control$elitism > 0 & control$elitism < 1) {
    nelite <- pmax(1, ceiling(n * control$elitism))
    elite <- seq_len(length.out = nelite)
    elite <- order(fitness, decreasing = T)[elite]
    evolved_population <- rbind(x[elite, ], evolved_population)
  }

  return(evolved_population)

}

#' @author David Senhora Navega
#' @noRd
#'
correlation_merit <- function(x, y) {

  if (NROW(x) != NROW(y))
    stop("\n(-) x and y do not have the same length.")

  k <- NROW(x)
  merit <- (k * mean(y)) / sqrt(k + (k * (k - 1)) * mean(x))
  return(merit)

}

#' @author David Senhora Navega
#' @noRd
#' @import stats
#'
gafs.test <- function(fittest, alpha = 0.01, digits = 3) {

  decision <- rep("Tentative", ncol(fittest))
  nloci <- colSums(fittest)
  nfittest <- nrow(fittest)

  # Accept Feature
  p.accept <- stats::pbinom(nloci - 1, nfittest, 0.5, F)
  p.accept <- stats::p.adjust(p = p.accept, method = "bonferroni")
  accept <- p.accept < alpha

  # Reject Feature
  p.reject <- stats::pbinom(nloci, nfittest, 0.5, T)
  p.reject <- stats::p.adjust(p = p.reject, method = "bonferroni")
  reject <- p.reject < alpha

  decision[accept] <- "Accept"
  decision[reject] <- "Reject"
  names(decision) <- colnames(fittest)
  select <- decision %in% c("Accept", "Tentative")
  names(select) <- colnames(fittest)

  object <- list(
    select = select,
    probability = round(x = colMeans(fittest), digits = digits),
    decision = decision
  )

  return(object)

}

