#' @title Illustrate.ONE.association.BY.STRATUM
#' @details Plots association blocked by stratum
#' @param data.for.plotting dataset to use for plotting
#' @param TYPES variable properties, ordinal, categorical, continous, etc.
#' @param stratum stratum to check
#' @param x.limits x-axis limits for plotting
#' @param y.limits y-axis limits for plotting
#' @param log_scale use log-scale to plot abundances
#' @author Sudarshan A. Shetty
#'
#' @references
#' Ferreira JA, Fuentes S. (2020). Some comments on certain statistical aspects of
#' the study of the microbiome.
#' \emph{Briefings in bioinformatics} 21(4), pp.1487-1494.
#' \url{https://doi.org/10.1093/bib/bbz077}
#'
#' @export

illustrate.ONE.association.BY.STRATUM <- function(data.for.plotting, TYPES,
                                                  stratum, x.limits, y.limits,
                                                  log_scale) {

  x <- y <- columns <- NULL
  columns <- c(1, 2)
  plotted <- FALSE
  colour.palette <- c(
    "lavender", "lightblue", "cornflowerblue", "aquamarine2",
    "lightgreen", "olivedrab2", "palegreen1", "dodgerblue", "cyan",
    "bisque2", "burlywood2", "darkgoldenrod3",
    "bisque4", "black", "blueviolet", "darkorchid4", "blue",
    "cyan2", "aquamarine4", "chartreuse3", "darkolivegreen2", "darkseagreen3", "grey",
    "darkkhaki", "darkorange4", "brown3", "red", "darksalmon", "orange1", "yellow", "yellow3",
    "springgreen2", "lavender", "sienna2", "pink3", "khaki1",
    "darkmagenta", "azure2"
  )

  data.to.plot <- na.omit(data.for.plotting[, columns])
  type.1 <- TYPES[columns[1]]
  type.2 <- TYPES[columns[2]]
  TYPES[columns]
  if ((type.2 == "binary" & type.1 != "binary") | (type.2 == "categorical" & type.1 != "binary")) {
    data.to.plot <- data.to.plot[, c(2, 1)]
    type.2 <- TYPES[columns[1]]
    type.1 <- TYPES[columns[2]]
    aux.x.limits <- x.limits
    x.limits <- y.limits
    y.limits <- aux.x.limits
  }
  good.names <- create.names(data.to.plot)
  names(data.to.plot) <- c("x", "y")

  if (type.1 == "dirac.and.continuous") {
    if (sum(data.to.plot$x == 0) <= 1) {
      type.1 <- "continuous"
    }
  }
  if (type.2 == "dirac.and.continuous") {
    if (sum(data.to.plot$y == 0) <= 1) {
      type.2 <- "continuous"
    }
  }

  if (is.element(type.1, c("categorical", "binary")) & type.2 == "dirac.and.continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 2])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
      y = ifelse(y > atom, paste(">", atom, sep = ""),
        ifelse(y < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
      )
    )
    contingency.table <- table(aux.data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
      beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
      ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)),
      xlab = good.names[1], main = stratum
    )
    legend(
      x = "top", legend = colnames(frequency.table), horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:nrow(frequency.table)], title = good.names[2]
    )
    aux.data.to.plot <- subset(data.to.plot, y != atom)
    if (length(unique(aux.data.to.plot$x)) == 1) {
      x.label <- paste(good.names[1], "=", unique(aux.data.to.plot$x), collapse = "")
    } else {
      x.label <- good.names[1]
    }
    if (log_scale) {
      boxplot(log(y) ~ x,
        data = aux.data.to.plot, col = "lavender", xlab = x.label, main = stratum,
        ylab = paste("log(", good.names[2], ")", collapse = ""), ylim = y.limits
      )
      stripchart(log(y) ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(y ~ x,
        data = aux.data.to.plot, col = "lavender", xlab = x.label, main = stratum,
        ylab = good.names[2], ylim = y.limits
      )
      stripchart(y ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }

  if (is.element(type.1, c("categorical", "binary")) &
    is.element(type.2, c("categorical", "binary"))) {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    contingency.table <- table(data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
      beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
      ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)), main = stratum,
      xlab = good.names[1]
    )
    legend(
      x = "top", legend = dimnames(frequency.table)$y, horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[2]
    )
    plotted <- TRUE
  }

  if (is.element(type.1, c("categorical", "binary", "ordinal")) & type.2 == "continuous") {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    if (log_scale) {
      boxplot(log(y) ~ x,
        data = data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
        ylab = paste("log(", good.names[2], ")", collapse = ""), ylim = y.limits
      )
      stripchart(log(y) ~ x, data = data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(y ~ x,
        data = data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
        ylab = good.names[2], ylim = y.limits
      )
    }
    plotted <- TRUE
  }

  if (is.element(type.2, c("categorical", "binary", "ordinal")) & type.1 == "continuous") {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    if (log_scale) {
      boxplot(log(x) ~ y,
        data = data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
        ylab = paste("log(", good.names[1], ")", collapse = ""), ylim = x.limits
      )
      stripchart(log(x) ~ y, data = data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(x ~ y,
        data = data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
        ylab = good.names[1], ylim = x.limits
      )
      stripchart(x ~ y, data = data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }

  if ((type.1 == "ordinal" & type.2 == "ordinal") |
    (type.1 == "binary" & type.2 == "ordinal") | (type.1 == "categorical" & type.2 == "ordinal")) {
    aux.table.1 <- table(data.to.plot[, 1])
    aux.table.2 <- table(data.to.plot[, 2])
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    if (length(aux.table.2) > length(aux.table.1)) {
      aux.means <- aggregate(y ~ x, data = data.to.plot, mean)
      boxplot(y ~ x,
        data = data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
        ylab = good.names[2], main = stratum, ylim = y.limits
      )
      points(y ~ as.factor(x),
        data = aux.means, col = "red", xlab = good.names[1], ylab = good.names[2],
        cex = 1.25, pch = 19
      )
    } else {
      aux.means <- aggregate(x ~ y, data = data.to.plot, mean)
      boxplot(x ~ y,
        data = data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
        ylab = good.names[1], main = stratum, ylim = x.limits
      )
      points(x ~ as.factor(y),
        data = aux.means, col = "red", pch = 19, cex = 1.25,
        xlab = good.names[2], ylab = good.names[1]
      )
    }
    plotted <- TRUE
  }

  if (type.1 == "continuous" & type.2 == "continuous") {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    plot(y ~ x,
      data = data.to.plot, col = "blue", xlab = good.names[1], ylab = good.names[2],
      main = stratum, ylim = y.limits, xlim = x.limits
    )
    plotted <- TRUE
  }

  if (type.1 == "ordinal" & type.2 == "dirac.and.continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 2])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
      y = ifelse(y > atom, paste(">", atom, sep = ""),
        ifelse(y < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
      )
    )
    contingency.table <- t(table(aux.data.to.plot))
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
      beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
      ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)),
      xlab = good.names[2]
    )
    legend(
      x = "top", legend = colnames(frequency.table), horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[1]
    )
    aux.data.to.plot <- subset(data.to.plot, y != atom)
    if (log_scale) {
      boxplot(log(y) ~ x,
        data = aux.data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
        ylab = paste("log(", good.names[2], ")", collapse = ""), ylim = y.limits
      )
      stripchart(log(y) ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(y ~ x,
        data = aux.data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
        ylab = good.names[2], ylim = y.limits
      )
      stripchart(y ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }

  if (type.2 == "ordinal" & type.1 == "dirac.and.continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 1])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
      x = ifelse(x > atom, paste(">", atom, sep = ""),
        ifelse(x < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
      )
    )
    contingency.table <- table(aux.data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
      beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
      ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)), main = stratum,
      xlab = good.names[1]
    )
    legend(
      x = "top", legend = colnames(frequency.table), horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[2]
    )
    aux.data.to.plot <- subset(data.to.plot, x != atom)
    if (log_scale) {
      boxplot(log(x) ~ y,
        data = aux.data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
        ylab = paste("log(", good.names[1], ")", collapse = ""), ylim = x.limits
      )
      stripchart(log(x) ~ y, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(x ~ y,
        data = aux.data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
        ylab = good.names[1], ylim = x.limits
      )
      stripchart(x ~ y, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }

  if (type.1 == "dirac.and.continuous" & type.2 == "dirac.and.continuous") {
    aux.table.1 <- table(data.to.plot[, 1])
    aux.table.2 <- table(data.to.plot[, 2])
    atom.1 <- as.numeric(dimnames(aux.table.1)[[1]][which.max(as.vector(aux.table.1))])
    atom.2 <- as.numeric(dimnames(aux.table.2)[[1]][which.max(as.vector(aux.table.2))])
    aux.data.to.plot <- base::transform(data.to.plot,
      x = ifelse(x > atom.1, paste(">", atom.1, sep = ""),
        ifelse(x < atom.1, paste("<", atom.1, sep = ""), paste("=", atom.1, sep = ""))
      ),
      y = ifelse(y > atom.2, paste(">", atom.2, sep = ""),
        ifelse(y < atom.2, paste("<", atom.2, sep = ""), paste("=", atom.2, sep = ""))
      )
    )
    contingency.table <- table(aux.data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    aux.data.to.plot <- subset(data.to.plot, x != atom.1 & y != atom.2)
    if (nrow(aux.data.to.plot) >= 1) {
      par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    }
    barplot(t(frequency.table),
      beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
      ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)), main = stratum,
      xlab = good.names[1]
    )
    legend(
      x = "top", legend = colnames(frequency.table), horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[2]
    )
    if (nrow(aux.data.to.plot) >= 1) {
      plot(y ~ x,
        data = aux.data.to.plot, col = "blue", xlab = good.names[1],
        ylab = good.names[2], main = stratum, ylim = y.limits, xlim = x.limits
      )
      stripchart(y ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    }
    plotted <- TRUE
  }

  if (type.1 == "dirac.and.continuous" & type.2 == "continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 1])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
      x = ifelse(x > atom, paste(">", atom, sep = ""),
        ifelse(x < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
      )
    )
    if (log_scale) {
      boxplot(log(y) ~ x,
        data = aux.data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
        ylab = paste("log(", good.names[2], ")", collapse = ""), ylim = y.limits
      )
      stripchart(log(y) ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(y ~ x,
        data = aux.data.to.plot, col = "lavender", xlab = good.names[1], main = stratum,
        ylab = good.names[2], ylim = y.limits
      )
      stripchart(y ~ x, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
      aux.data.to.plot <- subset(data.to.plot, x != atom)
      plot(y ~ x,
        data = aux.data.to.plot, col = "blue", xlab = good.names[1], ylab = good.names[2],
        main = stratum, ylim = y.limits, xlim = x.limits
      )
    }
    plotted <- TRUE
  }

  if (type.2 == "dirac.and.continuous" & type.1 == "continuous") {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
    aux.table <- table(data.to.plot[, 2])
    atom <- as.numeric(dimnames(aux.table)[[1]][which.max(as.vector(aux.table))])
    aux.data.to.plot <- base::transform(data.to.plot,
      y = ifelse(y > atom, paste(">", atom, sep = ""),
        ifelse(y < atom, paste("<", atom, sep = ""), paste("=", atom, sep = ""))
      )
    )
    if (log_scale) {
      boxplot(log(x) ~ y,
        data = aux.data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
        ylab = paste("log(", good.names[1], ")", collapse = ""), ylim = x.limits
      )
      stripchart(log(x) ~ y, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
    } else {
      boxplot(x ~ y,
        data = aux.data.to.plot, col = "lavender", xlab = good.names[2], main = stratum,
        ylab = good.names[1], ylim = y.limits
      )
      stripchart(x ~ y, data = aux.data.to.plot, col = "red", add = TRUE, vertical = TRUE, cex = 0.5)
      aux.data.to.plot <- subset(data.to.plot, y != atom)
      plot(x ~ y,
        data = aux.data.to.plot, col = "blue", xlab = good.names[2], ylab = good.names[1],
        main = stratum, xlim = x.limits, ylim = y.limits
      )
    }
    plotted <- TRUE
  }

  if (((type.1 == "binary" & type.2 == "ordinal") | (type.1 == "categorical" & type.2 == "ordinal")) &
    (plotted == FALSE)) {
    par(mfrow = c(1, 1), mgp = c(2.5, 1, 0))
    contingency.table <- table(data.to.plot)
    frequency.table <- contingency.table / rowSums(contingency.table)
    barplot(t(frequency.table),
      beside = TRUE, col = colour.palette[1:ncol(frequency.table)],
      ylab = "Proportion", ylim = c(0, 1.15 * max(frequency.table)), main = stratum,
      xlab = good.names[1]
    )
    legend(
      x = "top", legend = dimnames(frequency.table)$y, horiz = TRUE, cex = 1.1, bty = "n",
      fill = colour.palette[1:ncol(frequency.table)], title = good.names[2]
    )
    plotted <- TRUE
  }

  if (plotted == FALSE) {
    print(paste("Did not plot ", paste(columns, collapse = ","), "!", sep = ""))
  }
}
