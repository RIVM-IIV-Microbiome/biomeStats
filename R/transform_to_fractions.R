#' @title Transform to Fractions
#' @details Prepares the taxa abundance table for testing.
#' @param x \code{\link{phyloseq-class}} object
#' @param filter.abund remove taxa with less than this abundance in total dataset
#' @export

transform_to_fractions <- function(x, filter.abund=0.1){

  microbiome.data <- microbiome.data.1 <- codes.bacteria <- ID <- asv_tab <- NULL
  Shannon <- aux.data.i <- sum.i <- average <- SD <- prop.positive <- CV <- entropy <- NULL
  asv_tab.2 <- other <- test.taxa <- asv_tab.test <- NULL

  microbiome.data <- as.data.frame(microbiome::abundances(x))
  microbiome.data <- dplyr::as_tibble(microbiome.data, rownames = "sequence")
  codes.bacteria <- microbiome.data$sequence

  microbiome.data.1 <- t(as.matrix(microbiome.data[,-1]))
  microbiome.data.1 <- data.frame(microbiome.data.1)

  ID <- row.names(microbiome.data.1)
  row.names(microbiome.data.1) <- NULL
  names(microbiome.data.1) <- as.character(codes.bacteria)

  asv_tab <- cbind(ID, microbiome.data.1)
  #names(asv_tab)[1] <- "ID"
  # get column of bacterial abundances
  message("calculating ....")
  columns.of.bacteria <- (2:ncol(asv_tab))
  Shannon <- function(p){entropy <- ifelse(p>0,-p*log(p),0); return(entropy)}
  asv_tab <- base::transform(asv_tab,
                             average=rep(0,nrow(asv_tab)),
                             SD=rep(0,nrow(asv_tab)),
                             prop.positive=rep(0,nrow(asv_tab)),
                             entropy=rep(0,nrow(asv_tab)))

  for(i in (1:nrow(asv_tab))){
    aux.data.i <- asv_tab[i,columns.of.bacteria]
    aux.data.i <- as.numeric(aux.data.i)
    asv_tab[i,]$average <- mean(aux.data.i,na.rm=TRUE)
    asv_tab[i,]$SD <- sd(aux.data.i,na.rm=TRUE)
    asv_tab[i,]$prop.positive <- mean(aux.data.i>0,na.rm=TRUE)
    sum.i <- sum(aux.data.i,na.rm=TRUE)
    asv_tab[i,columns.of.bacteria] <- asv_tab[i,columns.of.bacteria]/sum.i
    asv_tab[i,]$entropy <- sum(sapply(asv_tab[i,columns.of.bacteria],Shannon))
  }
  # names(data.set); str(data.set), edit(data.set)
  asv_tab <- base::transform(asv_tab,CV=SD/average)

  # Filter our data
  # get col number with last taxa value
  subcols <- ncol(asv_tab) - 5

  asv_tab.2 <- asv_tab[,c(2:subcols)]
  asv_tab.2 <- asv_tab.2[,colSums(asv_tab.2) >= filter.abund]
  #colnames(asv_tab.2)

  test.taxa <- names(asv_tab.2)

  # last 5 cols
  other <- c("average","SD","prop.positive","entropy","CV")

  asv_tab.test <- subset(asv_tab,select=c("ID",test.taxa,other))

  message("done ....")

  return(asv_tab.test)

}
