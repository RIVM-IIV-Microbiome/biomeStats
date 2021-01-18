#' @title Transform to Fractions
#' @details Prepares the taxa abundance table for testing.
#' @param x \code{\link{phyloseq-class}} object
#' @param det.thres remove taxa with less than this abundance in total dataset
#' @param prev.thres remove taxa with less than this prevalence in total dataset
#'
#' @examples
#' # data("ili_data")
#' # ps <- ili_data
#' # gen_tab <- transform_to_fractions(ps, det.thres= 0.0001, prev.thres=5/100)
#' @importFrom tidyr pivot_wider
#' @importFrom tibble rownames_to_column column_to_rownames
#' @export

transform_to_fractions <- function(x,
                                   det.thres = 0.001,
                                   prev.thres=5/100){

  microbiome.data <- microbiome.data.1 <- codes.bacteria <- ID <- asv_tab <- NULL
  Shannon <- value.reabund <- value <- variable <- average <- SD <- NULL
  prop.positive <- CV <- entropy <- NULL
  asv_tab.2 <- other <- test.taxa <- asv_tab.test <- NULL

  #x <- ps

  microbiome.data <- t(as.data.frame(microbiome::abundances(x)))
  asv_tab <- dplyr::as_tibble(microbiome.data, rownames = "ID")

  #names(asv_tab)[1] <- "ID"
  # get column of bacterial abundances
  #message("calculating ....")
  Shannon <- function(p){
    entropy <- ifelse(p>0,-p*log(p),0);
    return(entropy)
    }

  asv_tab2 <- reshape2::melt(asv_tab)
  asv_tab2 <- asv_tab2 %>%
    group_by(ID) %>%
    mutate(average = mean(value, na.rm = TRUE),
           SD = sd(value, na.rm = TRUE),
           prop.positive = mean(value > 0, na.rm = TRUE),
           sum.i = sum(value,na.rm=TRUE),
           value.reabund = value/sum.i,
           entropy = sum(Shannon(value.reabund)),
           CV=SD/average) %>% ungroup()

  #head(subset(asv_tab2, ID=="052_056_S56"))

  #other <- c("average","SD","prop.positive","entropy","CV")
  asv_tab_sub <-  asv_tab2 %>%
    select(ID,average,SD,prop.positive,entropy,CV) %>%
    distinct(ID, .keep_all = T)

  asv_tabr <-  asv_tab2 %>%
    tidyr::pivot_wider(id_cols = ID,
                       names_from = variable,
                       values_from = value.reabund) %>%
    as.data.frame() %>% tibble::column_to_rownames("ID")

  asv_tab.core <- asv_tabr[,core_members(t(asv_tabr),
                                       detection = det.thres,
                                       prevalence =  prev.thres)]

  asv_tab.core <- asv_tab.core %>%
    tibble::rownames_to_column("ID") %>%
    left_join(asv_tab_sub, by="ID")
  # Filter our data

  #message("done ....")

  return(asv_tab.core)

}
