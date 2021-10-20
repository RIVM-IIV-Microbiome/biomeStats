#' @title Check for taxa-metadata associations
#'
#' @name doAssociationAnalysis
#'
#' @details The associations between taxa and metadata variables are tested.
#'          These are done using the method/approach described by
#'          JA Ferreira and S Fuentes
#'          \url{https://doi.org/10.1093/bib/bbz077}
#'
#' @param data_for_testing A data frame combined from prepTaxaAbundanceData and prepSampleData
#'
#' @param B Resampling default to 100000
#'
#' @param no_of_factors Number of factor from `data_for_testing` to test.
#'
#' @param aux_columns specific column numbers with taxa abundances.
#'
#' @param log_scale Logical. Plot scale. Default and suggested TRUE
#'
#' @param nominal_bound_on_FDR cut-off for multiple testing. Default=0.25
#'
#' @param phen_data data frame output from prepSampleData()
#'
#' @param compare_label label to save files.
#'
#' @param aux_TYPES How to treat each of the variables in the data_for_testing
#'
#' @param cat_ypes Categorical to plot proportions positive
#'
#' @param path_loc Location to store/save output
#'
#' @param name_of_stratification label of stratum e.g."Sex and age-group"
#'
#' @param verbose TRUE
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Ferreira JA, Fuentes S. (2020). Some comments on certain statistical aspects of
#' the study of the microbiome.
#' \emph{Briefings in bioinformatics} 21(4), pp.1487-1494.
#' \url{https://doi.org/10.1093/bib/bbz077}
#'
#' @importFrom readr write_tsv
#' @export

doAssociationAnalysis <- function(data_for_testing,
                                  B = 100000,
                                  no_of_factors= NULL ,
                                  aux_columns = NULL,
                                  log_scale = TRUE,
                                  nominal_bound_on_FDR = 0.25,
                                  phen_data = NULL,
                                  compare_label="case_control",
                                  aux_TYPES=NULL,
                                  cat_ypes = NULL,
                                  path_loc=".",
                                  name_of_stratification= "Sex and age-group",
                                  verbose=TRUE){

  if(is.null(no_of_factors) ||
     is.null(aux_columns) ||
     is.null(phen_data) ||
     is.null(aux_TYPES) ||
     is.null(cat_ypes)){
    warning("Please specify no_of_factors, aux_columns, phen_data, aux_TYPES, cat_ypes")
    stop("These arguments cannot be NULL check on eor more of the arguments")

  }

  stratum <-NULL

  frequency.of.strata <- table(phen_data$stratum)
  #variables.of.interest <- names(data_for_testing)
  variables.of.interest <- names(data_for_testing)
  name_of_stratification <- name_of_stratification

  for(f in (1:no_of_factors)){

    #if(verbose==TRUE){
    # message(paste0("Processing.. ", test.i[1], " vs ", test.i[2]))
    #}

    start.time <- Sys.time()

    columns.f <- cbind(f,aux_columns)
    #head(columns.f)
    tests <- NULL
    cat_ypes <- cat_ypes
    aux_TYPES <- aux_TYPES
    # STEP 1
    for(i in (1:nrow(columns.f))){

      test.i <- test.ONE.association.0(as.vector(columns.f[i,]),
                                       data_for_testing,
                                       aux_TYPES,variables.of.interest)

      #if(verbose==TRUE){
      # message(paste0("Processing.. ", test.i[1], " vs ", test.i[2]))
      #}

      p.value.i <- as.numeric(test.i[3])
      #p.value.i
      if(is.na(p.value.i)){p.value.i <- 2; B.i <- 2*B}
      if(p.value.i<1 & p.value.i>0){
        B.i <- min(max(B,ceiling(400/((1-p.value.i)*p.value.i))),1000000*10*10*10)
      }else{
        B.i <- B
      }

      if(B.i<=B){
        test.i <- test.ONE.association(as.vector(columns.f[i,]),
                                       data_for_testing,aux_TYPES,
                                       variables.of.interest,
                                       B.i)
        p.value.i <- as.numeric(test.i[3])
        no.of.successes <- p.value.i*B.i
        CI.p.value.i <- Wilson.interval(no.of.successes,B.i,0.95)
        CI.p.value.i <- paste("[",CI.p.value.i[1],",",CI.p.value.i[2],"]",collapse="")
      }else{
        #print(c(i,B.i))
        tests.i <- test.ONE.association(as.vector(columns.f[i,]),
                                        data_for_testing,aux_TYPES,
                                        variables.of.interest,
                                        ceiling(B.i/10))

        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),
                                              data_for_testing,
                                              aux_TYPES,
                                              variables.of.interest,
                                              ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),
                                              data_for_testing,
                                              aux_TYPES,
                                              variables.of.interest,
                                              ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),
                                              data_for_testing,
                                              aux_TYPES,
                                              variables.of.interest,
                                              ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),
                                              data_for_testing,
                                              aux_TYPES,
                                              variables.of.interest,
                                              ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),
                                              data_for_testing,aux_TYPES,
                                              variables.of.interest,
                                              ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),
                                              data_for_testing,aux_TYPES,
                                              variables.of.interest,
                                              ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),
                                              data_for_testing,
                                              aux_TYPES,
                                              variables.of.interest,
                                              ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),
                                              data_for_testing,
                                              aux_TYPES,
                                              variables.of.interest,
                                              ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),
                                              data_for_testing,
                                              aux_TYPES,
                                              variables.of.interest,
                                              ceiling(B.i/10)))

        p.values.i <- as.numeric(tests.i[,3])
        no.of.successes <- sum(p.values.i*ceiling(B.i/10))
        CI.p.value.i <- Wilson.interval(no.of.successes,
                                        10*ceiling(B.i/10),0.95)
        p.value.i <- mean(p.values.i)
        CI.p.value.i <- paste("[",CI.p.value.i[1],",",CI.p.value.i[2],"]",collapse="")

      }

      output.i <- data.frame(i=as.vector(columns.f[i,])[1],
                             j=as.vector(columns.f[i,])[2],
                             variable.1=as.character(test.i[1]),
                             variable.2=as.character(test.i[2]),
                             p.value=p.value.i,
                             sign=as.numeric(test.i[4]),
                             test=as.character(test.i[5]),
                             sampling.characteristics=as.character(test.i[6]),
                             CI.p.value=CI.p.value.i)
      #print(output.i)
      tests <- rbind(tests,output.i)
      #tests

    }

    aux.tests <- tests; rownames(aux.tests) <- NULL
    aux.tests <- as.data.frame(aux.tests)
    names(aux.tests) <- c("i","j","variable.1","variable.2","p.value","sign","test","sampling.characteristics","CI.95.pc.p.value")

    aux.tests$i <- as.numeric(as.character(aux.tests$i))
    aux.tests$j <- as.numeric(as.character(aux.tests$j))
    aux.tests$variable.1 <- as.character(aux.tests$variable.1)
    aux.tests$variable.2 <- as.character(aux.tests$variable.2)
    aux.tests$p.value <- as.numeric(as.character(aux.tests$p.value))
    aux.tests$sign <- as.numeric(as.character(aux.tests$sign))
    aux.tests$test <- as.character(aux.tests$test)
    aux.tests$sampling.characteristics <- as.character(aux.tests$sampling.characteristics)
    #head(aux.tests); str(aux.tests)

    path.save.assoc <- paste0(path_loc, "associations.between.",compare_label,".",names(data_for_testing)[f],".and.bacteria.tsv")

    #file.name <- paste("ms/taxa_test/associations.between.", compare_label,".", names(data_for_testing)[f],"and.bacteria.CSV",sep=".")
    readr::write_tsv(aux.tests, file=path.save.assoc)
    #write.csv2(aux.tests,file=path.save.assoc,row.names=FALSE,quote=TRUE)

    if(verbose){
      message(paste0("Processing.. multiple.tests"))
    }
    multiple.tests <- na.omit(aux.tests[,-ncol(aux.tests)])
    multiple.tests$variable.1 <- as.character(multiple.tests$variable.1)
    multiple.tests$variable.2 <- as.character(multiple.tests$variable.2)
    multiple.tests$test <- as.character(multiple.tests$test)
    multiple.tests <- multiple.tests[order(multiple.tests$p.value),]
    row.names(multiple.tests) <- NULL
    aux.histogram <- hist(multiple.tests$p.value, plot=FALSE)
    gamma.hat <- min(1,aux.histogram$density[length(aux.histogram$density)]); gamma.hat
    multiple.tests <- base::transform(multiple.tests,rank=(1:nrow(multiple.tests)))
    bound.FDR <- formatC(nrow(multiple.tests)*multiple.tests$p.value/(1:nrow(multiple.tests)),
                         digits=4,format="f")
    better.bound.FDR <- formatC(gamma.hat*nrow(multiple.tests)*multiple.tests$p.value/(1:nrow(multiple.tests)),
                                digits=4,
                                format="f")
    multiple.tests <- base::transform(multiple.tests,
                                      bound.FDR=as.numeric(bound.FDR),
                                      better.bound.FDR=as.numeric(better.bound.FDR))
    #head(multiple.tests); str(multiple.tests)


    if(any(multiple.tests$better.bound.FDR<=nominal_bound_on_FDR)){

      cut.off <- max((1:nrow(multiple.tests))[multiple.tests$better.bound.FDR<=nominal_bound_on_FDR])
      selected.tests <- multiple.tests[1:cut.off,]
      rownames(selected.tests) <- NULL
      path.save.select <- paste0(path_loc, "SELECTED.associations.between.",compare_label,".",names(data_for_testing)[f],".and.bacteria.tsv")

      if(verbose==TRUE){
        message(paste0("Plotting.. multiple.tests"))
      }
      #path.save <- paste0("ms/taxa_test/iliv1v2/", "SELECTED.associations.between.",compare_label,".",names(data_for_testing)[f],".and.bacteria.CSV",sep=".")
      #file.name <- paste("ms/taxa_test/SELECTED.associations.between", compare_label,".",names(data_for_testing)[f],"and.bacteria.CSV",sep=".")
      readr::write_tsv(selected.tests, file=path.save.select)

      #write.csv2(selected.tests,file=path.save.select,row.names=FALSE,quote=TRUE)

      path.save.select2 <- paste0(path_loc, "SELECTED.associations.between.",compare_label,".",names(data_for_testing)[f],".and.bacteria.pdf")
      #file.name <- paste("ms/taxa_test/SELECTED.associations.between", compare_label,".",names(data_for_testing)[f],"and.bacteria.pdf",sep=".")
      pdf(file=path.save.select2,width=9,height=7)

      par(mfrow=c(1,1))
      hist(multiple.tests$p.value,col="lavender",xlab="P-value",
           main=paste("Associations between ",names(data_for_testing)[f],
                      " and bacterial abundance;\n estimated proportion of potential associations: ",
                      formatC(1-gamma.hat,digits=2,format="f"),sep=""),
           prob=TRUE,cex.lab=1,cex.axis=1,cex.main=1)
      abline(h=1,col="red")
      abline(h=gamma.hat,col="blue",lty=2)

      for(r in (1:nrow(selected.tests))){
        columns.r <- as.numeric(selected.tests[r,1:2])
        illustrate.ONE.association(columns=columns.r,data_for_testing,cat_ypes,log_scale)
      }

      aux.data_for_testing <- subset(data_for_testing,
                                     select=c(selected.tests[1,3],selected.tests[,4],"stratum"))
      aux_TYPES <- cat_ypes[c(selected.tests[1,1],selected.tests[,2])]
      frequency.of.strata <- table(aux.data_for_testing$stratum)
      if(is.factor(aux.data_for_testing[,1]) | is.character(aux.data_for_testing[,1])){
        min.by.stratum <- aggregate(as.formula(paste(names(aux.data_for_testing)[1],
                                                     "stratum",sep="~")),
                                    data=aux.data_for_testing,min)
        max.by.stratum <- aggregate(as.formula(paste(names(aux.data_for_testing)[1],
                                                     "stratum",sep="~")),
                                    data=aux.data_for_testing,max)
        good.strata <- min.by.stratum[max.by.stratum[,2]!=min.by.stratum[,2],1]
      }else{
        sd.by.stratum <- aggregate(as.formula(paste(names(aux.data_for_testing)[1],
                                                    "stratum",sep="~")),
                                   data=aux.data_for_testing,sd)

        good.strata <- sd.by.stratum[sd.by.stratum[,2]>0,1]
      }

      frequency.of.strata <- frequency.of.strata[is.element(dimnames(frequency.of.strata)[[1]],
                                                            good.strata)]

      minimal.sample.size <- 25
      strata.of.interest <- dimnames(frequency.of.strata)[[1]][frequency.of.strata>=minimal.sample.size]
      strata.of.interest

      for(r in (2:(ncol(aux.data_for_testing)-1))){

        if(log_scale){
          x.limits <- range(aux.data_for_testing[,1],na.rm=TRUE)
          y.limits <- range(log(aux.data_for_testing[aux.data_for_testing[,r]>0,r]),na.rm=TRUE)
        }else{
          x.limits <- range(aux.data_for_testing[,1],na.rm=TRUE)
          y.limits <- range(aux.data_for_testing[,r],na.rm=TRUE)
        }

        for(s in strata.of.interest){

          data_for_testing.s <- na.omit(subset(aux.data_for_testing,stratum==s)[c(1,r)])
          title.s <- paste(name_of_stratification," stratum: ",s,sep="")
          if(nrow(data_for_testing.s)>=minimal.sample.size & sd(data_for_testing.s[,2])>0){
            illustrate.ONE.association.BY.STRATUM(data_for_testing.s,aux_TYPES[c(1,r)],stratum=title.s,
                                                  x.limits,y.limits,log_scale=log_scale)
          }
        }
      }

      dev.off()

    }

    end.time <- Sys.time()
    time.taken <- end.time-start.time
    print(time.taken)

  }

}

