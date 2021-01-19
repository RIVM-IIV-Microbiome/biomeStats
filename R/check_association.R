#' @title Check for taxa-metadata associations
#' @details The associations between taxa and metadata variables are tested. These are done using
#' the method/appraoched described by JA Ferreira and S Fuentes https://doi.org/10.1093/bib/bbz077
#' @param data.for.testing A data frame combined from transform_to_fractions and prep_metadata.
#' @param B Resampling default to 100000
#' @param no.of.factors Number of factor from `data.for.testing` to test.
#' @param aux.columns specific column numbers with taxa abundances.
#' @param log.scale Logical. Plot scale. Default and suggested TRUE
#' @param nominal.bound.on.FDR cut-off for multiple testing. Default=0.25
#' @param phen.data data frame output from prep_metadata()
#' @param compare.label label to save files.
#' @param aux.TYPES How to treat each of the variables in the data.for.testing
#' @param cat.types Categorical to plot proportions positive
#' @param path_loc Location to store/save output
#' @param name.of.stratification label of stratum e.g."Sex and age-group"
#' @param verbose TRUE
#' @export

check_association <- function(data.for.testing,
                              B = 100000,
                              no.of.factors= NULL ,
                              aux.columns = NULL,
                              log.scale = TRUE,
                              nominal.bound.on.FDR = 0.25,
                              phen.data = NULL,
                              compare.label="case_control",
                              aux.TYPES=NULL,
                              cat.types = NULL,
                              path_loc=".",
                              name.of.stratification= "Sex and age-group",
                              verbose=TRUE){

  if(is.null(no.of.factors) | is.null(aux.columns) | is.null(phen.data) | is.null(aux.TYPES) | is.null(cat.types)){
    warning("Please specify no.of.factors, aux.columns, phen.data, aux.TYPES, cat.types")
    stop("These arguments cannot be NULL check on eor more of the arguments")

  }



  stratum <-
    frequency.of.strata <- table(phen.data$stratum)
  variables.of.interest <- names(data.for.testing)
  name.of.stratification <- name.of.stratification

  # f <- 1
  for(f in (1:no.of.factors)){

    #if(verbose==TRUE){
    # message(paste0("Processing.. ", test.i[1], " vs ", test.i[2]))
    #}

    start.time <- Sys.time()

    columns.f <- cbind(f,aux.columns)
    #head(columns.f)
    tests <- NULL
    cat.types <- cat.types
    aux.TYPES <- aux.TYPES

    if(verbose==TRUE){
      print(paste0(variables.of.interest[as.vector(columns.f[i,])][1], " vs Taxa"))
    }


    # STEP 1
    for(i in (1:nrow(columns.f))){

      test.i <- test.ONE.association.0(as.vector(columns.f[i,]),data.for.testing,
                                       aux.TYPES,variables.of.interest)

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
        test.i <- test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,B.i)
        p.value.i <- as.numeric(test.i[3])
        no.of.successes <- p.value.i*B.i
        CI.p.value.i <- Wilson.interval(no.of.successes,B.i,0.95)
        CI.p.value.i <- paste("[",CI.p.value.i[1],",",CI.p.value.i[2],"]",collapse="")
      }else{
        #print(c(i,B.i))
        tests.i <- test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10))

        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))
        tests.i <- rbind(tests.i,
                         test.ONE.association(as.vector(columns.f[i,]),data.for.testing,aux.TYPES,variables.of.interest,ceiling(B.i/10)))

        p.values.i <- as.numeric(tests.i[,3])
        no.of.successes <- sum(p.values.i*ceiling(B.i/10))
        CI.p.value.i <- Wilson.interval(no.of.successes,10*ceiling(B.i/10),0.95)
        p.value.i <- mean(p.values.i)
        CI.p.value.i <- paste("[",CI.p.value.i[1],",",CI.p.value.i[2],"]",collapse="")

      }

      output.i <- data.frame(i=as.vector(columns.f[i,])[1],j=as.vector(columns.f[i,])[2],
                             variable.1=as.character(test.i[1]),variable.2=as.character(test.i[2]),
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

    path.save.assoc <- paste0(path_loc, "associations.between.",compare.label,".",names(data.for.testing)[f],".and.bacteria.CSV")

    #file.name <- paste("ms/taxa_test/associations.between.", compare.label,".", names(data.for.testing)[f],"and.bacteria.CSV",sep=".")
    write.csv2(aux.tests,file=path.save.assoc,row.names=FALSE,quote=TRUE)
    if(verbose==TRUE){
      message(paste0("Processing.. multiple.tests"))
    }
    multiple.tests <- na.omit(aux.tests[,-ncol(aux.tests)])
    multiple.tests$variable.1 <- as.character(multiple.tests$variable.1)
    multiple.tests$variable.2 <- as.character(multiple.tests$variable.2)
    multiple.tests$test <- as.character(multiple.tests$test)
    multiple.tests <- multiple.tests[order(multiple.tests$p.value),]; row.names(multiple.tests) <- NULL
    aux.histogram <- hist(multiple.tests$p.value,plot=FALSE)
    gamma.hat <- min(1,aux.histogram$density[length(aux.histogram$density)]); gamma.hat
    multiple.tests <- base::transform(multiple.tests,rank=(1:nrow(multiple.tests)))
    bound.FDR <- formatC(nrow(multiple.tests)*multiple.tests$p.value/(1:nrow(multiple.tests)),
                         digits=4,format="f")
    better.bound.FDR <- formatC(gamma.hat*nrow(multiple.tests)*multiple.tests$p.value/(1:nrow(multiple.tests)),
                                digits=4,format="f")
    multiple.tests <- base::transform(multiple.tests,
                                      bound.FDR=as.numeric(bound.FDR),
                                      better.bound.FDR=as.numeric(better.bound.FDR))
    #head(multiple.tests); str(multiple.tests)


    if(any(multiple.tests$better.bound.FDR<=nominal.bound.on.FDR)){

      cut.off <- max((1:nrow(multiple.tests))[multiple.tests$better.bound.FDR<=nominal.bound.on.FDR])
      selected.tests <- multiple.tests[1:cut.off,]
      rownames(selected.tests) <- NULL
      path.save.select <- paste0(path_loc, "SELECTED.associations.between.",compare.label,".",names(data.for.testing)[f],".and.bacteria.CSV")

      if(verbose==TRUE){
        message(paste0("Plotting.. multiple.tests"))
      }
      #path.save <- paste0("ms/taxa_test/iliv1v2/", "SELECTED.associations.between.",compare.label,".",names(data.for.testing)[f],".and.bacteria.CSV",sep=".")
      #file.name <- paste("ms/taxa_test/SELECTED.associations.between", compare.label,".",names(data.for.testing)[f],"and.bacteria.CSV",sep=".")
      write.csv2(selected.tests,file=path.save.select,row.names=FALSE,quote=TRUE)
      path.save.select2 <- paste0(path_loc, "SELECTED.associations.between.",compare.label,".",names(data.for.testing)[f],".and.bacteria.pdf")
      #file.name <- paste("ms/taxa_test/SELECTED.associations.between", compare.label,".",names(data.for.testing)[f],"and.bacteria.pdf",sep=".")
      pdf(file=path.save.select2,width=9,height=7)

      par(mfrow=c(1,1))
      hist(multiple.tests$p.value,col="lavender",xlab="P-value",
           main=paste("Associations between ",names(data.for.testing)[f],
                      " and bacterial abundance;\n estimated proportion of potential associations: ",
                      formatC(1-gamma.hat,digits=2,format="f"),sep=""),
           prob=TRUE,cex.lab=1,cex.axis=1,cex.main=1)
      abline(h=1,col="red")
      abline(h=gamma.hat,col="blue",lty=2)

      for(r in (1:nrow(selected.tests))){
        columns.r <- as.numeric(selected.tests[r,1:2])
        illustrate.ONE.association(columns=columns.r,data.for.testing,cat.types,log.scale)
      }

      aux.data.for.testing <- subset(data.for.testing,
                                     select=c(selected.tests[1,3],selected.tests[,4],"stratum"))
      aux.types <- cat.types[c(selected.tests[1,1],selected.tests[,2])]
      frequency.of.strata <- table(aux.data.for.testing$stratum)
      if(is.factor(aux.data.for.testing[,1]) | is.character(aux.data.for.testing[,1])){
        min.by.stratum <- aggregate(as.formula(paste(names(aux.data.for.testing)[1],"stratum",sep="~")),
                                    data=aux.data.for.testing,min)
        max.by.stratum <- aggregate(as.formula(paste(names(aux.data.for.testing)[1],"stratum",sep="~")),
                                    data=aux.data.for.testing,max)
        good.strata <- min.by.stratum[max.by.stratum[,2]!=min.by.stratum[,2],1]
      }else{
        sd.by.stratum <- aggregate(as.formula(paste(names(aux.data.for.testing)[1],"stratum",sep="~")),
                                   data=aux.data.for.testing,sd)
        good.strata <- sd.by.stratum[sd.by.stratum[,2]>0,1]
      }

      frequency.of.strata <- frequency.of.strata[is.element(dimnames(frequency.of.strata)[[1]],good.strata)]

      minimal.sample.size <- 25
      strata.of.interest <- dimnames(frequency.of.strata)[[1]][frequency.of.strata>=minimal.sample.size]
      strata.of.interest

      for(r in (2:(ncol(aux.data.for.testing)-1))){

        if(log.scale){
          x.limits <- range(aux.data.for.testing[,1],na.rm=TRUE)
          y.limits <- range(log(aux.data.for.testing[aux.data.for.testing[,r]>0,r]),na.rm=TRUE)
        }else{
          x.limits <- range(aux.data.for.testing[,1],na.rm=TRUE)
          y.limits <- range(aux.data.for.testing[,r],na.rm=TRUE)
        }

        for(s in strata.of.interest){

          data.for.testing.s <- na.omit(subset(aux.data.for.testing,stratum==s)[c(1,r)])
          title.s <- paste(name.of.stratification," stratum: ",s,sep="")
          if(nrow(data.for.testing.s)>=minimal.sample.size & sd(data.for.testing.s[,2])>0){
            illustrate.ONE.association.BY.STRATUM(data.for.testing.s,aux.types[c(1,r)],stratum=title.s,
                                                  x.limits,y.limits,log.scale=TRUE)
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
