##' This tool allows for the visualization and analysis of the harvest output from an
##'     Atlantis model. This tool includes analysis of catch, bycatch, and discarding
##'     either at the level of the entire population or by age class. In addition to
##'     graphical output, it is also possible to perform a skill assessment of the
##'     model performance (currently only for the time series of catch and
##'     bycatch). The quantitative metrics used to analyze the performance of the
##'     model is based on the approach described by Olsen \emph{et al.} (2016) and Stow \emph{et
##'     al.} (2009) using: the correlation coefficient, root mean squared error,
##'     reliability index, average error, average absolute error and the modeling
##'     efficiency routines.
##' @title Analysis of the Harvest in Atlantis
##' @param grp.csv Character string with the path to the Groups \emph{*.csv} file (Atlantis input file).
##' @param fish.csv Character string with the path to the fisheries \emph{*.csv} file (Atlantis input file)
##' @param catch.nc Character string with the path to the catch netcdf output
##'     file from Atlantis. Usually this file has the name of \emph{[Your_Model]CATCH.nc},
##'     where [Your_Model] is the name of your Atlantis model
##' @param ext.catch.f (Default = NULL) Character string with the path to the
##'     external file with the shiny::observed catches and discards by year. This helps to
##'     calibrate the harvest section of Atlantis and it is required to perform the
##'     skill assessment of the model.
##' @return This function provides 3 different sets of analyzes of the catches in the following tabs:
##' \itemize{
##' \item \bold{Biomass}: This function helps to analyze change through time of all
##'     the variables (i.e. catch and discard) for the selected Fishery. This
##'     information can be displayed by each recording time step (\emph{toutfinc}
##'     parameter) or by year.
##' \item \bold{Numbers}:  This function allows the user to analyze the change by age and
##'     through time of the variable catch and discards for the selected functional group
##' \item \bold{Compare}: This function allows for several skill assessment of the model
##'     (based on the analysis of the simulated and shiny::observed time series of catch and
##'     bycatch) to be performed. The analysis is based on the approach described by
##'     Olsen \emph{et al.} (2016) and Stow \emph{et al.} 2009, which is composed of the following
##'     quantitative metrics:
##' \itemize{
##' \item \bold{Correlation coefficient (\eqn{r})}: measures the tendency of the
##'     predicted \eqn{P} and shiny::observed \eqn{O} values to vary together. The values of
##'     correlation can range from -1 to 1, with negative values of correlation for
##'     time series that vary inversely.
##' \deqn{r = \frac{\displaystyle\sum_{i=1}^{n} (O_{i} - \bar{O})(P_{i} - \bar{P})}{\sqrt{\displaystyle\sum_{i=1}^{n} (O_{i} - \bar{O})^2 \displaystyle\sum_{i=1}^{n} (P_{i} - \bar{P})^2}}}
##' \item \bold{Average Error (\eqn{AE}; or model bias)}: The average error is a
##'     measure of aggregated model bias.
##' \deqn{AE = \frac{\displaystyle\sum_{i=1}^{n} (P_{i} - O_{i})}{n}  = \bar{P} - \bar{O}}
##' \item \bold{average absolute error (\eqn{AAE}) and root mean squared error
##'     (\eqn{RMSE})} : Both equations calculate the bias of the
##'     model, considering the magnitude rather than the direction of each
##'     discrepancy. That is because values near zero in the \eqn{AE} can be
##'     misleading,  as negative and positive discrepancies can cancel each other.
##' \deqn{AAE = \frac{\displaystyle\sum_{i=1}^{n} |P_{i} - O_{i}|}{n}} \cr
##' \deqn{RMSE = \sqrt{\frac{\displaystyle\sum_{i=1}^{n} (P_{i} - O_{i})}{n}}}
##' \item \bold{Reliability index \eqn{RI}}: This quantifies the average factor by
##'     which the predicted values differ from observations. If the model predictions
##'     do not differ too much from the shiny::observed then the value of \eqn{RI} should be close to 1. But if the
##'     value of \eqn{RI} is 2 it means that a model predicts the observations
##'     within a multiplicative factor of two, on average.
##' \deqn{RI = exp\sqrt{\frac{1}{n} \displaystyle\sum_{i=1}^{n} (log\frac{O_{i}}{P_{i}})^2}}
##' \item \bold{Modeling efficiency \eqn{MEF}}: The modeling efficiency measures how
##'     well a model predicts relative to the average of the observations. A value on
##'     \eqn{MEF} close to 1 indicates that the model matches the
##'     observations. A value of \eqn{MEF} less than zero means that the average of
##'     the observations is a better predictor than the model.
##' \deqn{MEF = \frac{\displaystyle\sum_{i=1}^{n} (O_{i} - \bar{O})^2 - \displaystyle\sum_{i=1}^{n} (P_{i} - O_{i})^2}{ \displaystyle\sum_{i=1}^{n} (O_{i} - \bar{O})^2 }}
##' }
##' where \eqn{n} is the number of observations;  \eqn{O_{i}} the \eqn{ith} of \eqn{n} observations;
##'     \eqn{P_{i}} ith of \eqn{n} predictions, and \eqn{\bar{O}} and \eqn{\bar{P}} are the
##'     observation and prediction averaged, respectively.
##' }
##' @import stats utils grDevices ggplot2 graphics
##' @author Demiurgo
##' @export
catch <- function(grp.csv, fish.csv, catch.nc, ext.catch.f = NULL){
    ## General setting
    mycol    <- c(RColorBrewer::brewer.pal(8, "Dark2"), c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
    mycol    <- grDevices::colorRampPalette(mycol)
    col.bi   <- mycol(15)[c(14, 13, 12, 10)]
    nc.data  <- ncdf4::nc_open(catch.nc)
    grp      <- utils::read.csv(grp.csv)
    if(!is.null(ext.catch.f)){
        ext.f  <- scan(ext.catch.f, nlines = 1, what = character(), sep = ',')
        ext.c  <- utils::read.table(ext.catch.f, skip = 1, header = TRUE, sep = ',', check.names = FALSE)
    }
    grp      <- grp[grp$IsImpacted == 1, ]
    fsh      <- utils::read.csv(fish.csv)
    nam.var  <- names(nc.data$var)
    orign    <- unlist(strsplit(ncdf4::ncatt_get(nc.data, 't')$units, ' ', fixed = TRUE))
    if(orign[1] == 'seconds') {
        Time <- ncdf4::ncvar_get(nc.data, 't') / 86400
    } else {
        Time <- ncdf4::ncvar_get(nc.data, 't')
    }
    Time     <- as.Date(Time, origin = orign[3])
    shiny::shinyApp(
        ## Create the different tabs
        ui <- shiny::navbarPage("Catch",
                         shiny::tabPanel('Biomass',
                                  shiny::fluidRow(
                                      shiny::column(2,
                                             shiny::wellPanel(
                                                 shiny::selectInput('FISHB', 'Fishery :', as.character(fsh$Code)),
                                                 shiny::selectInput('is.CB', label = shiny::strong("Data type"), c('Catch', 'Discard')),
                                                 shiny::checkboxInput('rem.cb', label = shiny::strong("Remove Values"), value = FALSE),
                                                 shiny::numericInput("lag.cb", "Observations:", 1, min = 1, max = 100),
                                                 shiny::checkboxInput('b.year', label = shiny::strong("By year"), value = FALSE),
                                                 shiny::downloadButton("DL_Biomass", "Download")
                                                 )),
                                      shiny::column(10,
                                             shiny::plotOutput('plotB', width = "100%", height = "800px")
                                             ))),
                         shiny::tabPanel('Numbers',
                                  shiny::fluidRow(
                                      shiny::column(2,
                                             shiny::wellPanel(
                                                 shiny::selectInput('FG', 'Functional Group :', as.character(grp$Code)),
                                                 shiny::selectInput('is.catch', label = shiny::strong("Data type"), c('Catch', 'Discard')),
                                                 shiny::checkboxInput('rem.ab', label = shiny::strong("Remove Values"), value = FALSE),
                                                 shiny::numericInput("lag.ab", "Observations:", 1, min = 1, max = 100),
                                                 shiny::checkboxInput('gen', label = shiny::strong("limit-axis"), value = TRUE),
                                                 shiny::downloadButton("DL_Abun", "Download")
                                             )),
                                      shiny::column(10,
                                             shiny::plotOutput('plot1', width = "100%", height = "800px")
                                             ))),
                         shiny::tabPanel('Compare',
                                  shiny::fluidRow(
                                      shiny::column(2,
                                             shiny::wellPanel(
                                                 shiny::selectInput('FISHC', 'Fishery :', as.character(fsh$Code)),
                                                 shiny::selectInput('is.CC', label = shiny::strong("Data type"), c('Catch', 'Discard')),
                                                 shiny::checkboxInput('rem.CC', label = shiny::strong("Remove Values"), value = FALSE),
                                                 shiny::numericInput("lag.CC", "Observations:", 1, min = 1, max = 100),
                                                 shiny::downloadButton("DL.cp.dat", "Download")
                                             )),
                                      shiny::column(10,
                                             shiny::plotOutput('plotC', width = "100%", height = "800px"),
                                             shiny::p(shiny::strong("\nModel skill assessment (quantitative metrics)")),
                                             DT::dataTableOutput('TabStat'),
                                             shiny::downloadButton("DL.cp.stat", "Download")
                                             ))),
                         ## -- Exit --
                         shiny::tabPanel(
                             shiny::actionButton("exitButton", "Exit")
                         )
                         ),
        ## Link the input for the different tabs with your original data
        ## Create the plots
        function(input, output, session){
            num <- shiny::reactive({
                num <- read.var(input$FG, nc.data, input$is.catch, grp)
                if(input$rem.ab){
                    rep <- rep(1, length(num[[1]]))
                    rep[c(1 : input$lag.ab)] <- NA
                    num <- lapply(num, function(x) x * rep)
                }
                num
            })
            bio <- shiny::reactive({
                cat.obs <- paste0(grp$Code, '_', input$is.CB, '_FC', which(fsh$Code == input$FISHB))
                cat.obs <- which(cat.obs %in% nam.var)
                arry    <- array(NA, dim = c(length(Time), length(cat.obs)))
                for(da in 1 : length(cat.obs)){
                    arry[, da] <- var.fish(grp$Code[cat.obs[da]], input$FISHB, nc.data, fsh = fsh, is.C = input$is.CB, grp = grp, by.box = FALSE)
                }
                colnames(arry) <- grp$Code[cat.obs]
                rem <- which(colSums(arry) == 0)
                if(length(rem) > 0){
                    arry <- arry[, -rem]
                    }
                if(input$rem.cb){
                    arry[c(1 : input$lag.cb), ] <- NA
                }
                if(input$b.year){
                    arry <- rowsum(arry, format(Time, '%Y'))
                }
                arry
            })
            ## External files
            ext <- shiny::reactive({
                cat.obs <- paste0(grp$Code, '_', input$is.CC, '_FC', which(fsh$Code == input$FISHC))
                cat.obs <- which(cat.obs %in% nam.var)
                C.arry  <- array(NA, dim = c(length(Time), length(cat.obs)))
                for(da in 1 : length(cat.obs)){
                    C.arry[, da] <- var.fish(grp$Code[cat.obs[da]], input$FISHC, nc.data, fsh = fsh, is.C = 'Catch', grp = grp, by.box = FALSE)
                }
                colnames(C.arry) <- grp$Code[cat.obs]

                if(input$rem.CC){
                    C.arry[c(1 : input$lag.CC), ] <- NA
                }

                rem <- which(colSums(C.arry, na.rm = TRUE) == 0)
                if(length(rem) > 0){
                    C.arry <- C.arry[, -rem]
                }
                C.arry <- rowsum(C.arry, format(Time, '%Y'), na.rm = TRUE)
                ext    <- ext.c[, c(1, which(ext.f %in% input$FISHC))]
                ext    <- list(A.catch = C.arry, external.c = ext)
                ext$Stats <- NULL
                for(i in  2 : ncol(ext$external.c)){
                   pos <- which(colnames(ext$A.catch) %in% colnames(ext$external.c)[i])
                   if(length(pos) == 0) next()
                   FG <- colnames(ext$external.c)[i]
                   na.rmv    <- which(!is.na(ext$external.c[, i]))
                   t.match   <- which(unique(format(Time, '%Y')) %in% ext$external.c$Time[na.rmv])
                   t.mat.ext <- which(ext$external.c$Time %in% unique(format(Time, '%Y'))[t.match])
                   if(sum(ext$A.catch[t.match, pos]) == 0 | sum(ext$external.c[t.mat.ext, i]) == 0) next()
                   ext$Stats <- rbind(ext$Stats, stats(as.vector(ext$A.catch[t.match, pos]), as.vector(ext$external.c[t.mat.ext, i]), FG))
                }
                ext
            })
            ## exit
            shiny::observeEvent(input$exitButton, {
                shiny::stopApp()
            })
            ## Plots
            output$plot1 <- shiny::renderPlot({
                rp           <- length(num())
                ylm          <- NULL
                if(input$gen) ylm <- range(sapply(num(), range, na.rm = TRUE))
                graphics::par(mfrow = c(t(grDevices::n2mfrow(rp))), cex = 1.2, oma = c(3, 3, 1, 1), cex = 1.1)
                for( i in 1 : rp){
                    plot.catch(num()[[i]], Time, ylm, coh = i, col.bi)
                }
                mtext("Time (days)", side = 1, outer = TRUE, cex = 2)
                mtext("Numbers", side = 2, outer = TRUE, line = 2, cex = 2)
            })
            output$plotB <- shiny::renderPlot({
                rp       <- ncol(bio())
                col.cur  <- col.bi[2]
                ylm      <- NULL
                if(input$is.CB == 'Discards'){
                    col.cur <- col.bi[1]
                }
                graphics::par(mfrow = c(t(grDevices::n2mfrow(rp))), cex = 1.2, oma = c(3, 3, 1, 1), cex = 1.1)
                for( i in 1 : rp){
                    plot.catch(bio()[, i], Time, ylm = NULL, coh = NULL, col.cur, bio.n = colnames(bio())[i], by.year = input$b.year)
                }
                mtext(2, text = 'Biomass (t)', outer = TRUE, cex = 2)
                mtext("Time (days)", side=1, line = 2, outer = TRUE, cex = 2)

            })
            output$plotC <- shiny::renderPlot({
                rp       <- ncol(ext()$A.catch)
                col.cur  <- col.bi[2]
                graphics::par(mfrow = c(t(grDevices::n2mfrow(rp))), cex = 1.2, oma = c(3, 3, 1, 1), cex = 1.1)
                for( i in 1 : rp){
                    plot.catch(ext()$A.catch[, i], Time, ylm = NULL, coh = NULL, col.cur, bio.n = colnames(ext()$A.catch)[i], by.year = TRUE, external = ext()$external.c)
                }
                mtext("Biomass", side=2, outer = TRUE, cex = 2)
                mtext("Time (days)", side=1, line = 2, outer = TRUE, cex = 2)
            })
            output$TabStat <- DT::renderDataTable(ext()$Stats)
            ## Save data
            output$DL_Biomass <- shiny::downloadHandler(
                filename = function(){
                    paste0(input$dataset, ".csv")
                },
                content = function(file) {
                    if(input$b.year){
                        Tim <- as.Date(unique(format(Time, '%Y')), format = '%Y')
                    }
                    write.csv(data.frame(Date = Tim, bio()), file, row.names = FALSE)
                }
            )
            output$DL_Abun <- shiny::downloadHandler(
                filename = function(){
                    paste0(input$dataset, ".csv")
                },
                content = function(file) {
                    out.abu <- matrix(unlist(num()), ncol = length(num()), byrow = TRUE)
                    colnames(out.abu) <- paste0(input$FG, '-Age', seq(1 : length(num())))
                    write.csv(data.frame(Date = Time, out.abu), file, row.names = FALSE)
                }
            )
             output$DL.cp.stat <- shiny::downloadHandler(
                filename = function(){
                    paste0(input$dataset, ".csv")
                },
                content = function(file) {
                    write.csv(data.frame(ext()$Stats), file, row.names = FALSE)
                }
            )
            output$DL.cp.dat <- shiny::downloadHandler(
                filename = function(){
                    paste0(input$dataset, ".csv")
                },
                content = function(file) {
                    out.ext <- ext()$external
                    out.int <- data.frame(Time = unique(format(Time, '%Y')), Simulated = ext()$A.catch)
                    diff    <- nrow(out.int) - nrow(out.ext)
                    if(diff > 0){
                        out.ext[c((nrow(out.ext) + 1) : (nrow(out.ext) + diff)), ] <- NA
                    } else if (diff < 0) {
                        out.int[c((nrow(out.int) + 1) : (nrow(out.int) + (diff *  - 1))), ] <- NA
                    }
                    out.f           <- data.frame(out.ext, out.int)
                    colnames(out.f) <- c(paste0('Shiny::Observed-', colnames(out.ext)), paste0('Simulated-', colnames(out.int)))
                    write.csv(out.f, file, row.names = FALSE)
                }
            )
        }
    )
}
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~          Internal Functions      ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##' @title Total catch and discard in biomass
##' @param FG Functional groups
##' @param FISH Fishery
##' @param nc.data NetCDF file with the information of catch
##' @param fsh Fisheries csv file (Input from Atlantis)
##' @param grp Groups csv file (Input from Atlantis)
##' @param is.C Bolean,  is Catch variable
##' @param by.box Default FALSE. If the information is needed by box (TRUE) or not (FALSE)
##' @return The information of the catch by biomass for all the FG
##' @author Demiurgo
var.fish <- function(FG, FISH, nc.data, fsh, grp, is.C = NULL, by.box = FALSE){
    pos     <- which(grp$Code == FG)
    f.pos   <- which(fsh$Code == FISH)
    name.fg <- paste0(grp$Code[pos], '_', is.C, '_FC', f.pos)
    B.catch <- ncdf4::ncvar_get(nc.data, name.fg)
    if(!by.box){
        B.catch <- colSums(B.catch, na.rm = TRUE)
    }
    return(B.catch)
}

##' @title Reading Discard and Catch by number and box
##' @param FG Csv with the Functional groups
##' @param nc.data NetCDF file with the information of catch
##' @param is.C Default TRUE,  if the data needed is Catch (TRUE) or Discard (FALSE)
##' @param grp Specific group to get from the csv file
##' @param by.box Default FALSE. If the information is needed by box (TRUE) or not (FALSE)
##' @return A list witht the information for all the functional groups
##' @author Demiurgo
read.var <- function(FG, nc.data, is.C = NULL, grp, by.box = FALSE){
    pos    <- which(grp$Code == FG)
    n.coh  <- grp$NumCohorts[pos]
    Tcatch <-list()
    if(n.coh > 1){
        for(coh in 1 : n.coh){
            name.fg <- paste0(grp$Name[pos], coh,'_', is.C)
            tmp    <- ncdf4::ncvar_get(nc.data, name.fg)
            if(!by.box){
                tmp <- colSums(tmp, na.rm=TRUE)
            }
            Tcatch[[coh]] <- tmp
        }
    } else {
        name.fg <- paste0(grp$Name[pos], ifelse(is.C, '_Catch', '_Discard'))
        tmp     <- ncdf4::ncvar_get(nc.data, name.fg)
        if(!by.box){
            tmp <- colSums(tmp, na.rm=TRUE)
        }
        Tcatch[[coh]] <- tmp
    }
    return(Tcatch)
}
##' @title Plot catch function
##' @param ctch Catch data.frame
##' @param Time Time on Date format
##' @param ylm Ylimit option for comparison
##' @param coh cohort number
##' @param col.bi Color of the lines
##' @param bio.n Name of the FG
##' @param by.year Plot and information is by year
##' @param external External file
##' @param ... Other options for plotting
##' @return plot of the catch
##' @author Demiurgo
plot.catch <- function(ctch, Time, ylm = NULL, coh = NULL, col.bi, bio.n = NULL, by.year = FALSE, external = NULL, ...){
    graphics::par(mar = c(1, 4, 3, 1) + 0.1)
    if(is.null(ylm)) ylm <- range(ctch, na.rm = TRUE)
    if(!is.null(external)){
        plp <- which(colnames(external) %in% bio.n)
        if(length(plp) == 0){  ## creating a fake variable to avoid problems
            external[, ncol(external) + 1] <- 0
            colnames(external)[ncol(external)] <- bio.n
            plp <- ncol(external)
            }
        ylm <- c(0, max(c(range(ctch, na.rm = TRUE), range(external[, plp], na.rm = TRUE))))

    }
    if(by.year){
        Time <- as.Date(unique(format(Time, '%Y')), format = '%Y')
    }
    if(is.null(bio.n)){
        main.lab <- paste0("AgeClass - ", ifelse(coh > 10, '', '0'), coh)
    } else{
        main.lab <- as.character(bio.n)
    }
#    yticks <- seq(ylm[1], round(ylm[2], 3), by = round(ylm[2], 3) / 5)
    tickx  <- seq(from = 1, to = length(Time), length = 5)
    plot(Time, ctch, type = 'l', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n',
         ylim = ylm, bty = 'n', lwd = 2, col = col.bi[1], main = main.lab)
    if(!is.null(external)){
        plp   <- which(colnames(external) %in% bio.n)
        Time2 <- as.Date(as.character(external$Time), '%Y')
        points(Time2, external[, plp], pch = 21, col = 'black', bg = 'red')
    }
    if(by.year){
        axis(1, at = Time[tickx], labels = format.Date(Time[tickx], '%Y'), lwd = 2,  col = 1)
    } else {
        axis(1, at = Time[tickx], labels = format.Date(Time[tickx], '%Y.%m'), lwd = 2,  col = 1)
    }
    axis(side = 2, at = pretty(ylm), lwd = 2)
}

##' @title Skill assessment of the model
##' @param obs Shiny::Observed values
##' @param mod modeled values
##' @param FG Functional group
##' @return metrics  =  AAE; AE; MEF; RMSE; COR
##' @author Demiurgo
stats <- function(obs, mod, FG){
    ## Stimation of Correlation
    COR  <- cor.test(obs, mod, method = 'spearman', use = "pairwise.complete.obs",  exact = FALSE)
    ## Average error
    AE   <- mean(obs, na.rm = TRUE) - mean(mod, na.rm = TRUE)
    ## Average absolute error
    diff <- mod - obs
    AAE  <- mean(abs(diff), na.rm = TRUE)
    ## Mean squared error
    RMSE <- sqrt(mean((diff) ^ 2, na.rm = TRUE))
    ## Reliability index
    ## Avoiding (inf values)
    tmp                   <- log(obs / mod) ^ 2
    tmp[is.infinite(tmp)] <- NA
    RI                    <-  exp(sqrt(mean(tmp, na.rm = TRUE)))
    ## Modeling efficiency
    ME <- 1 - (RMSE ^ 2) / var(obs, na.rm = TRUE) #option 1
    ## ME <- (sum((obs - mean(obs, na.rm = TRUE)) ^ 2, na.rm = TRUE) - sum((mod - obs) ^ 2, na.rm = TRUE)) / sum((obs - mean(obs, na.rm = TRUE)) ^ 2, na.rm = TRUE)
    if(COR$p.value == 0) COR$p.value <- '< 2.2e-16'
    out <- data.frame( FGroup = c(FG, NA, NA, NA, NA, NA),
                      Metrics = c('Correlation (Spearman)', 'Average Error (AE)',
                                  'Average Absolute Error (AAE)', 'Mean Squared Error (RMSE)', 'Reliability index',
                                  'Model Efficiency (ME)'),
                      Results = c(COR$estimate, AE, AAE, RMSE, RI, ME),
                      p.val = c(COR$p.value, NA, NA, NA, NA, NA))
    return(out)
}
