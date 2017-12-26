##' @title Analysis of the Harvest in Atlatnis
##' @param grp.csv Character string with the connection to the Groups \code{*.csv} file (Atlantis input file).
##' @param fish.csv Character string with the connection to the fisheries \code{*.csv} file (Atlantis input file)
##' @param catch.nc Character string with the connection to the catch netcdf output file from Atlantis. Usually with the name of \code{[Your_Model]CATCH.nc}, where [Your_Model] is the name of your Atlantis model
##' @param ext.catch.f (Default = NULL) Character string with the connection to the External File with the Observed catches and discards by year. This helps to calibrate the harvest section of atlatnis.
##' @return A shiny output (reactive html)
##' @author Demiurgo
##' @export
catch <- function(grp.csv, fish.csv, catch.nc, ext.catch.f = NULL){
## Libraries
    if (!require('shiny', quietly = TRUE)) {
        stop('The package shiny was not installed')
    }
    if (!require('ncdf4', quietly = TRUE)) {
        stop('The package ncdf4 was not installed')
    }
    if (!require('reshape', quietly = TRUE)) {
        stop('The package reshape was not installed')
    }
    if (!require('tidyverse', quietly = TRUE)) {
        stop('The package tidyverse was not installed')
    }
    if (!require('stringr', quietly = TRUE)) {
        stop('The package stringr was not installed')
    }
    if (!require('data.table', quietly = TRUE)) {
        stop('The package data.table was not installed')
    }
    if (!require('RColorBrewer', quietly = TRUE)) {
        stop('The package RColorBrewer was not installed')
    }
    ## General setting
    mycol    <- c(brewer.pal(8, "Dark2"), c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
    mycol    <- colorRampPalette(mycol)
    col.bi   <- mycol(15)[c(14, 13, 12, 10)]
    nc.data  <- nc_open(catch.nc)
    #nc.out   <- nc_open(out.nc)
    grp      <- read.csv(grp.csv)
    if(!is.null(ext.catch.f)){
        ext.f  <- scan(ext.catch.f, nlines = 1, what = character(), sep = ',')
        ext.c  <- read.table(ext.catch.f, skip = 1, header = TRUE, sep = ',')
    }
    grp      <- grp[grp$IsImpacted == 1, ]
    fsh      <- read.csv(fsh.csv)
    nam.var  <- names(nc.data$var)
    orign    <- unlist(strsplit(ncatt_get(nc.data, 't')$units, ' ', fixed = TRUE))
    if(orign[1] == 'seconds') {
        Time <- ncvar_get(nc.data, 't') / 86400
    } else {
        Time <- ncvar_get(nc.data, 't')
    }
    Time     <- as.Date(Time, origin = orign[3])
    shinyApp(
        ## Create the different tabs
        ui <- navbarPage("Catch",
                         tabPanel('Biomass',
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 selectInput('FISHB', 'Fishery :', as.character(fsh$Code)),
                                                 selectInput('is.CB', label = strong("Data type"), c('Catch', 'Discard')),
                                                 checkboxInput('b.year', label = strong("By year"), value = FALSE)
                                                 )),
                                      column(10,
                                             plotOutput('plotB', width = "100%", height = "700px")
                                             ))),
                         tabPanel('Numbers',
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 selectInput('FG', 'Functional Group :', as.character(grp$Code)),
                                                 selectInput('is.catch', label = strong("Data type"), c('Catch', 'Discard')),
                                                 checkboxInput('gen', label = strong("limit-axis"), value = TRUE)
                                             )),
                                      column(10,
                                             plotOutput('plot1', width = "100%", height = "700px")
                                             ))),
                         tabPanel('Compare',
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 selectInput('FISHC', 'Fishery :', as.character(fsh$Code)),
                                                 selectInput('is.CC', label = strong("Data type"), c('Catch', 'Discard'))
                                                 )),
                                      column(10,
                                             plotOutput('plotC', width = "100%", height = "700px"),
                                             p(strong("\nModel skill assessment (quantitative metrics)")),
                                             dataTableOutput('TabStat')
                                             ))),
                         ## -- Exit --
                         tabPanel(
                             actionButton("exitButton", "Exit")
                         )
                         ),
        ## Link the input for the different tabs with your original data
        ## Create the plots
        function(input, output, session){
            num <- reactive({
                num <- read.var(input$FG, nc.data, input$is.catch, grp) #, input$is.box)
            })
            bio <- reactive({
                cat.obs <- paste0(grp$Code, '_', input$is.CB, '_FC', which(fsh$Code == input$FISHB))
                cat.obs <- which(cat.obs %in% nam.var)
                arry    <- array(NA, dim = c(length(Time), length(cat.obs)))
                for(da in 1 : length(cat.obs)){
                    arry[, da] <- var.fish(grp$Code[cat.obs[da]], input$FISHB, nc.data, fsh = fsh, is.C = 'Catch', grp = grp, by.box = FALSE)
                }
                colnames(arry) <- grp$Code[cat.obs]
                if(input$b.year){
                    arry <- rowsum(arry, format(Time, '%Y'))
                }
                arry
            })
            ## External files
            ext <- reactive({
                cat.obs <- paste0(grp$Code, '_', input$is.CC, '_FC', which(fsh$Code == input$FISHC))
                cat.obs <- which(cat.obs %in% nam.var)
                C.arry  <- array(NA, dim = c(length(Time), length(cat.obs)))
                for(da in 1 : length(cat.obs)){
                    C.arry[, da] <- var.fish(grp$Code[cat.obs[da]], input$FISHC, nc.data, fsh = fsh, is.C = 'Catch', grp = grp, by.box = FALSE)
                }
                colnames(C.arry) <- grp$Code[cat.obs]
                C.arry <- rowsum(C.arry, format(Time, '%Y'))
                ext    <- ext.c[, c(1, which(ext.f %in% input$FISHC))]
                ext    <- list(A.catch = C.arry, external.c = ext)
                for(i in  2 : ncol(ext$external.c)){
                   pos <- which(colnames(ext$A.catch) %in% colnames(ext$external.c)[i])
                   if(length(pos) == 0) next()
                   na.rmv    <- which(!is.na(ext$external.c[, i]))
                   t.match   <- which(unique(format(Time, '%Y')) %in% ext$external.c$Time[na.rmv])
                   t.mat.ext <- which(ext$external.c$Time %in% unique(format(Time, '%Y'))[t.match])
                   ext$Stats <- stats(as.vector(ext$A.catch[t.match, pos]), as.vector(ext$external.c[t.mat.ext, i]))
                }
                ext
            })
            ## exit
            observeEvent(input$exitButton, {
                stopApp()
            })
            output$plot1 <- renderPlot({
                rp           <- length(num())
                ylm          <- NULL
                if(input$gen) ylm <- range(sapply(num(), range, na.rm = TRUE))
                par(mfrow = c(t(n2mfrow(rp))), cex = 1.2, oma = c(1, 1, 1, 1))
                for( i in 1 : rp){
                    plot.catch(num()[[i]], Time, ylm, coh = i, col.bi)
                }
            })
            output$plotB <- renderPlot({
                rp       <- ncol(bio())
                col.cur  <- col.bi[2]
                ylm      <- NULL
                if(input$is.CB == 'Discards'){
                    col.cur <- col.bi[1]
                }
                par(mfrow = c(t(n2mfrow(rp))), cex = 1.2, oma = c(1, 1, 1, 1))
                for( i in 1 : rp){
                    plot.catch(bio()[, i], Time, ylm = NULL, coh = NULL, col.cur, bio.n = colnames(bio())[i], by.year = input$b.year)
                }
            })
            output$plotC <- renderPlot({
                rp       <- ncol(ext()$A.catch)
                col.cur  <- col.bi[2]
                par(mfrow = c(t(n2mfrow(rp))), cex = 1.2, oma = c(1, 1, 1, 1))
                for( i in 1 : rp){
                    plot.catch(ext()$A.catch[, i], Time, ylm = NULL, coh = NULL, col.cur, bio.n = colnames(ext()$A.catch)[i], by.year = TRUE, external = ext()$external.c)
                }
            })
            output$TabStat <- renderDataTable(ext()$Stats)
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
##' @param catch Default TRUE,  if the data needed is Catch (
##' @param fsh Fisheries csv file (Input from Atlantis)
##' @param grp Groups csv file (Input from Atlantis)
##' @param catch Default TRUE,  if the data needed is Catch (TRUE) or Discard (FALSE)
##' @param by.box Default FALSE. If the information is needed by box (TRUE) or not (FALSE)
##' @return The information of the catch by biomass for all the FG
##' @author Demiurgo
var.fish <- function(FG, FISH, nc.data, fsh, grp, is.C = NULL, by.box = FALSE){
    pos     <- which(grp$Code == FG)
    f.pos   <- which(fsh$Code == FISH)
    name.fg <- paste0(grp$Code[pos], '_', is.C, '_FC', f.pos)
    B.catch <- ncvar_get(nc.data, name.fg)
    if(!by.box){
        B.catch <- colSums(B.catch, na.rm = TRUE)
    }
    return(B.catch)
}

##' @title Reading Discard and Catch by number and box
##' @param FG Csv with the Functional groups
##' @param nc.data NetCDF file with the information of catch
##' @param is.C Default TRUE,  if the data needed is Catch (TRUE) or Discard (FALSE)
##' @param grp.csv Groups csv file (Input from Atlantis)
##' @param by.box Default FALSE. If the information is needed by box (TRUE) or not (FALSE)
##' @return A list witht the information for all the functional groups
##' @author Demiurgo
read.var <- function(FG, nc.data, is.C = NULL, grp, by.box = FALSE){
    pos    <- which(grp$Code == FG)
    n.coh  <- grp$NumCohorts[pos]
    Tcatch <-list()
    if(n.coh > 1){
        for(coh in 1 : n.coh){
#            browser()
            name.fg <- paste0(grp$Name[pos], coh,'_', is.C)
            tmp    <- ncvar_get(nc.data, name.fg)
            if(!by.box){
                tmp <- colSums(tmp, na.rm=TRUE)
            }
            Tcatch[[coh]] <- tmp
        }
    } else {
        name.fg <- paste0(grp$Name[pos], ifelse(is.C, '_Catch', '_Discard'))
        tmp     <- ncvar_get(nc.data, name.fg)
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
    par(mar = c(1, 4, 1, 1) + 0.1)
    if(is.null(ylm)) ylm <- range(ctch)
    if(!is.null(external)){
        plp <- which(colnames(external) %in% bio.n)
        ylm <- c(0, max(c(range(ctch, na.rm = TRUE), range(external[, plp], na.rm = TRUE))))
    }
    if(by.year){
        Time <- as.Date(unique(format(Time, '%Y')), format = '%Y')
    }
    yticks <- seq(ylm[1], round(ylm[2], 3), by = round(ylm[2], 3) / 5)
    tickx  <- seq(from = 1, to = length(Time), length = 5)
    plot(Time, ctch, type = 'l', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n',
         ylim = ylm, bty = 'n', lwd = 2, col = col.bi[1])
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
    axis(2, at = yticks , labels = format(yticks, digits = 2), lwd = 2)
    mtext("Time (days)", side=1, line = 3)
    if(is.null(bio.n)){
        mtext(2, text = paste0('(Number)'), line = 3)
        mtext(paste0("AgeClass - ", ifelse(coh > 10, '', '0'), coh))
    } else {
        mtext(bio.n)
        mtext(2, text = 'Biomass (t)', line = 3)
    }
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Skill assessment of the model
##' @param obs Observed values
##' @param mod modeled values
##' @return metrics  =  AAE; AE; MEF; RMSE; COR
##' @author Demiurgo
stats <- function(obs, mod){
    ## Stimation of Correlation
    COR  <- cor(obs, mod, method = 'spearman')
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
    ME <- 1 - (RMSE ^ 2) / (var(obs, na.rm = TRUE) ^ 2)
    out <- data.frame(Metrics = c('Correlation (Spearman)', 'Average Error (AE)',
                                  'Average Absolute Error (AAE)', 'Mean Squared Error (RMSE)', 'Reliability index',
                                  'Model Efficiency (ME)'),
                      Results = c(COR, AE, AAE, RMSE, RI, ME))
    return(out)
}
