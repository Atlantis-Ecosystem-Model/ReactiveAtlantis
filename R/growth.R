##' This function helps you visualize what are the possible factors that limit the
##'     growth of primary producers. Primary producers in Atlantis are limited in
##'     their biomass by variables such as: nutrients, light, space, and by processes
##'     such as acidification and the presence of eddies (Audijyontze et al. 2017). In this
##'     function, most of these factors are analyzed, except for the acidification
##'     process and the spatial limitation.
##' \itemize{
##'    \item \bold{Nutrient limitation}: The nutrients limitation scalar
##'     (\eqn{\delta_{nutrients}}) in the case of species limited only by nitrogen is
##'     calculated as:
##'     \deqn{\delta_{nutrients}  = \frac{DIN}{KN + DIN}} \cr
##'     were \eqn{DIN} is the concentration of \eqn{NH_{3}} + \eqn{NO_{3}} and
##'     \eqn{KN} is the half-saturation constant of nutrient uptake. \cr
##'     For species that are limited for multiple nutrients like diatoms that can be
##'     limited by silicate \eqn{Si} and nitrogen you have 3 different options to choose on
##'     Atlantis (\emph{flagnut}) and the nutrient limitation can have different options:
##'     \itemize{
##'       \item \emph{Leibig limitation \bold{flagnut = 0}}
##'       \deqn{\delta_{nutrients} = min(\frac{DIN}{KN + DIN}, \frac{Si}{KS + Si})}
##'       \item \emph{multiplicative limitation \bold{flagnut = 1}}
##'       \deqn{\delta_{nutrients} = \sqrt{\frac{DIN}{KN + DIN} * \frac{Si}{KS + Si}}}
##'       \item \emph{ERSEM WQI limitation \bold{flagnut = 2}}
##'       \deqn{\delta_{nutrients} = \frac{2}{\frac{DIN}{KN + DIN} * \frac{Si}{KS + Si}}}
##'     were \eqn{KS} is the growth half-saturation constant for the functional group.
##'        }
##'   \item \bold{Light limitation}: In Atlantis photosynthesis in  primary producers
##'     can be limited by light. The light limitation factor \eqn{\delta_{light}} is calculated by:
##'  \deqn{\delta_{light} = min(\frac{IRR}{KI}, 1)} \cr
##'  were \eqn{IRR} is irradiance or available light  and\eqn{KI} is the light
##'     saturation coefficient. \cr
##' In Atlantis there is a second option \emph{flaglight = 1} that allows the calculation of the light
##'     limitation factor allowing for light adaptation. This is intended to capture
##'     the ability of the primary producers to rapidly adapt to different light
##'     conditions (Audijyontze et al. 2017). This option was not considered in this version
##'     of the tool but would be implemented in future updates.
##'  \item \bold{Effect of eddies on primary production}: The eddy effect on primary
##'     produces \eqn{\delta_{eddy}}  is calculated by multiplying the scale parameter \eqn{eddy_{scale}}
##'     (\emph{eddy_scale} on Atlantis) by the eddy strength \eqn{eddy_{strength}}.
##'    \deqn{\delta_{eddy} = eddy_{scale} * eddy_{strength} }
##' }
##' @title Limitation Growth primary producers
##' @param ini.nc.file Character string with the path to the netcdf file to
##'     read in. This netcdf file contains the initial conditions for the Atlantis
##'     model usually ends in \emph{.nc}.
##' @param grp.file Character string with the path to the Groups \emph{*.csv}
##'     file (Atlantis input file).
##' @param prm.file Character string with the path to the biology parameter
##'     file \emph{*.prm}.
##' @param out.nc.file Character string with the path to the netcdf file to
##'     read in. This netcdf file is a generic output from an Atlantis run and
##'     usually starts with output and ends in \emph{.nc}.
##' @return  A reactive HTML with graphical output by functional group, for layer and box
##'     of the following variables:
##' \itemize{
##'   \item \bold{Growth}:  This variable is the  Primary producer growth (\eqn{G_{pp}}) which is
##'     determined by the biomass for the primary producer (\eqn{Biom_{pp}}),  the maximum effective
##'     growth (\eqn{mum}) and the limiting factors: Light limitation (\eqn{\delta_{light}}),
##'     Nutrients limitation \eqn{\delta_{nutrients}}, and Eddy scalar (\eqn{\delta_{eddy}}).
##'     \deqn{G_{pp} = Biom_{pp} * mum * \delta_{light} * \delta_{nutrients} *
##'     \delta_{eddy}}
##' \item \bold{Light limitation}: The available light by box and layer. This option
##'     uses only the basic formula (\eqn{\delta_{light}}, described before) and currently does not use the light adaptation
##'     option. The light attenuation values used for the calculations are the values
##'     estimated by Atlantis and reported in the \bold{\emph{out.nc.file}}.
##' \item \bold{Nutrient limitation}: Nutrient limitation is calculated by layer and
##'     by box, using the option chosen by the user in the Atlantis run
##'     (\emph{flagnut = 0, 1 or 2}).
##' \item \bold{Eddy scalar}: This value is calculated by functional group and by
##'     box. The values of eddies strength (\eqn{eddy_{strength}}) used are from the
##'     Atlantis output and the value of eddy scale (\eqn{eddy_{scale}}) is read from
##'     the Atlantis configuration file.
##'   }
##' @author Demiurgo
##' @export
growth.pp <- function(ini.nc.file, grp.file, prm.file, out.nc.file){
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
    library(RColorBrewer)
    color <- palette(brewer.pal(n = 12, name = "Paired"))

    nc.ini    <- nc_open(ini.nc.file)
    group.csv <- read.csv(grp.file)
    is.off    <- which(group.csv$IsTurnedOn == 0)
    nc.out    <- nc_open(out.nc.file)
    prm       <- readLines(prm.file, warn = FALSE)
    ## Selecting primary producers
    pp.grp    <- with(group.csv, which(GroupType %in% c('PHYTOBEN', 'SM_PHY', 'LG_PHY', 'SEAGRASS', 'DINOFLAG'))) ## Just primary producer
    pp.grp    <- pp.grp[ - which(pp.grp %in% is.off)]
    cod.fg    <- with(group.csv, as.character(Code[pp.grp])) ## cod
    nam.fg    <- with(group.csv, as.character(Name[pp.grp])) ## name
    options(warn =  - 1)
    flagnut   <- text2num(prm, 'flagnut ', FG = 'look')
    flaglight <- text2num(prm, 'flaglight ', FG = 'look')
    ##Parameters needed
    mum <- sp.dep <- KI <- KS <- KN <- NULL
    Kiop.min <- Kiop.shift <- Ki.avail <- K.depth <- NULL
    for( i in 1: length(pp.grp)) {
        KN     <- rbind(KN, text2num(prm, paste0('KN_', cod.fg[i]), FG = cod.fg[i]))
        mum    <- rbind(mum, text2num(prm, paste0('mum_', cod.fg[i], '_T15'), FG = cod.fg[i]))
        KS     <- rbind(KS, text2num(prm, paste0('KS_', cod.fg[i]), FG = cod.fg[i]))
        KI     <- rbind(KI, text2num(prm, paste0('KI_', cod.fg[i]), FG = cod.fg[i]))
        sp.dep <- rbind(sp.dep, text2num(prm, paste0('habitat_', cod.fg[i]), FG = cod.fg[i], Vector = TRUE))
        if(flaglight$Value == 1){
            Kiop.min   <-rbind(Kiop.min, text2num(prm, paste0('KIOP_min', cod.fg[i]), FG = cod.fg[i]))
            Kiop.shift <-rbind(Kiop.shift, text2num(prm, paste0('KIOP_shift', cod.fg[i]), FG = cod.fg[i]))
            Ki.avail   <-rbind(Ki.avail, text2num(prm, paste0('KIOP_avail', cod.fg[i]), FG = cod.fg[i]))
            K.depth    <-rbind(K.depth, text2num(prm, paste0('KIOP_addepth', cod.fg[i]), FG = cod.fg[i]))
        }
    }
    ed.scl   <- text2num(prm, 'eddy_scale', FG = 'look')[, 2]
    options(warn =  0)
    ## nutrient
    DIN      <- ncvar_get(nc.out, 'NO3')  * ncvar_get(nc.out, 'NH3')
    Si       <- ncvar_get(nc.out, 'Si')
    l.eddy   <- ncvar_get(nc.out, 'eddy') * ed.scl
    nlayers  <- ncvar_get(nc.out, 'numlayers')[, 1]
    IRR      <- ncvar_get(nc.out, 'Light')
    l.space  <- l.light <- l.nut <- biom<- list()
    ## type of sediments
    for(fg in 1 : length(nam.fg)){
        biom[[fg]] <- ncvar_get(nc.out, paste0(nam.fg[fg], '_N'))
        ## nutrients
        if(flagnut$Value == 0){
            l.nut[[fg]] <- pmin((DIN / (KN$Value[fg] + DIN)) ,  (Si /(KS$Value[fg] + Si)), na.rm = TRUE)
        } else if(flagnut$Value == 1){
            l.nut[[fg]] <- sqrt((DIN / (KN$Value[fg] + DIN)) * (Si /(KS$Value[fg] + Si)))
        } else if(flagnut$Value == 2){
            l.nut[[fg]] <- 2 / ((DIN / (KN$Value[fg] + DIN)) + (Si /(KS$Value[fg] + Si)))
        } else if(is.null(flagnut$Value )){
            l.nut[[fg]] <- (DIN / (KN$Value[fg] + DIN))
        } else {
            stop('No such nut_case defined: ', KN$Value[fg],  '- value must be between 0 and 3 currently\n')
        }
        ## Light
        if(flaglight$Value == 0){
            l.light[[fg]] <- pmin((IRR / KI$Value[fg]) , 1, na.rm = TRUE)
        } else if(flagnut$Value > 0){
            ## not implemented
            stop('Option: ', flagnut$Value, 'for light limitation, not implemented yet\n')
        }

        ## space limitation for Phytobentos and seagrass
        if(group.csv$GroupType[pp.grp[fg]] %in% c('PHYTOBEN', 'SEAGRASS')){
            ratio <- 1 # I used this value just to avoid the calculation of area
            SPmax   <- text2num(prm, paste0(cod.fg[fg], 'max'), FG = 'look')[, 2]
            tmp     <- biom[[fg]] / (SPmax * ratio)
            ## removing nan and non - finite numbers
            tmp[!is.finite(tmp)] <- 0
            l.space[[fg]] <- tmp
        }
    }
    shinyApp(
        ui <- navbarPage('Growth primary producers',
                         tabPanel('Growth  - Limiting factors',
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 tags$h3('Functional Group'),
                                                 selectInput('sp', 'Functional Group', as.character(cod.fg)),
                                                 selectInput('box', 'Box', seq(0, (length(nlayers) - 1))),
                                                 checkboxInput("log", "Logarithmic scale", FALSE)
                                             )
                                             ),
                                      column(10,
                                             plotOutput('plot1', width = "100%", height = "700px")
                                             )
                                  )
                                  ),
                         ## -- Exit --
                         tabPanel(
                             actionButton("exitButton", "Exit")
                         )
                         ),
        function(input, output, session){
            p.fg <- reactive({which(cod.fg %in% input$sp)})
            lay  <- reactive({
                if(is.null(l.space[[p.fg()]])){
                    ((max(nlayers) + 1) - nlayers[as.numeric(input$box) + 1]): (max(nlayers) + 1)
                } else {
                    (((max(nlayers) + 1) - nlayers[as.numeric(input$box) + 1]) : (max(nlayers) + 1))[1]
                }
            })
            name.lay <- reactive({
                if(length(lay()) == 1){
                    c('Sediment')
                } else {
                    c(paste0('Layer_', seq(from = 0, to = (length(lay()) - 2))), 'Sediment')
                }
            })

            lim.nut   <- reactive({l.nut[[p.fg()]][lay(), as.numeric(input$box) + 1, ]})
            lim.light <- reactive({l.light[[p.fg()]][lay(), as.numeric(input$box) + 1, ]})
            lim.eddy  <- reactive({l.eddy[as.numeric(input$box) + 1, ]})
            step1     <- reactive({mum$Value[p.fg()] * lim.nut() * lim.light()})
            step2     <- reactive({
                if(!is.null(l.space[[p.fg()]])){
                    temp1 <- unlist(step1() * unlist(biom[[p.fg()]][as.numeric(input$box) + 1, ]) * lim.eddy())
                    data.frame(temp1)
                } else {
                    temp <- step1() * biom[[p.fg()]][lay(), as.numeric(input$box) + 1, ]
                    apply(temp, 1, FUN = function(x) x * lim.eddy())
                }
            })
            growth.fin    <- reactive({
                if(isTRUE(input$log) == TRUE){
                    if(all(step2() == 0)){
                        step2()
                    } else {
                        log(step2())
                    }
                }else{
                    step2()}
            })
            lim.nut.fin   <- reactive({
                if(isTRUE(input$log) == TRUE){
                    log(lim.nut())
                } else {
                    lim.nut()}
            })
            lim.light.fin <- reactive({
                if(isTRUE(input$log) == TRUE){
                    log(lim.light())
                } else {
                    lim.light()}
            })
            lim.eddy.fin  <- reactive({
                if(isTRUE(input$log) == TRUE){
                    log(lim.eddy())
                } else {
                    lim.eddy()}
            })
            ## STOP
            observeEvent(input$exitButton, {stopApp()})
            output$plot1 <- renderPlot({
               par(mfrow = c(2, 2), mar = c(0, 3, 1, 5.1), oma = c(4, 4, 0.5, 2), xpd = TRUE)
                ## growth
                ran <- range(unlist(growth.fin()), finite = TRUE)
                plot(growth.fin()[, 1], axes = FALSE, ylim = ran, bty = 'n', type = 'l',
                     lty = 1, pch = 19, col = color[1], ylab = '', main = 'Growth')
                axis(2, at = round(seq(ran[1], ran[2], length.out = 5), 3), las = 1, line =  - .7)
                if(ncol(growth.fin()) > 1){
                    for( j in 2 : ncol(growth.fin())){
                        lines(growth.fin()[, j], type = 'l', pch = 19, lty = 1, col = color[j])
                    }
                }
                ## Nutrients
                ran <-range(lim.nut.fin(), finite = TRUE)
                if(is.null(dim(lim.nut.fin()))){
                    plt.nut <- matrix(lim.nut.fin(), nrow = 1)
                } else {
                    plt.nut <- lim.nut.fin()
                }
                plot(plt.nut[1, ], axes = FALSE, ylim = ran, bty = 'n', type = 'l',
                     lty = 1, pch = 19, col = color[1], ylab = '', main = 'Nutrients limitation')
                axis(2, at = round(seq(ran[1], ran[2], length.out = 5), 3), las = 1, line =  - .7)
                if(nrow(plt.nut) > 1){
                    for( j in 2 : nrow(plt.nut)){
                        lines(plt.nut[j, ], type = 'l', pch = 19, lty = 1, col = color[j])
                    }
                }
                ## light
                ran <- range(lim.light.fin(), finite = TRUE)
                if(is.null(dim(lim.light.fin()))){
                    plt.light <- matrix(lim.light.fin(), nrow = 1)
                } else {
                    plt.light <- lim.light.fin()
                }
                plot(plt.light[1, ], yaxt = 'n', ylim = ran, bty = 'n', type = 'l',
                     lty = 1, pch = 19, col = color[1], ylab = '', main = 'Light limitation')
                axis(2, at = round(seq(ran[1], ran[2], length.out = 5), 3), las = 1, line =  - .7)
                if(nrow(plt.light) > 1){
                    for( j in 2 : nrow(plt.light)){
                        lines(plt.light[j, ], type = 'l', pch = 19, lty = 1, col = color[j])
                    }
                }
                ## eddyes
                ran <- range(lim.eddy.fin(), finite = TRUE)
                plot(lim.eddy.fin(), yaxt = 'n', ylim = ran, bty = 'n', type = 'l',
                          lty = 1, pch = 19, col = 'orangered2', ylab = '', main = 'Eddy scalar')
                axis(2, at = round(seq(ran[1], ran[2], length.out = 5), 3), las = 1, line =  - .7)
                legend("topright", inset = c(-0.1, 0), legend= c(paste0('Layer - ', 1 : (ncol(growth.fin())-1)), 'Sediment', 'all-Layers'),
                       fill = c(color[1 : ncol(growth.fin())], 'orangered2'), title = 'Water Layers', bty = 'n')
            })
        }
    )
}

## ## functions
## ##' @title Parameter file reader
## ##' @param text Biological parametar file for Atlantis
## ##' @param pattern Text that you are looking
## ##' @param FG Name of the functional groups
## ##' @param Vector Logic argument, if the data is on vectors or not
## ##' @return A matrix with the values from the .prm file
## ##' @author Demiurgo
## text2num <- function(text, pattern, FG = NULL, Vector = FALSE){
##     if(!isTRUE(Vector)){
##         text <- text[grep(pattern = pattern, text)]
##         txt  <- gsub(pattern = '[[:space:]]+' ,  '|',  text)
##         col1 <- col2 <- vector()
##         for( i in 1 : length(txt)){
##             tmp     <- unlist(strsplit(txt[i], split = '|', fixed = TRUE))
##             if(grepl('#', tmp[1])) next
##             tmp2    <- unlist(strsplit(tmp[1], split = '_'))
##             if(FG[1] == 'look') {
##                 col1[i] <- tmp2[1]
##             } else {
##                 id.co   <- which(tmp2 %in% FG )
##                 if(sum(id.co) == 0) next
##                 col1[i] <- tmp2[id.co]
##             }
##             col2[i] <- as.numeric(tmp[2])
##         }
##         if(is.null(FG)) col1 <- rep('FG', length(col2))
##         out.t <- data.frame(FG = col1, Value = col2)
##         if(any(is.na(out.t[, 1]))){
##             out.t <- out.t[-which(is.na(out.t[, 1])), ]
##         }
##         return(out.t)
##     } else {
##         l.pat <- grep(pattern = pattern, text)
##         if(length(l.pat) == 0){
##             warning(paste('you dont have values for the parameter', pattern, sep = ' '))

##         } else {
##             nam   <- gsub(pattern = '[ ]+' ,  '|',  text[l.pat])
##             fg    <- vector()
##             pos   <- 1
##             for( i in 1 : length(nam)){
##                 tmp     <- unlist(strsplit(nam[i], split = '|', fixed = TRUE))
##                 if(grepl('#', tmp[1])) next
##                 fg[pos] <- tmp[1]
##                 if(pos == 1) {
##                     t.text <- gsub('"[[:space:]]"', ' ',  text[l.pat[i] + 1])
##                     pp.mat <- matrix(as.numeric(unlist(strsplit(t.text, split = ' +', fixed = FALSE))), nrow = 1)
##                     pos    <- pos + 1
##                 } else {
##                     t.text <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", text[l.pat[i] + 1], perl=TRUE)
##                     pp.tmp <- matrix(as.numeric(unlist(strsplit(t.text, split = ' ', fixed = TRUE))), nrow = 1)
##                     if(ncol(pp.mat) != ncol(pp.tmp)) stop('\nError: The pPrey vector for ', tmp[1], ' has ', ncol(pp.tmp), 'columns and should have ', ncol(pp.mat))
##                     pp.mat <- rbind(pp.mat, pp.tmp)
##                     pos    <- pos + 1
##                 }
##             }
##             if(all(is.na(pp.mat[, 1]))) pp.mat <- pp.mat[, - 1]
##             row.names(pp.mat)                  <- fg
##             return(pp.mat)
##         }
##         return()
##     }
## }

##' @title Parameter file reader
##' @param text Biological parametar file for Atlatnis
##' @param pattern Text that you are looking
##' @param FG Name of the functional groups
##' @param Vector Logic argument, if the data is on vectors or not
##' @return A matrix with the values from the .prm file
##' @author Demiurgo
text2num <- function(text, pattern, FG = NULL, Vector = FALSE){
    if(!isTRUE(Vector)){
        text <- text[grep(pattern = pattern, text)]
        txt  <- gsub(pattern = '[[:space:]]+' ,  '|',  text)
        col1 <- col2 <- vector()
        for( i in 1 : length(txt)){
            tmp     <- unlist(strsplit(txt[i], split = '|', fixed = TRUE))
            if(grepl('#', tmp[1])) next
            tmp2    <- unlist(strsplit(tmp[1], split = '_'))
            if(FG[1] == 'look') {
                col1[i] <- tmp2[1]
            } else {
                id.co   <- which(tmp2 %in% FG )
                if(sum(id.co) == 0) next
                col1[i] <- tmp2[id.co]
            }
            col2[i] <- as.numeric(tmp[2])
        }
        if(is.null(FG)) col1 <- rep('FG', length(col2))
        out.t <- data.frame(FG = col1, Value = col2)
        if(any(is.na(out.t[, 1]))){
            out.t <- out.t[-which(is.na(out.t[, 1])), ]
        }
        return(out.t)
    } else {
        l.pat <- grep(pattern = pattern, text)
        nam   <- gsub(pattern = '[ ]+' ,  '|',  text[l.pat])
        fg    <- vector()
        pos   <- 1
        for( i in 1 : length(nam)){
            tmp     <- unlist(strsplit(nam[i], split = '|', fixed = TRUE))
            if(grepl('#', tmp[1]) || !grepl('^pPREY', tmp[1])) next
            fg[pos] <- tmp[1]
            if(pos == 1) {
                t.text <- gsub('"[[:space:]]"', ' ',  text[l.pat[i] + 1])
                pp.mat <- matrix(as.numeric(unlist(strsplit(t.text, split = ' +', fixed = FALSE))), nrow = 1)
                pos    <- pos + 1
            } else {
                t.text <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", text[l.pat[i] + 1], perl=TRUE)
                pp.tmp <- matrix(as.numeric(unlist(strsplit(t.text, split = ' ', fixed = TRUE))), nrow = 1)
                if(ncol(pp.mat) != ncol(pp.tmp)) stop('\nError: The pPrey vector for ', tmp[1], ' has ', ncol(pp.tmp), 'columns and should have ', ncol(pp.mat))
                pp.mat <- rbind(pp.mat, pp.tmp)
                pos    <- pos + 1
            }
        }
        if(all(is.na(pp.mat[, 1]))) pp.mat <- pp.mat[, - 1]
        row.names(pp.mat)                  <- fg
        return(pp.mat)
    }
}
