##' This function allows you to explore the changes in the food web structure and
##'     trophic level for an specific functional group through time. The estimation
##'     of the trophic position is based on the approach of Pauly \emph{et al.}(1998)
##'     \deqn{TL_{i}  = 1  +  \frac{\sum_{j = 1}^{n}(DC_{i, j}  *  TL_{j})}{\sum_{j = 1}^{n} DC_{i, j}}}
##' were \eqn{TL_{i}} in the trophic position of the predator, \eqn{DC_{i, j}} is the diet composition; \eqn{TL_{j}} is the trophic level of the
##'     prey,  and n in the number of preys. These initial values for the trophic position were based on Pauly \emph{et
##'     al.} (1998) and Tucker and Roger (2014) (Table 1).
##' \cr
##' Table 1 : Initial trophic position Assigned for each functional groups.
##'   \tabular{lc}{
##'     \bold{Functional Group} \tab \bold{Trophic level} \cr
##'     Small Phytoplankton   \tab 1.2   \cr
##'     Large Phytoplankton   \tab 2   \cr
##'     Coral    \tab 2.2 \cr
##'     Small zooplankton   \tab 2.4 \cr
##'     Medium zooplankton   \tab 2.5 \cr
##'     Large zooplankton   \tab 2.7 \cr
##'     Benthic invertebrates \tab 2.52 \cr
##'     Cephalopods   \tab 3.2 \cr
##'     Fishes   \tab 3.24 \cr
##'     Birds   \tab 3.5 \cr
##'     Shark and Mammals   \tab 4.0 \cr
##'  }
##'     The flexibility of this function allows you to change: a) \bold{the \emph{focus}
##'     functional group} over which all calculations would be made; b) \bold{Maximum
##'     trophic connection}, this allows for simplification of the food-web and the
##'     setting of the maximum number of connections from the focus functional group;
##'     c) \bold{the minimum proportion}, which sets the minimum proportion of a
##'     functional group in a diet that is considered as prey for the analysis; and
##'     d) \bold{time step} which sets the time step for the calculation of the food-web.
##' @title Trophic level
##' @param diet.file Character string with the path to the Diet output file. This
##'     file contains the diets of each functional group at each (recorded) time
##'     step. If the Atlantis simulation is for several years, it is highly
##'     recommended that a low frequency of recording of this output is used (set
##'     using toutinc) otherwise file sizes can become substantial (which can be
##'     very difficult to handle in R).
##' @param diet.file.bypol Detailed diet check file, this can be obtained as an extra output from Atlantis \emph{"DetailedDietCheck.txt"}. To get this file from Atlatnis turn on the option \emph{flagdietcheck} on the \emph{Run.prm} file on Atlantis.
##' @param grp.file Character string with the path to the groups \emph{*.csv} file (Atlantis input file).
##' @param quiet (Default = TRUE) this parameter helps during the process of debugging.
##' @return Reactive output with the plot of the food web for the specific time step,
##'     and a table with the trophic level of each prey and predator in the food web
##' @author Demiurgo
##' @importFrom Rdpack reprompt
##' @export
food.web <- function(diet.file, grp.file,  diet.file.bypol = NULL, quiet = TRUE){
    txtHelp <- "<h2>Summary</h2>"
    txtHelp <- paste(txtHelp, "<p>This bit of code help to visualize the food web change during the simulation output from <b>Atlantis</b> run. </p>")
    txtHelp <- paste(txtHelp, "<p>It calculate the trophic position of each functional group at each time step</p>")
    ## Libraries
    if(!quiet) cat('\n\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n # -     Step 1    -   #')
    if(!quiet) cat('\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n\n Loading libraries')
    if (!require('shiny', quietly = TRUE)) {
        stop('The package shiny was not installed')
    }
    if (!require('ncdf4', quietly = TRUE)) {
        stop('The package ncdf4 was not installed')
    }
    if (!require('plotrix', quietly = TRUE)) {
        stop('The package plotrix was not installed')
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
    if(!quiet) cat('  ...Done!')
    color.p <- brewer.pal(8, 'RdYlBu')[c(7, 2)]
    ## Reading files
    if(!quiet) cat('\n Reading files')
    dat        <- data.frame(fread(diet.file, header = TRUE, sep = ' ', showProgress = FALSE))
    dat.bp     <- pol <- NA
    if(!is.null(diet.file.bypol)){
        dat.bp <- data.frame(fread(diet.file.bypol, header = TRUE, sep = ' ', showProgress = FALSE))
        pol    <- unique(dat.bp$Box)
    }
    if(any(names(dat) == 'Group')) names(dat)[which(names(dat) == 'Group')] <- 'Predator'
    time.stp <- round(range(dat$Time), 0)
    lag      <- diff(unique(dat$Time))[1]
    grp.dat  <- read.csv(grp.file)
    stk      <- unique(dat$Stock)
    if(any(names(grp.dat) == 'InvertType')){
        names(grp.dat)[which(names(grp.dat) == 'InvertType')] <- 'GroupType'
        warning('You should change the name of the column \'InvertType\' for \'GroupType\' on your group csv file')
    }
    if(any(names(grp.dat) == 'isPredator')){
        names(grp.dat)[which(names(grp.dat) == 'isPredator')] <- 'IsPredator'
        warning('You should change the name of the column \'isPredator\' for \'IsPredator\' on your group csv file')
    }
    code.fg  <- grp.dat$Code[grp.dat$IsPredator > 0]
    if(!quiet) cat('      ...Done!')
    if(!quiet) cat('\n\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n # -     Step 2    -   #')
    if(!quiet) cat('\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n Processing and plotting')
    shinyApp(
        ui <- navbarPage('Atlantis Food Web Tool',
                         tabPanel('Food Web',
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 tags$h3('Functional Group'),
                                                 selectInput('foc.fg', 'Functional Group', c('All', as.character(code.fg))),
                                                 numericInput("m.stg", "Max. Trophic Connections:", min = 1,  max = 10, value = 4, step = 1),
                                                 numericInput("min", "Min. Proportion :", min = 0.001,  max = 1, value = 0.01, step = 0.001),
                                                 numericInput("time", "Time Step :", min = time.stp[1],  max = time.stp[2], value = time.stp[1], step = lag),
                                                 selectInput('stock', 'Stock', stk),
                                                 downloadButton("Dwnl", "Download-data"))
                                             ),
                                      column(10,
                                             plotOutput('plot1', width = "100%", height = "800px"),
                                             downloadButton("Dwnl.tl", "Download-table"),
                                             tableOutput('table')
                                             )
                                  )
                                  ),
                         tabPanel('Food Web by polygon',
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 tags$h3('Functional Group'),
                                                 selectInput('foc.fg.bp', 'Functional Group', c('All', as.character(code.fg))),
                                                 numericInput("m.stg.bp", "Max. Trophic Connections:", min = 1,  max = 10, value = 4, step = 1),
                                                 numericInput("min.bp", "Min. Proportion :", min = 0.001,  max = 1, value = 0.01, step = 0.001),
                                                 numericInput("time.bp", "Time Step :", min = time.stp[1],  max = time.stp[2], value = time.stp[1], step = lag),
                                                 selectInput('poly.bp', 'Polygon', pol),
                                                 downloadButton("Dwnl.bp", "Download-data"))
                                             ),
                                      column(10,
                                             plotOutput('plot.bp', width = "100%", height = "800px"),
                                             downloadButton("Dwnl.tl.bp", "Download-table"),
                                             tableOutput('table.bp')
                                             )
                                  )
                                  ),
                         ## -- Exit --
                         tabPanel(
                             actionButton("exitButton", "Exit")
                         )
                         ),
        function(input, output, session) {
            time.prey <- reactive({
                time.prey <- dat[dat$Time == input$time & dat$Stock  == input$stock, c(2, 6 : ncol(dat))]
            })
            time.prey.bp <- reactive({
                if(any(!is.na(dat.bp)))
                    time.prey.bp <- dat.bp[dat.bp$Time == input$time.bp & dat.bp$Box  == input$poly.bp, c(2, 6 : ncol(dat.bp))]
            })
            ## Original Recruitment
            t.prey <- reactive({
                t.prey <- tot.prey(input$foc.fg, input$m.stg, code.fg, time.prey(), input$min, grp.dat)
            })
            t.prey.bp <- reactive({
                if(any(!is.na(dat.bp)))
                    t.prey.bp <- tot.prey(input$foc.fg.bp, input$m.stg.bp, code.fg, time.prey.bp(), input$min.bp, grp.dat)
            })
            ## assing the value of trophic level for the prey
            TL <- reactive({
                TL <- NULL
                if(!is.null(t.prey())){
                    TL <- trophic.lvl(t.prey(), grp.dat)
                }
            })
            TL.bp <- reactive({
                if(any(!is.na(dat.bp)))
                TL.bp <- NULL
                if(!is.null(t.prey.bp())){
                    TL.bp <- trophic.lvl(t.prey.bp(), grp.dat)
                }
            })
            output$plot1 <- renderPlot({
                validate(
                    need((length(TL()) != 0),  'Apparently there is no interaction between predators and prey.')
                )
                plot.tlvl(TL(), t.prey(), input$foc.fg, pol = NULL, color.p)
            })
            output$plot.bp <- renderPlot({
                validate(
                    need((length(TL.bp()$Tlevel) != 0),  paste('Apparently there is no interaction between predators and prey in the box ', input$poly.bp, '. \n\nOr the Detailed diet check file is not provided.'))
                )
                plot.tlvl(TL.bp(), t.prey.bp(), input$foc.fg.bp, pol = input$poly.bp, color.p)
            })
            output$table <- renderTable({
                validate(
                    need((length(TL()) != 0),  'Apparently there is no interaction between predators and prey.')
                )
                table    <- with(TL(), data.frame(Functional.group = FG, Trophic.Level = Tlevel))
            })
            output$table.bp <- renderTable({
                validate(
                    need((length(TL.bp()$Tlevel) != 0),  paste('Apparently there is no interaction between predators and prey in the box ', input$poly.bp, '. \n\nOr the Detailed diet check file is not provided.'))
                )
                table    <- with(TL.bp(), data.frame(Functional.group = FG, Trophic.Level = Tlevel))
            })
            output$Dwnl.tl <- downloadHandler(
                filename = function(){
                    paste0(input$dataset, ".csv")
                },
                content = function(file){
                    write.csv(data.frame(TL()), file, row.names = FALSE)
                }
            )
            output$Dwnl.tl.bp <- downloadHandler(
                filename = function(){
                    paste0(input$dataset, ".csv")
                },
                content = function(file){
                    write.csv(data.frame(TL.bp()), file, row.names = FALSE)
                }
            )
            output$Dwnl <- downloadHandler(
                filename = function(){
                    paste0(input$dataset, ".csv")
                },
                content = function(file){
                    t.out <- as.data.frame(t.prey()[, 1 : 3])
                    write.csv(t.out, file, row.names = FALSE)
                }
            )
            output$Dwnl.bp <- downloadHandler(
                filename = function(){
                    paste0(input$dataset, ".csv")
                },
                content = function(file){
                    t.out <- as.data.frame(t.prey.bp()[, 1 : 3])
                    write.csv(t.out, file, row.names = FALSE)
                }
            )
            observeEvent(input$exitButton, {
                stopApp()
            })
        }
    )
}

## functions
##' @title Assing trophic level of the prey the values are based on Pauly et al 1998 [Diet composition and trophic levels of marine mammals]
##'        And Tucket and Rogers 2014 [Examining predator–prey body size, trophic level and body mass across marine and terrestrial mammals]
##' @param FGs Preys
##' @param grp.dat csv groups file
##' @return The trophic level of the Functional group
##' @author Demiurgo
prey.pos <- function(FGs, grp.dat){
    ctg <- vector('numeric')
    for( i in 1 : length(FGs)){
        ctg[i] <- which(grp.dat$Code %in% FGs[i], arr.ind = TRUE)
    }
    typ  <-  grp.dat$GroupType[ctg]
    TL.v <- ifelse(typ %in% c('SHARK', 'MAMMAL'), 4,
            ifelse(typ == 'BIRD', 3.5,
            ifelse(typ == 'FISH', 3.24,
            ifelse(typ == 'CEP', 3.2,
            ifelse(typ == 'LG_ZOO', 2.7,
            ifelse(typ == 'FISH_INVERT', 2.52,
            ifelse(typ %in% c('PWN', 'REPTILE', 'MOB_EP_OTHER', 'MED_ZOO'), 2.4,
            ifelse(typ == 'SM_ZOO', 2.4,
            ifelse(typ == 'SED_EP_FF', 2.3,
            ifelse(typ == 'CORAL', 2.2,
            ifelse(typ == 'SED_EP_OTHER', 2.1,
            ifelse(typ == 'LG_PHY', 1.5,
            ifelse(typ %in% c('MICROPHYTOBENTHOS', 'DINOFLAG'),  1.2, 1)))))))))))))
    return(TL.v)
}

##' @title Throphic level estimation. Equation based on Pauly et al 1998 [Diet composition and trophic levels of marine mammals]
##' @param DC Diet composition
##' @param TLp Trophic level of the prey
##' @return The trophic level of the predator
##' @author Demiurgo
Tlevel <- function(DC, TLp){
    tl <- 1 + (sum(DC * TLp, na.rm = TRUE) / sum(DC, na.rm = TRUE))
    return(tl)
}

##' @title Piramid order
##' @param vector Vector of numbers or character that needs to be order with the bigest value in the center
##' @return A vector with the higest value in the center
##' @author Demiurgo
ord <- function(vector){
    vec     <- vector[order(vector)]
    or.v    <- order(vector)
    vec.out <- NA * vec
    rev <- length(vec)
    fwd <- 1
    for(v in 1  :  length(vec)){
        if(v%%2 == 0){
            vec.out[rev] <- or.v[v]
            rev <- rev - 1
            next
        }
        vec.out[fwd] <- or.v[v]
        fwd <- fwd + 1
    }
    return(vec.out)
}

##' @title Total predation
##' @param foc.fg Focus functional group
##' @param connect Total connection
##' @param Fg.code Original code of the functional groups from the csv file
##' @param real.pred Realized predation
##' @param min Minimum proportion of the prey in the diet
##' @param grp.dat Functional group basic information
##' @return Total consumed prey by predator
##' @author Demiurgo
tot.prey <- function(foc.fg, connect, Fg.code, real.pred, min, grp.dat){
    if(foc.fg == 'All') foc.fg <- as.character(Fg.code)
    t.prey <- NULL
    stg    <- 1
    while(stg <= connect){
        for(fg in foc.fg){
            ## approach remove not common prey but keeps the proportion
            diet     <- colSums(real.pred[real.pred$Predator == fg, 2 : ncol(real.pred)])
            diet     <- diet[diet > 0] / sum(diet, na.rm = TRUE)
            if(length(diet) == 0) next
            diet     <- diet[diet >= min]
            n.prey   <- names(diet)
            m.prey   <- data.frame(Pred = fg, Prey = n.prey, value = diet, Stage = stg)
            t.prey   <- rbind(t.prey, m.prey)
        }
        foc.fgp <- foc.fg
        foc.fg  <- as.character(unique(t.prey$Prey[t.prey$Stage == stg]))
        if(all(length(foc.fg) == length(foc.fgp), all(foc.fg %in% foc.fgp))) break
        stg <- stg + 1
    }
    if(!is.null(t.prey$Prey)){
        t.prey$TLprey <- prey.pos(t.prey$Prey, grp.dat)
        t.prey <- t.prey[!duplicated(t.prey[, c(1, 2)]), ]  ## removing duplicates
        t.prey <- as.data.frame(t.prey)
        return(t.prey)
    } else{
        return(NULL)
    }
}

##' @title Trophic level
##' @param rel.prey realized diet by predator
##' @param grp.dat Functional group basic information
##' @return The trophic level by predator
##' @author Javier Porobic
trophic.lvl <- function(rel.prey, grp.dat){
    npred         <- unique(rel.prey$Pred)
    TL            <- NULL
    for(pred in npred){
        pospred <- which(rel.prey$Pred %in% pred)
        TL      <- c(TL, Tlevel(rel.prey$value[pospred], rel.prey$TLprey[pospred]))
    }
    pp.prey  <- unique(rel.prey$Prey[ - which(rel.prey$Prey %in% npred)])
    pp.prey  <- data.frame(FG = pp.prey, Tlevel = prey.pos(pp.prey, grp.dat))
    TL       <- data.frame(FG = npred, Tlevel = TL)
    TL       <- rbind(TL, pp.prey)
    ## Updating the trophic level based on the bsic TL and the total realized diet
    prey.f   <- as.data.frame(rel.prey)
    for(fg in 1 : nrow(TL)){
        pntl <- which(prey.f$Prey %in% TL$FG[fg])
        if(length(pntl) == 0) next
        prey.f$TLprey[pntl] <- TL$Tlevel[fg]
    }
    npred         <- unique(prey.f$Pred)
    TL            <- NULL
    for(pred in npred){
        pospred <- which(prey.f$Pred %in% pred)
        TL      <- c(TL, Tlevel(prey.f$value[pospred], prey.f$TLprey[pospred]))
    }
    pp.prey  <- unique(prey.f$Prey[ - which(prey.f$Prey %in% npred)])
    pp.prey  <- data.frame(FG = pp.prey, Tlevel = prey.pos(pp.prey, grp.dat))
    TL       <- data.frame(FG = npred, Tlevel = TL)
    TL       <- rbind(TL, pp.prey)
    ## location on the plot
    his     <- hist(TL$Tlevel, breaks = c(0.5 : 6.5), plot = FALSE)
    v.lev   <- max(his$counts)
    brk     <- his$breaks
    h.lev   <- max(his$mids[his$counts > 0]) + 1 ## getting the top of the plot
    TL$v.lev <- v.lev
    TL$h.lev <- h.lev
    TL$vpos <- NA
    for(i in 1 : (length(brk) - 1)){
        nfg.ly          <- which(TL$Tlevel >= brk[i] & TL$Tlevel < brk[i + 1] )
        tot.fg          <- length(nfg.ly)
        vpos            <-  cumsum(rep(v.lev / tot.fg, tot.fg))  - (v.lev / tot.fg) * 0.5
        pos.or          <- vector('numeric', length(TL$FG[nfg.ly]))
        ## it's necesary to order the fg to get a nice plot
        for(i in 1 : length(TL$FG[nfg.ly])){
            pos.or[i] <- sum(c(prey.f$Pred %in% TL$FG[nfg.ly][i], prey.f$Prey %in% TL$FG[nfg.ly][i]))
        }
        TL$vpos[nfg.ly] <- vpos[ord(pos.or)]
    }
    TL <- as.data.frame(TL)
    return(TL)
}

##' @title Plot trophic levels
##' @param T.lvl Trophic level by predator
##' @param rel.prey Realized prey
##' @param foc.fg Focus funciontal group
##' @param pol Polygon
##' @param color.p Colors used
##' @return A food web plot
##' @author Javier Porobic
plot.tlvl <- function(T.lvl, rel.prey, foc.fg, pol = NULL, color.p){
    rad     <- 0.25
    y.lab <- "Trophic-level"
    if(!is.null(pol)) y.lab <- paste0("Trophic-level at polygon ",  as.numeric(pol))
    plot(1, type = "n", xlab = '', ylab = y.lab,
         xlim = c(0, unique(T.lvl$v.lev)), ylim = c(0.5, unique(T.lvl$h.lev)), axes = FALSE)
    axis(2, at = c(0.5 : unique(T.lvl$h.lev)), las = 1)
    for( i in 1 : nrow(rel.prey)){
        p.pred <- which(T.lvl$FG %in% rel.prey$Pred[i])
        p.prey <- which(T.lvl$FG %in% rel.prey$Prey[i])
        ## Direction of predation
        col    <- ifelse(T.lvl$Tlevel[p.pred] < T.lvl$Tlevel[p.prey], color.p[1], color.p[2])
        y      <- c(T.lvl$Tlevel[p.pred], T.lvl$Tlevel[p.prey] )
        x      <- c(T.lvl$vpos[p.pred], T.lvl$vpos[p.prey] )
        lines(x, y , lty = 1, lwd = (rel.prey$value[i] * 3), col = col)
    }
    ## circles
    for( i in 1 : nrow(T.lvl)){
        draw.circle(T.lvl$vpos[i], T.lvl$Tlevel[i], radius = rad, nv = 1500, border = NULL,
                    col = ifelse(i == 1 && foc.fg != 'All', 'steelblue', 'gray91'), lty = 1, lwd = 1)
        text(T.lvl$vpos[i], T.lvl$Tlevel[i], T.lvl$FG[i], cex = .8, font = 2, col = ifelse(i == 1 && foc.fg != 'All', 'white', 1))
    }
    legend('topright', c('DownTop', 'TopDown'), lty = 1, col = c(color.p[1], color.p[2]), lwd = 1.5, bty = 'n')
}
