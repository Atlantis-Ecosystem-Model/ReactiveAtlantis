##' This function helps to estimate the recruitment for age structured functional
##' groups and the primary and secondary production.
##' \itemize{
##'     \item \bold{Recruits and YOY}: This option shows the parameter values and the
##' recruitment equation for the chosen functional groups (i.e. Ricker or
##' Beverton-holt). Based on the larval survival and the biomass at each reproductive
##' time step the function allows the user to re-calculate
##'     the recruitment based on a new set of parameters provided by the user.
##'      \item \bold{Primary production and light limitation}:
##'     This option provides a view of the
##'     primary production (phytoplankton, seagrass and macroalgae),  light and
##'     secondary production (zooplankton) by box and layer. This helps to  calibrate
##'     the growth of primary producer and the consumption of light.}
##' @title Estimation of Recruitment and primary producer growth
##' @param ini.nc.file Character string with the path to the netcdf file to read
##' in. This netcdf file contains the initial conditions for the Atlantis model and
##' usually ends in \emph{.nc}".
##' @param out.nc.file Character string with the path to the netcdf file to read
##' in. This netcdf file contains is a generic output from an Atlantis run usually
##' starts with \emph{output} and ends in \emph{.nc}".
##' @param yoy.file Character string with the path to Young of the Year standard
##' output file. Usually the name for this file is : \emph{output[YOUR_MODEL]YOY.txt}
##' where [YOUR_MODEL] is the name of your Atlantis model.
##' @param grp.file Character string with the path to the Groups \emph{*.csv} file
##' (Atlantis input file).
##' @param prm.file Character string with the path to the biology parameter file
##' \emph{*.prm}.
##' @param quiet (Default = TRUE) this parameter helps during the process of
##' debugging.
##' @return The outputs of this function are divided into 3 tabs:
##' \itemize{
##'         \item \bold{recruitment and YOY}: Shows the recruitment and YOY curves
##' from the Atlantis output for each functional group and provides the option to
##' test different parameter values to obtain new recruitment and YOY curves.
##'         \item \bold{Growth Zoo and PP's}: Allows the user to check the light
##' levels and biomass values for primary producers and zooplankton by box and
##' layer. The user has the option to set those values as a proportion by box,
##' proportion by layer, or their logarithmic value.
##'     \item \bold{Help}: Shows information about the inputs, parameter values and
##' output. It also, provides an overview of the different options for the user.}
##' @author Demiurgo
##' @export
recruitment.cal <- function(ini.nc.file, out.nc.file, yoy.file, grp.file, prm.file,  quiet = TRUE){
    txtHelp <- "<h2>Summary Recruit and YOY page</h2>"
    txtHelp <- paste(txtHelp, "<p>This code Helps to calibrate the recruitment for <b>Atlantis</b> based on the recruitment and allows you to test new values</p>")
    txtHelp <- paste(txtHelp, "<p><b>Rec_model</b> Recruitment model used for the functional group</p>")
    txtHelp <- paste(txtHelp, "<p><b>Alpha</b> Current value of Alpha ( on the .prm file)</p>")
    txtHelp <- paste(txtHelp, "<p><b>Beta</b> Current value of Alpha ( on the .prm file)</p>")
    txtHelp <- paste(txtHelp, "<p><b>Initial YOY</b> Initial value of Young of the Year (YOY) for this functional groups reported in the 'yoy.file' Atlantis output</p>")
    txtHelp <- paste(txtHelp, "<p><b>New Alpha</b> You can put a new value for Alpha and calculate a new YOY for that new value</p>")
    txtHelp <- paste(txtHelp, "<p><b>New Beta</b> You can put a new value for Beta and calculate a new YOY for that new value</p>")
    txtHelp <- paste(txtHelp, "<h4>Outputs</h4>")
    txtHelp <- paste(txtHelp, "<p><b>Plot 1</b> YOY and Larvaes calculated for ATLANTIS; new larvae and new YOY values calculated using <b>New Alpha</b> and <b>New Beta</b> values </p>")
    txtHelp <- paste(txtHelp, "<p><b>Plot 2</b> Proportion of the YOY compared with the initial value (YOY<sub>0</sub></))p>")
    txtHelp <- paste(txtHelp, "<p><b>Table</b> Output from both plots by each reproduction period</p>")
    txtHelp <- paste(txtHelp, "<h2>Summary Growth Zooplankton and Primary producer</h2>")
    txtHelp <- paste(txtHelp, "<p><b>Functional group 1 </b> and <b>Functional group 2</b> allow you to highlight the FGs on the plots</p>")
    txtHelp <- paste(txtHelp, "<p><b>Box</b> Select the box that you are interested</p>")
    txtHelp <- paste(txtHelp, "<p><b>Layer-Proportion</b> Each time series on the plot is scaled between 0 and 1 in each layer (divided by the highest value of the time series on that specific layer and box)  (<i>Default = TRUE</i>)</p>")
    txtHelp <- paste(txtHelp, "<p><b>Box-Proportion</b> Each time series on the plot is scaled between 0 and 1 (divided by the highest value of the time series on that box)  (<i>Default = FALSE</i>)</p>")
    txtHelp <- paste(txtHelp, "<p><b>Logarithm</b> Logarithm of the time series(<i>Default = FALSE</i>)</p>")
    txtHelp <- paste(txtHelp, "<h4>Outputs</h4>")
    txtHelp <- paste(txtHelp, "<p>One Plot for layer with the time series of <b><i>Zooplankton</i></b>, <b><i>Primary Producers</i></b>, <b><i>Light</i></b> and <b><i>Eddies</i></b></p>")
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
    if (!require('reshape', quietly = TRUE)) {
        stop('The package reshape was not installed')
    }
    if (!require('tidyverse', quietly = TRUE)) {
        stop('The package tidyverse was not installed')
    }
    if (!require('stringr', quietly = TRUE)) {
        stop('The package stringr was not installed')
    }
    if (!require('RColorBrewer', quietly = TRUE)) {
        stop('The package stringr was not installed')
    }
    ## general settings
    colors    <- brewer.pal(n = 8, name = "Set1")
    if(!quiet) cat('  ...Done!')
    ## Reading files
    if(!quiet) cat('\n Reading files')
    nc.ini    <- nc_open(ini.nc.file)
    yoy       <- read.csv(yoy.file, sep = ' ')
    group.csv <- read.csv(grp.file)
    ## remove here those that are turned off!!
    nc.out    <- nc_open(out.nc.file)
    prm       <- readLines(prm.file, warn = FALSE)
    mg2t      <-  0.00000002 ## mg C converted to wet weight in tonnes == 20 / 1000000000
    if(!quiet) cat('      ...Done!')
    if(!quiet) cat('\n\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n # -     Step 2    -   #')
    if(!quiet) cat('\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n Processing')
    sp.dat    <- with(group.csv, which(IsTurnedOn == 1 & NumCohorts > 1)) ## Age structure groups
    options(warn =  - 1)

    rec       <- text2num(prm, '^flagrecruit', FG = 'look') ## Avoinding the annoying warnings
    options(warn =  0)
    rec       <- rec[complete.cases(rec), ] ## avoiding NAs
    rec       <- cbind(rec, Alpha = NA, Beta = NA, KSPA = NA, FSP = NA, Time.sp  = NA, Period.sp = NA,
                       Time.rec = NA, Period.rec = NA, XCN = text2num(prm, 'X_CN', FG = 'look')[1, 2],
                       XRS = text2num(prm, 'X_RS', FG = 'look')[1, 2], Rec.SNW = NA, Rec.RNW = NA)
    sps    <- gsub(pattern = '^flagrecruit', '', rec$FG)
    rec$FG <- sps
    max.gr <- with(group.csv, max(NumCohorts[IsTurnedOn == 1]))
    FSPB   <- NULL
    #m.spw  <- with(group.csv, Code[which(NumSpawns > 1)]) ## special option for multiple spawn
    if(!quiet) cat('\n Reading parameters')
    t = 1 ## counter
    for(fg.r in 1 : length(sps)){
        if(!sps[fg.r] %in% group.csv$Code[sp.dat]) next()
        if(rec$Value[fg.r] == 3){  ## Beverton Holt Recruitment
            rec$Alpha[fg.r]     <- text2num(prm, paste0('BHalpha_', sps[fg.r]), FG = 'look')[1, 2]
            rec$Beta[fg.r]      <- text2num(prm, paste0('BHbeta_', sps[fg.r]), FG = 'look')[1, 2]
            rec$KSPA[fg.r]      <- text2num(prm, paste0('KSPA_', sps[fg.r]), FG = 'look')[1, 2]
            rec$FSP[fg.r]       <- text2num(prm, paste0('FSP_', sps[fg.r]), FG = 'look')[1, 2]
            ## Proportion of mature at age
            fspb.tmp            <- text2num(prm, paste0('FSPB_', sps[fg.r]), FG = sps[fg.r], Vector = TRUE)
            fspb.tmp            <- fspb.tmp[!is.na(fspb.tmp)]
            fspb.tmp            <- c(fspb.tmp, rep(0, (max.gr - length(fspb.tmp))))
            FSPB                <- rbind(FSPB, fspb.tmp)
            rownames(FSPB)[t]   <- as.character(sps[fg.r])
            t                   <- t + 1
        } else if(rec$Value[fg.r] == 1 || rec$Value[fg.r] == 12){ ## constant recruitment
            rec$Alpha[fg.r]     <- text2num(prm, paste0('KDENR_', sps[fg.r]), FG = 'look', Vector = TRUE)[1]
            fspb.tmp            <- text2num(prm, paste0('FSPB_', sps[fg.r]), FG = sps[fg.r], Vector = TRUE)
            fspb.tmp            <- fspb.tmp[!is.na(fspb.tmp)]
            fspb.tmp            <- c(fspb.tmp, rep(0, (max.gr - length(fspb.tmp))))
            FSPB                <- rbind(FSPB, fspb.tmp)
            rownames(FSPB)[t]   <- as.character(sps[fg.r])
            t                   <- t + 1
        }
        rec$Time.sp[fg.r]   <- text2num(prm, paste0(sps[fg.r], '_Time_Spawn'), FG = 'look')[1, 2]
        rec$Period.sp[fg.r] <- text2num(prm, paste0(sps[fg.r], '_spawn_period'), FG = 'look')[1, 2]
        rec$Time.rec[fg.r]  <- text2num(prm, paste0(sps[fg.r], '_Recruit_Time'), FG = 'look')[1, 2]
        rec$Period.rec[fg.r]<- text2num(prm, paste0('Recruit_Period_', sps[fg.r]), FG = 'look')[1, 2]
        if(rec$Value[fg.r] != 1 && rec$Value[fg.r] != 12){
            rec$Rec.SNW[fg.r]   <- text2num(prm, paste0('KWSR_', sps[fg.r]), FG = 'look')[1, 2]
            rec$Rec.RNW[fg.r]   <- text2num(prm, paste0('KWRR_', sps[fg.r]), FG = 'look')[1, 2]
        }
    }
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Primary producers section
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~
    pp.pos  <- with(group.csv, which(GroupType %in% c('MED_ZOO', 'LG_ZOO', 'LG_PHY', 'SM_PHY', 'PHYTOBEN', 'DINOFLAG') & IsTurnedOn == 1))
    pp.fg   <- group.csv$Name[pp.pos]
    pp.cod  <- as.character(group.csv$Code[pp.pos])
    pp.list <- list()
    for(l.pp in 1 : length(pp.fg)){
        pp.list[[l.pp]]      <- ncvar_get(nc.out, paste0(pp.fg[l.pp], '_N'))
        names(pp.list)[l.pp] <- pp.cod[l.pp]
    }
    pp.list[['Light']] <- ncvar_get(nc.out, 'Light')
    pp.list[['Eddy']]  <- ncvar_get(nc.out, 'eddy')
    numlay             <- ncvar_get(nc.out, 'numlayers')
    n.box              <- dim(pp.list[['Light']])[2]
    if(!quiet) cat('          ...Done!')
    if(!quiet) cat('\n Reading YOY from Atlantis')
    ##~~~~~~~~~~~~~~~~~~~~~~~~~##
    ##    YOY file array       ##
    ##~~~~~~~~~~~~~~~~~~~~~~~~~##
    pwn.op   <- group.csv$Name[which(group.csv$GroupType == 'PWN')]
    tmp.code <- paste0(group.csv$Code[sp.dat], '.0')
    tmp.code <- tmp.code[tmp.code %in% names(yoy)]
    cod.yoy  <- data.frame(Code = tmp.code, Initial = NA)
    nam.fg   <- group.csv$Name[sp.dat]
    yoy.tmp  <- yoy * 0
    for(c in cod.yoy$Code){
        yoy.tmp[, c]   <- (yoy[, c] / yoy[1, c])
    }
    f.yoy <- cbind(Time = yoy$Time, yoy.tmp[, which(names(yoy.tmp) %in% cod.yoy$Code)])
    if(!quiet) cat('   ...Done!')
    if(!quiet) cat('\n Calculating recruits')
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ## Number and weight of individual at age in each reproduction perior  ##
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    coh.fg  <- group.csv$NumCohorts[sp.dat]
    cod.fg  <- group.csv$Code[sp.dat]
    time    <- ncvar_get(nc.out, 't') / 86400  ## to have the time step in days
    spw     <- NULL
    nam     <- NULL
    SSB.tot <- NULL
    num.tot <- NULL
    fg.spw  <- list()
    fg.ssb  <- list()
    loo     <- 1
    ## estimation of Sp by reproduction event and By FG
    for( fg in 1 : length(nam.fg)){ ## loop through functional groups
        pos.fspb <- which(rownames(FSPB) %in% cod.fg[fg])
        if(length(pos.fspb) == 0) next()
        fg.row   <- which(rec$FG %in% cod.fg[fg])
        xrs      <- rec$XRS[fg.row]
        FSP      <- rec$FSP[fg.row]
        KSPA     <- rec$KSPA[fg.row]
        sp.tmp   <- NULL
        SSB.fg   <- NULL
        num.fg   <- NULL
        ## time of spawning
        time.stp <- seq(from = 0, by = 365, to = tail(time, 1))  + rec$Time.sp[fg.row]
        time.stp <- time.stp[time.stp < tail(time, 1)]
        time.stp <- sapply(time.stp, function(x) which.min(abs(x - time)))
        spw.coh  <- list()
        ssb.coh  <- list()
        for(coh in 1 : coh.fg[fg]){
            if(nam.fg[fg] %in% pwn.op){
                name.fg <- paste0(nam.fg[fg], '_N', coh)
                nums    <- ncvar_get(nc.out, name.fg)[, , time.stp]
                SSB.tmp <- nums
                spawn   <- nums *  FSPB[pos.fspb, coh]
                num.sp  <- spawn
            } else {
                name.fg <- paste0(nam.fg[fg], coh)
                nums    <- ncvar_get(nc.out, paste0(name.fg, '_Nums'))[, , time.stp]
                mask    <- ifelse(nums > 1.e-16, 1, 0)
                nums    <- nums * mask
                rn      <- ncvar_get(nc.out, paste0(name.fg, '_ResN'))[, , time.stp] * mask
                sn      <- ncvar_get(nc.out, paste0(name.fg, '_StructN'))[, , time.stp] * mask
                wspi    <- sn * (1 + xrs)       ## minimum weigth for spawning
                rat     <- ((rn  + sn ) - wspi) ## weight deficit
                rat[(rat < 0)]   <- 0
                spawn            <- ((wspi * FSP - KSPA)  -  rat)
                SSB.tmp <- (ncvar_get(nc.out, paste0(name.fg, '_ResN'))[, , time.stp]  +
                            ncvar_get(nc.out, paste0(name.fg, '_StructN'))[, , time.stp])  *
                    ncvar_get(nc.out, paste0(name.fg, '_Nums'))[, , time.stp]
                spawn[spawn < 0] <- 0
                spawn  <- spawn *  FSPB[pos.fspb, coh] ## individual spawn
                spawn  <- spawn * nums      ## total spawn
                num.sp <- nums  * FSPB[pos.fspb, coh]
            }
            SSB.tmp <- SSB.tmp * mask
            SSB.tmp <- SSB.tmp
            spw.coh[[coh]] <- spawn    ## spawning by time step and cohort
            ssb.coh[[coh]] <- SSB.tmp ## Biomass by time step and cohort
            if(length(dim(spawn)) > 2){
                spawn   <- apply(spawn, 3, sum, na.rm = TRUE)
                SSB.tmp <- apply(SSB.tmp, 3, sum, na.rm = TRUE)
                num.tmp <- apply(num.sp, 3, sum, na.rm = TRUE)
            }
            num.fg      <- rbind(num.fg, num.tmp)
            sp.tmp      <- rbind(sp.tmp, spawn)   ## Spawning by functional group and Age class
            SSB.fg      <- rbind(SSB.fg, SSB.tmp) ## Spawning Stock by functional group and Age class
        }
        fg.spw[[loo]] <- spw.coh
        fg.ssb[[loo]] <- ssb.coh
        loo     <- loo + 1
        nam     <- c(nam, as.character(cod.fg[fg]))
        fin.sp  <- colSums(sp.tmp)
        SSB.fg  <- colSums(SSB.fg)
        num.fg  <- colSums(num.fg)
        ## all the estimation will have the same length
        if(length(num.tot) > 0 && length(num.fg) != ncol(num.tot)){
            ln <- max(length(num.fg), ncol(num.tot))
            length(num.fg)                <- ln
            tmp.num                       <- matrix(NA, ncol = ln, nrow = nrow(num.tot))
            tmp.num[, seq(ncol(num.tot))] <- num.tot
            num.tot                       <- tmp.num
        }
        if(length(spw) > 0 && length(fin.sp) != ncol(spw)){
            ln <- max(length(fin.sp), ncol(spw))
            length(fin.sp)            <- ln
            tmp.spw                   <- matrix(NA, ncol = ln, nrow = nrow(spw))
            tmp.spw[, seq(ncol(spw))] <- spw
            spw                       <- tmp.spw
        }
        if(length(SSB.tot) > 0 && length(SSB.fg) != ncol(SSB.tot)){
            ## increase the number of columns
            ln <- max(length(SSB.fg), ncol(SSB.tot))
            length(SSB.fg) <- ln
            SSB.2.tot      <- matrix(NA, ncol = ln, nrow = nrow(SSB.tot))
            SSB.2.tot[, seq(ncol(SSB.tot))] <- SSB.tot
            SSB.tot        <- SSB.2.tot
        }
        num.tot <- rbind(num.tot, num.fg)
        SSB.tot <- rbind(SSB.tot, SSB.fg)
        spw     <- rbind(spw, fin.sp)
    }
    ## FG Names
    rownames(num.tot) <- nam
    rownames(spw)     <- nam
    rownames(SSB.tot) <- nam
    if(!quiet) cat('        ...Done!')
    if(!quiet) cat('\n\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n # -     Step 3    -   #')
    if(!quiet) cat('\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n\n  -  - Plotting -  -  \n\n')
    shinyApp(
        ui <- navbarPage('Atlantis Recruitment Tool',
                         tabPanel('Recruits and YOY',
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 tags$h3('Functional Group'),
                                                 selectInput('sp', 'Functional Group', as.character(cod.fg)),
                                                 mainPanel(strong("Rec model:"), textOutput("Rec.mod")),
                                                 mainPanel(strong("Alpha: "), verbatimTextOutput("Alpha.mod", placeholder = TRUE)),
                                                 mainPanel(strong("Beta: "), verbatimTextOutput("Beta.mod", placeholder = TRUE)),
                                                 mainPanel("Initial YOY: ", textOutput("Ini.YOY")),
                                                 br(),
                                                 numericInput("new.alpha", label = "New Alpha", value = 0),
                                                 br(),
                                                 numericInput("new.beta", label = "New Beta", value = 0),
                                                 mainPanel("Recomened Alpha: ", textOutput("ratio"))
                                             )
                                             ),
                                      column(10,
                                             plotOutput('plot1', width = "100%", height = "400px"),
                                             plotOutput('plot2', width = "100%", height = "400px"),
                                             tableOutput('table')
                                             )
                                  )
                                  ),
                         tabPanel('Growth Zoo and PPs',
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 tags$h3('Functional Group'),
                                                 selectInput('sp.pp', 'Functional Group 1', as.character(c(pp.cod, 'Eddy', 'Light'))),
                                                 selectInput('sp2.pp', 'Functional Group 2', as.character(c('Light', 'Eddy', pp.cod))),
                                                 selectInput('s.box', 'Box', 0 : (n.box - 1)),
                                                 checkboxInput('l.prop', 'Layer-Proportion', TRUE),
                                                 checkboxInput('b.prop', 'Box-Proportion', FALSE),
                                                 checkboxInput('log.v', 'Logarithm', FALSE)
                                             )
                                             ),
                                      column(10,
                                             plotOutput('plot3', width = "100%", height = "800px")
                                             )
                                  )
                                  ),
                         ## -  -  Help
                         tabPanel("Help",
                                  fluidPage(HTML(txtHelp)
                                            )
                                  ),
                         ## -  - Exit
                         tabPanel(
                             actionButton("exitButton", "Exit")
                         )
                         ),
        function(input, output, session) {
            time.stp <- reactive({
                time.stp <- seq(from = 0, by = 365, to = tail(time, 1))  + rec$Time.sp[rec$FG == input$sp]
                time.stp <- time.stp[time.stp < tail(time, 1)]
            })
            ## recruitment model
            output$Rec.mod <- renderText({
                mod   <- rec$Value[rec$FG == input$sp]
                model <- ifelse(mod == 1, 'Constant recruitment',
                         ifelse(mod == 2, 'Determined by chlA',
                         ifelse(mod == 3, 'Beverton-Holt',
                         ifelse(mod == 4, 'Random Lognormal',
                         ifelse(mod == 12, 'Fixed offspring', 'Other')))))
            })
            ## Original Recruitment
            rec.bio <- reactive({
                mod      <- rec$Value[rec$FG == input$sp]
                spawn.fg <- spw[input$sp, ]
                biom.fg  <- SSB.tot[input$sp, ]
                num.fg   <- num.tot[input$sp, ]
                sp.plt   <- paste0(input$sp, '.0')
                if(mod == 3){
                    recruit  <- unlist(BH.rec(spawn.fg, rec$Alpha[rec$FG == input$sp], rec$Beta[rec$FG == input$sp], biom.fg))
                    new.rec  <- unlist(BH.rec(spawn.fg, input$new.alpha, input$new.beta, biom.fg))
                } else if(mod == 12){
                    recruit <- rec$Alpha[rec$FG == input$sp] * num.fg
                    new.rec <- input$new.alpha  * num.fg
                }else if(mod == 1){
                    recruit <- rec$Alpha[rec$FG == input$sp]  * rep(1, length(num.fg))
                    new.rec <- input$new.alpha * rep(1, length(num.fg))
                }
                ## Avoiding NaNs
                spawn.fg[is.na(spawn.fg)] <- 0
                biom.fg[is.na(biom.fg)] <- 0
                recruit[is.na(recruit)] <- 0
                new.rec[is.na(new.rec)] <- 0
                if(mod == 1 || mod == 12){
                    rec.bio  <- recruit
                    new.bio  <- new.rec
                } else {
                    rec.bio  <- recruit * (rec$Rec.SNW[rec$FG == input$sp] + rec$Rec.RNW[rec$FG == input$sp]) * rec$XCN[rec$FG == input$sp] * mg2t
                    new.bio  <- new.rec * (rec$Rec.SNW[rec$FG == input$sp] + rec$Rec.RNW[rec$FG == input$sp]) * rec$XCN[rec$FG == input$sp] * mg2t
                }
                yoy.fg   <- data.frame(Time = yoy$Time, Rec = yoy[, which(names(yoy) == sp.plt)])
                dif      <- sapply(time.stp(), function(x){ which.min(abs(x - yoy.fg[, 1]))})
                n.yoy    <- yoy.fg[dif, 2]
                Time.yoy <- yoy.fg[dif, 1]
                n.f.yoy  <- f.yoy[dif, c(1, which(names(f.yoy) == sp.plt))]
                rec.bio  <- rec.bio[1 : length(time.stp())]
                new.bio  <- new.bio[1 : length(time.stp())]
                prop.dif <- yoy.fg$Rec[dif] / rec.bio
                fst.val  <- (new.bio * prop.dif) /  yoy.fg$Rec[1]
                N.YOY    <- (new.bio * prop.dif)
                df.end   <- data.frame(Rec = rec.bio, N.YOY = N.YOY, N.Rec = new.bio, BYOY = n.yoy, TYOY = Time.yoy, PrpYOY = n.f.yoy[, 2], Prp.st = fst.val, P.diff = prop.dif, Model = mod)
                if(mod %in% c(3, 10, 19)) df.end$Ratio <-  biom.fg[1 : length(time.stp())] / recruit[1 : length(time.stp())]
                df.end
            })
            ## Out Primary producers List
            o.pp <- reactive({
                box         <- as.numeric(input$s.box) + 1
                ly.box      <- numlay[box]
                out.pp.list <- list()
                nam.plot    <- paste0('layer ', c(ly.box : 1))
                nam.plot    <- c('Time', 'FG', paste0(nam.plot[1], ' [Deepest]'), nam.plot[c( - 1,  - ly.box)],
                                 paste0(nam.plot[ly.box], ' [Surface]'))
                for( i in 1 : length(pp.list)){
                    if(length(dim(pp.list[[i]])) == 3){
                        nr               <- nrow(pp.list[[i]])
                        out.pp.list[[i]] <- pp.list[[i]][c((nr - ly.box) : (nr - 1)), box, ]
                    } else {
                        out.pp.list[[i]]       <- array(0, dim = c(ly.box, length(pp.list[[i]][box, ])))
                        out.pp.list[[i]][1, ]  <- pp.list[[i]][box, ]
                        if(names(pp.list)[i]  == 'Eddy'){
                            out.pp.list[[i]] <- matrix(rep(pp.list[[i]][box, ], ly.box), nrow = ly.box, byrow = TRUE)
                        }
                    }
                    out.pp.list[[i]][out.pp.list[[i]] <= 1e-8] <- 0 ## removing ceros from atlatnis
                    if(input$l.prop == TRUE & input$b.prop == FALSE){
                        out.pp.list[[i]] <- out.pp.list[[i]] / apply(out.pp.list[[i]] , 1, max, na.rm = TRUE)
                    } else if (input$l.prop == FALSE & input$b.prop == TRUE){
                        out.pp.list[[i]] <- out.pp.list[[i]] / max(out.pp.list[[i]] , na.rm = TRUE)
                    } else {
                        warning('\nYou need to choose between proportion by box or proportion by layer, but you cannot use both\n')
                    }
                    if(input$log.v == TRUE){
                        out.pp.list[[i]] <- log(out.pp.list[[i]])
                    }
                    out.pp.list[[i]]           <- t(out.pp.list[[i]])
                    out.pp.list[[i]]           <- data.frame(1 : nrow(out.pp.list[[i]]), names(pp.list)[i], out.pp.list[[i]])
                    colnames(out.pp.list[[i]]) <- nam.plot
                    names(out.pp.list)[i]      <- names(pp.list)[i]
                }
                ## To get the proper plot in the right order
                ordn        <- c(names(pp.list)[names(pp.list) != input$sp.pp & names(pp.list) != input$sp2.pp], input$sp.pp, input$sp2.pp)
                out.pp.list <- out.pp.list[ordn]
                ## getting ready for ggplot
                sel.data <- do.call(rbind.data.frame, out.pp.list)
                sel.data <- melt(sel.data, id = c('Time', 'FG'))
            })
            observeEvent(input$exitButton, {
                stopApp()
            })
            output$Alpha.mod <- renderText({
                Alpha   <- rec$Alpha[rec$FG == input$sp]
            })
            output$Beta.mod <- renderText({
                Beta   <- rec$Beta[rec$FG == input$sp]
            })
            output$Ini.YOY <- renderText({
                sp.plt   <- paste0(input$sp, '.0')
                Ini.YOY  <- yoy[1, which(names(yoy) == sp.plt)]
            })
            output$ratio <- renderText({
                if(rec.bio()$Model[1] %in% c(3, 10, 19)) {
                    ratio <- rec.bio()$Ratio[1] * as.numeric(rec$Alpha[rec$FG == input$sp])
                } else {
                    ratio <- 'Not Available for this model'
                }
                ratio
            })
            output$plot1 <- renderPlot({
                par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
                plot(rec.bio()$TYOY, rec.bio()$BYOY , xlab = 'Time (days)', ylab = ifelse(rec.bio()$Model[1] %in% c(1, 12), 'Numbers', 'Biomass [Tonnes]'), las = 1, bty = 'n', pch = 20,type = 'b',
                     col = 'royalblue', ylim = c(0, max(rec.bio()$Rec, rec.bio()$BYOY, rec.bio()$N.Rec)),
                     xlim = range(c(rec.bio()$TYOY, time.stp())))
                lines(time.stp(), rec.bio()$Rec, type = 'b', pch = 20, col = 'red4')
                lines(time.stp(), rec.bio()$N.Rec, type = 'b', pch = 20, col = 'green4')
                lines(time.stp(), rec.bio()$N.YOY, type = 'b', pch = 20, col = 'yellowgreen')
                legend("topright", inset = c(-0.1, 0), legend = c('Atlantis YOY', 'Larvaes', 'New Larvaes', 'New YOY'),
                       lty=c(1, 1, 1, 1), col=c('royalblue', 'red4', 'green4', 'yellowgreen'), pch = c(20, 20, 20, 20), bty = 'n', lwd = 2)
            })
            output$plot2 <- renderPlot({
                par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
                plot(rec.bio()$TYOY, rec.bio()$PrpYOY, ylim = c(0, ifelse(max(rec.bio()$PrpYOY) < 1, 1, max(rec.bio()$PrpYOY))), bty = 'n', type = 'b',
                     lty = 2, pch = 19, col = 'olivedrab4', ylab = 'proportion from initial yoy (t)', xlab = 'Time (days)', las = 1,
                     xlim = range(c(rec.bio()$TYOY, time.stp())))
                lines(time.stp(), rec.bio()$Prp.st, type = 'b', pch = 19, lty = 2, col = 'yellow3')
                legend('topright', inset=c(-0.1, 0), legend = c('YOY prop', 'New prop'), lty = 2, col=c('olivedrab4', 'yellow3' ), pch = c(19, 19), bty = 'n', lwd = 2)
            })
            output$plot3 <- renderPlot({
                ## colors
                colo  <- c(rep('grey', (length(pp.list) - 2)), colors[c(1, 2)])
                ggplot(o.pp(), aes(x = Time, y = value, colour = FG)) + geom_line(na.rm = TRUE) +
                    facet_wrap(~ variable, ncol = 2) + ylim(ifelse(input$log.v == TRUE, NA, 0), max(o.pp()$value, na.rm = TRUE)) +
                    scale_colour_manual(values = colo)
            })
            output$table <- renderTable({
                table <- with(rec.bio(), data.frame(Time.Larv = time.stp(), TimeYOY = TYOY, Larvaes.Atlantis = Rec, YOY.Atlantis = BYOY, Diff.Prop = P.diff * 100, Est.Larvaes = N.Rec, Est.YOY = N.YOY))
                if(rec.bio()$Model[1] %in% c(3, 10, 19)) table$Ratio <- rec.bio()$Ratio
                table <- as.data.frame(table)
            })
        }
    )
}

##' @title Beverton Equation
##' @param sp spawning power
##' @param bha Alpha parameter
##' @param bhb Beta parameter
##' @param bio Biomass
##' @return The amout of recruit (Larvaes)
##' @author Demiurgo
BH.rec <- function(sp, bha, bhb, bio){
    num <- lapply(sp,  '*',  bha)
    den <- lapply(bio,  '+',  bhb)
    recruit <- mapply('/', num,  den, SIMPLIFY = FALSE)
    return(recruit)
}
## functions
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
            tmp2    <- unlist(strsplit(tmp[1], split = '_'))
            if(FG[1] == 'look') {
                col1[i] <- tmp2[1]
            } else {
                id.co   <- which(tmp2 %in% FG )
                col1[i] <- tmp2[id.co]
            }
            col2[i] <- as.numeric(tmp[2])
        }
        if(is.null(FG)) col1 <- rep('FG', length(col2))
        return(data.frame(FG = col1, Value = col2))
    } else {
        l.pat <- grep(pattern = pattern, text)
        nam   <- gsub(pattern = '[ ]+' ,  '|',  text[l.pat])
        fg    <- vector()
        pos   <- 1
        for( i in 1 : length(nam)){
            tmp     <- unlist(strsplit(nam[i], split = '|', fixed = TRUE))
            if(tmp[1] %in% c('#','##', '###')) next  ## check this part!!
            fg[pos] <- tmp[1]
            if(pos == 1) {
                pp.mat <- matrix(as.numeric(unlist(strsplit(text[l.pat[i] + 1], split = ' ', fixed = TRUE))), nrow = 1)
                pos    <- pos + 1
            } else {
                pp.tmp <- matrix(as.numeric(unlist(strsplit(text[l.pat[i] + 1], split = ' ', fixed = TRUE))), nrow = 1)
                if(ncol(pp.mat) != ncol(pp.tmp)) stop('\nError: The pPrey vector for', tmp[1], ' has ', ncol(pp.tmp))
                pp.mat <- rbind(pp.mat, pp.tmp)
                pos    <- pos + 1
            }
        }
        if(all(is.na(pp.mat[, 1]))) pp.mat <- pp.mat[, - 1]
        row.names(pp.mat)                  <- fg
        return(pp.mat)
    }
}
