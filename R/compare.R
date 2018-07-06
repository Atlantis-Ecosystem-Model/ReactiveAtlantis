##' This function allows visualization of either the output from a simulation or
##'     comparison between two different Atlantis simulations. The user can compare
##'     biomasses, abundance, structural or reserve nitrogen between two different
##'     outputs. The user is also able to check the changes is these variables by
##'     age for the different functional groups.
##' @title Compare Outputs
##' @param nc.out.current Current netcdf output file. Character string with the path
##'     to the netcdf file to read in. This netcdf file contains a generic output
##'     from an Atlantis run and usually starts with \emph{output} and ends in
##'     \emph{.nc}.
##' @param nc.out.old (Default = NULL) This is the character string with the path to
##'     the second netcdf output file to read in (used only for comparisons). This
##'     netcdf file contains a generic output from an Atlantis run and usually
##'     starts with \emph{output} and ends in \emph{.nc}.
##' @param grp.csv Character string with the path to the Groups \emph{*.csv}
##'     file (Atlantis input file).
##' @param bgm.file Character string with the path to the XY coordinates
##'     Atlantis input file \emph{*.bgm} with the information in metres.
##' @param cum.depths Cumulative depth of the Atlantis model
##' @return A reactive HTML that is divided in 3 different panels:
##' \itemize{
##' \item \bold{Biomass}: This panel shows the results for:
##'     \itemize{
##'       \item \emph{Total Biomass} which displays the changes in total biomass for
##'     each functional group.
##'       \item \emph{Relative Biomass} which displays the variation of the
##'     total biomass relative to the initial biomass (\eqn{B_{0}}).
##'      }
##'  If only one netcdf file is read in it provides a single output, which expands
##'     two sub-panels in the case of comparisons.
##' \item \bold{Total}: This function shows the changes in the variables
##'     (i.e. Biomass, Numbers and Structural and reserve nitrogen) for the entire
##'     population over time. If a second output is provided the function allows
##'     comparison of the outputs. In this function the user is able to:
##'   \itemize{
##'     \item \emph{Visualize any combination of variables}, that allows the user to activate and
##'     deactivate variables from the visualization. This provides increased flexibility
##'     and comparing between variables.
##'     \item \emph{Scaled}, this option allows for visualization of the values
##'     normalized by the value in the initial conditions (this option is very
##'     helpful for checking the main changes in variables during the calibration
##'     process).
##'     \item \emph{Limit-axis}, indicating whether the user desires total or normalized
##'     values, this option allows for the y-axis to be unscaled or to set the upper
##'     limit for the scaled version(from 0 to 3).
##'    }
##' \item \bold{By AgeClass}: This tab visualizes the changes through time of the
##'     variables by age. In this function the user is able to:
##'     \itemize{
##'     \item \emph{Visualize any combination of variables}, by activating and
##'     deactivating variables from the visualization (this provides increased
##'     flexibility to explore the changes in the variables and compare between
##'     them).
##'     \item \emph{Scaled}, this option allows for visualization of the values
##'     normalised by the value in the initial conditions (this option is very
##'     helpful for checking the main changes in variables during the calibration
##'     process).
##'     \item \emph{Limit-axis} which allows the user to set the visualization to
##'     either total or scaled values (this option allows the y-axis to be unscaled
##'     or to set the limit from 0 to 3 for the normalized version.
##' }
##' }
##' @author Demiurgo
##' @export
compare <- function(nc.out.current, nc.out.old = NULL, grp.csv, bgm.file, cum.depths){
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
    ## General configuration
    mycol  <- c(brewer.pal(8, "Dark2"), c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
    mycol  <- colorRampPalette(mycol)
    col.bi <- mycol(15)[c(14, 13, 12, 10)]
    mg2t   <- 0.00000002  ## mgC converted to wet weight in tonnes = 20/1000000000
    x.cn   <- 5.7   	  ## Redfield ratio of C:N 5.7
    ## Reading data
    grp     <- read.csv(grp.csv)
    inf.box <- boxes.prop(bgm.file,  cum.depths)
    nc.cur  <- nc_open(nc.out.current)
    ## Time Vector
    orign   <- unlist(strsplit(ncatt_get(nc.cur, 't')$units, ' ', fixed = TRUE))
    if(orign[1] == 'seconds') {
        Time <- ncvar_get(nc.cur, 't') / 86400
    } else {
        Time <- ncvar_get(nc.cur, 't')
    }
    Time    <- as.Date(Time, origin = orign[3])
    if(!is.null(nc.out.old)) nc.old <- nc_open(nc.out.old)
    grp     <- grp[grp$IsTurnedOn == 1, c('Code', 'Name', 'LongName', 'GroupType', 'NumCohorts')]
    ## Getting the total biomass
    age.grp <- grp[grp$NumCohorts > 1 & !grp$GroupType %in% c('PWN', 'CEP', 'MOB_EP_OTHER', 'SEAGRASS', 'CORAL', 'MANGROVE'), ]
    pol.grp <- grp[grp$NumCohorts == 1, ]
    ## Some model use Agestructured biomass pools
    pwn.grp <- grp[grp$NumCohorts > 1 & grp$GroupType %in% c('PWN', 'CEP', 'MOB_EP_OTHER', 'SEAGRASS', 'CORAL', 'MANGROVE'), ]
    ## Reading biomass outputs
    ## this approach allows to use only the current output file
    if(nrow(age.grp) > 0){
        age.cur.bio  <- bio.age(age.grp, nc.cur, 'Current', mg2t, x.cn)
    } else {
        age.cur.bio <- NULL
        warning('You do not have any Age structure Active in this simulation.\n Check your group.csv file\n')
    }
    pool.cur.bio <- bio.pool(pol.grp, nc.cur, 'Current', mg2t, x.cn, inf.box)
    old.bio      <- NULL
    if(!is.null(nc.out.old)){
        if(nrow(age.grp) > 0){
            age.old.bio  <- bio.age(age.grp, nc.old, 'Previous', mg2t, x.cn)
        } else {
            age.old.bio  <- NULL
        }
        pool.old.bio <- bio.pool(pol.grp, nc.old, 'Previous', mg2t, x.cn, inf.box)
        old.bio      <- rbind(pool.old.bio, age.old.bio)
    }
    pwn.bio      <- NULL
    pwn.old.bio  <- NULL
    if(nrow(pwn.grp) > 0){
        pwn.cur.bio <- bio.pwn(pwn.grp, nc.cur, 'Current', mg2t, x.cn, inf.box)
        if(!is.null(nc.out.old)){
            pwn.old.bio <- bio.pwn(pwn.grp, nc.old, 'Current', mg2t, x.cn, inf.box)
        }
        pwn.bio <- rbind(pwn.cur.bio, pwn.old.bio)
    }
    if(is.null(pwn.bio)){
        only.cur <- rbind(age.cur.bio, pool.cur.bio)
    } else {
        only.cur <- rbind(age.cur.bio, pool.cur.bio, pwn.cur.bio)
    }
    ## binding all the outputs for the biomass comparision
    t.biomass <- rbind(age.cur.bio, pool.cur.bio, old.bio, pwn.bio)
    ## relative
    rel.bio <- by(t.biomass, t.biomass$FG, relative)
    rel.bio <- do.call(rbind.data.frame, rel.bio)
    ## Start the Shiny application
    shinyApp(
        ## Create the different tabs
        ui <- navbarPage("Compare outputs",
                         tabPanel('Biomass',
                                  tabsetPanel(
                                      tabPanel('Total Biomass',
                                               plotOutput('plot1', width = "100%", height = "1000px"),
                                               downloadButton("dwn.bio", "Download")
                                               ),
                                      tabPanel('Relative Biomass',
                                               plotOutput('plot1B', width = "100%", height = "1000px"),
                                               downloadButton("dwn.rel", "Download")
                                               ))),
                         tabPanel('Total',
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 selectInput('FG2', 'Functional Group :', as.character(grp$Code)),
                                                 checkboxInput('rn2', label = strong("Reserve Nitrogen"), value = FALSE),
                                                 checkboxInput('sn2', label = strong("Structural Nitrogen"), value = FALSE),
                                                 checkboxInput('num2', label = strong("Numbers"), value = FALSE),
                                                 checkboxInput('bio2', label = strong("Biomass"), value = TRUE),
                                                 checkboxInput('scl2', label = strong("Scaled"), value = TRUE),
                                                 checkboxInput('limit2', label = strong("limit-axis"), value = TRUE),
                                                 selectInput('right2', label = strong("Legend position"), c('Right', 'Left'))
                                             ),
                                             wellPanel(
                                                 checkboxInput('cur2', label = strong("Compare current output"), value = FALSE),
                                                 selectInput('FG2b', 'Functional Group :', as.character(grp$Code)),
                                                 checkboxInput('rn2b', label = strong("Reserve Nitrogen"), value = FALSE),
                                                 checkboxInput('sn2b', label = strong("Structural Nitrogen"), value = FALSE),
                                                 checkboxInput('num2b', label = strong("Numbers"), value = FALSE),
                                                 checkboxInput('bio2b', label = strong("Biomass"), value = TRUE),
                                                 checkboxInput('scl2b', label = strong("Scaled"), value = TRUE),
                                                 checkboxInput('limit2b', label = strong("limit-axis"), value = TRUE),
                                                 selectInput('right2b', label = strong("Legend position"), c('Right', 'Left'))
                                             )                                             ),
                                      column(10,
                                             plotOutput('plot2a', width = "100%", height = "450px"),
                                             plotOutput('plot2b', width = "100%", height = "450px")
                                             ))),
                         tabPanel('By AgeClass',
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 selectInput('FG3a', 'Functional Group :', as.character(grp$Code[grp$NumCohorts > 1 & !grp$GroupType %in% c('PWN', 'CEP', 'MOB_EP_OTHER', 'SEAGRASS', 'CORAL', 'MANGROVE')])),
                                                 checkboxInput('rn3a', label = strong("Reserve Nitrogen"), value = FALSE),
                                                 checkboxInput('sn3a', label = strong("Structural Nitrogen"), value = FALSE),
                                                 checkboxInput('num3a', label = strong("Numbers"), value = FALSE),
                                                 checkboxInput('bio3a', label = strong("Biomass"), value = TRUE),
                                                 checkboxInput('scl3a', label = strong("Scaled"), value = TRUE),
                                                 checkboxInput('limit3a', label = strong("limit-axis"), value = TRUE),
                                                 selectInput('right3a', label = strong("Legend position"), c('Right', 'Left')),
                                                 downloadButton("Dwn.age", "Download")
                                             )),
                                      column(10,
                                             plotOutput('plot3a', width = "100%", height = "1000px")
                                             ))),
                         ## -- Exit --
                         tabPanel(
                             actionButton("exitButton", "Exit")
                         )
                         ),
        ## Link the input for the different tabs with your original data
        ## Create the plots
        function(input, output, session){
            total <- reactive({
                total <- nitro.weight(nc.cur, grp, input$FG2, By = 'Total', inf.box, mg2t, x.cn)
                if(input$scl2){
                    rmv         <- which(sapply(total, function(x) length(x) == 0 || is.null(x) || is.character(x)))
                    total.tmp   <- total[-rmv]
                    total.tmp   <- lapply(total.tmp, function(x) x / x[1])
                    total[-rmv] <- total.tmp
                }
                total
            })
            total2 <- reactive({
                if(!is.null(nc.out.old) | input$cur2){
                    if(input$cur2){
                        total2 <- nitro.weight(nc.cur, grp, input$FG2b, By = 'Total', inf.box, mg2t, x.cn)
                    } else {
                        total2 <- nitro.weight(nc.old, grp, input$FG2b, By = 'Total', inf.box, mg2t, x.cn)
                    }
                    if(input$scl2b){
                        rmv         <- which(sapply(total2, function(x) length(x) == 0 || is.null(x) || is.character(x)))
                        total2.tmp   <- total2[-rmv]
                        total2.tmp   <- lapply(total2.tmp, function(x) x / x[1])
                        total2[-rmv] <- total2.tmp
                    }
                } else {
                    total2 <- NULL
                }
                total2
            })
            coho <- reactive({
                coho <- nitro.weight(nc.cur, grp, FG = input$FG3a, By = 'Cohort', inf.box, mg2t, x.cn)
                if(input$scl3a){
                    rmv         <- which(sapply(coho, function(x) length(x) == 0 || is.null(x) || is.character(x)))
                    coho.tmp   <- coho[-rmv]
                    coho.tmp   <- lapply(coho.tmp, function(x) lapply(x, function(x) x / x[1]))
                    coho[-rmv] <- coho.tmp
                }
                coho
            })
            observeEvent(input$exitButton, {
                stopApp()
            })
            output$plot1 <- renderPlot({
                plot <- ggplot(t.biomass, aes(x = Time, y = Biomass, colour = Simulation)) +
                    geom_line() + facet_wrap(~ FG, ncol = 4,  scale = 'free_y') + theme_minimal() +
                    scale_color_manual(values = c('firebrick3', 'darkolivegreen'))
                plot <- update_labels(plot, list(x = 'Time step', y = 'Biomass (tons)'))
                plot
            })
            output$plot1B <- renderPlot({
                plot <- ggplot(rel.bio, aes(x = Time, y = Relative, colour = Simulation)) +
                    geom_line() + facet_wrap( ~ FG, ncol = 4) + ylim(0, 2) + theme_minimal()  +
                    annotate('rect', xmin =  - Inf, xmax = Inf, ymax = 1.5, ymin = 0.5, alpha = .1, colour = 'royalblue', fill = 'royalblue') +
                    scale_color_manual(values = c('firebrick3', 'darkolivegreen'))
                plot <- update_labels(plot, list(x = 'Time step', y = 'Relative Biomass (Bt/B0)'))
                plot
            })
            output$plot2a <- renderPlot({
                plot.age.total(total(), Time = Time, input$rn2, input$sn2, input$num2, input$bio2, input$scl2, input$limit2, input$right2, colors = col.bi)
            })
            output$plot2b <- renderPlot({
                validate(
                    need(total2() != '',  'To Display this plot you need to provide an old .nc output file or activate the box compare current (Default = FALSE)')
                )
                plot.age.total(total2(), Time = Time, input$rn2b, input$sn2b, input$num2b, input$bio2b, input$scl2b, input$limit2b, input$right2b, colors = col.bi)
            })
            output$plot3a <- renderPlot({
                n.coh <- grp$NumCohorts[grp$Code == input$FG3a]
                par(mfrow = n2mfrow(n.coh), cex = 1.2, oma = c(1, 1, 1, 1))
                for( i in 1 : n.coh){
                    plot.age.total(coho(), Time, input$rn3a, input$sn3a, input$num3a, input$bio3a, input$scl3a, input$limit3a, input$right3a, colors = col.bi, coh = i, max.coh = n.coh)
                }
            })
            output$Dwn.age <- downloadHandler(
                filename = function(){
                    paste0(input$dataset, ".csv")
                },
                content = function(file) {
                    out <- to.save(coho(), input$sn3a, input$rn3a, input$num3a, input$bio3a, Time)
                    write.csv(out, file, row.names = FALSE)
                }
            )
            output$dwn.bio <- downloadHandler(
                filename = function(){
                    paste0(input$dataset, ".csv")
                },
                content = function(file) {
                    biom.save <- list(Biomass = lapply(split(t.biomass, paste(t.biomass$FG, t.biomass$Simulation, sep = '_')),
                                                       function(x){x$Biomass}))
                    out       <- to.save(list = biom.save, b = TRUE, Time = Time, tot = TRUE)
                    write.csv(out, file, row.names = FALSE)
                }
             )
            output$dwn.rel <- downloadHandler(
                filename = function(){
                    paste0(input$dataset, ".csv")
                },
                content = function(file) {
                    biom.save <- list(Biomass = lapply(split(rel.bio, paste(rel.bio$FG, rel.bio$Simulation, sep = '_')),
                                                       function(x){x$Relative}))
                    out       <- to.save(list = biom.save, b = TRUE, Time = Time, tot = TRUE)
                    write.csv(out, file, row.names = FALSE)
                }
             )
        }
    )
}

##' @title Biomass for age groups
##' @param age.grp Age groups
##' @param nc.out ncdf atlantis' output file
##' @param ctg category
##' @return a dataframe with the biomass for all the functional groups with age classes
##' @author Demiurgo
bio.age <- function(age.grp, nc.out, ctg, mg2t, x.cn){
    grp.bio <- NULL
    for(age in 1 : nrow(age.grp)){
        cohort <- NULL
        for(coh in 1 : age.grp[age, 'NumCohorts']){
            name.fg <- paste0(age.grp$Name[age], coh)
            b.coh   <- (ncvar_get(nc.out, paste0(name.fg, '_ResN'))  +
                        ncvar_get(nc.out, paste0(name.fg, '_StructN')))  *
                ncvar_get(nc.out, paste0(name.fg, '_Nums')) * mg2t * x.cn
            b.coh   <- apply(b.coh, 3, sum, na.rm = TRUE)
            cohort  <- cbind(cohort, b.coh)
        }
        grp.bio <- rbind(grp.bio, data.frame(Time = seq(1, nrow(cohort)), FG  = as.character(age.grp$Code[age]), Biomass  = rowSums(cohort, na.rm  = TRUE), Simulation = ctg))
    }
    return(grp.bio)
}
##' @title Biomass for age groups
##' @param pol.grp Biomass pool groups
##' @param nc.out ncdf atlantis' output file
##' @param ctg category
##' @param mg2t Miligrams to tons scalar
##' @param x.cn Ration from nitrogen to carbon
##' @param area Box area
##' @param age.grp Age groups
##' @return a dataframe with the biomass for all the functional groups with age classes
##' @author Demiurgo
bio.pool <- function(pol.grp, nc.out, ctg, mg2t, x.cn, box.info){
    grp.bio <- NULL
    for(pool in 1 : nrow(pol.grp)){
        name.fg <- paste0(pol.grp$Name[pool], '_N')
        biom    <- ncvar_get(nc.out, name.fg)
        if(length(dim(biom))== 3){
           if(pol.grp$GroupType[pool] == 'LG_INF'){
               biom <- apply(biom, 3, '*', box.info$VolInf)
           }else{
               biom <- apply(biom, 3, '*', box.info$Vol)
           }
            biom <- apply(biom, 2, sum, na.rm = TRUE)
        } else {
            biom <- apply(biom, 2, function(x) x * box.info$info$Area)
            biom <- apply(biom, 2, sum, na.rm = TRUE)
        }
        biom    <- biom * mg2t * x.cn
        biom[1] <- biom[2]
        grp.bio <- rbind(grp.bio, data.frame(Time = seq(1, length(biom)), FG  = as.character(pol.grp$Code[pool]), Biomass  = biom, Simulation = ctg))
    }
    return(grp.bio)
}
##' @title Biomass for age groups
##' @param pwn.grp Biomass pools with age classes
##' @param nc.out ncdf atlantis' output file
##' @param ctg category
##' @param mg2t mg C converted to wet weight in tonnes == 20 / 1000000000
##' @param x.cn Redfield ratio of C:N 5.7
##' @return a dataframe with the biomass for all the functional groups with Biomass pool age classes
##' @author Demiurgo
bio.pwn <- function(pwn.grp, nc.out, ctg, mg2t, x.cn, box.info){
    grp.bio <- NULL
    for(pwn in 1 : nrow(pwn.grp)){
        cohort <- NULL
        for(coh in 1 : pwn.grp[pwn, 'NumCohorts']){
            name.fg <- paste0(pwn.grp$Name[pwn], '_N', coh)
            b.coh   <- ncvar_get(nc.out, name.fg) * mg2t * x.cn
            if(length(dim(b.coh)) > 2){
                ## for Volumen
                b.coh   <- apply(b.coh, 3, '*', box.info$Vol)
            } else {
                ## for area
                b.coh   <- apply(b.coh, 2, '*', box.info$info$Area)
            }
            b.coh   <- apply(b.coh, 2, sum, na.rm = TRUE)
            cohort  <- cbind(cohort, b.coh)
        }
        grp.bio <- rbind(grp.bio, data.frame(Time = seq(1, nrow(cohort)), FG  = as.character(pwn.grp$Code[pwn]), Biomass  = rowSums(cohort, na.rm  = TRUE), Simulation = ctg))
    }
    return(grp.bio)
}
##' @title Calculate the value of a time serie based on the first value
##' @param df Data frame of biomass by FG
##' @return A vector of relative values
##' @author Demiurgo
relative <- function(df, biomass = TRUE, Vec = NULL){
    if(biomass)       vector <- df$Biomass
    if(!is.null(Vec)) vector <- df[Vec]
    df$Relative  <- vector / vector[1]
    return(df)
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
##' @title Box information
##' @param bgm.file BGM file,  Atlantis input
##' @param depths Cummulative depths (de max depth of each layer)
##' @return a dataframe with information (box-id, area (m), )
##' @author Demiurgo
boxes.prop <- function(bgm.file, depths){
    bgm       <- readLines(bgm.file, warn = FALSE)
    boxes     <- text2num(bgm, 'nbox', FG = 'look')
    out       <- NULL
    depths    <- depths[ - which(depths == 0)]
    max.nlyrs <- length(depths)              ## maximum number of water layers
    vol       <- array(NA, dim = c(boxes$Value, length(depths)))
    for(b in 1 : boxes$Value){
        area      <- text2num(bgm, paste0('box', b - 1,'.area'), FG = 'look')
        z         <- text2num(bgm, paste0('box', b - 1,'.botz'), FG = 'look')
        box.lyrs  <- sum(depths <  - z$Value)
        box.lyrs  <- pmin(box.lyrs, max.nlyrs) # bound by maximum depth
        out       <- rbind(out, data.frame(Boxid = b - 1, Area  = area$Value, Volumen = area$Value *  -z$Value,
                                           Depth =  -z$Value, Layers = box.lyrs))
        vol[b, 1 : box.lyrs] <- area$Value * depths[1 : box.lyrs]
    }
    vol                <- cbind(out$Area, vol) # to include the sediment layer for later calculations
    vol                <- t(vol[, ncol(vol) : 1])
    vol2               <- vol ## arrange not for infauna
    vol2[nrow(vol2), ] <- 0
    if(any(out$Area < 1)) warning('\nOne (or more) of the boxes areas is less than 1m2,  Check if the right BGM file in xy coordinates')
    out[out$Depth <= 0, 2 : ncol(out)] <- 0
    out <- list(info   = out,
                Vol    = vol2,
                VolInf = vol)
    return(out)
}
##' @title Calculate the different weigt (estructural reserve) and number of indivduals
##' @param nc.out Netcdf out file. this is the traditional .nc file from atlantis
##' @param grp groups csv file
##' @param FG Functional group (selected)
##' @param By option to aggregate the time series in just one or leave the values by cohort
##' @param box.info Internal function
##' @return List witht he weight (reserve and structural), total biomass (tons) and number
##' @author Demiurgo
nitro.weight <- function(nc.out, grp, FG, By = 'Total', box.info, mg2t, x.cn){
    ## Age classes
    pos.fg <- which(grp$Code == FG)
    Bio <- Num <- SN  <- RN  <- list()
    if(grp[pos.fg, 'NumCohorts'] > 1 & grp[pos.fg, 'GroupType'] != 'PWN'){
        n.coh <- grp[pos.fg, 'NumCohorts']
        for(coh in 1 : n.coh){
            name.fg <- paste0(grp$Name[pos.fg], coh)
            resN    <- ncvar_get(nc.out, paste0(name.fg, '_ResN'))
            resN[which(resN == 0)] <- NA
            strN    <- ncvar_get(nc.out, paste0(name.fg, '_StructN'))
            strN[which(strN == 0)] <- NA
            nums    <- ncvar_get(nc.out, paste0(name.fg, '_Nums'))
            b.coh   <- (resN  + strN)  * nums * mg2t * x.cn
            if(By %in% c('Total', 'Cohort')){
                nums    <- apply(nums, 3, sum, na.rm = TRUE)
                b.coh   <- apply(b.coh, 3, sum, na.rm = TRUE)
                strN    <- apply(strN, 3, mean, na.rm = TRUE)
                resN    <- apply(resN, 3, mean, na.rm = TRUE)
            }
                RN[[coh]]  <- resN
                SN[[coh]]  <- strN
                Bio[[coh]] <- b.coh
                Num[[coh]] <- nums
        }
        if(By ==  'Total'){
            RN  <- rowSums(matrix(unlist(RN), ncol = n.coh))
            SN  <- rowSums(matrix(unlist(SN), ncol = n.coh))
            Bio <- rowSums(matrix(unlist(Bio), ncol = n.coh))
            Num <- rowSums(matrix(unlist(Num), ncol = n.coh))
        }
        type <- 'AgeClass'
    } else if (grp[pos.fg, 'NumCohorts'] == 1){ ## Biomass pool
        name.fg <- paste0(grp$Name[pos.fg], '_N')
        biom    <- ncvar_get(nc.out, name.fg)
        if(length(dim(biom)) == 3){
            if(grp$GroupType[pos.fg] == 'LG_INF'){
                biom <- apply(biom, 3, '*', box.info$VolInf)
            }else{
                biom <- apply(biom, 3, '*', box.info$Vol)
            }
            if(By == 'Total'){
                biom <- apply(biom, 2, sum, na.rm = TRUE)
            }
        } else {
            biom <- apply(biom, 2, function(x) x * box.info$info$Area)
            if(By == 'Total'){
                biom <- apply(biom, 2, sum, na.rm = TRUE)
            }
        }
        Bio <- biom
        type <- 'BioPool'
    } else if(grp[pos.fg, 'NumCohorts'] > 1 & grp[pos.fg, 'GroupType']  == 'PWN'){
        ## Some model use Agestructured biomass pools
        for(coh in 1 : pwn.grp[pwn, 'NumCohorts']){
            name.fg <- paste0(grp$Name[pos.fg],'_N', coh)
            biom    <- ncvar_get(nc.out, name.fg)
            biom    <- apply(biom, 3, '*', box.info$Vol)
            if(By == 'Total'){
                biom <- apply(biom, 2, sum, na.rm = TRUE)
            }
            Bio[[coh]] <- biom
        }
        type <- 'AgeBioPool'
    }
    return <- list(Biomass    = Bio,
                   Numbers    = Num,
                   Structural = SN,
                   Reserve    = RN,
                   Type       = type)
}
##' @title Interactive Plot for biomass (structural and reserve nitrogen)
##' @param total Dataframe with the weight output from the nitro.weight function
##' @param rn2 True if active the plot of reserve nitrogen
##' @param sn2 True if active the plot of structural nitrogen
##' @param num2 True if active the plot of numbers
##' @param bio2 True if active the plot of biomass
##' @param scl2 True if active the scaling base on the initial value of the time serie
##' @param limit True if active limiting the plot area from 0 to maximum 3
##' @param right choose left of right for the legend of the plor
##' @return A interactive plot
##' @author Demiurgo
plot.age.total <- function(total, Time, rn2, sn2, num2, bio2, scl2, limit, right, colors, coh = NULL, max.coh = NULL){
    ## Definig limits and general configuration of plots
    l.lab  <- 4 * ifelse(scl2, 1, sum(rn2, sn2, num2, bio2))
    yli    <- cbind(c(0, 3), sapply(total[-5], range, na.rm = TRUE))
    yli    <- yli[, c(TRUE, bio2, num2, sn2, rn2)]
    yli    <- cbind(range(yli), yli)
    tickx <- seq(from = 1, to = length(Time), length = 5)
    par(mar = c(5, l.lab, 4, 4) + 0.1)
    if(!is.null(coh)){
        xtlab  <- max.coh - rev(n2mfrow(max.coh))[1]
        nr.lab <- rev(n2mfrow(max.coh))[1]
        par(mar = c(1, l.lab, 1, 1) + 0.1)
    }
    if(bio2){
        ylim <- yli[, ifelse(limit, 2, ifelse(scl2, 1, which(colnames(yli) == 'Biomass')))]
        if(is.null(coh)){
            bms <- unlist(total$Biomass)
        } else {
            bms <- unlist(total$Biomass[[coh]])
        }
        plot(Time, bms, type = 'l', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n',
             ylim = ylim, bty = 'n', lwd = 2, col = colors[1])
        axis(1, at = tickx, labels = FALSE, lwd = 2,  col = 1)
        if(scl2){
            yticks <- seq(ylim[1], round(ylim[2]), by = round(ylim[2]) / 5)
                                        #axis(2, at = yticks , labels = yticks, lwd = 2,  col = 1)
            if(is.null(coh) || coh%%nr.lab == 1){
                axis(2, at = yticks , labels = yticks, lwd = 2,  col = 1)
                mtext(2, text = 'Relative Values (X_t/X_0)', line = 2)
            }
        } else {
            yticks <- seq(ylim[1], round(ylim[2]), by = round(ylim[2]) / 5)
            axis(2, at = yticks , labels = format(yticks, scientific = TRUE, digits = 2), lwd = 2,  col = colors[1])
            mtext(2, text = 'Biomass (tons)', line = 2)
        }
    }
    if(num2 & total$Type == 'AgeClass'){
        if(bio2) par(new = TRUE)
        sum.l <- sum(bio2) * 4
        ylim  <- yli[, ifelse(limit, 2, ifelse(scl2, 1, which(colnames(yli) == 'Numbers')))]
        if(is.null(coh)) {
            nms <- unlist(total$Numbers)
        } else {
            nms <- unlist(total$Numbers[[coh]])
        }
        plot(Time, nms, type = 'l', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n',
             ylim = ylim, bty = 'n', lwd = 2, col = colors[2])
        axis(1, at = tickx, labels = FALSE, lwd = 2,  col = 1)
        if(scl2 & sum.l < 4) {
            yticks <- seq(ylim[1], round(ylim[2]), by = round(ylim[2]) / 5)
            axis(2, at = yticks , labels = yticks, lwd = 2,  col = 1)
            if(is.null(coh) || coh%%nr.lab == 1)
                mtext(2, text = 'Relative Values (X_t/X_0)', line = 2)
        } else if(!scl2){
            yticks <- seq(ylim[1], round(ylim[2]), by = round(ylim[2]) / 5)
            axis(2, at = yticks , labels = format(yticks, scientific = TRUE, digits = 2), lwd = 2, line = sum.l, col.axis = colors[2], col = colors[2])
            mtext(2, text = 'Abundance (Numbers)', line = sum.l + 2, col = colors[2])
        }
    }
    if(rn2 & total$Type == 'AgeClass'){
        if(bio2 || num2) par(new = TRUE)
        sum.l <- sum(bio2, num2) * 4
        ylim <- yli[, ifelse(limit, 2, ifelse(scl2, 1, which(colnames(yli) == 'Reserve')))]
        if(is.null(coh)){
            rns <- unlist(total$Reserve)
        } else {
            rns <- unlist(total$Reserve[[coh]])
        }
        plot(Time, rns, type = 'l', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n',
             ylim = ylim, bty = 'n', lwd = 2, col = colors[4])
        axis(1, at = tickx, labels = FALSE, lwd = 2,  col = 1)
        if(scl2 & sum.l < 4) {
            yticks <- seq(ylim[1], round(ylim[2]), by = round(ylim[2]) / 5)
            axis(2, at = yticks , labels = yticks, lwd = 2,  col = 1)
            if(is.null(coh) || coh%%nr.lab == 1)
                mtext(2, text = 'Relative Values (X_t/X_0)', line = 2)
        } else if(!scl2){
            yticks <- seq(ylim[1], round(ylim[2]), by = round(ylim[2]) / 5)
            axis(2, at = yticks , labels = format(yticks, scientific = TRUE, digits = 2), lwd = 2, line = sum.l, col.axis = colors[4], col = colors[4])
            mtext(2, text = 'Reserve Nitrogen (mgNm3)', line = sum.l + 2, col = colors[4])
        }
    }
    if(sn2 & total$Type == 'AgeClass'){
        if(bio2 || num2 || rn2) par(new = TRUE)
        sum.l <- sum(bio2, num2, rn2) * 4
        ylim  <- yli[, ifelse(limit, 2, ifelse(scl2, 1, which(colnames(yli) == 'Structural')))]
        if(is.null(coh)){
            sts <- unlist(total$Structural)
        } else {
            sts <- unlist(total$Structural[[coh]])
        }
        plot(Time, sts, type = 'l', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n',
             ylim = ylim, bty = 'n', lwd = 2, col = colors[3])
        axis(1, at = tickx, labels = FALSE, lwd = 2,  col = 1)
        if(scl2 & sum.l < 4) {
            yticks <- seq(ylim[1], round(ylim[2]), by = round(ylim[2]) / 5)
            axis(2, at = yticks , labels = yticks, lwd = 2,  col = 1)
            if(is.null(coh) || coh%%nr.lab == 1)
                mtext(2, text = 'Relative Values (X_t/X_0)', line = 2)
        } else if(!scl2){
            yticks <- seq(ylim[1], round(ylim[2]), by = round(ylim[2]) / 5)
            axis(2, at = yticks , labels = format(yticks, scientific = TRUE, digits = 2), lwd = 2, line = sum.l, col.axis = colors[3], col = colors[3])
            mtext(2, text = 'Structural Nitrogen (mgNm3)', line = sum.l + 2, col = colors[3])
        }
    }
    if(is.null(coh) || coh > xtlab){
        axis(1, at = Time[tickx], labels = format.Date(Time[tickx], '%Y.%m'), lwd = 2,  col = 1)
        mtext("Time (days)", side=1, line = 3)
    }
    if(scl2){
        on <- seq(1:4)[c(bio2, num2, sn2, rn2)]
        legend(ifelse(right == 'Right', 'topright', 'topleft'), names(total)[on], lty=1, col = colors[on], lwd=2, bty='n')
    }
    if(!is.null(coh)){
        mtext(paste0("AgeClass - ", ifelse(coh>10, '', '0'), coh))}
}
##' @title Function to save data
##' @param list object from the function
##' @param sn Structural weigth Object
##' @param rn Reserve weight Object
##' @param n Abundance Object
##' @param b Biomass Object
##' @param Time Time Object
##' @param tot Total or partial
##' @return csv file
##' @author Demiurgo
to.save <- function(list, sn = FALSE, rn = FALSE, n = FALSE, b = FALSE, Time = Time, tot = FALSE){
    out <- NULL
    nam <- NULL
    if(rn){
        n.c    <- length(list$Reserve)
        reserv <- matrix(unlist(list$Reserve, use.names = FALSE), ncol = n.c, byrow = FALSE)
        out    <- cbind(out, reserv)
        nam    <- c(nam, paste0('RN.Age', 1 : n.c))
    }
    if(sn){
        n.c    <- length(list$Structural)
        struct <- matrix(unlist(list$Structural, use.names = FALSE), ncol = n.c, byrow = FALSE)
        out    <- cbind(out, Struct)
        nam    <- c(nam, paste0('SN.Age', 1 : n.c))
    }
    if(b){
        n.c    <- length(list$Biomass)
        biom   <- matrix(unlist(list$Biomass, use.names = FALSE), ncol = n.c, byrow = FALSE)
        out    <- cbind(out, biom)
        if(!tot){
            nam    <- c(nam, paste0('Biom.Age', 1 : n.c))
        } else {
            nam    <- names(list$Biomass)
        }
    }
    if(n){
        n.c    <- length(list$Numbers)
        num    <- matrix(unlist(list$Numbers, use.names = FALSE), ncol = n.c, byrow = FALSE)
        out    <- cbind(out, num)
        nam    <- c(nam, paste0('Num.Age', 1 : n.c))
    }
    colnames(out) <- c(nam)
    out <- data.frame(Date = Time, out)
    return(data.frame(out))
}
