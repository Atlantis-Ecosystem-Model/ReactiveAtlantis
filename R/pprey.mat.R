##' This function helps the user to calibrate and explore different aspect of the
##'     predator-prey interaction. In Atlantis, interaction between the predator and
##'     the prey is mainly defined by the predator-prey matrix which represents the
##'     maximum availability of each prey (age) group to a specific predator (age
##'     group). This consumption can be shiny::strongly affected by different processes such
##'     as the spatial overlap, the biomass of the prey and the gape limitation of
##'     the predator. With this tool you can explore the availability matrix and the
##'     processes that affect consumption. This tool also allows you to explore new
##'     values for the predator-prey matrix and shiny::observe online how these new values
##'     affect the predator-prey interaction.
##' @title Atlantis feeding tool
##' @param prm.file Character string with the path to the biological parameter \emph{*.prm} file (Atlantis input file).
##' @param grp.file Character string with the path to the Groups \emph{*.csv} file (Atlantis input file).
##' @param nc.file Character string with the path to the netcdf file to read
##'     in. This netcdf file contains the initial conditions for the Atlantis model
##'     usually ends in \emph{.nc}.
##' @param bgm.file Character string with the path to the XY coordinates Atlantis input file \emph{*.bgm} with the information in metres.
##' @param cum.depths Vector with the cumulative depths of the different layers
##'     (e.g. \emph{cum.depths <- c(0, 20, 100, 200, 500)})
##' @param quiet (Default = TRUE) this parameter helps during the process of debugging.
##' @return Shiny::Reactive shiny::HTML that displays predator-prey information such as: the ppPREY
##'     matrix, the initial abundance of prey, the overlap matrix based on gape size,
##'     predator preference and predator prey spatial overlap. All this information
##'     is divided into two different tabs:
##'  \itemize{
##'   \item \bold{Non-spatial}: This tab allows the user to check and modify
##'     availability values from the pprey matrix. This Function provides 5 outputs:
##'    \itemize{
##'      \item \bold{Availability matrix (\eqn{\lambda})}: Matrix of prey
##'     availability values (pprey) for each adult, young or biomass pool prey and
##'     predator.
##'      \item \bold{Overlap matrix}: Calculates if a predator is able to eat a prey
##'     given its gape size limitations. This function uses the knife-edge size
##'     selectivity, where availability of the prey is either available (1) or not
##'     available to the predator (0).
##'  \deqn{\omega = \left\{
##'     \begin{array}{ll}
##'     1 & SN_{predator} * KLP < SN_{prey} < SN_{predator} * KUP \\
##'     0 & Otherwise
##'     \end{array}
##'     \right. }{ (non-Latex version) }
##'  were \eqn{SN} is the structural weight; \eqn{KLP} Minimum gape limit of the
##'     predator (age structured or biomass pool); \eqn{KUP} Maximum gape limit of
##'     the predator (age structured or biomass pool).
##' \item \bold{Effective Predation}: this provides an approximation of the total
##'     biomass of prey consumed by a predator. It assumes a perfect spatial match
##'     between prey and predator.
##'     \deqn{\gamma  =  B * \lambda * \omega}
##' were \eqn{\gamma} is the effective predation; \eqn{B} is the biomass;
##'     \eqn{\lambda} is the predator-prey availability (pprey matrix); and
##'     \eqn{\omega} is the predator-prey gape overlap.
##' \item \bold{Predation pressure}: this function calculates the percentage of the
##'     predator diet corresponding to a specific prey.
##' \item \bold{Total biomass prey}: this is the total biomass per functional group
##'     and maturity stage (adult and juvenile).
##'}
##' \item \bold{Spatial Overlap}- this tab allows the user to check the spatial
##'     overlap between the prey and the predator in all the boxes and layers based
##'     on the values in the initial conditions file and the parameterized gape
##'     limitation.
##'}
##' @import stats utils grDevices ggplot2 graphics
##' @importFrom ggplot2 ggplot aes geom_bar coord_flip scale_color_manual geom_line facet_wrap theme_minimal update_labels geom_hline
##' @author Demiurgo
##' @export
feeding.mat <- function(prm.file, grp.file, nc.file, bgm.file, cum.depths, quiet = TRUE){
    txtHelp <- "<h2>Summary</h2>"
    txtHelp <- paste(txtHelp, "<p>This program displays data for the predator prey relationship for  <b>Atlantis</b> run. Also,  provide help for the tuning of the pprey matrix</p>")
    txtHelp <- paste(txtHelp, "<h3>Details</h3>")
    txtHelp <- paste(txtHelp, "<p>Plots have a zoom feature. Draw a box and double click to zoom into the box. Double click to reset zoom.</p>")
    txtHelp <- paste(txtHelp, "<p>Using the input panel,  the user can select the focal prey (<b>Prey panel</b>) and predator <b>Predator panel</b> .</p>")
    txtHelp <- paste(txtHelp, "<p>The <b>Original value :</b> box,  display the value that is used on the pprey matrix on the Biological parameter file.</p>")
    txtHelp <- paste(txtHelp, "<p>The <b>Current value :</b> box,  display the value that the user set in the current run. If the user don't save the run,  the value will be lost.</p>")
    txtHelp <- paste(txtHelp, "<p>On the <b>New value :</b> box,  the user can enter the new value for the pprey matrix for the selected predator and prey. The result of the application would be reflected on the current run in all the plots after the used click the box <b>Change value</b>.</p>")
    txtHelp <- paste(txtHelp, "<p>The <b>Write pPREY Matrix</b> bottom create a txt file that contain the new pPrey matrix created by the user</p>")
    txtHelp <- paste(txtHelp, "<p><b>Effective predation</b> Represent the total predation based on the Predator\'s biomass,  is the same approach that Beth used on the spreadsheet to calibrate predation. Values are on logarithmic scale</p>")
    txtHelp <- paste(txtHelp, "<p><b>Availability matrix</b> The raw pPREY matrix from the biological parameter file.</p>")
    txtHelp <- paste(txtHelp, "<p><b>Overlap matrix</b> shows if the predator is able to eat the prey based on the gape limitations.</p>")
    txtHelp <- paste(txtHelp, "<p><b>% of predation pressure</b> Which percentage of each prey corresponds to the total consumed by the predator.</p>")
    txtHelp <- paste(txtHelp, "<p><b>Total biomass prey</b> Total biomass of each functional group on logarithmic scale.</p>")
    if(!quiet) cat('\n\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n # -     Step 1    -   #')
    if(!quiet) cat('\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n\n Loading libraries')
    if(!quiet) cat('      ...Done!')
    ## Reading files
    if(!quiet) cat('\n Reading files')
    groups.csv        <- utils::read.csv(grp.file)
    names(groups.csv) <- tolower(names(groups.csv))
    if(any(grepl('invertype', names(groups.csv)))){
        names(groups.csv)[which(grepl('invertype', names(groups.csv)))] <- 'grouptype'
    }
    prm        <- readLines(prm.file, warn = FALSE)
    numlayers  <- find.z(bgm.file, cum.depths)
    min.depth  <- text2num(prm, '_mindepth', FG = 'look')
    max.depth  <- text2num(prm, '_maxdepth', FG = 'look')
    depth.dst  <- data.frame(FG = min.depth[, 1], Min = min.depth[, 2], Max = max.depth[which(max.depth[, 1] %in% min.depth[,1]), 2])
    ## availability matrix
    Ava.mat            <- text2num(prm, 'pPREY', Vector=TRUE, pprey = TRUE)
    colnames(Ava.mat)  <- c(as.character(groups.csv$code), 'DLsed', 'DRsed', 'DCsed')
    if(!quiet) cat('          ...Done!')
    ## Biomass,  age and Gape size
    if(!quiet) cat('\n\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n # -     Step 2    -   #')
    if(!quiet) cat('\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n\n Calculating Biomass and spatial distribution')
    out.Bio  <- Bio.func(nc.file, groups.csv, numlayers)
    Struct   <- out.Bio[[1]]
    Biom.N   <- out.Bio[[2]]
    if(!quiet) cat('       ...Done!')
    if(!quiet) cat('\n Calculating gape limitation and prey size')
    age_transition  <- text2num(prm, '_age_mat', FG = as.character(groups.csv$code))
    is.off   <- which(groups.csv$isturnedon == 0)
    if(length(is.off) > 0){ ## removing the groups that are turned off
        out <- which(age_transition$FG %in% groups.csv$code[is.off])
        if(sum(out) != 0) age_transition      <- age_transition[- out, ]
    }
    Ages      <- data.frame(FG = groups.csv$code, Adul = groups.csv$numcohorts)
    Ages$Juv  <- apply(Ages, 1, function(x){if(x[1] %in% age_transition$FG){ age_transition$Value[which(age_transition$FG %in% x[1])]} else {NA}})
    Gape     <- gape.func(groups.csv, Struct, Biom.N, prm)
    if(!quiet) cat('          ...Done!')
    if(!quiet) cat('\n Calculating size and spatial overlap')
    browser()
    Over.mat <- Over.mat.func(Ava.mat, Gape[[1]])
    bio.a    <- Bio.ages(Biom.N, age = Gape[[2]], Over.mat)
    bio.juv  <- bio.a[[1]]
    bio.adl  <- bio.a[[2]]
    bio.juv  <- data.frame(FG = bio.a[[1]][, 1], Biomass = as.numeric(bio.a[[1]][, 2]))
    bio.adl  <- data.frame(FG = bio.a[[2]][, 1], Biomass = as.numeric(bio.a[[2]][, 2]))
    b.juv    <- bio.juv[complete.cases(bio.juv), ]
    b.adl    <- bio.adl[complete.cases(bio.adl), ]
    if(!quiet) cat('               ...Done!')
    if(!quiet) cat('\n Calculating feeding pressure')
    ## Total feeding
    real.feed  <- Over.mat * NA
    pred       <- row.names(Over.mat)
    for( pd in 1 : nrow(Over.mat)){
        ## Getting the number of biomass needed by each functional group
        c.pred      <- unlist(strsplit(pred[pd],'pPREY'))[2]
        predator    <- gsub(pattern = "[[:digit:]]+", '\\1', c.pred)
        a.pred.prey <- as.numeric(unlist(strsplit(c.pred, predator)))
        pry.loc     <- which(bio.adl[, 1] %in% predator)
        if(length(a.pred.prey) == 0 || is.na(a.pred.prey)) a.pred.prey[2] <- 2
        ## Young Predator
        if(a.pred.prey[2] == 1){
            ## Young Prey
            real.feed[pd, ] <- (Over.mat[pd, ] * as.numeric(bio.juv[, 2]))
        } else {
            ## Adult Prey
            real.feed[pd, ] <- (Over.mat[pd, ] * as.numeric(bio.adl[, 2]))
        }
    }
    if(!quiet) cat('                       ...Done!')
    ## Real Overlap matrix Including pPREY and Overlap matrix
    t.o.mat <- t(Over.mat * Ava.mat)
    t.o.mat[which(t.o.mat > 0)] <- 1
    ## Plot output
    real.feed <- real.feed * Ava.mat
    ## Gape overlap for spatial output
    ntrans   <- reshape::melt(t.o.mat)
    over.tmp <- do.call(rbind.data.frame, apply(ntrans, 1, sepText))
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Spatial Overlap functions and procedures
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(!quiet) cat('\n Reading and preparing the spatial data for plotting')
    juv.sp.ov <- as.character(unlist(apply(Ages, 1, function(x) if(!is.na(x[3])) paste(rep(x[1], x[3]), 1 : x[3], sep = '_'))))
    ad.sp.ov  <- as.character(unlist(apply(Ages, 1, function(x) if(!is.na(x[3])){paste(rep(x[1], x[2]), x[3] : x[2], sep = '_')} else {paste(x[1], 1, sep = '_')} )))
    ad.sp.ov  <- setdiff(ad.sp.ov,  juv.sp.ov)
    juv.sp.ov <- out.Bio[[3]][juv.sp.ov]
    juv.sp.ov <- sapply(as.character(Ages$FG), function(x) rowSums(juv.sp.ov[, grep(x, names(juv.sp.ov)), drop = FALSE]))
    ## Adults
    ## removing groups that are not in the matrix
    ad.sp.ov  <- ad.sp.ov[which(ad.sp.ov %in% colnames(out.Bio[[3]]))]
    ad.sp.ov  <- out.Bio[[3]][ad.sp.ov]
    ad.sp.ov  <- sapply(as.character(groups.csv$code), function(x) rowSums(ad.sp.ov[, grep(x, names(ad.sp.ov)), drop = FALSE]))
    ## Binary values
    ad.sp.ov[ad.sp.ov > 1]   <- 1
    juv.sp.ov[juv.sp.ov > 1] <- 1
    land      <- rep(numlayers[, 3], each = length(cum.depths))
    ad.sp.ov  <- data.frame (Stage = 'Adult', Land = land, ad.sp.ov)
    juv.sp.ov <- data.frame (Stage = 'Juvenile', Land = land, juv.sp.ov)
    ad.sp.ov  <- cbind(out.Bio[[3]][, 1 : 2], ad.sp.ov)
    juv.sp.ov <- cbind(out.Bio[[3]][, 1 : 2], juv.sp.ov)
    sp.ov     <- rbind(reshape::melt(ad.sp.ov, id = c('Layer', 'Box', 'Stage', 'Land')),
                       reshape::melt(juv.sp.ov, id = c('Layer', 'Box', 'Stage', 'Land')))
    sediment  <- max(sp.ov$Layer, na.rm = TRUE)
    ## Geting the map from the bgm file
    map <- make.map(bgm.file)
    ## resolution table
    plotHeigh <- paste(as.character(40 * (ncol(t.o.mat) %/% 3 + 1)), "px", sep = "")
    plotWidth <- paste(as.character(45 * (nrow(t.o.mat) %/% 3 + 1)), "px", sep = "")
    if(!quiet) cat(' ...Done!')
    if(!quiet) cat('\n\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n # -     Step 3    -   #')
    if(!quiet) cat('\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n\n Plotting \n\n')
    ## Shiny Application
    shiny::shinyApp(
        ui <- shiny::navbarPage('Atlantis Diet Tool',
                         shiny::tabPanel('Non-Spatial',
                                  shiny::fluidRow(
                                      shiny::column(2,
                                             shiny::wellPanel(
                                                 shiny::tags$h3('Predator'),
                                                 shiny::selectInput('ycol', 'Functional Group', as.character(groups.csv$code[which(groups.csv$ispredator == 1)])),
                                                 shiny::selectInput('y1col', 'Stage', c('Juvenile', 'Adult', 'Biomass pool')),
                                                 shiny::tags$h3('Prey'),
                                                 shiny::selectInput('xcol', 'Functional Group', colnames(Ava.mat)),
                                                 shiny::selectInput('x1col', 'Stage', c('Juvenile', 'Adult')),
                                                 shiny::mainPanel("Original Value: ", shiny::verbatimTextOutput("numPoints")),
                                                 shiny::mainPanel("Current Value: ", shiny::verbatimTextOutput("CurPoints")),
                                                 shiny::numericInput("num", label = "New Value", value = 0, min = 0, max = 1, step = 0.00001)
                                             ),
                                             shiny::actionButton("do", label = "Change Value"),
                                             shiny::br(),
                                             shiny::br(),
                                             shiny::actionButton("save", label = "Write pPPREY Matrix")
                                             ),
                                      shiny::column(10,
                                                    shiny::wellPanel(
                                                               shiny::tabsetPanel(
                                                                          shiny::tabPanel('Availability matrix',#
                                                                                          shiny::plotOutput('Ava_plot', width = plotWidth, height = plotHeigh)
                                                                                          ),
                                                                          shiny::tabPanel('Effective Predation',#
                                                                                          shiny::plotOutput('Eff_pred', width = plotWidth, height = plotHeigh)
                                                                                          ),
                                                                          shiny::tabPanel('Overlap matrix',#
                                                                                          shiny::plotOutput('Overlap_plot', width = plotWidth, height = plotHeigh)
                                                                                          ),
                                                                          shiny::tabPanel('% of predation pressure',#
                                                                                          shiny::plotOutput('Press_plot', width = plotWidth, height = plotHeigh)
                                                                                          ),
                                                                          shiny::tabPanel('Total biomass prey',#
                                                                                          shiny::h3('Juvenile biomass'),
                                                                                          shiny::plotOutput('plot5', width = "100%", height = "400px"),
                                                                                          shiny::br(),
                                                                                          shiny::h3('Adult biomass'),
                                                                                          shiny::plotOutput('plot6', width = "100%", height = "400px")
                                                                                          )
                                                                      )
                                                           )
                                                    )
                                  )
                                  ),
                         shiny::tabPanel('Spatial Overlap',
                                  shiny::fluidRow(
                                      shiny::column(2,
                                             shiny::wellPanel(
                                                 shiny::selectInput('pred', 'Predator',  groups.csv$code),
                                                 shiny::selectInput('prey', 'Prey', groups.csv$code),
                                                 shiny::selectInput("layer", "Layer:", c(1 : (sediment - 1), 'Sediment')),
                                                 shiny::br(),
                                                 shiny::br(),
                                                 shiny::plotOutput('plot11', height = "300px")
                                             )
                                             ),
                                      shiny::column(10,
                                             shiny::wellPanel(
                                                 shiny::h3(shiny::textOutput("text1")),
                                                 shiny::plotOutput('plot10', height = "700px", dblclick = "plot1_dblclick",
                                                            brush = shiny::brushOpts( id = "plot1_brush", resetOnNew = TRUE
                                                                              ))
                                             )
                                             )
                                  )
                                  ),
                         shiny::tabPanel("Help",
                                  shiny::fluidPage(shiny::HTML(txtHelp)
                                            )
                                  ),
                         ## -- Exit --
                         shiny::tabPanel(
                             shiny::actionButton("exitButton", "Exit")
                         )
                         ),
        function(input, output, session) {
            ## Combine the selected variables into a new data frame
            N.mat     <- shiny::reactiveValues()
            N.mat$Ava <- Ava.mat
            pred.name <- shiny::reactive({
                nprey <- ifelse(input$x1col ==  'Juvenile', 1, 2)
                npred <- ifelse(input$y1col ==  'Juvenile', 1, ifelse(input$y1col ==  'Adult', 2, 0))
                if(npred > 0){
                    paste('pPREY', nprey, input$ycol, npred, sep = '')
                } else {
                    paste('pPREY', input$ycol, sep = '')
                }
            })
            newEntry  <- shiny::observe({
                if(input$do > 0) {
                    newval <- shiny::isolate(input$num)
                    col.ch <- shiny::isolate(which(colnames(Ava.mat) == input$xcol))
                    row.ch <- shiny::isolate(which(row.names(Ava.mat) == pred.name()))
                    shiny::isolate(N.mat$Ava[row.ch, col.ch] <- newval)
                }
            })
            shiny::observeEvent(input$save, {
                saveData(N.mat$Ava)
            })
            linex <- shiny::reactive( {
                which(sort(colnames(Ava.mat)) == input$xcol)
            })
            liney <- shiny::reactive({
                which(sort(row.names(Ava.mat)) == pred.name())
            })
            rff <- shiny::reactive({
                rff <- log(real.feed * N.mat$Ava)
                rff[!is.finite(rff)] <- 0
                t(rff)
            })
            rff2 <- shiny::reactive({
                t((real.feed * N.mat$Ava) / rowSums(real.feed * N.mat$Ava, na.rm=TRUE)) * 100
            })
            ## zoom plot
            ranges <- shiny::reactiveValues(x = NULL, y = NULL)
            shiny::observeEvent(input$plot1_dblclick, {
                brush <- input$plot1_brush
                if (!is.null(brush)) {
                    ranges$x <- c(brush$xmin, brush$xmax)
                    ranges$y <- c(brush$ymin, brush$ymax)
                } else {
                    ranges$x <- NULL
                    ranges$y <- NULL
                }
            })
            spatial <- shiny::reactive({
                gp.pred  <- with(sp.ov, sp.ov[which(Predator == input$pred && Prey == input$prey, arr.ind = TRUE)])
                gp.pred  <- dplyr::filter(data = over.tmp, .data$Predator == input$pred, .data$Prey == input$prey)
                input.layer <- as.character(ifelse(input$layer == 'Sediment', max(sp.ov$Layer, na.rm = TRUE), input$layer))
                pred.ad  <- dplyr::filter(data = sp.ov, .data$variable == input$pred, .data$Stage == 'Adult', .data$Layer == input.layer)
                pred.juv <- dplyr::filter(data = sp.ov, .data$variable == input$pred, .data$Stage == 'Juvenile', .data$Layer == input.layer)
                prey.ad  <- dplyr::filter(data = sp.ov, .data$variable == input$prey, .data$Stage == 'Adult', .data$Layer == input.layer)
                prey.juv <- dplyr::filter(data = sp.ov, .data$variable == input$prey, .data$Stage == 'Juvenile', .data$Layer == input.layer)
                ## Checking for Juveniles on the biomass pools
                ## avoiding inf problems
                if(length(prey.juv[, 6]) == 0){
                    juv.pry <- NA
                } else {
                    juv.pry <- prey.juv[, 6]
                }
                if(length(pred.juv[, 6]) == 0){
                    juv.prd <- NA
                } else {
                    juv.prd <- pred.juv[, 6]
                }
                ## Checking the gape overlap
                AoA <- ifelse(length(dplyr::filter(data = gp.pred, .data$Stg.predator == 'Adult', .data$Stg.prey == 'Adult')[, 5]) == 0, 0, as.numeric(as.character(dplyr::filter(data = gp.pred, .data$Stg.predator == 'Adult', .data$Stg.prey == 'Adult')[, 5])))
                AoJ <- ifelse(length(dplyr::filter(data = gp.pred, .data$Stg.predator == 'Adult', .data$Stg.prey == 'Juvenile')[, 5]) == 0, 0, as.numeric(as.character(dplyr::filter(data = gp.pred, .data$Stg.predator == 'Adult', .data$Stg.prey == 'Juvenile')[, 5])))
                JoA <- ifelse(length(dplyr::filter(data = gp.pred, .data$Stg.predator == 'Juvenile', .data$Stg.prey == 'Adult')[, 5]) == 0, 0, as.numeric(as.character(dplyr::filter(data = gp.pred, .data$Stg.predator == 'Juvenile', .data$Stg.prey == 'Adult')[, 5])))
                JoJ <- ifelse(length(dplyr::filter(data = gp.pred, .data$Stg.predator == 'Juvenile', .data$Stg.prey == 'Juvenile')[, 5]) == 0, 0, as.numeric(as.character(dplyr::filter(data = gp.pred, .data$Stg.predator == 'Juvenile', .data$Stg.prey == 'Juvenile')[, 5])))
                ## Merging
                ad.on.juv    <- data.frame(Land = prey.ad[, 4], Layer = prey.ad[, 1], Box = prey.ad[, 2], overlap = juv.pry  * pred.ad[, 6], Stage.prey = 'Juvenile - PREY',
                                           Stage.pred = 'Adult - PREDATOR', Gape.Overlap = AoJ)
                ad.on.ad     <- data.frame(Land = prey.ad[, 4], Layer = prey.ad[, 1], Box = prey.ad[, 2], overlap = prey.ad[, 6] * pred.ad[, 6], Stage.prey = 'Adult - PREY',
                                           Stage.pred = 'Adult - PREDATOR', Gape.Overlap = AoA)
                juv.on.juv   <- data.frame(Land = prey.ad[, 4], Layer = prey.ad[, 1], Box = prey.ad[, 2], overlap = juv.pry * juv.prd, Stage.prey = 'Juvenile - PREY',
                                           Stage.pred = 'Juvenile - PREDATOR', Gape.Overlap = JoJ)
                juv.on.ad    <- data.frame(Land = prey.ad[, 4], Layer = prey.ad[, 1], Box = prey.ad[, 2], overlap = prey.ad[, 6] * juv.prd, Stage.prey = 'Adult - PREY',
                                           Stage.pred = 'Juvenile - PREDATOR', Gape.Overlap = JoA)
                pred.tot     <- rbind(ad.on.juv, ad.on.ad, juv.on.juv, juv.on.ad)
                pred.tot$Box <- pred.tot$Box - 1
                pred.tot$overlap[is.na(pred.tot$overlap)] <- 0
                pred.tot$Rel.overlap <- with(pred.tot, ifelse(sediment == Layer & Gape.Overlap == 1 & overlap == 1, 'Sediment Layer - Gape - Spatial',
                                                       ifelse(sediment == Layer & overlap == 1 & Gape.Overlap == 0, 'Sediment Layer - Spatial - No Gape',
                                                       ifelse(Land < Layer, 'Land or Sediment',
                                                       ifelse(overlap == 0 & Gape.Overlap == 0, 'No Gape - No Spatial',
                                                       ifelse(overlap == 0 & Gape.Overlap == 1, 'Gape - No Spatial',
                                                       ifelse(overlap == 1 & Gape.Overlap == 0, 'No Gape - Spatial', 'Gape - Spatial')))))))
                overlap.pred.prey <- suppressMessages(dplyr::left_join(map, pred.tot))
            })
            dpt <- shiny::reactive({
                dpt <- rbind(dplyr::filter(depth.dst, .data$FG == input$pred), dplyr::filter(depth.dst, .data$FG == input$prey))
            })
            title <- shiny::reactive({
                paste('Realized Spatial overlap between the predator', groups.csv$long.name[which(groups.csv$code %in% input$pred)] ,'and the prey',
                      groups.csv$long.name[which(groups.csv$code %in% input$prey)], sep = ' ')
            })
            output$text1 <- shiny::renderText({
                title()
            })
            shiny::observeEvent(input$exitButton, {
                shiny::stopApp()
            })
            ## ~~~~~~~~~~~~~~~~~~~~~~ ##
            ## ~      Build Plots   ~ ##
            ## ~~~~~~~~~~~~~~~~~~~~~~ ##
            ## Effective Predation
            Eff_pred_out <- shiny::reactive({
                dat      <- reshape::melt(rff())
                dat[which(dat == 0,  arr.ind = TRUE)] <- NA
                p <- ggplot2::ggplot(data = dat, aes(x = .data$X1, y = .data$X2, fill = .data$value)) + geom_tile(colour="grey45",size=0.2)
                p <- p + scale_fill_distiller(palette = "YlGnBu", limits=c(0, max(rff(), na.rm = TRUE)), name = 'Predation \nLn()',  na.value = 'white', direction = 1)
                #p <- p + scale_fill_distiller()
                p <- p + theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
                p <- p + labs(x = 'Prey', y = 'Predator') + scale_x_discrete(position = "top")
                p <- p + annotate("rect", xmin = linex() -.5, xmax = linex() +.5, ymin = 0, ymax = ncol(rff()) + 1, alpha = .1, colour = 'goldenrod')
                p <- p + annotate("rect", xmin =  - .5, xmax = nrow(rff()) + .5, ymin = liney() - .5, ymax = liney() + .5, alpha = .1, colour = 'goldenrod')
                p
            })
            ## Overlap Matrix
            Overlap_out <- shiny::reactive({
                col.tmp <- RColorBrewer::brewer.pal(11, 'RdBu')[c(6, 11)]
                dat2 <- reshape::melt(t.o.mat)
                p <- ggplot2::ggplot(data = dat2, aes(x = .data$X1, y = .data$X2, fill = .data$value)) + geom_tile(colour= 'grey45', aes( fill = factor(.data$value)))
                p <- p + theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = 'Prey', y = 'Predator')
                p <- p + scale_x_discrete(position = "top")  + scale_fill_manual(values = col.tmp, name = 'Gape overlap', labels = c('No', 'Yes')) #+ scale_fill_manual(name = 'Gape overlap', labels = c('No', 'Yes'))
                p <- p + annotate("rect", xmin = linex() -.5, xmax = linex() +.5, ymin = 0, ymax = ncol(t.o.mat) + 1, alpha = .1, colour = 'darkorange')
                p <- p + annotate("rect", xmin = -.5, xmax = nrow(t.o.mat) + .5, ymin = liney() - .5, ymax = liney() + .5, alpha = .1, colour = 'darkorange')
                p
            })
            ## Availability matrix
            Avail_out <- shiny::reactive({
                dat.A <- reshape::melt(t(N.mat$Ava))
                dat.A[which(dat.A == 0,  arr.ind = TRUE)] <- NA
                p <- ggplot2::ggplot(data = dat.A, aes(x = .data$X1, y = .data$X2, fill = .data$value)) + geom_tile(colour="grey45",size=0.2)
                p <- p  + scale_fill_distiller(palette = "YlOrRd", limits=c(0, max(N.mat$Ava, na.rm = TRUE)), name = 'Availability', na.value = 'white', direction = 1)
                p <- p + theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = 'Prey', y = 'Predator') + scale_x_discrete(position = "top")
                p <- p +annotate("rect", xmin = linex() -.5, xmax = linex() +.5, ymin = 0, ymax = ncol(t(N.mat$Ava)) + 1, alpha = .1, colour = 'royalblue')
                p <- p + annotate("rect", xmin = -.5, xmax = nrow(t(N.mat$Ava)) + .5, ymin = liney() - .5, ymax = liney() + .5, alpha = .1, colour = 'royalblue')
                p
            })
            ## % of pressure
            Press_out <- shiny::reactive({
                dat.P <- reshape::melt(rff2())
                dat.P[which(dat.P == 0,  arr.ind = TRUE)] <- NA
                p <- ggplot2::ggplot(data = dat.P, aes(x = .data$X1, y = .data$X2, fill = .data$value)) + geom_tile(colour="grey45", size = 0.2)
                p <- p + scale_fill_distiller(palette = "RdPu",, limits=c(0, 100), name = 'Precentage of pressure', na.value = 'white', direction = 1)
                p <- p + theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = 'Prey', y = 'Predator') + scale_x_discrete(position = "top")
                p <- p + annotate("rect", xmin = linex() -.5, xmax = linex() +.5, ymin = 0, ymax = ncol(rff2()) + 1, alpha = .1, colour = 'goldenrod')
                p <- p + annotate("rect", xmin =  - .5, xmax = nrow(rff2()) + .5, ymin = liney() - .5, ymax = liney() + .5,alpha = .1, colour = 'goldenrod')
                p
            })

            ## ~~~~~~~~~~~~~~~~~~~~ ##
            ## ~      Call Plots  ~ ##
            ## ~~~~~~~~~~~~~~~~~~~~ ##
            output$Eff_pred <- shiny::renderPlot({
                Eff_pred_out()
            })
            output$Overlap_plot <- shiny::renderPlot({
                Overlap_out()
            })
            output$Ava_plot <- shiny::renderPlot({
                Avail_out()
            })
            output$Press_plot <- shiny::renderPlot({
                Press_out()
            })


            output$plot5 <- shiny::renderPlot({
                ggplot2::ggplot(data = b.juv, aes(x = .data$FG, y = log(.data$Biomass), fill=.data$FG)) +
                    geom_bar(colour="black", stat="identity") +
                    guides(fill = FALSE)+
                    xlab("Functional Groups") + ylab("Biomass [MgN] or Density [MgNm-3]")
            })
            output$plot6 <- shiny::renderPlot({
                ggplot2::ggplot(data = b.adl, aes(x = .data$FG, y = log(.data$Biomass), fill=.data$FG)) +
                    geom_bar(colour="black", stat="identity") +
                    guides(fill = FALSE)+
                    xlab("Functional Groups") + ylab("Biomass [MgN] or Density [MgNm-3]")
            })
            output$plot10 <- shiny::renderPlot({
                ggplot2::ggplot(data = spatial(), aes(x = .data$lon, y = .data$lat, group = .data$Box, fill = .data$Rel.overlap)) +
                    geom_polygon(colour = "black", size = 0.25, na.rm = TRUE) +
                    scale_fill_manual('Realized overlap\n', values = c('No Gape - No Spatial'              = 'azure1',
                                                                       'No Gape - Spatial'                 = 'royalblue',
                                                                       'Gape - No Spatial'                 = 'mistyrose',
                                                                       'Gape - Spatial'                     = 'firebrick3',
                                                                       'Sediment Layer - Gape - Spatial'    = 'goldenrod2',
                                                                       'Sediment Layer - Spatial - No Gape' = 'darkolivegreen1',
                                                                       'Land or Sediment'                   = 'grey50'
                                                                       )) +
                    facet_grid(.data$Stage.prey~.data$Stage.pred) +
                    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
                          strip.text = element_text(size = 14)) +
                    labs(x = 'Longitude', y = 'Latitude') +
                    coord_cartesian(xlim = ranges$x, ylim = ranges$y)
            })
            output$plot11 <- shiny::renderPlot({
                plot(c(1,1), -dpt()[1,2:3], xaxt = 'n', bty = 'n', ylab = 'Depth (m)', xlab = 'Functional group', type = 'l', las = 1,
                     ylim = -c(max(dpt()$Max),min(dpt()$Min)), xlim = c(0.5, 2.5), lwd = 15, col = 'royalblue',main='Depth Distribution')
                lines(c(2,2), -dpt()[2,2:3], col='goldenrod', lwd = 15)
                axis(side = 1, at= c(1,2), labels=as.character(dpt()[1:2,1]))
            })
            output$numPoints <- shiny::renderText({
                Ava.mat[which(row.names(Ava.mat) == pred.name()), which(colnames(Ava.mat) == input$xcol)]
            })
            output$CurPoints <- shiny::renderText({
                N.mat$Ava[which(row.names(Ava.mat) == pred.name()), which(colnames(Ava.mat) == input$xcol)]
            })
        }
    )
}

## ~~~~~~~~~~~~~~~~~~~~~~ ##
## ~      FUNCTIONS!!   ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~ ##
## getting the RN and SN from the NC file
##' @title Biomass from nc file
##' @param nc.file Atlantis initial condition file
##' @param groups.csv Atlantis groups file
##' @param numlayers Number of layers by box
##' @return The biomass for age class and the sturctural nitrogen by age class
##' @author Demiurgo
Bio.func <- function(nc.file, groups.csv, numlayers){
    nc.out     <- ncdf4::nc_open(nc.file)
    m.depth    <- nc.out$dim$z$len ## get the max depth to avoid problem with the bgm file
    active     <- which(groups.csv$isturnedon == 1)
    is.off     <- which(groups.csv$isturnedon == 0)
    FG         <- as.character(groups.csv$name)
    FGN        <- as.character(groups.csv$code)
    TY         <- as.character(groups.csv$grouptype)
    Biom.N     <- array(data = NA, dim = c(length(FG), max(groups.csv$numcohorts)))
    Struct     <- Biom.N
    over.sp    <-reshape::melt(matrix(NA, m.depth, nrow(numlayers)))[,1: 2]
    names.temp <- c('Layer', 'Box')
    #over.sp <- NULL
    special <- c('PWN', 'PRAWNS', 'PRAWN', 'CEP', 'MOB_EP_OTHER', 'SEAGRASS', 'CORAL', 'MANGROVE', 'MANGROVES', 'SPONGE')
    for(code in active){
        #if(code %in% is.off) next
        if(TY[code] %in% special && groups.csv$numcohorts[code] > 1){
            ## This bit is for Aged structured Biomass pools
            sed     <- ncdf4::ncatt_get(nc.out, varid = paste(FG[code], "_N1", sep = ""), attname = "insed")$value
            unit    <- ncdf4::ncatt_get(nc.out, varid = paste(FG[code], "_N1", sep = ""), attname = "units")$value
            epi     <- ncdf4::ncatt_get(nc.out, varid = paste(FG[code], "_N1", sep = ""), attname = "epibenthos")$value == 'epibenthos'
        } else {
            sed     <- ncdf4::ncatt_get(nc.out, varid = paste(FG[code], "_N", sep = ""), attname = "insed")$value
            unit    <- ncdf4::ncatt_get(nc.out, varid = paste(FG[code], "_N", sep = ""), attname = "units")$value
            epi     <- ncdf4::ncatt_get(nc.out, varid = paste(FG[code], "_N", sep = ""), attname = "bmtype")$value == 'epibenthos'
        }
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## ~                       Age structured biomass pools and biomass pool                    ~ ##
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if(groups.csv$numcohorts[code] == 1 || TY[code] %in% special){
            for(coh in 1 : groups.csv$numcohorts[code]){
                if(TY[code] %in% special && groups.csv$numcohorts[code] > 1){
                    N.tot <- ncdf4::ncvar_get(nc.out, paste(FG[code], "_N", coh, sep = ""))
                } else {
                    N.tot <- ncdf4::ncvar_get(nc.out, paste(FG[code], "_N", sep = ""))
                }
                if(all(is.na(N.tot)) || all(N.tot == 0) || sum(N.tot, na.rm = TRUE) == 0){
                    ## Getting the total volumen
                    water.t <- ncdf4::ncvar_get(nc.out, 'volume')
                    w.depth <- ncdf4::ncvar_get(nc.out, 'nominal_dz')
                    w.depth[is.na(w.depth)] <- 0
                    w.m2    <- colSums(water.t,na.rm=TRUE) / apply(w.depth, 2, function(x) max(x, na.rm = TRUE))
                    w.m2[is.infinite(w.m2)] <- NA
                    w.m2    <- sum(w.m2, na.rm=TRUE)
                    w.m3    <- sum(water.t, na.rm = TRUE)
                    if(TY[code] %in% special){
                        Biom.N[code, 1] <- ncdf4::ncatt_get(nc.out, varid = paste(FG[code], "_N", coh, sep = ""), attname = "_FillValue")$value
                    } else {
                        Biom.N[code, 1] <- ncdf4::ncatt_get(nc.out, varid = paste(FG[code], "_N", sep = ""), attname = "_FillValue")$value
                    }
                    Biom.N[code, 1] <- ifelse(sed == 1 || epi, Biom.N[code, 1] * w.m2, Biom.N[code, 1] * w.m3)
                } else {
                    if(length(dim(N.tot)) >= 3){
                        N.tot <- N.tot[, , 1]
                    } else if(unit == "mg N m-3"){
                        water.t <- ncdf4::ncvar_get(nc.out, 'volume')
                        N.tot   <- N.tot * water.t
                        Temp.N  <- N.tot
                    } else if (unit == 'mg N m-2'){
                        water.t <- ncdf4::ncvar_get(nc.out, 'volume')
                        w.depth <- ncdf4::ncvar_get(nc.out, 'nominal_dz')
                        w.depth[is.na(w.depth)] <- 0
                        w.m2    <- colSums(water.t,na.rm=TRUE) / apply(w.depth, 2, function(x) max(x, na.rm = TRUE))
                        w.m2[is.infinite(w.m2)] <- NA
                        Temp.N  <- N.tot
                        N.tot   <- sum(N.tot * w.m2, na.rm = TRUE)
                    }
                    Biom.N[code, coh] <- sum(N.tot, na.rm = TRUE)
                    Numb.tmp <- matrix(NA, m.depth, nrow(numlayers))
                    if(code == 1 && coh == 1 && !(code %in% is.off)){
                        if(length(dim(Temp.N)) == 1){
                            Numb.tmp[nrow(Numb.tmp), ] <- Temp.N
                        } else {
                            for(box in 1 : ncol(Temp.N)){
                                if(numlayers[box, 2]  == 1) next()
                                arreg           <- c(numlayers[box, 3] : 1, nrow(Numb.tmp), (numlayers[box, 3]  +  1) :(nrow(Numb.tmp) - 1))[1 : nrow(Numb.tmp)]
                                Numb.tmp[, box] <- Temp.N[arreg, box]
                            }
                        }
                        new.sp     <- as.vector(reshape::melt(Numb.tmp)[, 3])
                        over.sp    <- cbind(over.sp, new.sp)
                        names.temp <- c(names.temp, paste(FGN[code], coh, sep = '_'))
                    } else if(!(code %in% is.off)){
                        if(length(dim(Temp.N)) == 1){
                            Numb.tmp[nrow(Numb.tmp), ] <- Temp.N
                        } else {
                            for(box in 1 : ncol(Temp.N)){
                                if(numlayers[box, 2]  == 1) next()
                                arreg           <- c(numlayers[box, 3] : 1, nrow(Numb.tmp), (numlayers[box, 3]  +  1) :(nrow(Numb.tmp) - 1))[1 : nrow(Numb.tmp)]
                                Numb.tmp[, box] <- Temp.N[arreg, box]
                            }
                        }
                        new.sp     <- as.vector(reshape::melt(Numb.tmp)[, 3])
                        over.sp    <- cbind(over.sp, new.sp)
                        names.temp <- c(names.temp, paste(FGN[code], coh, sep = '_'))
                    }
                }
            }
        } else if(groups.csv$numcohorts[code] > 1 && groups.csv$isturnedon[code] == 1 && !(TY[code] %in% special)) {
            for(cohort in 1 : groups.csv$numcohorts[code]){
                StructN <- ncdf4::ncvar_get(nc.out, paste(FG[code], as.character(cohort), "_StructN", sep = ""))
                if(all(is.na(StructN))){
                    ## Some model don't have the values by box and layer,  they use the FillValue attribute
                    StructN <- ncdf4::ncatt_get(nc.out, paste(FG[code], as.character(cohort), "_StructN", sep = ""), attname = "_FillValue")$value
                }
                ReservN <- ncdf4::ncvar_get(nc.out, paste(FG[code], as.character(cohort), "_ResN", sep = ""))
                if(all(is.na(ReservN))){
                    ## Some model don't have the values by box and layer,  they use the FillValue attribute
                    ReservN <- ncdf4::ncatt_get(nc.out, paste(FG[code], as.character(cohort), "_ResN", sep = ""), attname = "_FillValue")$value
                }
                Numb    <- ncdf4::ncvar_get(nc.out, paste(FG[code], as.character(cohort), "_Nums", sep = ""))
                if(length(dim(ReservN)) > 2){
                    StructN <- StructN[, , 1]
                    ReservN <- ReservN[, , 1]
                }
                if(length(dim(Numb)) > 2){
                    Numb    <- Numb[, , 1]
                }
                if(code == 1 && cohort == 1 && (code %in% active)){
                    Numb.tmp <- matrix(NA, m.depth, nrow(numlayers))
                    for(box in 1 : ncol(Numb.tmp)){
                        if(numlayers[box, 2]  == 1) next()
                        arreg <- c(numlayers[box, 3] : 1, m.depth, (numlayers[box, 3]  +  1) : (m.depth - 1))[1 : m.depth]
                        Numb.tmp[, box] <- Numb[arreg, box]
                    }
                    new.sp     <- as.vector(reshape::melt(Numb.tmp)[, 3])
                    over.sp    <- cbind(over.sp, new.sp)
                    names.temp <- c(names.temp, paste(FGN[code], cohort, sep = '_'))
                }else if((code %in% active)){
                    Numb.tmp <- matrix(NA, m.depth, nrow(numlayers))
                    for(box in 1 : ncol(Numb)){
                        if(numlayers[box, 2]  == 1) next()
                        arreg <- c(numlayers[box, 3] : 1, m.depth, (numlayers[box, 3]  +  1) :(m.depth - 1))[1 : m.depth]
                        Numb.tmp[, box] <- Numb[arreg, box]
                    }
                    new.sp     <- as.vector(reshape::melt(Numb.tmp)[, 3])
                    over.sp    <- cbind(over.sp, new.sp)
                    names.temp <- c(names.temp, paste(FGN[code], cohort, sep = '_'))
                }
                Biom.N[code, cohort] <- sum(StructN + ReservN * Numb, na.rm = TRUE)
                Struct[code, cohort] <- max(StructN, na.rm = TRUE)
            }
        }
    }
    ncdf4::nc_close(nc.out)
    names(over.sp)    <- names.temp
    row.names(Biom.N) <- as.character(groups.csv$code)
    row.names(Struct) <- as.character(groups.csv$code)
    if(length(is.off) > 0){
        ## Remove groups that are not ON in the model
        Struct            <- Struct[ - is.off, ]
        Biom.N            <- Biom.N[ - is.off, ]
    }
    mom.t            <- over.sp[, 3 : ncol(over.sp)]
    mom.t[mom.t > 0] <- 1
    over.sp[, 3 : ncol(over.sp)] <- mom.t
    return(list(Struct, Biom.N, over.sp))
}

## ##' @title Parameter file reader
## ##' @param text Biological parametar file for Atlatnis
## ##' @param pattern Text that you are looking
## ##' @param FG Name of the functional groups
## ##' @param Vector Logic argument, if the data is on vectors or not
## ##' @param pprey Logic argument, if the data is a pprey matrix or not
## ##' @return A matrix with the values from the .prm file
## ##' @author Demiurgo
## text2num <- function(text, pattern, FG = NULL, Vector = FALSE, pprey = FALSE){
##     if(!isTRUE(Vector)){
##         text <- text[grep(pattern = pattern, text)]
##         if(length(text) == 0) warning(paste0('\n\nThere is no ', pattern, ' parameter in your file.'))
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
##         nam   <- gsub(pattern = '[[:space:]]+' ,  '|',  text[l.pat])
##         fg    <- vector()
##         pos   <- 1
##         for( i in 1 : length(nam)){
##             tmp     <- unlist(strsplit(nam[i], split = '|', fixed = TRUE))
##             if(grepl('#', tmp[1]) || (!grepl('^pPREY', tmp[1]) && pprey  == TRUE)) next
##             fg[pos] <- tmp[1]
##             if(pos == 1) {
##                 #t.text  <- gsub('"[[:space:]]"', ' ',  text[l.pat[i] + 1])
##                 t.text <- gsub('+[[:space:]]+', ' ',  text[l.pat[i] + 1])
##                 pp.mat <- matrix(as.numeric(unlist(strsplit(t.text, split = ' +', fixed = FALSE))), nrow = 1)
##                 pos    <- pos + 1
##             } else {
##                 #t.text <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", text[l.pat[i] + 1], perl=TRUE)
##                 t.text <- gsub('+[[:space:]]+', ' ',  text[l.pat[i] + 1])
##                 pp.tmp <- matrix(as.numeric(unlist(strsplit(t.text, split = ' ', fixed = TRUE))), nrow = 1)
##                 if(ncol(pp.mat) != ncol(pp.tmp)) stop('\nError: The pPrey vector for ', tmp[1], ' has ', ncol(pp.tmp), ' columns and should have ', ncol(pp.mat))
##                 pp.mat <- rbind(pp.mat, pp.tmp)
##                 pos    <- pos + 1
##             }
##         }
##         if(all(is.na(pp.mat[, 1]))) pp.mat <- pp.mat[, - 1]
##         row.names(pp.mat)                  <- fg
##         return(pp.mat)
##     }
## }

##' @title Calculate the gape limitation for each functional group
##' @param groups.csv Atlantis group file
##' @param Struct Structural weight by age group
##' @param Biom.N Biomass by age groups
##' @param prm Atlantis paramter file
##' @return the limits for the prey based on the predator gape size and prey size
##' @author Demiurgo
gape.func <- function(groups.csv, Struct, Biom.N, prm){
    ## Gape size and adult and young age
    KLP                     <- text2num(prm, 'KLP', FG = as.character(groups.csv$code))
    KUP                     <- text2num(prm, 'KUP',  FG = as.character(groups.csv$code))

    if(length(unique(KLP[,1])) != nrow(KLP)){
        stop('You have repeated values of KLP for ', KLP[which(duplicated(KLP[, 1])), 1], ', fix that and run again the tool')
    }
    if(length(unique(KUP[,1])) != nrow(KUP)){
        stop('You have repeated values of KUP for ', KLP[which(duplicated(KLP[, 1])), 1], ', fix that and run again the tool')
    }

    if(nrow(KUP) != nrow(KLP)) {
        warning(cat('You have', nrow(KUP),  ' values of KUP_ and ', nrow(KLP),  'values of KLP_. Its a good idea to fix that. For now the extra values will not be considered'))
        if(nrow(KUP) > nrow(KLP)){
            KUP <- KUP[which(KLP[, 1] %in% KUP[, 1]), ]
        } else {
            KLP <- KLP[which(KUP[, 1] %in% KLP[, 1]), ]
        }
    }
    age             <- text2num(prm, '_age_mat', FG = as.character(groups.csv$code))
    names(age)      <- c('FG', 'Age.Adult')
    Gape            <- data.frame(FG = as.factor(KLP$FG), KLP = KLP$Value, KUP = KUP$Value, adult.Min = NA, adult.Max = NA, juv.Min = NA,
                                  juv.Max = NA, JminS = NA, JmaxS = NA, AminS = NA, AmaxS = NA)
    Gape            <- dplyr::left_join(Gape, age, by = 'FG')
    Gape$Age.Young  <- Gape$Age.Adult - 1
    Gape$Age.Young  <- ifelse(Gape$Age.Young == 0,  1, Gape$Age.Young)
    ## Pre-Calculations
    Biom.N        <- Biom.N[order(row.names(Biom.N)), ]
    Struct        <- data.frame(as.factor(row.names(Struct)), Struct)
    names(Struct) <- c('FG',  paste0('Age_', 1 : (ncol(Struct) - 1)))
    Gape          <- dplyr::left_join(Gape, Struct, by = 'FG')
    ages.pos      <- grep('Age_', colnames(Gape))
    Gape$juv.Min  <- Gape[, ages.pos[1]] * Gape$KLP
    for( i in 1 : nrow(Gape)){
        if(all(is.na(Gape[i, ages.pos]))) next()
        ## Gape limitation as a predator
        Gape$adult.Min[i]  <- Gape[i, ages.pos[Gape$Age.Adult[i]]] * Gape$KLP[i]
        Gape$adult.Max[i]  <- Gape[i, ages.pos[sum(!is.na(Gape[i, ages.pos]))]] * Gape$KUP[i]
        Gape$juv.Max[i]    <- Gape[i, ages.pos[Gape$Age.Young[i]]] * Gape$KLP[i]
        ## Gape limitation as prey
        Gape$JminS[i]      <- Gape[i, ages.pos[1]]
        Gape$AminS[i]      <- Gape[i, ages.pos[Gape$Age.Adult[i]]]
        Gape$JmaxS[i]      <- Gape[i, ages.pos[Gape$Age.Young[i]]]
        Gape$AmaxS[i]      <- Gape[i, ages.pos[sum(!is.na(Gape[i, ages.pos]))]]
    }
    return(list(Gape, age))
}
##' @title Overlap matrix
##' @param Ava.mat Availavility matrix (or pPREY matrix from the parameter file
##' @param Gape Limit of prey by Gape size
##' @return Overlap Matrix based only in the gape limitation
##' @author Demiurgo
Over.mat.func <- function(Ava.mat, Gape){
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~        Overlap-Matrix    ~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    Over.mat <- Ava.mat * 0
    Prey     <- colnames(Ava.mat)
    Pred     <- row.names(Ava.mat)
    for( py in 1: length(Prey)){
        for( pd in 1: length(Pred)){
            c.pred      <- unlist(strsplit(Pred[pd],'pPREY'))[2]
            predator    <- gsub(pattern = "[[:digit:]]+", '\\1', c.pred)
            a.pred.prey <- as.numeric(unlist(strsplit(c.pred, predator)))
            pry.loc     <- which(Gape$FG %in% Prey[py])
            prd.loc     <- which(Gape$FG %in% predator)
            if(length(pry.loc) == 0 || is.na(a.pred.prey)){
                Over.mat [pd, py] <- 1
            } else {
                if(a.pred.prey[1] == 1){
                    ## Young Predator
                    if(a.pred.prey[2] == 1){
                        ## Young Prey
                        Over.mat [pd, py]  <- ifelse(Gape$JminS[pry.loc] >= Gape$juv.Min[prd.loc],
                                              ifelse(Gape$JminS[pry.loc] <= Gape$juv.Max[prd.loc], 1, 0),
                                              ifelse(Gape$JmaxS[pry.loc] >= Gape$juv.Min[prd.loc],
                                              ifelse(Gape$JmaxS[pry.loc] <= Gape$juv.Max[prd.loc], 1, 1), 0))
                        if(is.na(Over.mat[pd, py])) Over.mat[pd, py] <- 1
                    } else {
                        ## Adult Prey
                        Over.mat [pd, py]  <- ifelse(Gape$AminS[pry.loc] >= Gape$juv.Min[prd.loc],
                                              ifelse(Gape$AminS[pry.loc] <= Gape$juv.Max[prd.loc], 1, 0),
                                              ifelse(Gape$AmaxS[pry.loc] >= Gape$juv.Min[prd.loc],
                                              ifelse(Gape$AmaxS[pry.loc] <= Gape$juv.Max[prd.loc], 1, 1), 0))
                        if(is.na(Over.mat[pd, py])) Over.mat[pd, py] <- 1
                    }
                } else {
                    ## Adult Predator
                    if(a.pred.prey[2] == 1){
                        ## Young Prey
                        Over.mat [pd, py]  <- ifelse(Gape$JminS[pry.loc] >= Gape$adult.Min[prd.loc],
                                              ifelse(Gape$JminS[pry.loc] <= Gape$adult.Max[prd.loc], 1, 0),
                                              ifelse(Gape$JmaxS[pry.loc] >= Gape$adult.Min[prd.loc],
                                              ifelse(Gape$JmaxS[pry.loc] <= Gape$adult.Max[prd.loc], 1, 1), 0))
                        if(is.na(Over.mat[pd, py])) Over.mat[pd, py] <- 1
                    } else {
                        ## Adult Prey
                        Over.mat [pd, py]  <- ifelse(Gape$AminS[pry.loc] >= Gape$adult.Min[prd.loc],
                                              ifelse(Gape$AminS[pry.loc] <= Gape$adult.Max[prd.loc], 1, 0),
                                              ifelse(Gape$AmaxS[pry.loc] >= Gape$adult.Min[prd.loc],
                                              ifelse(Gape$AmaxS[pry.loc] <= Gape$adult.Max[prd.loc], 1, 1), 0))
                        if(is.na(Over.mat[pd, py])) Over.mat[pd, py] <- 1
                    }
                }
            }
        }
    }
    return(Over.mat)
}
##' @title Biomass by age class
##' @param Biom.N Biomass by Cohort
##' @param age Age of maturation
##' @param Over.mat Overlap matrix
##' @return Biomass by adult and juveniles
##' @author Demiurgo
Bio.ages <- function(Biom.N, age, Over.mat){
    ## total biomasss by Juv and Adults
    browser()
    Biom.N  <- Biom.N[order(row.names(Biom.N)), ]
    fg      <- row.names(Biom.N)
    bio.juv <- bio.adl <- matrix(NA, ncol = 2, nrow = nrow(Biom.N))
    for( i in 1 : nrow(Biom.N)){
        l.age <- which(age$FG  ==  fg[i])
        if(length(l.age) !=  0){
            bio.juv[i, ] <- c(fg[i], sum(Biom.N[i, 1 : (age$Age.Adult[l.age] - 1)], na.rm = TRUE))
            bio.adl[i, ] <- c(fg[i], sum(Biom.N[i, age$Age.Adult[l.age] : ncol(Biom.N)], na.rm = TRUE))
        } else {
            ## For Biomass pool,  only adult
            bio.juv[i, ] <- c(fg[i], sum(Biom.N[i, 1], na.rm = TRUE))
            bio.adl[i, ] <- c(fg[i], sum(Biom.N[i, 1], na.rm = TRUE))
        }
    }
    ## Sort based on the order of the prey in the Availavility matrix   ##
    or.prey <- match(colnames(Over.mat), bio.juv[, 1])
    bio.juv <- bio.juv[or.prey, ]
    bio.adl <- bio.adl[or.prey, ]
    return(list(bio.juv, bio.adl))
}
##' @title Save pPREY matrix
##' @param data matrix of pPrey modified
##' @return a txt file of the pPREY matrix
##' @author Demiurgo
saveData <- function(data) {
    fileName <- sprintf("%s_NewpPrey.txt", as.integer(Sys.time()))
    rows     <- row.names(data)
    cols     <- c('##    ', colnames(data))
    nprey    <- ncol(data)
    sink(fileName)
    cat(cols)
    for( i in 1 : length(rows)){
        cat(paste('\n', rows[i], '  ', nprey, '\n', sep = ''))
        cat(data[i, ])
    }
    sink()
}
##' @title Function to get the vertices of the bgm map for Atlantis
##' @param bgm.file BGM file for atlantis
##' @return A dataframe con witht 3 columns,  first the box_id,  then the latitude and the longitude
##' @author Demiurgo,  based on Shane function \href{https://github.com/shanearichards/shinyrAtlantis}{shinyrAtlantis}
make.map <- function(bgm.file){
    bgm          <- readLines(bgm.file)
    numboxes     <- as.numeric(gsub('nbox', '', grep("nbox", bgm, value = TRUE)))
    proj         <- gsub("projection[[:space:]]+", "", grep("projection", bgm, value = TRUE))
    if(!grepl('\\+', proj)){
        proj <- gsub(' ', ' \\+', proj)
        proj <- gsub('proj', '\\+proj', proj)
    }
    proj <- gsub('#', '', proj, ignore.case = TRUE) ## some bgm file have the '#' symbol at the begining
    map.vertices <- data.frame()
    for(i in 1 : numboxes){
        txt.find <- paste("box", i - 1, ".vert", sep = "")
        j <- grep(txt.find, bgm)
        for (jj in 1 : length(j)) {
            text.split <- unlist(stringr::str_split(
                gsub(pattern = "[\t ]+", x = bgm[j[jj]], replacement = " "), " "))
            if (text.split[1] == txt.find) {
                map.vertices <- rbind(map.vertices, cbind(i - 1, as.numeric(text.split[2]),
                                                          as.numeric(text.split[3])))
            }
        }
    }
    ## some datum are not included on the proj4 poject (NZ and AU)
    proj = gsub('\\+datum=NZGD2000', '', proj)
    ## in case you use alb for albers equal area (aes)
    proj = gsub('alb', 'aea', proj)
    ## in case someone is using grs and no GRS. case sensitive for R - proj4
    proj = gsub('grs', 'GRS', proj)
    ## In case the file have some spaces after or before any variable setting
    proj = gsub(' = ', '=', proj); proj = gsub('= ', '=', proj); proj = gsub(' =', '=', proj)
    ## Convert latlon coordinates!
    latlon    <- proj4::project(map.vertices[, 2 : 3], proj = proj, inverse = T)
    map       <- data.frame(Box = map.vertices[, 1],
                            lat   = latlon$y,
                            lon   = latlon$x)
    return(map)
}
##' @title text separatror
##' @param ortext Original text
##' @return data.frame with 3 columns with names and stage of the predator and the prey
##' @author Demiurgo
sepText <- function(ortext){
    ortext      <- as.vector(ortext)
    text        <- as.character(ortext[2])
    prey        <- as.character(ortext[1])
    c.pred      <- unlist(strsplit(text, 'pPREY'))[2]
    predator    <- gsub(pattern = "[[:digit:]]+", '\\1', c.pred)
    stage.prey  <- as.numeric(unlist(strsplit(c.pred, predator)))
    st.pred     <- ifelse(stage.prey[1] != 1 || is.na(stage.prey[1]), "Adult", "Juvenile")
    st.prey     <- ifelse(stage.prey[2] != 1 || is.na(stage.prey[2]), "Adult", "Juvenile")
    out         <- data.frame(Predator = predator,    Prey = prey,
                              Stg.predator = st.pred, Stg.prey = st.prey, Overlap = ortext[3])
    return(out)
}
##' @title Num layer  per box
##' @param bgm.file BGM file for the atlantis model
##' @param cum.depths Cummulativce depths of the model
##' @return dataframe with the number of layer per box including the box that are island
##' @author Demiurgo
find.z <- function(bgm.file, cum.depths){
    bgm         <- readLines(bgm.file)
    numboxes    <- as.numeric(gsub('nbox', '', grep("nbox", bgm, value = TRUE)))
    box.indices <- vector()
    for(i in 1 : numboxes){ # box depth
        box.indices[i] <- grep(paste("box", i - 1, ".botz", sep = ""), bgm)
    }
    z.tmp         <- strsplit(bgm[box.indices], "\t")
    z             <- as.numeric(sapply(z.tmp,`[`,2)) # - depth of water column
    depth.inf     <- data.frame(depth = z, is.Islan = ifelse(z >= 0, 1, 0))
    max.numlayers <- length(cum.depths) - 1 # maximum number of water layers
    ## calculate the number of water layers
    box.numlayers <- rep(0, numboxes) # vector containing number of water layers
    for (i in 1: numboxes) {
        box.numlayers[i] <- sum(depth.inf$depth[i] < - cum.depths)
    }
    max.numlayers    <- length(cum.depths)              ## maximum number of water layers
    box.numlayers    <- pmin(box.numlayers, max.numlayers) # bound by maximum depth
    depth.inf$numlay <- box.numlayers
    return(depth.inf)
}
