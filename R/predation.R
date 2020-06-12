##' This function is useful for dynamically evaluating and visualizing predator-prey
##'     interaction from Atlantis. Providing information of the predation trough
##'     time for biomass pool and for age class functional groups. This tool can helps to
##'     analyze the predatory pressure from the predator and from the prey
##'     perspective.
##' @title Predation analysis
##' @param biom.file Character string with the path to \emph{biomass output}
##'     file from the Atlantis simulation. Usually this
##'     output has the name "output\bold{\emph{_Model_}}BiomIndx.txt", were
##'     \bold{\emph{Model}} is the name of your Atlantis model.
##' @param groups.csv Character string with the path to the Groups \emph{*.csv}
##'     file (Atlantis input file).
##' @param diet.file Character string with the path to the \emph{diet output}
##'     file. This file contain the diets of each functional group at each recorded
##'     time step. If the Atlantis simulation is for several years, it is highly
##'     recommended a low frequency recording periodicity of this output file
##'     (toutinc). General high frequency engravings very large files and are very difficult
##'     to handle in R.
##' @param age.biomass Character string with the path to the output file
##'     \bold{\emph{biomass by age}}. This is a non-standard Atlantis output file,
##'     make sure to turn on the flag \bold{\emph{"flag_age_output"}} in the
##'     \emph{run.prm} file if you want Atlantis to write this this output file.
##' @return The output of this function is a shiny::reactive html, useful for dynamically evaluating and visualizing
##'     predator-prey relationships from an Atlantis output. The outputs of the predation analysis
##'     is divided in three parts:
##'    \itemize{
##'    \item \bold{Biomass}: This panel is divided in two sub-panels:
##'     \itemize{
##'      \item \emph{Total Biomass} which displays the changes in total biomass for
##'     each functional group.
##'      \item \emph{Relative Biomass} displays the variation of the
##'     total biomass relative to the initial biomass (\eqn{B_{0}}).
##'     }
##'    \item \bold{Predation through time}: This panel
##'     displays in the upper section the proportion of prey ingested by a specific
##'     predator during the
##'     simulation. The lower section is the proportion of
##'     biomass that is predated from the selected functional group by different
##'     predators at each time step during the entire simulation. This function has
##'     different option that includes:
##'  \itemize{
##'     \item \bold{Threshold} this option filters the
##'     predation values, only the values of proportion of predation greater than
##'     this limit are displayed;
##'     \item \bold{Scale} this option
##'     allows you to normalize the predation values between 0 and 1 for each time
##'     step. Otherwise, predation values are weighted by the predator's biomass. This
##'     allows to evaluate the impact of the predator on each prey at each time
##'     step.
##'      }
##'    \item \bold{Predation by age group}: This panel
##'     visualize the predation by an age group and it is composed by two sub-panels:
##'    \itemize{
##'     \item \bold{Total Biomass} shows in the upper panel the time series of biomass of
##'     the functional group selected. The lower section shows the proportion of each prey in the diet of the selected
##'     predator at a given time step. The time step can be changed allowing the user
##'     to explore the diets of the functional at any point in the simulation; and
##'     \item \bold{Biomass by age}: This sub-panel shows
##'     the biomass of each prey consumed by each age class. This allows visualizing
##'     amount of biomass consumed from each prey by the predator at each time step. The upper panel
##'     shows the biomass of the prey throughout the simulation, indicating the
##'     selected time step.
##'      }
##'}
##' @import stats utils grDevices ggplot2 graphics
##' @author Demiurgo
##' @export
predation <- function(biom.file, groups.csv, diet.file, age.biomass = NULL ){
    txtHelp <- "<h2>Summary</h2>"
    txtHelp <- paste(txtHelp, "<p>This program is useful for dynamically evaluating and visualizing predator-prey relationships from Atlantis</p>")
    txtHelp <- paste(txtHelp, "<h3>Details</h3>")
    txtHelp <- paste(txtHelp, "<p>The outputs of this program are divided into three parts: 1) <b>Biomass</b>, 2) <b>Predation through time</b>,  and 3) <b>Predation by age group</b>.</p>")
    txtHelp <- paste(txtHelp, "<p><b>1) Biomass:</b> This panel displays the changes in total biomass for each functional group. Also, the sub-panel displays the variation of the total biomass relative to B0.</p>")
    txtHelp <- paste(txtHelp, "<p><b>2) Predation through time:</b> This panel displays in the upper section the proportion of prey ingested by a specific predator (<i>Functional Group</i>) and stock (<i>Stocks</i>) during the entire Atlantis simulation. The lower plot section is the proportion of biomass that is predated from the selected functional group by different predators at each time step</p>")
    txtHelp <- paste(txtHelp, "<p>The <i>Threshold</i> option filters the predation values, only the values of proportion of predation greater than this limit are displayed (values from 0 to 1).</p>")
    txtHelp <- paste(txtHelp, "<p>The option scale (default <b>scale = TRUE</b>) allows you to integrate predation values between 0 and 1 for each time step. Otherwise, predation values are weighed by the predator's biomass. This allows to evaluate the impact of the predator on each prey at each time step.</p>")
    txtHelp <- paste(txtHelp, "<p><b>2) Predation by Age group</b> This panel visualize the predation by age group. The first sub-panel (<b>Total Biomass</b>) shows the proportion of each prey in the diet of the selected predator at a given time step (lower plot section)</p>")
    txtHelp <- paste(txtHelp, "<p>The upper plot is the variation of the total biomass of the functional group during the Atlantis simulation. The second sub-panel (<b>Biomass by Age</b>) shows the biomass of each prey consumed by each age class. This allows visualizing the impact (in biomass) of the predator in each of the prey. The upper panel shows the biomass of the prey throughout the simulation, indicating the selected time step.</p>")
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


    ## Reading Data
    cur.dat   <- data.frame(data.table::fread(biom.file, header = TRUE, sep = ' ', showProgress = FALSE))
    diet.data <- data.frame(data.table::fread(diet.file, header=TRUE, sep = ' ', showProgress = FALSE))
    grp       <- utils::read.csv(groups.csv)
    grp.prd   <- grp[grp$IsTurnedOn == 1 & grp$NumCohorts > 1, ]$Code
    grp       <- grp[grp$IsTurnedOn == 1, ]$Code
    sub.cur   <- cbind(cur.dat[c('Time', as.character(grp))])
    sp.name   <- c("Time", paste0('Rel', grp),"PelDemRatio", "PiscivPlankRatio")
    sp.name2  <- c("Time", as.character(grp))
    biom.tot  <- reshape::melt(cur.dat[, sp.name2], id.vars = 'Time')
    rel.bio   <- reshape::melt(cur.dat[, sp.name], id.vars = 'Time')
    if(any('Updated' == colnames(diet.data))){
        rem       <- which(names(diet.data) == 'Updated')
        diet.data <- diet.data[,  - rem]
    }
    ## Preparing the data
    ## checking for predator or groups
    if(any(colnames(diet.data) == 'Group')){
        colnames(diet.data)[which(colnames(diet.data) == 'Group')] <- 'Predator'
    }
    predators         <- as.character(unique(diet.data$Predator))
    time              <- unique(diet.data$Time)
    stocks            <- unique(diet.data$Stock)
    sub.cur           <- cbind(Time = sub.cur$Time, sub.cur[, names(sub.cur) %in% predators])
    new.bio           <- reshape::melt(sub.cur, id = c('Time'))    ## current Biomass of the functional groups
    colnames(new.bio) <- c('Time', 'Predator', 'Biomass') ## This is the current Biomass
    new.bio$Predator  <- as.character(new.bio$Predator)   ## Removing factors
    ## Thes bit was done by Sieme,  she was trying to look for the impact of each age class on predation
    if(!is.null(age.biomass)){
        age.gr.pred   <- data.frame(data.table::fread(age.biomass, header=TRUE, sep = ' ', showProgress = FALSE))
        col.nam       <- stringr::str_extract(names(age.gr.pred), "[aA-zZ]+")
    }

    ## Colors
    g.col    <- c(RColorBrewer::brewer.pal(9, "BuPu") [2 : 9], RColorBrewer::brewer.pal(9, "BrBG")[1 : 3],  RColorBrewer::brewer.pal(9, "OrRd") [2 : 9])
    g.col    <- data.frame(grp, col = grDevices::colorRampPalette(g.col)(length(grp)))

    ## Start the Shiny application
    shiny::shinyApp(
        ## Create the different tabs
        ui <- shiny::navbarPage("Predation",
                         shiny::tabPanel('Biomass',
                                  shiny::tabsetPanel(
                                      shiny::tabPanel('Total Biomass',
                                               shiny::plotOutput('plot1', width = "100%", height = "1000px")
                                               ),
                                      shiny::tabPanel('Relative Biomass',
                                               shiny::plotOutput('plot1B', width = "100%", height = "1000px")
                                               )
                                  )
                                  ),
                         shiny::tabPanel('Predation through time',
                                  shiny::fluidRow(
                                      shiny::column(2,
                                             shiny::wellPanel(
                                                 shiny::selectInput('FG', 'Functional Group :', as.character(grp)),
                                                 shiny::selectInput('Stocks', 'Stocks :', stocks),
                                                 shiny::numericInput("Thr", "Threshold :", min = 1e-16,  max = 1, value = 1e-3, step = 0.001),
                                                 shiny::checkboxInput(inputId = "scl", label = shiny::strong("scaled to 1"), value = TRUE),
                                                 shiny::checkboxInput(inputId = "merT", label = shiny::strong("Melt time step"), value = FALSE)
                                             )
                                             ),
                                      shiny::column(10,
                                             shiny::plotOutput('plot2', width = "100%", height = "400px"),
                                             shiny::plotOutput('plot3', width = "100%", height = "400px")
                                             )
                                  )
                                  ),
                         shiny::tabPanel('Predation by Age group',
                                  shiny::tabsetPanel(
                                      shiny::tabPanel('Total Biomass',
                                               shiny::fluidRow(
                                                   shiny::column(2,
                                                          shiny::wellPanel(
                                                              shiny::selectInput('FG2', 'Functional Group :', as.character(grp)),
                                                              shiny::sliderInput("Time", "Simulation Time :", min = min(time),  max = max(time), value = min(time), step = diff(time)[1]),
                                                              shiny::selectInput('Stocks2', 'Stocks :', stocks),
                                                              shiny::numericInput("Thr2", "Threshold :", min = 1e-16,  max = 1, value = 1e-3, step = 0.001)
                                                          )
                                                          ),
                                                   shiny::column(10,
                                                          shiny::plotOutput('plot4A', width = "100%", height = "300px"),
                                                          shiny::plotOutput('plot4B', width = "100%", height = "500px")
                                                          )
                                               )
                                               ),
                                      shiny::tabPanel('Biomass by Age',
                                               shiny::fluidRow(
                                                   shiny::column(2,
                                                          shiny::wellPanel(
                                                              shiny::selectInput('FG3', 'Functional Group :', as.character(grp.prd)),
                                                              shiny::sliderInput("Time3", "Simulation Time :", min = min(time),  max = max(time), value = min(time), step = diff(time)[1]),
                                                              shiny::selectInput('Stock3', 'Stocks :', stocks),
                                                              shiny::numericInput("Thr3", "Threshold :", min = 1e-16,  max = 1, value = 1e-3, step = 0.001)
                                                          )
                                                          ),
                                                   shiny::column(10,
                                                          shiny::plotOutput('plot5', width = "100%", height = "400px"),
                                                          shiny::plotOutput('plot6', width = "100%", height = "400px")
                                                          )
                                               )
                                               )
                                  )),
                         shiny::tabPanel("Help",
                                  shiny::fluidPage(shiny::HTML(txtHelp)
                                            )
                                  ),
                         ## -- Exit --
                         shiny::tabPanel(
                             shiny::actionButton("exitButton", "Exit")
                         )
                         ),
        ## Link the input for the different tabs with your original data
        ## Create the plots
        function(input, output, session){
            age.diet <- shiny::reactive({
                age.diet <- diet.data[diet.data$Predator ==  input$FG2 & diet.data$Stock == as.numeric(input$Stocks2) & diet.data$Time == as.numeric(input$Time), ]
                age.diet <- age.diet[,  -which(names(age.diet) %in% c('Predator','Time', 'Stock'))]
                age.diet <- age.diet[, (colSums(age.diet, na.rm = TRUE) > input$Thr2)]
                age.diet <- reshape::melt(age.diet, id = 'Cohort') ## Biomass of the predator, so we have an idea of the total pressure
                age.diet[age.diet$value > input$Thr2, ]
            })

            age.diet.bio <- shiny::reactive({
                new.bio[new.bio$Predator == input$FG2, ]
            })

            if(!is.null(age.biomass)){
                age.bio.gr <- shiny::reactive({
                    age.gr          <- age.gr.pred[, col.nam %in% c('Time', input$FG3)]
                    age.gr          <- age.gr[age.gr$Time == input$Time3, ]
                    age.gr          <- reshape::melt(age.gr, id.vars = c("Time"))
                    age.gr$variable <- as.numeric(stringr::str_extract(age.gr$variable, "[0-9]+"))
                    names(age.gr)   <- c('Time', 'Cohort', 'Biomass')
                    ## diet
                    diet <- diet.data[diet.data$Predator == input$FG3 & diet.data$Stock == input$Stock3 & diet.data$Time == input$Time3,  ]
                    diet <- diet[, -which(names(diet) %in% c('Time', 'Predator', 'Stock'))]
                    diet <- diet[, (colSums(diet, na.rm = TRUE) > input$Thr3)]
                    diet <- reshape::melt(diet, id = 'Cohort')
                    diet <- diet[diet$value > input$Thr3, ]
                    out  <- dplyr::left_join(age.gr, diet,  by = 'Cohort')
                    out$eff.pred <- out$Biomass * out$value
                    out <- out[complete.cases(out), ]
                    out
                })
            }

            if(!is.null(age.biomass)){
                prey.bio <- shiny::reactive({
                    prey <- as.character(unique(age.bio.gr()$variable))
                    prey <- biom.tot[biom.tot$variable %in% prey, ]
                })
            }

            predator <- shiny::reactive({
                predator        <- diet.data[diet.data$Predator ==  input$FG & diet.data$Stock == as.numeric(input$Stocks), ]
                predator        <- predator[,  -which(names(predator) %in% c('Predator', 'Cohort', 'Stock'))]
                predator        <- predator[, (colSums(predator, na.rm = TRUE) > input$Thr)]
                predator[which(predator < input$Thr,  arr.ind = TRUE)] <- NA
                predator        <- reshape::melt(predator, id.vars = "Time", na.rm = TRUE)
                predator        <- dplyr::left_join(predator, new.bio[new.bio$Predator == input$FG, ],  by = 'Time')
                predator$consum <- with(predator, value * Biomass)
                if(input$merT) predator$Time <- predator$Time / diff(unique(predator$Time))[1]
                predator        <- as.data.frame(predator)

            })
            shiny::observeEvent(input$exitButton, {
                shiny::stopApp()
            })
            prey <- shiny::reactive({
                prey          <- diet.data[, names(diet.data) %in% c('Time', 'Predator', input$FG)]
                prey          <- prey[prey[, input$FG] > 0, ] ## removing zero values, faster now
                prey          <- dplyr::left_join(prey, new.bio, by = c('Predator', 'Time')) ## Biomass of the predator, so we have an idea of the total pressure
                prey$eff.pred <- prey[, input$FG] * prey$Biomass
                trh.max       <- max(prey$eff.pred) * input$Thr
                if(input$merT) prey$Time <- prey$Time / diff(unique(prey$Time))[1]
                prey          <- prey[prey$eff.pred > trh.max, ]
            })
            output$plot1 <- shiny::renderPlot({
                mycol  <- c(RColorBrewer::brewer.pal(8, "Dark2"), c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
                mycol  <- grDevices::colorRampPalette(mycol)
                plot   <- ggplot2::ggplot(data = biom.tot, aes(x = .data$Time, y = .data$value)) +
                    geom_line(colour = 'darkorange3') + ggplot2::facet_wrap(~ .data$variable, ncol = 4,  scale = 'free_y') + theme_bw()+
                    scale_color_manual(values = mycol(2))
                plot <- update_labels(plot, list(x = 'Time step', y = 'Biomass (tons)'))
                plot
            })
            output$plot1B <- shiny::renderPlot({
                plot <- ggplot2::ggplot(data = rel.bio, aes(x = .data$Time, y = .data$value)) +
                    geom_line(colour = 'darkorange3', na.rm = TRUE) + ggplot2::facet_wrap( ~ .data$variable, ncol = 4) +
                    theme_bw() + ylim(0, 2) +
                    annotate('rect', xmin =  - Inf, xmax = Inf, ymax = 1.5, ymin = 0.5, alpha = .1, colour = 'royalblue', fill = 'royalblue')
                plot <- update_labels(plot, list(x = 'Time step', y = 'Relative Biomass (Bt/B0)'))
                plot
            })
            output$plot2 <- shiny::renderPlot({
                ## Colors
                mycol      <- c(RColorBrewer::brewer.pal(9, "BrBG")[1 : 3],  RColorBrewer::brewer.pal(9, "OrRd") [2 : 9])
                colorpp    <- grDevices::colorRampPalette(mycol)(length(unique(predator()$variable)))
                ggplot2::ggplot(data = predator(), aes(x = .data$Time, y = .data$consum, fill = .data$variable, width = 1)) + geom_bar(stat = "identity", position = ifelse(input$scl, 'fill', 'stack'), na.rm = TRUE) + scale_fill_manual(values = colorpp, name = 'Prey') +
                    labs(list(title = paste('Predator -', input$FG), x = 'Time step', y = ifelse(input$scl, 'Proportion', 'Biomass [tons]'), colour = 'Prey'))
            })
            output$plot3 <- shiny::renderPlot({
                shiny::validate(
                    shiny::need(length(prey()$eff.pred) != 0,  'Apparently this functional group has no predators.')
                )
                mycol      <- c(RColorBrewer::brewer.pal(9, "YlGnBu")[1 : 3], RColorBrewer::brewer.pal(9, "BuPu") [2 : 9])
                colorpp    <- grDevices::colorRampPalette(mycol)(length(unique(prey()$Predator)))
                ggplot2::ggplot(data = prey(), aes(x = .data$Time, y = .data$eff.pred, fill = .data$Predator, width = 1)) + geom_bar(stat = "identity", position = ifelse(input$scl, 'fill', 'stack'), na.rm = TRUE) + scale_fill_manual(values = colorpp) +
                    labs(list(title = paste('Prey -', input$FG), x = 'Time step', y = ifelse(input$scl, 'Proportion', 'Biomass [tons]')))
            })

            output$plot4A <- shiny::renderPlot({
                with(age.diet.bio(), plot(Time, Biomass, ylab = 'Biomass (tons)', xlab = 'Time step', bty = 'n', type = 'l',
                                          ylim = range(Biomass), las = 1, main = paste0('Biomass  - ', input$FG2)))
                with(age.diet.bio(), points(Time[Time == input$Time], Biomass[Time == input$Time], pch = 19, col = 'firebrick3', cex = 1.3))
            })
            output$plot4B <- shiny::renderPlot({
                color.pp <- as.character(g.col$col[which(g.col$grp %in% levels(age.diet()$variable))])
                ggplot2::ggplot(data = age.diet(), aes(x = .data$Cohort, y = .data$value, fill = .data$variable, width = .75)) + geom_bar(stat = "identity", position = 'fill') + scale_fill_manual(values = color.pp) +
                    labs(list(title = paste('Predator  -', input$FG2, 'on Time step :', input$Time), x = 'AgeGroup', y = 'Proportion'))
            })

            output$plot5 <- shiny::renderPlot({
                shiny::validate(
                    shiny::need(age.biomass != '',  'To Display this plot, please provide the non-standard atlantis output \'Biomass by age\' (default : age.biomass = NULL). Make sure to put the flag \'flag_age_output\' in the run.prm file if you want to look at this output')
                )
                color.pp <- as.character(g.col$col[which(g.col$grp %in% levels(age.bio.gr()$variable))])
                plot <- ggplot2::ggplot(data = prey.bio(), aes(x = .data$Time, y = .data$value)) +
                    geom_line(colour = 'firebrick3') + ggplot2::facet_wrap(~ .data$variable, ncol = 5,  scale = 'free_y') + theme_bw()+
                    scale_color_manual(values = color.pp)  + labs(list(title = paste('Prey(s) for', input$FG3, 'on Time step :', input$Time3), x = 'Time step', y = 'Biomass (tons)'))
                plot + geom_vline(xintercept = input$Time3, linetype = "dashed",  color = 'royalblue')
            })
            output$plot6 <- shiny::renderPlot({
                shiny::validate(
                    shiny::need(age.biomass != '',  'To Display this plot, please provide the non-standard atlantis output \'Biomass by age\' (default : age.biomass = NULL). Make sure to put the flag \'flag_age_output\' in the run.prm file if you want to look at this output')
                )
                color.pp <- as.character(g.col$col[which(g.col$grp %in% levels(age.bio.gr()$variable))])
                ggplot2::ggplot(data = age.bio.gr(), aes(x = .data$Cohort, y = .data$eff.pred, fill = .data$variable, width = .75)) + geom_bar(stat = "identity", position = 'stack') + scale_fill_manual(values = color.pp)  +
                    labs(list(title = paste('Predator  -', input$FG3, 'on Time step :', input$Time3), x = 'AgeGroup', y = 'Effective predation (tons)')) + theme_bw()
            })


        }

    )
}
