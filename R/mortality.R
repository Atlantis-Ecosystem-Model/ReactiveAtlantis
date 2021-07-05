##' This function helps to explore the different types of mortalities in an Atlantis
##'     model. Technically, this function takes the results of mortality from
##'     Atlantis and compares them. This includes temporal components of natural, fishing and predation
##'     mortality. This function also explores the different mortalities in space and
##'     time.
##' @title Mortality analysis for Atlantis
##' @param grp.file Character string with the path to the groups \emph{*.csv} file (Atlantis input file).
##' @param prm.file Parameter file
##' @param SpeMort Character string with the path to the plain text Specific
##'     mortality file \emph{*SpecificMort.txt} (Atlantis output file).
##' @param PredMort Character string with the path to the plain text predation
##'     morality file \emph{*SpecificPredMort.txt} (Atlantis output file).
##' @return A reactive output with plots to analyse mortality in Atlantis
##' @author Javier Porobic
##' @import stats utils grDevices ggplot2 graphics shiny RColorBrewer
##' @export
mortality <- function(grp.file, prm.file, SpeMort, PredMort){
    SpeMort    <- read.csv(SpeMort, sep = ' ')
    PredMort   <- read.csv(PredMort, sep = ' ')
    grp        <- read.csv(grp.file)
    names(grp) <- tolower(names(grp))
    ## Start the Shiny application
    shiny::shinyApp(
        ## Create the different tabs
        ui <- shiny::navbarPage("Mortality",
                         shiny::tabPanel('Total',
                                  shiny::fluidRow(
                                      shiny::column(2,
                                             shiny::wellPanel(
                                                 shiny::selectInput('FG1', 'Functional Group :', as.character(grp$longname))
                                             )),
                                      shiny::column(10,
                                             shiny::plotOutput('plot2a', width = "100%", height = "600px")
                                             ))),
                         ## -- Exit --
                         shiny::tabPanel(
                             shiny::actionButton("exitButton", "Exit")
                         )
                         ),
        ## Link the input for the different tabs with your original data
        ## Create the plots
        function(input, output, session){
            mort.sp <- shiny::reactive({
                FG      <- grp$code[which(grp$longname %in% input$FG1)]
                mort.sp <- .calc.mort(SpeMort, PredMort, FG)
                return(mort.sp)
            })
            ## input plots
            mortInputPlot <- shiny::reactive({
                sps  <- grp$code[which(grp$longname %in% input$FG1)]
                mort <- reshape2::melt(mort.sp(), id.var = c('Time', 'Mtype'))
                getPalette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
                cols <- getPalette(length(levels(mort$variable)))
                Title <- paste0('Mortalities for ', as.character(input$FG1))
                p <-  ggplot2::ggplot(mort, ggplot2::aes(.data$Time, .data$value, fill = .data$variable))
                p <- p + ggplot2::geom_bar(stat="identity", position = ggplot2::position_dodge())
                p <- p + ggplot2::facet_wrap(~Mtype, ncol = 1, scales="free_y") + ggplot2::scale_fill_manual('Age', values = cols)
                p <- p + ggplot2::theme_linedraw() + ggplot2::labs(title = Title, x = 'Simulation Time', y = 'Mortality level')
                p
            })
            ## Plots
            shiny::observeEvent(input$exitButton, {
                shiny::stopApp()
            })
            output$plot2a <- shiny::renderPlot({
                print(mortInputPlot())
            })
        }
    )

}

##' @title Extranting mortality values from the txt files provided by Atlantis
##' @param SpeMort Especific mortality
##' @param PredMort Predation mortality
##' @param FG functional group
##' @return List of mortalities by functional group by time
##' @author Javier Porobic
.calc.mort <- function(SpeMort, PredMort, FG){
    F.mort <- paste(as.character(FG), c(0 : 9), 'S0', 'F', sep = '.')
    N.mort <- paste(as.character(FG), c(0 : 9), 'S0', 'M1', sep = '.')
    P.mort <- paste(as.character(FG), c(0 : 9), 'S0', 'M2', sep = '.')
    ## Calculating Fishing mortality
    Fmort <- SpeMort[, c(1, which(colnames(SpeMort) %in% F.mort))]
    colnames(Fmort) <- c('Time', paste0('Age0', 1 : (ncol(Fmort) - 1)))
    ## Calculating Natural mortality
    Nmort <- SpeMort[, c(1, which(colnames(SpeMort) %in% N.mort))]
    colnames(Nmort) <- c('Time', paste0('Age0', 1 : (ncol(Nmort) - 1)))
    ## Calculating Predation mortality
    Pmort <- SpeMort[, c(1, which(colnames(SpeMort) %in% P.mort))]
    colnames(Pmort) <- c('Time', paste0('Age0', 1 : (ncol(Pmort) - 1)))
    Fmort$Mtype <- 'Fishing Mortality'
    Nmort$Mtype <- 'Natural Mortality'
    Pmort$Mtype <- 'Predation Mortality'
    return(rbind(Fmort, Nmort, Pmort))
}
