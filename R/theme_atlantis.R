##' @title Main Atlatnis Theme
##' @return The ggplot theme
##' @author Javier Porobic
##' @export
theme_atlantis <- function(){
    ggplot2::theme(
                 plot.title = ggplot2::element_text(colour = "#844D14", face = 'bold', family = "Times New Roman", size = ggplot2::rel(2)),
                 ## add border 1)
                 panel.border = ggplot2::element_rect(colour = "#562F0E", fill = NA, linetype = 2),
                 ## color background 2)
                 panel.background = ggplot2::element_rect(fill = scales::alpha('#fcda96', 0.1)),
                 ## modify grid 3)
                 panel.grid.major.x = ggplot2::element_line(colour = "#844D14", linetype = 3, size = 0.5),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_line(colour = "#844D14", linetype = 3, size = 0.5),
                 panel.grid.minor.y = ggplot2::element_blank(),
                 ## modify text, axis and colour 4) and 5)
                 axis.text = ggplot2::element_text(colour = "#292B15", face = "italic", family = "Times New Roman", size = ggplot2::rel(1.2)),
                 axis.title = ggplot2::element_text(colour = "#292B15", face = 'bold', family = "Times New Roman", size = ggplot2::rel(1.5)),
                 axis.ticks = ggplot2::element_line(colour = "#292B15", size = ggplot2::rel(1.2)),
                 ## legend at the bottom 6)
                 legend.text = ggplot2::element_text(colour = "#7A6F42", face = 'bold', family = "Times New Roman", size = ggplot2::rel(1.5)),
                 legend.title = ggplot2::element_text(colour = "#292B15", face = 'bold', family = "Times New Roman", size = ggplot2::rel(1.5)),
                 ## Faceting
                 strip.background = ggplot2::element_rect(fill = scales::alpha("#7A6F42", 0.7), linetype = 2),
                 strip.background.x = ggplot2::element_rect(fill = scales::alpha("#7A6F42", 0.7), linetype = 2),
                 strip.background.y = ggplot2::element_rect(fill = scales::alpha("#7A6F42", 0.7), linetype = 2),
                 strip.text = ggplot2::element_text(colour = "#241309", face = 'bold', family = "Times New Roman", size = ggplot2::rel(1.5))
             )
}
