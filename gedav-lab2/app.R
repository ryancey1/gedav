# load necessary libraries
library(shiny)
library(shinythemes)
library(htmltools)

# load in the data with headers and row names specified
spellman <- read.table(file = "spellman.txt",
                       header = TRUE,
                       row.names = 1)

# subset the cdc15 experiment
expt <- grepl("^cdc15", colnames(spellman))
cdc15 <- spellman[, expt]

# Define UI logic
ui <- fluidPage(theme = "lumen",
                titlePanel(
                    h2(
                        tagList("Interactive correlation plot of CDC15 yeast microarray data"),
                        align = "center"
                    ), windowTitle = "CDC15 correlation"
                ),
                sidebarLayout(
                    sidebarPanel(
                        # choice 1
                        selectInput(
                            inputId = 'xcol',
                            label = "X Variable",
                            choices = colnames(cdc15),
                            selected = colnames(cdc15)[1]
                        ),
                        # choice 2
                        selectInput(
                            inputId = 'ycol',
                            label = "Y Variable",
                            choices = colnames(cdc15),
                            selected = colnames(cdc15)[2]
                        ),
                        # point type choice
                        selectInput(
                            inputId = "pch",
                            label = "Point",
                            choices = c(
                                "Circle" = 21,
                                "Square" = 22,
                                "Diamond" = 23,
                                "Triangle" = 24
                            )
                        ),
                        # point color choice
                        selectInput(
                            inputId = "col",
                            label = "Color",
                            choices = c(
                                "Red" = "red",
                                "Orange" = "orange",
                                "Yellow" = "yellow",
                                "Green" = "green",
                                "Blue" = "blue",
                                "Purple" = "purple"
                            )
                        )
                    ),
                    # main panel correlation plot
                    mainPanel(plotOutput(outputId = 'xyplot'),
                              tagList("Source: ", tags$a(href = "https://www.molbiolcell.org/doi/pdf/10.1091/mbc.9.12.3273", "Spellman et al 1998")))
                ))

# Define server logic
server <- function(input, output) {
    # selected columns
    selected_data <- reactive({
        cdc15[, c(input$xcol, input$ycol)]
    })
    
    # selected color
    selected_color <- reactive({
        rgb(t(col2rgb(input$col)),
            alpha = (255 / 2),
            maxColorValue = 255)
    })
    
    # selected point
    selected_pch <- reactive({
        as.numeric(input$pch)
    })
    
    # plot
    output$xyplot <- renderPlot({
        # xy plot
        plot(
            selected_data(),
            type = "p",
            col = "black",
            bg = selected_color(),
            pch = selected_pch(),
            lwd = 2,
            cex = 1.5
        )
        # title
        title(paste(
            "Correlation plot:",
            colnames(selected_data())[2],
            "vs",
            colnames(selected_data())[1]
        ))
        # linear model
        abline(
            lm(selected_data()[, 2] ~ selected_data()[, 1]),
            lwd = 3,
            lty = 3,
            col = "gray25"
        )
        # legend
        legend(
            "bottomright",
            legend = "linear model",
            lwd = 3,
            lty = 3,
            xjust = 0.5,
            col = "gray25",
            inset = 0.02
        )
    })
    
}

# compile and execute
shinyApp(ui = ui, server = server)