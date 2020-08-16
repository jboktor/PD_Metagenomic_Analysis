#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)

rootdir = substr(getwd(), 1, nchar(getwd()) -7)

source(paste0(rootdir, "src/miscellaneous_funcs.R"))
source(paste0(rootdir, "src/load_packages.R"))
# source(paste0(rootdir, "src/load_phyloseq_obj.R"))



ui <- fluidPage(
    titlePanel("Explore Dataset"),
    h4("Instructions: .."),
    sidebarLayout(
        sidebarPanel(
            selectInput("phyloseqInput", "Dataset",
                        choices = names(Phylo_Objects)
            ),
            sliderInput("binsInput","Number of bins:",
                        min = 1, max = 200,
                        value = 30),
            radioButtons("normalizationInput", "Normalization",
                         choices = c("TSS", "CLR", "None"),
                         selected = "None"),
            radioButtons("transformationInput", "Transformation",
                         choices = c("ArcSinSqrt", "None"),
                         selected = "None"),
            textInput("featureInput","Input feature to plot:")
        ),
        mainPanel(h4("Data Distribution"),
                  plotOutput("distribution_plot"),
                  br(), br(),
                  DT::dataTableOutput("table_results"),
                  br(),
                  h4("Feature of interest"),
                  plotOutput("boxplot"),
                  br(), br(), br(), br(), br(), br(), br(), br() )
    )
) 


# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distribution_plot <- renderPlot({
        # generate bins based on input$bins from ui.R
        
        print(Phylo_Objects[[input$phyloseqInput]])
        datObj <- Phylo_Objects[[input$phyloseqInput]]
        
        if (input$normalizationInput == "TSS") {
            datObj <- datObj %>% microbiome::transform("compositional")
        } else if (input$normalizationInput == "CLR") {
            datObj <- datObj %>% microbiome::transform("clr")
        }
        
        dataset <- datObj %>%
            microbiome::abundances()
        
        if (input$transformationInput == "ArcSinSqrt") {
            dataset <- asin(sqrt(dataset))
        }

        dataset %>% distribution_sanity2(binN = input$binsInput)


    })
    output$boxplot <- renderPlot({
        
        # generate bins based on input$bins from ui.R
        datObj <- Phylo_Objects[[input$phyloseqInput]]
        
        # This method doesn't include transfm. yet
        plot_feature(datObj, input$featureInput)
    })
    
    output$table_results <- DT::renderDataTable({
        Phylo_Objects[[input$phyloseqInput]] %>% 
            microbiome::abundances()
    })
    
} 

# Run the application 
shinyApp(ui = ui, server = server)
