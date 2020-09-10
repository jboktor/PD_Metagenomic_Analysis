#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#

library(shiny)
library(rsconnect)
library(DT)
library(BiocManager)
options(repos = BiocManager::repositories())

source("miscellaneous_funcs_EDA.R")
source("load_packages_EDA.R")

load("PhyloseqObj.RData")
Phylo_Objects

# rootdir = substr(getwd(), 1, nchar(getwd()) -7)
# source(paste0(rootdir, "src/miscellaneous_funcs.R"))
# source(paste0(rootdir, "src/load_packages.R"))
# source(paste0(rootdir, "src/load_phyloseq_obj.R"))



ui <- fluidPage(
    titlePanel("Explore Dataset"),
    h4("Instructions: "),
    h6("1) Select a dataset of interest from drop-down menu below.", 
       br(), br(),
       "2) When table loads, you may search for particular features and copy and paste them into the textbox below. ", 
       br(), br(),
       "Suggested Normalization: TSS ", br(),
       "Suggested Transformation: ArcsinSqrt"),
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
            textInput("featureInput","Manually Input feature to plot:")
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("Data Distribution", plotOutput("distribution_plot")), 
                tabPanel("Table Results",  DT::dataTableOutput("table_results")), 
                tabPanel("Feature of Interest", plotOutput("boxplot"))
            )
            
    #         h4("Data Distribution"),
    #               plotOutput("distribution_plot"),
    #               br(), br(),
    #               DT::dataTableOutput("table_results"),
    #               hr(),
    #               h4("Feature of interest"),
    #               plotOutput("boxplot"),
    #               br(), br(), br(), br(), br(), br(), br(), br()
        )
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
    
    output$table_results <- 
        DT::renderDataTable(server = FALSE, selection = 'single', {
    
        
        if (input$normalizationInput == "TSS") {
            Phylo_Objects[[input$phyloseqInput]] %>% 
                microbiome::transform("compositional") %>% 
                microbiome::abundances()
            
        } else if (input$normalizationInput == "CLR") {
            Phylo_Objects[[input$phyloseqInput]] %>% 
                microbiome::transform("clr") %>% 
                microbiome::abundances()
        } else {
            Phylo_Objects[[input$phyloseqInput]] %>% 
                microbiome::abundances()
        }
        
    })
    
    output$boxplot <- renderPlot({
        
        s = input$table_results_rows_selected
        
        if (length(s)) {
            
            
            
            # datObj <- Phylo_Objects[[input$phyloseqInput]]
            
            
            # plot_feature(datObj, input$featureInput)
            
            if (input$normalizationInput == "TSS") {
                abund <- Phylo_Objects[[input$phyloseqInput]] %>% 
                    microbiome::transform("compositional") %>% 
                    microbiome::abundances() %>% 
                    as.data.frame() %>% 
                    rownames_to_column() 
                
            } else if (input$normalizationInput == "CLR") {
                abund <- Phylo_Objects[[input$phyloseqInput]] %>% 
                    microbiome::transform("clr") %>% 
                    microbiome::abundances() %>% 
                    as.data.frame() %>% 
                    rownames_to_column() 
                
            } else {
                abund <- Phylo_Objects[[input$phyloseqInput]] %>% 
                    microbiome::abundances() %>% 
                    as.data.frame() %>% 
                    rownames_to_column() 
            }
            
            
            # Select Feature from abundance df
            feature <- abund$rowname[s]
            d <- abund %>% 
                filter(rowname == feature) %>% 
                melt()
            
            dm <-  group_col_from_ids(d, id=d$variable)
            dm$value <- asin(sqrt(dm$value))
            dm$group <- factor(dm$group, levels = c("PC", "PD", "HC"))
            boxplot_all(dm, x=dm$group, y=dm$value, 
                        cols=c("PC"= "#bfbfbf", 
                               "PD" = "#ed7d31", 
                               "HC" = "#5b9bd5"), 
                        title=" ", 
                        ylabel= paste(unique(dm$rowname), "Abundance"))
        }
        
    })
    
    
} 

# Run the application 
shinyApp(ui = ui, server = server)
