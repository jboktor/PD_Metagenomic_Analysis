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



ui <- fluidPage(
    titlePanel("Explore Dataset"),
    h4("Instructions: "),
    h6("1) Select a dataset of interest from drop-down menu below.", 
       br(), br(),
       "2) When table loads, click on any feature to generate a plot in the 'Feature of Interest' Tab. ", 
       br(), br(),
       "Suggested Normalization: TSS ", br(),
       "Suggested Transformation: ArcsinSqrt"),
    sidebarLayout(
        sidebarPanel(
            selectInput("phyloseqInput", "Dataset",
                        choices = names(Phylo_Objects)
            ),
            # sliderInput("binsInput","Number of bins:",
            #             min = 1, max = 200,
            #             value = 30),
            radioButtons("normalizationInput", "Normalization",
                         choices = c("TSS", "CLR", "None"),
                         selected = "None"),
            radioButtons("transformationInput", "Transformation",
                         choices = c("ArcSinSqrt", "Log10(x)", "None"),
                         selected = "None")
            ),
        mainPanel(
            tabsetPanel(
                tabPanel("Data Distribution", plotOutput("distribution_plot")), 
                tabPanel("Table Results",  DT::dataTableOutput("table_results")), 
                tabPanel("Feature of Interest", plotOutput("boxplot"))
            )
        )
    )
) 


# Define server logic required to draw a histogram
server <- function(input, output) {
    
    df_products_upload <- reactive({
        
        print(Phylo_Objects[[input$phyloseqInput]])
        datObj <- Phylo_Objects[[input$phyloseqInput]]
        
        # Save pre-normalized minimum value if needed 
        raw.input <- datObj %>% 
            abundances()
        impute <- min(raw.input[raw.input > 0])/2
        
        #-------------- Normalization -------------- 
        if (input$normalizationInput == "TSS") {
            datObj <- datObj %>% microbiome::transform("compositional")
        } else if (input$normalizationInput == "CLR") {
            datObj <- datObj %>% microbiome::transform("clr")
        }
        
        #-------------- Transformation -------------- 
        if (input$transformationInput == "ArcSinSqrt") {
            asin.input <- datObj %>%
                microbiome::abundances()
            df1 <- asin(sqrt(asin.input))
        } else if (input$transformationInput == "Log10(x)") {
            df1 <- datObj %>% 
                microbiome::transform('shift', shift = impute) %>%
                microbiome::transform("log10") %>%
                microbiome::abundances()
        } else {
            df1 <- datObj %>%
                microbiome::abundances()
        }
        
        return(df1)
    })
        
    
    output$distribution_plot <- renderPlot({
        
        df_products_upload() %>% 
            distribution_sanity2()

    })
    
    output$table_results <- 
        DT::renderDataTable(server = FALSE, selection = 'single', {
            df_products_upload()
    })
    
    output$boxplot <- renderPlot({
        
        s = input$table_results_rows_selected
        
        if (length(s)) {

            abund <- df_products_upload() %>% 
                as.data.frame() %>% 
                rownames_to_column() 
            
            # Select Feature from abundance df
            feature <- abund$rowname[s]
            d <- abund %>% 
                filter(rowname == feature) %>% 
                melt()
            
            dm <-  group_col_from_ids(d, id=d$variable)
            dm$group <- factor(dm$group, levels = c("PC", "PD", "HC"))
            boxplot_all(dm, x=dm$group, y=dm$value, 
                        title=" ", 
                        ylabel= paste(unique(dm$rowname), "Abundance"))
        }
        
    })
    
    
} 

# Run the application 
shinyApp(ui = ui, server = server)
