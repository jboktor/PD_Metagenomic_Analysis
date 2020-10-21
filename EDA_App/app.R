#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

library(shiny)
library(rsconnect)
library(DT)
library(BiocManager)
options(repos = BiocManager::repositories())
source("miscellaneous_funcs_EDA.R")
source("load_packages_EDA.R")

load("PhyloseqObj.RData")

ui <- fluidPage(
    titlePanel("PD Metagenomics Feature Exploration"),
    actionButton("instruct", "Instructions"),
    br(), br(),
    sidebarLayout(
        sidebarPanel(
            selectInput("phyloseqInput", "Select dataset",
                        choices = names(Phylo_Objects)
            ),
            radioButtons("normalizationInput", "Normalization",
                         choices = c("TSS", "CLR", "None"),
                         selected = "None"),
            radioButtons("transformationInput", "Transformation",
                         choices = c("ArcSinSqrt", "Log10(x)", "None"),
                         selected = "None"),
            radioButtons("fileType", label = "Select file type", choices = list("png", "pdf"))
            ),
        mainPanel(
            tabsetPanel(
                tabPanel("Data Distribution", plotOutput("distribution_plot")), 
                tabPanel("Table Results",  DT::dataTableOutput("table_results")), 
                tabPanel("Feature of Interest", plotOutput("boxplot"),
                         downloadButton('downloadPlot', 'Download Figure'))
            )
        )
    )
) 


# Define server logic required to draw a histogram
server <- function(input, output) {
    
    observeEvent(input$instruct, {
        showModal(modalDialog(
            title = "Instructions",
            "1) Select a dataset of interest from the drop-down menu", 
            br(),
            "2) Click on any feature to generate a plot under the 'Feature of Interest' tab. ",
            br(), br(),
            "Suggested Normalization: TSS ", 
            br(),
            "Suggested Transformation: ArcsinSqrt",
            br(), br(),
            "(TSS) = total sum scaling (compositional transformation)", br(),
            "(CLR) = center-log Ratio", br(), br(),
            "Valid normalization/transformation combinations include the following: ", br(),
            "TSS + Log10x", br(),
            "TSS + ArcsinSqrt", br(),
            "TSS + None", br(),
            "None + Log10x", br(),
            "CLR + None", br(),
            "CLR + Log10x", br(), br(),
            
            easyClose = TRUE
        ))
    })

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
    
    
    boxplot_prep <- reactive({
        
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
            
            return(dm)
        }
        
    })
    
    output$boxplot <- renderPlot({
        dm <- boxplot_prep()
        boxplot_all(dm, x=dm$group, y=dm$value, 
                    title=" ", ylabel= paste(unique(dm$rowname), "Abundance"))
    })
    

    # downloadHandler contains 2 arguments as functions, namely filename, content
    output$downloadPlot <- downloadHandler(
        filename =  function() {
            dm <- boxplot_prep()
            paste0(unique(dm$rowname), "_Abundance.", input$fileType)
        },
        # content is a function with argument file. content writes the plot to the device
        content = function(file) {
            if(input$fileType == "png")
                png(file, width = 1024, height = 1228, res =300) # open the png device
            else
                pdf(file) # open the pdf device
            dm <- boxplot_prep()
            print(boxplot_all(dm, x=dm$group, y=dm$value, 
                        title=" ", 
                        ylabel= paste(unique(dm$rowname), "Abundance")))
            dev.off()  # turn the device off
            
        } 
    )
    
    
} 

# Run the application 
shinyApp(ui = ui, server = server)
