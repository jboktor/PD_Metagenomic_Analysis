#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

library(shiny)
library(rsconnect)
library(DT)
library(BiocManager)
library(plotly)
options(repos = BiocManager::repositories())
source("miscellaneous_funcs_EDA.R")
source("load_packages_EDA.R")
load("PhyloseqObj.RData")

ui <- fluidPage(
    titlePanel("PD Metagenomics Feature Exploration"),
    actionButton("instruct", "Instructions"),
    br(),
    br(),
    sidebarLayout(
        sidebarPanel(
            selectInput("phyloseqInput", "Select dataset",
                        choices = names(Phylo_Objects)),
            radioButtons(
                "normalizationInput",
                "Normalization",
                choices = c("TSS", "CLR", "None"),
                selected = "None"
            ),
            radioButtons(
                "transformationInput",
                "Transformation",
                choices = c("ArcSinSqrt", "Log10(x)", "None"),
                selected = "None"
            ),
            radioButtons(
                "fileType",
                label = "Select file type",
                choices = list("png", "pdf")
            )
        ),
        mainPanel(tabsetPanel(
            tabPanel("Data Distribution", plotOutput("distribution_plot")),
            tabPanel("Table Results",  DT::dataTableOutput("table_results")),
            tabPanel(
                "Feature of Interest",
                plotOutput("boxplot"),
                downloadButton('downloadPlot', 'Download Figure'),
                plotlyOutput("strat_plot")
            )
        ))
    )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
    observeEvent(input$instruct, {
        showModal(
            modalDialog(
                title = "Instructions",
                "1) Select a dataset of interest from the drop-down menu",
                br(),
                "2) Click on any feature to generate a plot under the 'Feature of Interest' tab. ",
                br(),
                br(),
                "Suggested Normalization: TSS ",
                br(),
                "Suggested Transformation: ArcsinSqrt",
                br(),
                br(),
                "(TSS) = total sum scaling (compositional transformation)",
                br(),
                "(CLR) = center-log Ratio",
                br(),
                br(),
                "Valid normalization/transformation combinations include the following: ",
                br(),
                "TSS + Log10x",
                br(),
                "TSS + ArcsinSqrt",
                br(),
                "TSS + None",
                br(),
                "None + Log10x",
                br(),
                "CLR + None",
                br(),
                "CLR + Log10x",
                br(),
                br(),
                
                easyClose = TRUE
            )
        )
    })
    #-------------------------------------------------------------------------------
    #                             Reactive variables
    #-------------------------------------------------------------------------------
    
    df_products_upload <- reactive({
        print(Phylo_Objects[[input$phyloseqInput]])
        datObj <- Phylo_Objects[[input$phyloseqInput]]
        
        # Save pre-normalized minimum value if needed
        raw.input <- datObj %>%
            abundances()
        impute <- min(raw.input[raw.input > 0]) / 2
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
            dm <-  group_col_from_ids(d, id = d$variable)
            dm$group <-
                factor(dm$group, levels = c("PC", "PD", "HC"))
            
            return(dm)
        }
        
    })
    barplot_prep <- reactive({
        s = input$table_results_rows_selected
        strat <- c("Pathways", "Enzymes", "KOs.all", "KOs")
        non.strat <- paste0(strat, ".slim")
        if (length(s)) {
            abund <- df_products_upload() %>%
                as.data.frame() %>%
                rownames_to_column()
            feature <- abund$rowname[s]
            if (input$phyloseqInput %in% strat) {
                if (!grepl("\\|", feature)) {
                    df.barplot <- Phylo_Objects[[input$phyloseqInput]] %>%
                        microbiome::transform("compositional") %>%
                        microbiome::abundances() %>%
                        as.data.frame() %>%
                        rownames_to_column() %>%
                        filter(grepl(feature, rowname, fixed = TRUE)) %>%
                        filter(grepl("\\|", rowname)) %>%
                        column_to_rownames() %>%
                        t() %>%
                        melt()
                    if (nrow(df.barplot) > 0) {
                        df.barplot <- df.barplot %>%
                            group_col_from_ids(ids = df.barplot$Var1) %>%
                            mutate(group = factor(group, levels = c("PC", "PD", "HC")))
                        return(df.barplot)
                    } else {
                        return(NULL)
                    }
                }
            } else if (input$phyloseqInput %in% non.strat) {
                
                restrat <- sub(".slim", "", input$phyloseqInput)
                datObj <- Phylo_Objects[[restrat]]
                # Troubleshooting 
                
                
                if (!grepl("\\|", feature)) {
                    df.barplot <- datObj %>%
                        microbiome::transform("compositional") %>%
                        microbiome::abundances() %>%
                        as.data.frame() %>%
                        rownames_to_column() %>%
                        filter(grepl(feature, rowname, fixed = TRUE)) %>%
                        filter(grepl("\\|", rowname)) %>%
                        column_to_rownames() %>%
                        t() %>%
                        melt()
                    if (nrow(df.barplot) > 0) {
                        df.barplot <- df.barplot %>%
                            group_col_from_ids(ids = df.barplot$Var1) %>%
                            mutate(group = factor(group, levels = c("PC", "PD", "HC")))
                        return(df.barplot)
                    } else {
                        return(NULL)
                    }
                } else {
                    return(NULL)
                }
        } 
            }})
    
    #-------------------------------------------------------------------------------
    #                             Outputs
    #-------------------------------------------------------------------------------
    
    output$distribution_plot <- renderPlot({
        df_products_upload() %>%
            distribution_sanity2()
    })
    
    output$table_results <-
        DT::renderDataTable(server = FALSE, selection = 'single', {
            df_products_upload()
        })
    
    output$boxplot <- renderPlot({
        dm <- boxplot_prep()
        boxplot_all(
            dm,
            x = dm$group,
            y = dm$value,
            title = " ",
            ylabel = paste(unique(dm$rowname), "Abundance")
        )
    })
    
    
    # downloadHandler contains 2 arguments as functions, namely filename, content
    output$downloadPlot <- downloadHandler(
        filename =  function() {
            dm <- boxplot_prep()
            paste0(unique(dm$rowname), "_Abundance.", input$fileType)
        },
        # content is a function with argument file. content writes the plot to the device
        content = function(file) {
            if (input$fileType == "png")
                png(file,
                    width = 1024,
                    height = 1228,
                    res = 300) # open the png device
            else
                pdf(file) # open the pdf device
            dm <- boxplot_prep()
            print(boxplot_all(
                dm,
                x = dm$group,
                y = dm$value,
                title = " ",
                ylabel = paste(unique(dm$rowname), "Abundance")
            ))
            dev.off()  # turn the device off
            
        }
    )
    output$strat_plot <- renderPlotly({
        stratplot.df <- barplot_prep()
        
        if (!is.null(stratplot.df)) {
            stratplot <-
                stratplot.df %>%
                ggplot(aes(
                    x = reorder(Var1,-value),
                    y = value,
                    fill = Var2
                )) +
                geom_bar(stat = "identity") +
                theme_bw() +
                labs(x = "Donor", y = "Abundance") +
                facet_wrap( ~ group, scales = "free_x") +
                theme(
                    legend.position = "none",
                    axis.text.x = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major.x = element_blank(),
                )
            # stratplot
            ggplotly(stratplot, tooltip = "all")
        }
        
    })
    
}

# Run the application
shinyApp(ui = ui, server = server)
