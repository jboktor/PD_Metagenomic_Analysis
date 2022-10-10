# Joe Boktor 2020
# Caltech - Mazmanian Lab
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

source("load_packages_EDA.R")
source("functions_EDA.R")
cat("Current Working Directory: ", getwd(), "\n")
# conflict_prefer("box", "shinydashboard")

strat <- c("Pathways", "Enzymes", "KOs", "GOs", "Pfams", "eggNOGs")
non.strat <- paste0(strat, ".slim")
functionalObj <- c(non.strat, strat)

# _______________________________________________________________________________
#                               APP/UI
# _______________________________________________________________________________

ui <- dashboardPage(
  skin = "yellow",
  dashboardHeader(title = "PD Microbiome"),
  dashboardSidebar(
    selectInput(
      inputId = "phyloseqInput",
      "Select dataset",
      choices = names(Phylo_Objects)
    ),
    selectInput(
      "normalizationInput",
      "Normalization",
      choices = c("TSS", "CLR", "None"),
      selected = "None"
    ),
    selectInput(
      "transformationInput",
      "Transformation",
      choices = c("ArcSinSqrt", "Log10(x)", "None"),
      selected = "None"
    ),
    selectInput(
      "fileType",
      label = "Select file type",
      choices = list("png", "pdf")
    ),
    sidebarMenu(
      menuItem(
        "Dashboard",
        tabName = "dashboard",
        icon = icon("dashboard")
      ),
      menuItem("Funding/Support",
        icon = icon("th"),
        tabName = "creds"
      )
    ),
    fluidRow(
      column(
        width = 1,
        box(imageOutput("caltechLogo"),
          height = 30
        ),
        box(imageOutput("SKM"),
          height = 30
        )
      )
    )
  ),
  dashboardBody(tabItems(
    tabItem(
      tabName = "dashboard",
      h2("Parkinson's Disease Microbiome WGS Feature Exploration"),
      actionButton(inputId = "instruct", label = "Instructions"),
      br(),
      br(),
      selectizeInput(
        inputId = "featureSelection",
        label = "Select feature",
        choices = NULL,
        width = 800,
        list(maxOptions = 5)
      ),
      fluidRow(box(plotOutput("boxplot"))),
      fluidRow(
        downloadButton("downloadBoxPlot", "Download Boxplot"),
        # Conditional stratification download button
        uiOutput("downloadStratPlot.button")
      )
    ),
    tabItem(
      tabName = "creds",
      # fluidRow(imageOutput("TBC"), imageOutput("MJFF"), align="center"),
      # fluidRow(imageOutput("MJFF"), align="center")
    )
  ))
)


# _______________________________________________________________________________
#                               SERVER
# _______________________________________________________________________________

# Define server logic required to draw a histogram
server <- function(input, output) {
  # session_store <- reactiveValues()

  observeEvent(input$instruct, {
    showModal(
      modalDialog(
        title = "Instructions",
        "This application allows users to explore metagenomic profiles of",
        "of the PD microbiome. This includes features at the  taxonomic, metabolic, and gene-level.",
        "You may search for features of interest to generate boxplots. Additionally,",
        "for metabolic and gene level profiles we provide the option to view",
        "the assigned taxonomic stratification for a given feature as a stacked barplot.",
        "When viewing a feature of this class, a button titled 'Download Stratified'",
        "will appear in the bottom left-hand corner",
        br(),
        br(),
        "1) Select a dataset of interest from the sidebar menu on the left",
        br(),
        "2) Enter the desired normalization/transformation methods:
                (see suggested combination below)",
        br(),
        "3) Visualize the feature of interest select by the following:",
        br(),
        "A) Use the drop down menu in the Dashboard tab",
        br(),
        "B) Click on a feature in the same menu, hit return/backspace,
               and type in the name of your feature of interest",
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
  # _______________________________________________________________________________
  #                             Reactive variables
  # _______________________________________________________________________________

  df_upload <- reactive({
    datObj <- Phylo_Objects[[input$phyloseqInput]]
    raw.input <- datObj %>%
      subset_samples(donor_id %ni% low_qc[[1]]) %>%
      abundances() %>%
      as.data.frame() %>%
      rownames_to_column(var = "tempvars") %>%
      decode_rfriendly_rows(passed_column = "tempvars") %>%
      dplyr::select(-tempvars) %>%
      dplyr::rename(features = fullnames)
    return(raw.input)
  })

  df_products_upload <- reactive({
    datObj <- Phylo_Objects[[input$phyloseqInput]] %>%
      subset_samples(donor_id %ni% low_qc[[1]])

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
        microbiome::transform("shift", shift = impute) %>%
        microbiome::transform("log10") %>%
        microbiome::abundances()
    } else {
      df1 <- datObj %>%
        microbiome::abundances()
    }
    df_final <-
      df1 %>%
      as.data.frame() %>%
      rownames_to_column("tempvars") %>%
      decode_rfriendly_rows(passed_column = "tempvars") %>%
      dplyr::select(-tempvars) %>%
      dplyr::rename(features = fullnames)
    return(df_final)
  })

  boxplot_prep <- reactive({
    if (!is.null(input$featureSelection)) {
      cat(paste0("Feature selected: ", input$featureSelection, "\n"))
      abund <- df_products_upload()
      # Select Feature from abundance df
      d <- abund %>%
        dplyr::filter(features == input$featureSelection) %>%
        melt()
      dm <- group_col_from_ids(d, id = d$variable)
      dm$group <-
        factor(dm$group, levels = c("PC", "PD", "HC"))
      return(dm)
    }
  })

  barplot_prep <- reactive({
    if (!is.null(input$featureSelection)) {
      if (input$phyloseqInput %in% strat) {
        if (!grepl("\\|", input$featureSelection)) {
          datObj <- Phylo_Objects[[input$phyloseqInput]]
          df.barplot <- datObj %>%
            microbiome::transform("compositional") %>%
            microbiome::abundances() %>%
            as.data.frame() %>%
            rownames_to_column("tempvars") %>%
            decode_rfriendly_rows(passed_column = "tempvars") %>%
            dplyr::select(-tempvars) %>%
            dplyr::rename(features = fullnames) %>%
            dplyr::filter(grepl(input$featureSelection, features, fixed = TRUE)) %>%
            dplyr::filter(grepl("\\|", features)) %>%
            column_to_rownames(var = "features") %>%
            t() %>%
            melt()
          if (nrow(df.barplot) > 0) {
            df.barplot <- df.barplot %>%
              group_col_from_ids(ids = df.barplot$Var1) %>%
              dplyr::mutate(group = factor(group, levels = c("PC", "PD", "HC")))
            return(df.barplot)
          } else {
            return(NULL)
          }
        }
      } else if (input$phyloseqInput %in% non.strat) {
        restrat <- sub(".slim", "", input$phyloseqInput)
        datObj <- Phylo_Objects[[restrat]]
        if (!grepl("\\|", input$featureSelection)) {
          df.barplot <- datObj %>%
            microbiome::transform("compositional") %>%
            microbiome::abundances() %>%
            as.data.frame() %>%
            rownames_to_column("tempvars") %>%
            decode_rfriendly_rows(passed_column = "tempvars") %>%
            dplyr::select(-tempvars) %>%
            dplyr::rename(features = fullnames) %>%
            dplyr::filter(grepl(input$featureSelection, features, fixed = TRUE)) %>%
            dplyr::filter(grepl("\\|", features)) %>%
            column_to_rownames(var = "features") %>%
            t() %>%
            melt()
          if (nrow(df.barplot) > 0) {
            df.barplot <- df.barplot %>%
              group_col_from_ids(ids = df.barplot$Var1) %>%
              dplyr::mutate(group = factor(group, levels = c("PC", "PD", "HC")))
            return(df.barplot)
          } else {
            return(NULL)
          }
        } else {
          return(NULL)
        }
      }
    }
  })

  # _______________________________________________________________________________
  #                             Outputs
  # _______________________________________________________________________________

  observeEvent(df_upload(), {
    updateSelectizeInput(
      inputId = "featureSelection",
      choices = df_upload()$features,
      server = TRUE
    )
  })

  output$boxplot <- renderPlot({
    dm <- boxplot_prep()
    plot <- boxplot_all(
      dm,
      x = dm$group,
      y = dm$value,
      title = "",
      ylabel = paste(unique(dm$features), "Abundance")
    ) +
      stat_compare_means(method = "kruskal.test")
    print(plot)
    return(plot)
  })


  # Download Function:
  # 2 arguments as functions {filename, content}
  output$downloadBoxPlot <- downloadHandler(
    filename = function() {
      dm <- boxplot_prep()
      paste0(
        Sys.Date(),
        "_",
        unique(dm$features),
        "_Abundance.",
        input$fileType
      )
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if (input$fileType == "png") {
        png(file,
          width = 1024,
          height = 1228,
          res = 300
        )
      } else {
        pdf(file)
      }
      dm <- boxplot_prep()
      print(boxplot_all(
        dm,
        x = dm$group,
        y = dm$value,
        title = "",
        ylabel = paste(unique(dm$features), "abundance")
      ))
      dev.off()
    }
  )

  output$downloadStratPlot <- downloadHandler(
    filename = function() {
      dm <- boxplot_prep()
      paste0(
        Sys.Date(),
        "_",
        unique(dm$features),
        "_stratified_abundance.html"
      )
    },
    content = function(file) {
      stratplot.df <- barplot_prep() %>%
        mutate(Var2 = sub(".*\\|", "", Var2)) %>%
        group_by(Var2) %>%
        mutate(mean_size = mean(value, na.rm = TRUE)) %>%
        ungroup() %>%
        dplyr::mutate(Var2 = fct_reorder(Var2, mean_size))
      stratplot <-
        stratplot.df %>%
        ggplot(aes(
          x = reorder(Var1, -value),
          y = value,
          fill = Var2
        )) +
        geom_bar(stat = "identity") +
        theme_bw() +
        labs(x = "Donor", y = "Abundance", fill = NULL) +
        facet_grid(
          cols = vars(group), space = "free_x", scales = "free_x",
          labeller = labeller(group = c(
            "PC" = "Population Controls",
            "PD" = "Parkinson's Disease",
            "HC" = "Household Controls"
          ))
        ) +
        scale_fill_manual(
          values = color_loop_generator(levels(stratplot.df$Var2)),
          guide = guide_legend(reverse = TRUE)
        ) +
        theme(
          axis.text.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          plot.margin = unit(c(1, 1, 1, 3), "cm")
        )
      stratPlot <-
        ggplotly(stratplot, tooltip = c("y", "x", "fill"))
      saveWidget(as_widget(stratPlot), file, selfcontained = TRUE)
    }
  )

  output$downloadStratPlot.button <-
    renderUI(expr = if (input$phyloseqInput %in% functionalObj &
      !grepl("\\|", input$featureSelection)) {
      downloadButton("downloadStratPlot", "Download Stratified Barplot")
    } else {
      NULL
    })

  output$strat_plot <- renderPlotly({
    stratplot.df <- barplot_prep()

    if (!is.null(stratplot.df)) {
      stratplot.final.df <- stratplot.df %>%
        group_by(Var2) %>%
        mutate(mean_size = mean(value, na.rm = TRUE)) %>%
        ungroup() %>%
        dplyr::mutate(Var2 = fct_reorder(Var2, mean_size))
      stratplot <-
        stratplot.final %>%
        ggplot(aes(
          x = reorder(Var1, -value),
          y = value,
          fill = Var2
        )) +
        geom_bar(stat = "identity") +
        theme_bw() +
        labs(x = "Donor", y = "Abundance") +
        facet_grid(
          cols = vars(group),
          space = "free_x",
          scales = "free_x"
        ) +
        scale_fill_manual(
          values = color_loop_generator(levels(stratplot.final.df$Var2)),
          guide = guide_legend(reverse = TRUE)
        ) +
        theme(
          legend.position = "none",
          axis.text.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin = unit(c(1, 1, 1, 3), "cm")
        )
      stratPlot <-
        ggplotly(stratplot, tooltip = c("y", "x", "fill"))
    }
  })

  # LOGOS
  output$caltechLogo <- renderImage(
    {
      return(
        list(
          src = "Caltech_LOGO-Orange_RGB.png",
          contentType = "image/png",
          alt = "caltechLogo",
          height = 30
        )
      )
    },
    deleteFile = FALSE
  )

  output$SKM <- renderImage(
    {
      return(list(
        src = "SarkisMazmanianLabLogo.png",
        contentType = "image/png",
        alt = "Error",
        height = 30
      ))
    },
    deleteFile = FALSE
  )

  output$TBC <- renderImage(
    {
      return(list(
        src = "TBC.png",
        contentType = "image/png",
        alt = "Error",
        width = 400
      ))
    },
    deleteFile = FALSE
  )

  output$MJFF <- renderImage(
    {
      return(list(
        src = "MJFF.png",
        contentType = "image/png",
        alt = "Error",
        width = 600
      ))
    },
    deleteFile = FALSE
  )
}

# Run the application
shinyApp(ui = ui, server = server)
