library(shiny)
library(ggplot2)
library(dplyr)

load(url("https://share.genome.au.dk/DdpQDYTVMqg/te_ins_hier_te_mots.RData"))
te_hier <- te.hier.counts %>% mutate(TE_fam = sub("_", "-", Alias))
te_mots_kds_te <- merge(te_mots_kds_te, te_hier, by = "TE_fam")
te_mots_kds_te$KD <- factor(te_mots_kds_te$KD)
levels(te_mots_kds_te$KD) <- c("Treatment 1", "Treatment 2", "Treatment 3")
# Define UI for application
ui <- fluidPage(
  # Application title
  titlePanel("Motif hunting"),
  
  # Sidebar layout with a slider input for the threshold
  sidebarLayout(
    sidebarPanel(
      sliderInput("threshold",
                  "TPM threshold:",
                  min = 1,
                  max = 333,
                  value = c(1, 333),  # Range values
                  dragRange = TRUE),
      # Add another checkbox group for motif selection
      tabsetPanel(
        # First tab for TE selection
        tabPanel("TE Selection",
                 actionButton("select_all", "Select All"),
                 actionButton("deselect_all", "Deselect All"),
                 actionButton("select_LTR", "Select LTR"),
                 actionButton("select_LINE", "Select LINE"),
                 actionButton("select_DNA", "Select DNA"),
                 # Add more buttons for other classes
                 uiOutput("checkbox_ui_te")),
        # Second tab for Motif selection
        tabPanel("Motif Selection",
                 uiOutput("checkbox_ui_motif"))
      )
    ),
    
    # Main panel for displaying the plot
    mainPanel(
      plotOutput("distPlot", width = "100%", height = "700px")
    )
  )
)

# Define server logic 
server <- function(input, output, session) {
  # Select the data frame based on the input selection
  # Initialize reactiveValues for TE counts
  rv <- reactiveValues()
  
  
  # Reactive expression for the filtered data
  filtered_data <- reactive({
    # Filter the data based on the input threshold and TE_insert
    te_mots_kds_te %>% 
      filter(TSS_TPM > input$threshold[1] & TSS_TPM < input$threshold[2])
  })
  # Initialize reactive expression for motifs
  
  motifs <- reactive({
    unique(na.omit(filtered_data()$Motif))
  })
  
  observe({
    rv$counts <- filtered_data() %>%
      group_by(TE_insert) %>% summarize(n = 1, TE_fam = unique(TE_fam),
                                        TE_class = unique(class),
                                        counts = unique(count)) %>%
      ungroup() %>% group_by(TE_fam) %>% summarize(n = sum(n), 
                                                   counts = unique(counts),
                                                   TE_class = unique(TE_class)) %>% 
      mutate(TE_color = case_when(TE_class == "LTR" ~ "red",
                                  TE_class == "LINE" ~ "green",
                                  TE_class == "DNA" ~ "purple",
                                  TRUE ~ "black"))
  })
  
  # Generate checkboxes based on unique TE_insert values
  output$checkbox_ui_te <- renderUI({
    if(!is.null(rv$counts)) {
      # Get unique TE_fam values
      te_fams <- rv$counts$TE_fam
      te_classes <- rv$counts$TE_class  # assuming 'class' column exists in rv$counts
      
      # Create labels for each checkbox with counts
      labels <- lapply(seq_along(te_fams), function(i) {
        te_fam <- te_fams[i]
        te_class <- te_classes[i]
        
        # Create labels for each checkbox with counts
        paste(te_fam, " (", rv$counts$n[rv$counts$TE_fam == te_fam], "/",
              rv$counts$counts[rv$counts$TE_fam == te_fam], ")", sep = "")
      })
      names(labels) <- te_fams
      
      checkboxGroupInput("TE_checkboxes", 
                         "Select TEs:",
                         # choices = setNames(te_fams, labels), 
                         choiceNames = lapply(te_fams, function(tf){
                           span(labels[[tf]], style = paste0("color: ",
                                                             rv$counts$TE_color[rv$counts$TE_fam == tf]))
                         }),
                         choiceValues = as.list(te_fams),
                         selected = te_fams,
                         inline = FALSE)
    }
  })
  
  # Generate checkboxes based on unique Motif values
  output$checkbox_ui_motif <- renderUI({
    
    # Get unique Motif values
    
    checkboxGroupInput("motif_checkboxes", "Select Motifs:", 
                       choices = motifs(),
                       selected = c(
                         "TATA",
                         "DPE",
                         "dInr",
                         "hInr", 
                         "dTCT",
                         "hTCT",
                         "PB",
                         "GAGA",
                         "MTE"
                       ))
    
  })
  
  observeEvent(input$threshold, {
    # Your logic to update the checkbox selections based on slider value changes
    # For instance, you can add the selected values back if they still exist after filtering:
    
    # TE checkboxes
    valid_TE_fams <- rv$counts$TE_fam  # These are the valid TE_fams after filtering
    current_selected_TE_fams <- isolate(input$TE_checkboxes)  # Get the current selections
    valid_selected_TE_fams <- current_selected_TE_fams[current_selected_TE_fams %in% valid_TE_fams]
    
    # Update the TE checkboxes with valid selections
    updateCheckboxGroupInput(session, "TE_checkboxes", selected = valid_selected_TE_fams)
    
    # Repeat similar steps for the Motif checkboxes:
    valid_motifs <- motifs()  # These are the valid motifs after filtering
    current_selected_motifs <- isolate(input$motif_checkboxes)  # Get the current selections
    valid_selected_motifs <- current_selected_motifs[current_selected_motifs %in% valid_motifs]
    
    # Update the Motif checkboxes with valid selections
    updateCheckboxGroupInput(session, "motif_checkboxes", selected = valid_selected_motifs)
  })
  
  
  observeEvent(input$select_all, {
    if (!is.null(rv$counts)) {
      updateCheckboxGroupInput(session, "TE_checkboxes", 
                               selected = rv$counts$TE_fam)
      valid_motifs <- motifs()  # These are the valid motifs after filtering
      current_selected_motifs <- isolate(input$motif_checkboxes)  # Get the current selections
      valid_selected_motifs <- current_selected_motifs[current_selected_motifs %in% valid_motifs]
      
      # Update the Motif checkboxes with valid selections
      updateCheckboxGroupInput(session, "motif_checkboxes", selected = valid_selected_motifs)
    }
  })
  
  observeEvent(input$deselect_all, {
    if (!is.null(input$TE_checkboxes)) {
      updateCheckboxGroupInput(session, "TE_checkboxes", selected = "")
    }
  })
  
  observeEvent(input$select_LTR, {
    if (!is.null(rv$counts)) {
      selected <- rv$counts$TE_fam[rv$counts$TE_class == "LTR"]
      updateCheckboxGroupInput(session, "TE_checkboxes", selected = selected)
    }
  })
  
  observeEvent(input$select_LINE, {
    if (!is.null(rv$counts)) {
      selected <- rv$counts$TE_fam[rv$counts$TE_class == "LINE"]
      updateCheckboxGroupInput(session, "TE_checkboxes", selected = selected)
    }
  })
  
  observeEvent(input$select_DNA, {
    if (!is.null(rv$counts)) {
      selected <- rv$counts$TE_fam[rv$counts$TE_class == "DNA"]
      updateCheckboxGroupInput(session, "TE_checkboxes", selected = selected)
    }
  })
  output$distPlot <- renderPlot({
    if (!is.null(input$motif_checkboxes)){
      data <- filtered_data() %>% filter(TE_fam %in% input$TE_checkboxes,
                                         Motif %in% input$motif_checkboxes) 
    } else {
      data <- filtered_data() %>% filter(TE_fam %in% input$TE_checkboxes,
                                         Motif %in% c(
                                           "TATA",
                                           "DPE",
                                           "dInr",
                                           "hInr", 
                                           "dTCT",
                                           "hTCT",
                                           "PB",
                                           "GAGA",
                                           "MTE"
                                         )) 
    }
    
    # Generate the plot
    ggplot(data, aes(x = Motif_start, fill = KD)) +
      geom_histogram(binwidth = 1) +
      facet_wrap(~Motif) +
      xlim(c(-40, 40))+
      theme_bw()+
      theme(text = element_text(size = 20))+
      labs(fill = "Condition")
  })
}


# Run the application 
shinyApp(ui = ui, server = server)
