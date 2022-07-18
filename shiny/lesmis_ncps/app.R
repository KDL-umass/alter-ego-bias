library(shiny)

source('lesmis.R')

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Les Mis!"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "n_ncps",
                  label = "Number of NCPs:",
                  min = 0,
                  max = graph.properties$n,
                  value = 11)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "exposureGraphPlot")
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Plot of Les Mis dialog graph ----
  # NCPs marked in orange
  # unexposed nodes marked in light blue 
  # NCP exposure 
  # with requested number of NCPs
  
  # This expression that generates a graph plot is wrapped in a call
  # to renderPlot to indicate that:
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$n_ncps) change
  # 2. Its output type is a plot
  
  output$exposureGraphPlot <- renderPlot({
    
    #adversaries.deg <- unlist(list(determine.adversaries(graph.properties, ncp.params)))
    #adversaries <- matrix(0,1,graph.properties$n)
    #adversaries[which(adversaries.deg==1)] <- 1
    
    #random NCP selection
    n_adversaries <- input$n_ncps
    random.adversaries <- sample(1:graph.properties$n, n_adversaries, replace=FALSE)
    adversaries <- matrix(0,1,graph.properties$n)
    adversaries[random.adversaries] <- 1
    
    treatment <- treatment.assignment(graph.properties$g, clusters)
    treatment.assignments <- treatment[clusters]
    
    uncovered.vertices <- 1 - adversaries %*% graph.properties$adj - adversaries
    pt.uncovered <- sum(uncovered.vertices == 1)/graph.properties$n
    #prepare.for.plots(g, adversaries, ncp.params$ncp.exposure, treatment.assignments, labels=TRUE)
    exposure.params <- exposure.probs(ncp.params, graph.properties, treatment.assignments, adversaries)
    print(exposure.params$ncp.exposure.neighbors)
    #V(g)$color <- adversaries.deg

    bdr <- rep("black", length(V(g)))
    if(length(treatment.assignments) > 2) bdr <- ifelse(treatment.assignments, "green", "black")
    
    g$palette <- grey.colors(100)
    # ncp.params$ncp.exposure.neighbors
    V(g)$color <- 100 - exposure.params$ncp.exposure.neighbors * 100 # grey.colors runs dark to light
    V(g)$color[which(exposure.params$ncp.exposure.neighbors == 0)] <- "lightblue"
    V(g)$color[which(adversaries > 0)] <- "orange"
    
    
    lo <- layout_with_kk(g) # create a layout
    lo <- norm_coords(lo, ymin=-1, ymax=1, xmin=-1, xmax=1)
    plot(g, layout=lo*2, vertex.color=V(g)$color, vertex.frame.color=bdr)
    
  })
  
}

shinyApp(ui = ui, server = server)