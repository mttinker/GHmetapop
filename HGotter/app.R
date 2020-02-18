#
# This is a Shiny web application for running sea otter population simulations
#  on the Island Archipelego of Haida Gwaii
library(shiny)
library(DT)
library(shinyWidgets)
library(rhandsontable)
library(gtools)
library(stats)
library(mvtnorm)
library(boot)
library(ggplot2)
library(rgdal)
library(dplyr) 
library(ggrepel)
library(ggsn)
library(reshape2)

Bdat = read.csv("./data/GHBlockdata.csv", header = TRUE)  # Data for coastal blocks
NBlks <- nrow(Bdat)
arealist <- as.list(seq(1,NBlks)); names(arealist) <- paste0("CS_",seq(1,NBlks))
# User interface part of app
ui <- shinyUI(
    navbarPage(
        theme = shinythemes::shinytheme("cerulean"),
        "Haida Gwaii Sea Otter Population Model (HGotter V 1.0c)",
        tabPanel("Set up Model Simulations",
                 # Sidebar panel: user input
                 sidebarPanel(
                     actionButton("RunSim", "Run Simulations Now", class = "btn-primary"),
                     helpText("NOTE: Please review and adjust user-parameters BEFORE running simulations"),
                     uiOutput("Reps_slider"),  
                     uiOutput("Nyrs_slider"), 
                     uiOutput("Ephase_slider"),
                     uiOutput("ImmRt_slider"),
                     uiOutput("V_sp_slider"),
                     uiOutput("Kmean_slider"),
                     uiOutput("Ksig_slider"),
                     uiOutput("Rmax_slider"),                     
                     uiOutput("Estoch_slider"),
                     uiOutput("Theta_slider"),
                     uiOutput("InitN_slider"),
                     selectInput(
                         "InitSect",
                         label = ("Select Initially Occupied Coastal Section(s), see map"),
                         choices = arealist,selected = 1, multiple = TRUE, selectize = FALSE,
                         width='100%'),
                     width = 4
                 ),
                 # Main panel: map of coastal sections 
                 mainPanel(
                     img(src='HGmap.png', align = "left"),
                     width = 8
                 )
        ),
        tabPanel("Model Output: graphs",             
                 tabsetPanel(
                         tabPanel("Simulation Results, Population Trends",
                                  plotOutput(outputId="sumgraph1",height = 800,width = 1000)),
                         tabPanel("Simulation Results, Range Expansion",
                                  plotOutput(outputId="sumgraph2",
                                             height = 800,width = 1000)),
                         tabPanel("Simulation Results, Projected Density Map",
                                  plotOutput(outputId="sumgraph3",
                                             height = 800,width = 800))
                     )
                 ),
        tabPanel("Model Output: tables",             
                 tabsetPanel(
                     tabPanel("Table 1: Estimated total abundance by year", 
                              tableOutput('sumtable1')),
                     tabPanel("Table 2: Summary by coastal section, final year", 
                              tableOutput('sumtable2'))
                 )
        )
    )
)
# Server part of app
sv <- shinyServer(function(input, output, session){
    # Create reactive values object to store sim results data frame
    values <- reactiveValues(Pop_Overall = NULL,dfDens = NULL,
                             Hab_Blocks = NULL,Hab_Blocks_Fin = NULL)
    source('runsims.r', local = TRUE)
    # Run sims when button pushed
    observeEvent(input$RunSim, {
        req(input$InitSect)
        datasim <- runsims() 
        # HGblk1 = merge(HGblk, datasim, by.x='BlockID', by.y = 'Block')
        Pop_Overall <- values$Pop_Overall 
        dfDens <- values$dfDens        
        Hab_Blocks <- values$Hab_Blocks
        Hab_Blocks_Fin <- values$Hab_Blocks_Fin
        Nyrs <- nrow(Pop_Overall)
        titletxt <- paste0("Sea Otter Population Projection, ", Nyrs," Years")
        maxN <- ceiling(Pop_Overall$upper[Nyrs]/100)*100
        maxD <- ceiling(100*max(dfDens$Density))/100
        plt1 <- reactive({
            ggplot(Pop_Overall, aes(Year, Mean))+
                geom_line(data=Pop_Overall)+
                geom_ribbon(data=Pop_Overall,aes(ymin=lower,ymax=upper),alpha=0.2)+
                geom_ribbon(data=Pop_Overall,aes(ymin=CImeanLo,ymax=CImeanHi),alpha=0.3)+
                ylim(0,maxN) +  
                xlab("Year") +
                ylab("Expected Abundance") +
                ggtitle(titletxt) + 
                theme_classic(base_size = 13)
        })
        plt2 <- reactive({
            ggplot(dfDens, aes(Year, Block)) +
                geom_tile(aes(fill = Density), color = "white") +
                scale_fill_gradient(low = "white", high = "steelblue",limits=c(0, maxD)) +
                xlab("Year in Future") +
                ylab("Coastal Block #") +
                theme(legend.title = element_text(size = 13),
                      legend.text = element_text(size = 13),
                      plot.title = element_text(size=15),
                      axis.title=element_text(size=13),
                      axis.text.y = element_text(size=8),
                      axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(fill = "Mean Expected Density",
                     title=paste0("Projected Density by Block (otters/km2) by Year")) 
        })
        output$sumgraph1 <- renderPlot({
            plt1() 
        })
        output$sumgraph2 <- renderPlot({
            plt2() 
        })
        output$sumgraph3 <- renderImage({
            if (is.null(values$Pop_Overall)) {
                return(NULL)
            } else {
                return(list(
                    src = "./www/mapdens.png",
                    contentType = "image/png",
                    width = 800,
                    height = 800,
                    alt = "Map of Haida Gwaii showing projected otter density"
                ))
            }
            }, deleteFile = TRUE)    
    })
    # Output tables
    #
    # Summary table 1 (rangewide sums by year), updated when sims completed
    output$sumtable1 <- renderTable(values$Pop_Overall) 
    # Summary table 2 (section sums at Final year), updated when sims completed
    output$sumtable2 <- renderTable(values$Hab_Blocks_Fin) 
    #
    # User input sliders etc. 
    output$Reps_slider <- renderUI({
        sliderInput("Reps","Number of times to iterate this population simulation",
                     min = 50, max = 1000, value = 500, step = 10, round = TRUE)
    })
    output$Nyrs_slider <- renderUI({
        sliderInput("Nyrs","Number of years for each simulation", min = 10, max = 100,
                    value = 25, step = 5, round = TRUE)
    })    
    output$Ephase_slider <- renderUI({
        sliderInput("Ephase","Number of years expected population establishment phase", 
                    min = 1, max = 20, value = 5, step = 1, round = TRUE)
    })
    output$ImmRt_slider <- renderUI({
        sliderInput("ImmRt","Average annual immigration rate (from outside HG))", 
                    min = 0, max = 10, value = 0, step = .1, round = FALSE)
    })
    output$V_sp_slider <- renderUI({
        sliderInput("V_sp","Mean expected range exansion rate (km/yr)", 
                    min = .5, max = 8, value = 3, step = .1, round = FALSE)
    })    
    output$Kmean_slider <- renderUI({
        sliderInput("Kmean","Overall mean density at K (otters/km2)", 
                    min = .5, max = 10, value = 3.5, step = .1, round = FALSE)
    })   
    output$Ksig_slider <- renderUI({
        sliderInput("Ksig","Spatial variation (sd) in density at K (otters/km2)", 
                    min = .1, max = 3, value = 1, step = .1, round = FALSE)
    })  
    output$Rmax_slider <- renderUI({
        sliderInput("Rmax","Maximum annual growth rate (at low density)", 
                    min = .1, max = .25, value = .2, step = .01, round = FALSE)
    })      
    output$Estoch_slider <- renderUI({
        sliderInput("Estoch","Environmental stochasticity (sd in annual growth rate)", 
                    min = .01, max = .2, value = .05, step = .01, round = FALSE)
    })  
    output$Theta_slider <- renderUI({
        sliderInput("Theta","Value of theta (for theta-logistic growth)", 
                    min = .5, max = 2, value = 1, step = .1, round = FALSE)
    })  
    output$InitN_slider <- renderUI({
        sliderInput("InitN","Initial population size (total otters in Year 1)", 
                    min = 2, max = 50, value = 20, step = 1, round = TRUE)
    })  
})
# Create Shiny app ----
shinyApp(ui = ui, server = sv)
