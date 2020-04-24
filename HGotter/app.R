#
# This is a Shiny web application for running sea otter population simulations
#  on the Island Archipelego of Haida Gwaii
library(shiny)
library(DT)
library(shinyWidgets)
library(shinyBS)
library(shinyjs)
library(rhandsontable)
library(gtools)
library(stats)
library(mvtnorm)
library(boot)
library(ggplot2)
library(mapproj)
library(rgdal)
library(dplyr) 
library(ggrepel)
library(ggsn)
library(reshape2)
#
Bdat = read.csv("./data/GHBlockdata.csv", header = TRUE)  # Data for coastal blocks
NBlks <- nrow(Bdat)

instructfile = "./data/Instruct.txt"
Blurb1 = readChar(instructfile, file.info(instructfile)$size)
blurbfile = "./data/Blurb.txt"
Blurb2 = readChar(blurbfile, file.info(blurbfile)$size)

arealist <- as.list(seq(1,NBlks)); names(arealist) <- paste0("CS_",seq(1,NBlks))
# User interface part of app
ui <- shinyUI(
    navbarPage(
        theme = shinythemes::shinytheme("cerulean"),
        "Haida Gwaii Sea Otter Population Model (HGotter V 1.10)",
        tabPanel("Setup Model Simulations",
                 # Sidebar panel: user input
                 useShinyjs(),
                 sidebarPanel(
                     helpText(Blurb1),
                     div(style="display:inline-block; horizontal-align: left",
                         actionButton("RunSim", "Run Simulations Now", 
                         class = "btn-primary"),width=10),
                     div(style="display:inline-block; margin-left:100px"),
                     div(style="display:inline-block; horizontal-align: right",
                         downloadButton('GetManual', 'Download User Manual')),
                     div(style="margin-top:10px") ,
                     selectizeInput(
                         "InitSect",
                         label = ("Select initially occupied coastal section(s): refer to map"),
                         choices = arealist, multiple = TRUE, width='100%', # selected = 1, 
                         options = list(placeholder = 'click to select one or more sections',
                         'plugins' = list('remove_button'),'create' = TRUE,'persist' = FALSE)), 
                     uiOutput("Reps_slider"),
                     bsPopover("Reps_slider",
                               "Number of iterations: more information",
                                paste0("Increasing the number of replications of a simulation improves the precision ",
                                       "of model predictions, but will take longer to run. At least 100 reps is suggested"),
                               placement = "top", trigger = "hover"),
                     uiOutput("Nyrs_slider"), 
                     uiOutput("InitN_slider"),
                     uiOutput("ImmRt_slider"),
                     uiOutput("Ephase_slider"),
                     uiOutput("V_sp_slider"),
                     uiOutput("Kmean_slider"),
                     uiOutput("Ksig_slider"),
                     uiOutput("Rmax_slider"),                     
                     uiOutput("Estoch_slider"),
                     uiOutput("Theta_slider"),
                     width = 4
                 ),
                 # Main panel: map of coastal sections 
                 mainPanel(
                     helpText(Blurb2),
                     img(src='HGmap.png', width = 600,height = 800,
                         align = "center"),
                     width = 8
                 )
        ),
        tabPanel("Model Output GRAPHS",             
                 tabsetPanel(
                         tabPanel("Population Trend",
                                  plotOutput(outputId="sumgraph1",height = 800,width = 1000)),
                         tabPanel("Range Expansion",
                                  plotOutput(outputId="sumgraph2",
                                             height = 800,width = 1000)),
                         tabPanel("Density Map",
                                  plotOutput(outputId="sumgraph3",
                                             height = 800,width = 1200))
                     )
                 ),
        tabPanel("Model Output TABLES",             
                 tabsetPanel(
                     tabPanel("Table 1: Projected Sea Otter Abundance by Year",
                              downloadButton('download1',"Download Table 1"),
                              tableOutput('sumtable1')),
                     tabPanel("Table 2: Projected Abundance by Coastal Section in Final Year", 
                              downloadButton('download2',"Download Table 2"),
                              tableOutput('sumtable2'))
                 )
        )
    )
)
# Server part of app
sv <- shinyServer(function(input, output, session){
    # Create reactive values object to store sim results data frame
    values <- reactiveValues(Pop_Overall = NULL,dfDens = NULL,Tab1=NULL,
                             Hab_Blocks = NULL,Hab_Blocks_Fin = NULL)
    source('runsims.r', local = TRUE)
    # Allow user manual to be downloaded
    output$GetManual <- downloadHandler(
        filename = "HGotter_Manual.pdf",
        content = function(file) {
            file.copy("www/Manual.pdf", file)
        })
    # Run sims when button pushed (but only enable if coastal sections selected)
    observe({
        toggleState(id = "RunSim", condition = input$InitSect)
    })
    observeEvent(input$RunSim, {
        req(input$InitSect)
        datasim <- runsims() 
        # HGblk1 = merge(HGblk, datasim, by.x='BlockID', by.y = 'Block')
        Pop_Overall <- values$Pop_Overall 
        dfDens <- values$dfDens        
        Hab_Blocks <- values$Hab_Blocks
        Hab_Blocks_Fin <- values$Hab_Blocks_Fin
        Nyrs <- nrow(Pop_Overall)
        titletxt <- paste0("Projected Sea Otter Population, ", Nyrs," Years")
        subtxt <- paste0("Average expected abundance (line) ",
                            "with 95% CI for the mean (dark shaded band) ",
                            "and 95% CI for projection uncertainty (light shaded band)")
        maxN <- ceiling(Pop_Overall$upper[Nyrs]/100)*100
        maxD <- ceiling(100*max(dfDens$Density))/100
        plt1 <- reactive({
            ggplot(Pop_Overall, aes(Year, Mean))+
                geom_line(data=Pop_Overall)+
                geom_ribbon(data=Pop_Overall,aes(ymin=lower,ymax=upper),alpha=0.2)+
                geom_ribbon(data=Pop_Overall,aes(ymin=CImeanLo,ymax=CImeanHi),alpha=0.3)+
                ylim(0,maxN) +  
                xlab("Year") +
                ylab("Expected sea otter abundance") +
                ggtitle(titletxt, subtitle=subtxt) + 
                theme_classic(base_size = 13)
        })
        plt2 <- reactive({
            ggplot(dfDens, aes(Year, Block)) +
                geom_tile(aes(fill = Density), color = "white") +
                scale_fill_gradient(low = "white", high = "steelblue",limits=c(0, maxD)) +
                xlab("Year") +
                ylab("Coastal section") +
                theme(legend.title = element_text(size = 13),
                      legend.text = element_text(size = 13),
                      plot.title = element_text(size=15),
                      axis.title=element_text(size=13),
                      axis.text.y = element_text(size=8),
                      axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(fill = "Mean Expected Density",
                     title= "Projected Sea Otter Density by Coastal Section by Year (sea otters/km2)") 
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
                    alt = "Map of Haida Gwaii showing projected sea otter density"
                ))
            }
            }, deleteFile = TRUE)    
    })
    # Output tables
    #
    # Summary table 1 (rangewide sums by year), updated when sims completed
    output$sumtable1 <- renderTable(values$Tab1) 
    output$download1 <- downloadHandler(
        filename = function(){"Table1.csv"}, 
        content = function(fname){
            write.csv(values$Tab1, fname)
        }
    )
    # Summary table 2 (section sums at Final year), updated when sims completed
    output$sumtable2 <- renderTable(values$Hab_Blocks_Fin) 
    output$download2 <- downloadHandler(
        filename = function(){"Table2.csv"}, 
        content = function(fname){
            write.csv(values$Hab_Blocks_Fin, fname)
        }
    )    
    #
    # User input sliders etc. 
    output$Reps_slider <- renderUI({
        sliderInput("Reps","Number of times to iterate this population simulation",
                     min = 50, max = 1000, value = 100, step = 10, round = TRUE)
    })
    output$Nyrs_slider <- renderUI({
        sliderInput("Nyrs","Number of years for each simulation", min = 10, max = 100,
                    value = 25, step = 5, round = TRUE)
    })
    addPopover(session,"Nyrs_slider",
               "Number of years: more information",
               paste0("The population model will be run for 'N' years into the future. Increasing 'N' can provide ", 
                      "insights into conditions farther in the future, but results become less reliable the farther ",
                      "ahead in time the model is projected."),
               placement = "top", trigger = "hover")
    output$Ephase_slider <- renderUI({
        sliderInput("Ephase","Number of years expected population establishment phase", 
                    min = 1, max = 20, value = 5, step = 1, round = TRUE)
    })
    addPopover(session,"Ephase_slider",
               "Number of years of population establishment: more information",
               paste0("Newly established sea otter populations often experience an initial period ",
                      "of reduced growth and limited range expansion, as the population becomes established. ",
                      "This establishment period has varied from 5-15 years in previous re-introductions."),
               placement = "top", trigger = "hover")
    output$ImmRt_slider <- renderUI({
        sliderInput("ImmRt","Average annual immigration rate (from outside Haida Gwaii)", 
                    min = 0, max = 10, value = 0, step = .1, round = FALSE)
    })
    addPopover(session,"ImmRt_slider",
               "Immigration from outside Haida Gwaii: more information",
               paste0("The first sea otters at Haida Gwaii were immigrants that dispersed from otter populations in neighbouring regions. ",
                      "Population growth at Haida Gwaii can occur from just these initial colonizers, but it is possible that additional ",
                      "otters immigrating from outside Haida Gwaii will enhance population recovery. Here you can set the average ",
                      "number of new immigrants expected per year (this can be a decimal, or even 0 if one assumes no further immigration)."),
               placement = "top", trigger = "hover")
    output$V_sp_slider <- renderUI({
        sliderInput("V_sp","Average expected rate of range exansion (km/yr)", 
                    min = .5, max = 8, value = 3, step = .1, round = FALSE)
    })    
    addPopover(session,"V_sp_slider",
               "Rate of range expansion: more information",
               paste0("Habitat use by the initial sea otter population will likely be limited to a relatively small area of the coast. ",
                      "As this initial population grows, it's distribution (range of occupancy) will spread outwards along the coastline, ",
                      "encompassing more habitat. The rate of range expansion is measured as the speed at which the population front moves ",   
                      "along the coastline. In other populations, this expansion speed has varied from 1-5 km/year."),
               placement = "top", trigger = "hover")
    output$Kmean_slider <- renderUI({
        sliderInput("Kmean","Overall average density at 'K' (otters/km2)", 
                    min = .5, max = 10, value = 3.5, step = .1, round = FALSE)
    })   
    addPopover(session,"Kmean_slider",
               "Average density at 'K': more information",
               paste0("The long-term equilibrium density eventually reached by a sea otter population in a given area is called carrying capacity, or 'K'. ",
                      "The population density at 'K' is not constant, but varies as a function of local habitat quality and prey dynamics. The model accounts ",
                      "for variation in relative density at K due to differences in local habitat quality: this parameter allows the user to adjust the ", 
                      "archipelago-wide AVERAGE. In other sea otter populations, the regional average density at 'K' is 2-7 otters/km2 habitat."),
               placement = "top", trigger = "hover")    
    output$Ksig_slider <- renderUI({
        sliderInput("Ksig","Spatial variation (SD) in density at 'K' (otters/km2)", 
                    min = .1, max = 3, value = 1, step = .1, round = FALSE)
    })  
    addPopover(session,"Ksig_slider",
               "Spatial variation (SD) in density at 'K': more information",
               paste0("This parameter allows the user to specify the degree of variation in equilibrium densities across coastal sections. ",
                      "The higher the number, the greater the variation in equilibrium density among sections."),
               placement = "top", trigger = "hover")        
    output$Rmax_slider <- renderUI({
        sliderInput("Rmax","Maximum annual growth rate (at low density)", 
                    min = .1, max = .25, value = .2, step = .01, round = FALSE)
    })      
    addPopover(session,"Rmax_slider",
               "Maximum annual growth rate: more information",
               paste0("Sea otter populations tend to show the highest rate of growth at low densities: as local abundance increases, ", 
                      "the growth rate slows until it eventually reaches 0 when population abundance reaches 'K'. This parameter allows the user ",
                      "to adjust the maximum rate of growth (at low densities): in most sea otter populations this value is close to 0.2"),
               placement = "top", trigger = "hover")    
    output$Estoch_slider <- renderUI({
        sliderInput("Estoch","Environmental stochasticity (SD in annual growth rate)", 
                    min = .01, max = .2, value = .05, step = .01, round = FALSE)
    })  
    addPopover(session,"Estoch_slider",
               "Environmental stochasticity: more information",
               paste0("The average rate of growth for a recovering sea otter population in a given area can be predicted as a function of ",
                      "the local density with respect to carrying capacity, or 'K'. However, year-to-year variation in environmental ",
                      "conditions and prey population dynamics can lead to unpredictable deviations in growth rate, referred to as 'environmental ",
                      "stochasticity'. This parameter controls the degree of variation in growth rates: typical values are 0.02 - 0.08"),
               placement = "top", trigger = "hover")        
    output$Theta_slider <- renderUI({
        sliderInput("Theta","Value of 'theta' (for theta-logistic growth)", 
                    min = .5, max = 2, value = 1, step = .1, round = FALSE)
    })  
    addPopover(session,"Theta_slider",
               "Value of 'theta': more information",
               paste0("The average rate of growth for a recovering sea otter population in a given area can be predicted as a function of ",
                      "the local density with respect to carrying capacity, or 'K'. One of the parameters of this function is 'theta', ",
                      "which determines the onset of reduced growth rates at higher densities: values <1 lead to onset of reduced ",
                      "growth rates at fairly low densities, while values >1 mean that reductions in growth occur only at high densities."),
               placement = "top", trigger = "hover")        
    output$InitN_slider <- renderUI({
        sliderInput("InitN","Initial population size (total otters in Year 1)", 
                    min = 2, max = 50, value = 20, step = 1, round = TRUE)
    })  
    addPopover(session,"InitN_slider",
               "Initial population size: more information",
               paste0("This parameter allows the user to specify the number of sea otters believed to comprise the initial population, ", 
                      "from which future growth will occur. The initial population is founded by otters immigrating from neighbouring ",
                      "populations, but also includes any progeny produced by these initial colonizers."),
               placement = "top", trigger = "hover")     
})
# Create Shiny app ----
shinyApp(ui = ui, server = sv)
