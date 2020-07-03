########piRNA User Interface####
###Amhed Vargas
###amhed.velazquez@kaust.edu.sa
###UI

#Load libraries
library(shiny)
library(shinythemes)

# Define User interface
shinyUI(
    fluidPage(
        ##Costum extra styles: single sliders background and title of navbar 
        tags$style(type = 'text/css', 
                   ".js-irs-none .irs-single, .js-irs-none .irs-bar-edge, .js-irs-none .irs-bar {
                          background: transparent;
                          border-top-color: transparent;
                          border-bottom-color: transparent;
                          border-left-color: transparent;
                          border-right-color: transparent}
               .navbar-default .navbar-brand:hover {color: #ffffff;}
               "),
        #Main tab pages
        navbarPage(
            title="piRNAi design",
            ###Theme of shiny
            theme = shinytheme("sandstone"),
            
            ####Genome browser
            tabPanel("piRNAi",
                     mainPanel(
                         textAreaInput("genesearch", label = "Pick a gene", value = "", placeholder= "WormbaseID, transcript or common name", rows=1),
                         actionButton("actiongenesearch", label = "Search gene"),
                         hr(),
                         uiOutput("DesignControls"),
                         hr(),
                         tableOutput(otherPis),
                         verbatimTextOutput("PartialResult"),
                         uiOutput("downloadseq")
                         )
            ),
            ###About
            tabPanel("About",
                     mainPanel(
                         h3("The app"),
                         HTML("<p align=\"justify\">Still in development
                      </p>")
                     )
            )
        ),
        hr(),
        HTML("<a href=\"https://syngenbio.kaust.edu.sa/\">Syntetic genome biology laboratory @KAUST</a><br>"),
        HTML("<a href=\"http://www.wormbuilder.org/\">Wormbuilder</a><br>"),
        HTML("<a href=\"mailto:amhed.velazquez@kaust.edu.sa\">Contact us!</a>")
        
    )
)
