########piRNA User Interface####
###Amhed Vargas
###amhed.velazquez@kaust.edu.sa
###UI

#Load libraries
library(shiny)
library(shinythemes)
library(DT)

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
            
            ###Simple
            tabPanel("Simple",
                     mainPanel(
                         textAreaInput("geneinput", label = "Pick a gene", value = "", placeholder= "WormbaseID, transcript or common name", rows=1),
                         actionButton("actiongenesearch", label = "Search gene"),
                         hr(),
                         uiOutput("DesignControls"),
                         hr(),
                         #tableOutput(otherPis),
                         htmlOutput("SelPiTabSummary"),
                         tableOutput('SelPiTab'),
                         uiOutput("downloadseq")
                         )
            ),
            
            tabPanel("Advanced",
                     mainPanel(
                         h3("Make your own construct:"),
                         fluidRow(
                             column(width = 2,textAreaInput("piRNAseq1", label="", rows=1, placeholder = "First piRNAi")),
                             column(width = 2,textAreaInput("piRNAseq2", label="", rows=1, placeholder = "Second piRNAi")),
                             column(width = 2,textAreaInput("piRNAseq3", label="", rows=1, placeholder = "Third piRNAi")),
                             column(width = 2,textAreaInput("piRNAseq4", label="", rows=1, placeholder = "Fourth piRNAi")),
                             column(width = 2,textAreaInput("piRNAseq5", label="", rows=1, placeholder = "Fifth piRNAi")),
                             column(width = 2,textAreaInput("piRNAseq6", label="", rows=1, placeholder = "Sixth piRNAi"))
                         ),
                         actionButton("actionconstruct", label = "Produce piRNAi fragment"),
                         verbatimTextOutput("AdvancedErrorMessage"),
                         uiOutput("downloadconstruct"),
                         hr(),
                         textAreaInput("Advancedgeneinput", label = "Pick a gene", value = "", placeholder= "WormbaseID, transcript or common name", rows=1),
                         actionButton("actionAdvsearch", label = "Search piRNAi fragments"),
                         hr(),
                         uiOutput("AdvDesignControls"),
                         hr(),
                         DT::dataTableOutput('AllPiTab')
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
