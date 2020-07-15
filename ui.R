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
        tags$style(type='text/css', '#AdvancedFragment {white-space: pre-wrap;}'),
        tags$style(type='text/css', '#SimpleFragment {white-space: pre-wrap;}'),
        #Main tab pages
        navbarPage(
            title="piRNAi design",
            ###Theme of shiny
            theme = shinytheme("sandstone"),
            
            ###Simple
            tabPanel("Simple",
                     mainPanel(
                         textAreaInput("geneinput", label = "Pick a gene", value = "", resize="none", placeholder= "WormbaseID, transcript or common name", rows=1),
                         actionButton("actiongenesearch", label = "Search gene"),
                         hr(),
                         uiOutput("DesignControls"),
                         hr(),
                         #tableOutput(otherPis),
                         htmlOutput("SelPiTabSummary"),
                         tableOutput('SelPiTab'),
                         verbatimTextOutput("SimpleFragment"),
                         uiOutput("downloadseq")
                         )
            ),
            
            tabPanel("Advanced",
                     mainPanel(
                         h3("Make your own construct:"),
                    fluidRow(
                        column(4,
                        textAreaInput("piRNAseq1", label="", rows=1, cols=21, resize="none", placeholder = "First piRNAi"),
                        textAreaInput("piRNAseq2", label="", rows=1, cols=21, resize="none", placeholder = "Second piRNAi"),
                        textAreaInput("piRNAseq3", label="", rows=1, cols=21, resize="none", placeholder = "Third piRNAi"),
                        textAreaInput("piRNAseq4", label="", rows=1, cols=21, resize="none", placeholder = "Fourth piRNAi"),
                        textAreaInput("piRNAseq5", label="", rows=1, cols=21, resize="none", placeholder = "Fifth piRNAi"),
                        textAreaInput("piRNAseq6", label="", rows=1, cols=21, resize="none", placeholder = "Sixth piRNAi")),
                        column(8,
                        verbatimTextOutput("AdvancedFragment"),
                        uiOutput("downloadconstruct"))),
                        
                     fluidRow(
                         column(2,actionButton("actionclean", label = "Clean boxes")),
                         column(2,actionButton("actionconstruct", label = "Produce piRNAi fragment"))),
                         verbatimTextOutput("AdvancedErrorMessage"),
                         hr(),
                         textAreaInput("Advancedgeneinput", label = "Pick a gene", value = "", resize="none", placeholder= "WormbaseID, transcript or common name", rows=1),
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
