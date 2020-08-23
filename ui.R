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
        tags$head(
            tags$link(rel="stylesheet",type = "text/css", href="bootstrap.min.css")
        ),
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
        tags$style(type='text/css', '#piBoxes .form-group {margin-bottom: 0px; margin-top: 0px;}'),
        #Main tab pages
        navbarPage(
            title="piRNAi cluster designer",
            ###Theme of shiny
            #theme = shinytheme("sandstone"),
            
            ###Simple
            tabPanel("Simple",
                     mainPanel(
                         textAreaInput("geneinput", label = "Target gene", value = "", resize="none", placeholder= "WormbaseID, transcript or gene name", rows=1),
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
                        ##Chose cluster
                        
                        radioButtons("clustercon", label = HTML("Select piRNA cluster
                                                     [<a href=\"\" onclick=\"$('#explain_cluster').toggle(); return false;\">info</a>]
                                                     "),
                                     choices = list("21ur-1224 (6 sites)" = 1, "21ur-1692 (6 sites)" = 2, "21ur-8831 (7 sites)" = 3, "21ur-1610 (8 sites)" = 4), selected = 1, width='100%', inline= TRUE),
                        HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_cluster\">
            We recommend to use the cluster 21ur-1224.
                                     </div></p>
                     "),
                        conditionalPanel(condition = "input.clustercon==1",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 1"),
                        textAreaInput("piRNAseq2_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 2"),
                        textAreaInput("piRNAseq3_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 3"),
                        textAreaInput("piRNAseq4_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 4"),
                        textAreaInput("piRNAseq5_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 5"),
                        textAreaInput("piRNAseq6_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 6")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==2",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 1"),
                        textAreaInput("piRNAseq2_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 2"),
                        textAreaInput("piRNAseq3_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 3"),
                        textAreaInput("piRNAseq4_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 4"),
                        textAreaInput("piRNAseq5_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 5"),
                        textAreaInput("piRNAseq6_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 6")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==3",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 1"),
                        textAreaInput("piRNAseq2_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 2"),
                        textAreaInput("piRNAseq3_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 3"),
                        textAreaInput("piRNAseq4_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 4"),
                        textAreaInput("piRNAseq5_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 5"),
                        textAreaInput("piRNAseq6_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 6"),
                        textAreaInput("piRNAseq7_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 7")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==4",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 1"),
                        textAreaInput("piRNAseq2_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 2"),
                        textAreaInput("piRNAseq3_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 3"),
                        textAreaInput("piRNAseq4_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 4"),
                        textAreaInput("piRNAseq5_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 5"),
                        textAreaInput("piRNAseq6_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 6"),
                        textAreaInput("piRNAseq7_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 7"),
                        textAreaInput("piRNAseq8_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 8")),
                        br()
                        )),
                        column(8,
                        verbatimTextOutput("AdvancedFragment"),
                        uiOutput("downloadconstruct"))
                        ),
                        
                     fluidRow(
                         column(2,actionButton("actionclean", label = "Reset")),
                         column(2,actionButton("actionconstruct", label = "Generate piRNAi cluster"))),
                         verbatimTextOutput("AdvancedErrorMessage"),
                         hr(),
                         textAreaInput("Advancedgeneinput", label = "Target gene", value = "", resize="none", placeholder= "WormbaseID, transcript or gene name", rows=1),
                         actionButton("actionAdvsearch", label = "Search for specific piRNAs"),
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
