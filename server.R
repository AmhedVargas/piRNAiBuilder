########Wormtracks server####
###Amhed Vargas
###amhed.velazquez@kaust.edu.sa
###Server

##Required packages
#install.packages("shiny")
#install.packages("shinythemes")


#Load libraries
library(shiny)
library(shinythemes)

shinyServer(function(input, output, session) {
    
    ##Retrieve unique ID for the session
    session_id <- session$token
    
    ##Create temporary folder for unique user
    system(paste("mkdir -p WormT/users/",session_id,sep=""))
    
    ###On exit, force the remove of directory
    ##Attention: Interactive sessions create problems because the listening of the server stops in main directory and not in sub directory
    session$onSessionEnded(function(){
        system(paste("rm -rf WormT/users/",session_id,sep=""))
    }
    )
    
    ###
    ###Data explorer
    ###Format input table
    ChrisPAT=read.table("WormT/PATCsGenesChristianData.tsv",sep="\t", header= TRUE)
    rownames(ChrisPAT)=as.character(ChrisPAT[,2])
    ##Assign Yes, No, unknown insted of TRUE, FALSE, NA
    ChrisPAT$GermlineExpression[which(!(ChrisPAT$GermlineExpression))] <- "No"
    ChrisPAT$GermlineExpression[which(ChrisPAT$GermlineExpression == TRUE)] <- "Yes"
    ChrisPAT$GermlineExpression[which(is.na(ChrisPAT$GermlineExpression))] <- "unknown"
    
    ##Click action - inputs to reactive variable geneval that is used in a function to generate table
    geneweb_tooltip <- function(x) {
        if (is.null(x)) return(NULL)
        if (is.null(x$Wormbase.ID)) return(NULL)
        # Pick out the gene with this ID
        geneval$id <- as.character(x$Wormbase.ID)
        return(NULL)
    }
    
    geneval <- reactiveValues(id=c(""))
    
    ##Function that returns the gene info from main table and formats it into html
    geneinfo = function(x){
        if (is.null(x)) return(NULL)
        gidi=""
        if(x %in% as.character(ChrisPAT$Wormbase.ID)){gidi = x}
        if(x %in% as.character(ChrisPAT$Gene.name)){gidi = rownames(ChrisPAT[which(as.character(ChrisPAT$Gene.name) == x),])}
        if(gidi ==""){return(paste("Gene not found! try with Wormbase ID"))}
        
        gene <- ChrisPAT[as.character(gidi),]
        
        paste0("<b>Gene: ", gene$Gene.name, "</b><br>Wormbase ID: ",
               "<a href=\"https://www.wormbase.org/species/c_elegans/gene/", gene$Wormbase.ID,"\">", gene$Wormbase.ID, "</a><br>",
               "Transcript analyzed: ", gene$Transcript.name, "<br>",
               "Chromosome: ", gene$Chromosome, "<br>",
               "Relative position: ", gene$Relative.loc.chromosome, "<br>",
               "Gene size: ", gene$Bin.size, "<br>",
               "Number of phased bases: ", gene$Bases.with.PATC.55, "<br>",
               "Phased bases in gene: ", as.integer(gene$Frequency.of.PATCs.55 * 100), "%<br>",
               "Total PATC score: ", gene$Total.value.of.PATC.algorithm, "<br>",
               "PATC density: ", gene$PATC.density, "<br>",
               "Germline expression: ", gene$GermlineExpression, "<sup>1</sup><br>",
               "RPKM Oocyte expression<sup>2</sup>: ", gene$stoekius_oocyte_rpkm, "<br>",
               "RPKM Sperm expression<sup>3</sup>: ", gene$spermatogenic_gonad_fem.3_RPKM_Ortiz_2014,"<br>","<br>",
               "<p align=\"justify\">",
               "<font size=\"2\">",
               "<sup>1</sup>: Oocyte Reads Per Kilobase of transcript per Million mapped reads (RPKM) > 2<br>",
               "<sup>2</sup>: Stoeckius, <i>et al.</i> (2019). Large-scale sorting of <i>C. elegans</i> embryos reveals the dynamics of small RNA expression.
                                  <i>Nat. Methods 6(10)</i>: 745-751<br>",
               "<sup>3</sup>: Ortiz, <i>et al.</i> (2014). A New Dataset of Spermatogenic <i>vs.</i> Oogenic Transcriptomes in the Nematode <i>Caenorhabditis elegans</i>.
             <i>G3: GENES, GENOMES, GENETICS 4(9)</i>: 1765-1772<br>",
               "</font>",
               "</p>"
        )
    }
    
    geteinfo = function(x){
        if (is.null(x)) return(NULL)
        gidi=""
        if(x %in% as.character(ChrisPAT$Wormbase.ID)){gidi = x}
        if(x %in% as.character(ChrisPAT$Gene.name)){gidi = rownames(ChrisPAT[which(as.character(ChrisPAT$Gene.name) == x),])}
        if(gidi ==""){return(paste("Gene not found! try with Wormbase ID"))}
        
        gene <- ChrisPAT[as.character(gidi),]
        
        #runjs(paste("igv.browser.search(",gene$Chromosome,":",gene$Gene.start,"-",gene$Gene.end,")", sep=""))
        #paste0("<script>igv.browser.search(",gene$Chromosome,":",gene$Gene.start,"-",gene$Gene.end,")</script>")

        
        igvloc=paste(gene$Chromosome,":",gene$Gene.start,"-",gene$Gene.end,sep="")
        session$sendCustomMessage("gene-coordinates", igvloc)
    }
    
    ##Observers for action button search
    observeEvent(input$actionsearch, {
        mygene =input$genetext
        output$geneid <- renderUI({
            HTML(
                geneinfo(mygene)
            )
        })
    }, ignoreInit = T)
    
    ##Observers for action button browser
    observeEvent(input$actionbrow, {
        gege =input$genebrow
        output$genebrowsearch <- renderUI({
            HTML(
                geteinfo(gege)
            )
        })
        
        getPage<-function(dir) {
            
            return(includeHTML(paste(paste("WormT/tracks/",dir,"/igv.html", sep=""))))
        }
        
        output$igvcoord <-renderUI({getPage(input$selectTRACK)})
        
        
    }, ignoreInit = T)
    
    ##Observe Js value
    observeEvent(input$jsValue, {
        cat("\nMessage Received\n")
        cat(paste(input$jsValue),sep="\n")
        output$igv_id <- renderText({
            paste("Your sequence to sinthetize is:", input$jsValue)
        })
    })
    
    #output$tablepiRNA = DT::renderDataTable(DT::datatable({
        #tt=read.table(paste("WormT/users/",session_id,"/tab.file", sep=""),sep="\t",header=TRUE)
        #colnames(tt)=c("Sequence ID", "Number of bases in phase", "Phasing frequency", "Total PATC value", "Phasing density")
        #tt[,3]=round(tt[,3], 2)
        #tt[,4]=round(tt[,4], 2)
        #tt[,5]=round(tt[,5], 2)
        #tt
        
    #}, options = list(dom = 't')))
    
    ##External function to include the HTML output
    #It has to be now relocated within temporary directory
    getPage<-function(dir) {
        
        return(includeHTML(paste(paste("WormT/tracks/",dir,"/igv.html", sep=""))))
    }
    #output$igvcoord<-renderUI({getPage(session_id)})
    
    output$igvcoord <-renderUI({getPage(input$selectTRACK)})
    
    #Observer reactive to changes in geneval to create table
    observeEvent(geneval$id, {
        mygene = geneval$id
        output$geneid<- renderUI({
            HTML(
                geneinfo(mygene)
            )
        })
    }, ignoreInit = T)
    
    observeEvent(input$actionbrow,{
        ###Add things to show output
        output$igvcoord <-renderUI({getPage(input$selectTRACK)})
        
    }, ignoreInit = T)
    
})  
