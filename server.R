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
    system(paste("mkdir -p WorkingSpace/users/",session_id,sep=""))
    
    ###On exit, force the remove of directory
    ##Attention: Interactive sessions create problems because the listening of the server stops in main directory and not in sub directory
    session$onSessionEnded(function(){
        system(paste("rm -rf WorkingSpace/users/",session_id,sep=""))
    }
    )
    
    ###
    ###Data
    Genes=read.table("WorkingSpace/Genes.txt",sep="\t", header= FALSE)
    locus=unique(as.character(Genes[,6]))
    locus=locus[-c(which(locus == ""))]
    #file=paste(c("DataBase/",as.character(Genes[1,4]),"_",as.character(Genes[1,5]),"_",as.character(Genes[1,6]),".txt"),sep="",collapse="")
    
    
    ##Main search function
    observeEvent(input$actiongenesearch, {
        type="Not"
        wbid=""
        mygene =input$geneinput
        
        if(as.character(mygene) %in% as.character(Genes[,4])){type=4}
        if(as.character(mygene) %in% as.character(Genes[,5])){type=5}
        if(as.character(mygene) %in% locus){type=6}
        
        if(type == "Not"){
            output$DesignControls <- renderUI({
                HTML("<b>Gene not found</b>")
                })
            }else{
                wbid=as.character(unique(Genes[which(mygene==Genes[,type]),4]))
        output$DesignControls <- renderUI({
            fluidRow(
                
                selectInput("isoform", label = HTML("<b>Select Isoform
                                                           [<a href=\"\" onclick=\"$('#explain_isoform').toggle(); return false;\">info</a>]
                                                           </b>"), 
                            unique(Genes[which(as.character(Genes[,4])==wbid),5])),
                
                HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_isoform\">
            Isoforms based on the WormBase 270 annotations</div></p>
                     "),
                selectInput("selectMM", label = HTML("<b>Uniqueness of sequence
                                                           [<a href=\"\" onclick=\"$('#explain_uniqueness').toggle(); return false;\">info</a>]
                                                           </b>"), 
                            choices = list("Up to 5 Mismatches across genome" = 1, "Up to 4 Mismatches across genome" = 2, "Up to 3 Mismatches across genome" = 3, "Up to 2 Mismatches across genome" = 4,
                                           "Up to 1 Mismatches across genome" = 5,"Up to 0 Mismatches across genome" = 6),
                            selected = 1),
                HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_uniqueness\">
            Select the specificity of piRNAi fragments. Each piRNAi fragment was mapped across the C. elegans genome and verified its uniqueness up uo n Mismatches.
                                                 </div></p>
                     "),
                radioButtons("cluster", label = HTML("Select piRNA cluster
                                                     [<a href=\"\" onclick=\"$('#explain_cluster').toggle(); return false;\">info</a>]
                                                     "),
                             choices = list("21ur-1224" = 1), selected = 1, width='100%'),
                HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_cluster\">
            For the moment, we use the cluster 21ur-1224 as a template to express 6 piRNAis fragments that are antisente to the transcript being targeted
                                     </div></p>
                     "),
            checkboxInput("FlaControl", label = HTML("<b>Design control fragment
                                                               [<a href=\"\" onclick=\"$('#explain_control').toggle(); return false;\">info</a>]
                                                               </b>"), value = FALSE, width='100%'),
            HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_control\">
            Reverse complement sequences so they have the same orientation as the coding sequence
                                     </div></p>"),
            actionButton("actionPI", label = "Create piRNAi fragment")
            
            )
        })
        }
    }, ignoreInit = T)
    
    
    ##Main search function
    observeEvent(input$actionPI, {
        ErrorFLag=0
        
        matches=as.integer(input$selectMM)
        mm=c(5,4,3,2,1,0)[matches]
        isoform = input$isoform
        wbid = as.character(unique(Genes[which(isoform==Genes[,5]),4]))
        loc = as.character(unique(Genes[which(isoform==Genes[,5]),6]))
        
        file=paste(c("DataBase/",as.character(wbid),"_",as.character(isoform),"_",as.character(loc),".txt"),sep="",collapse="")
        
        tab=read.table(file,sep="\t",header=F)
        tab[,4]=as.integer((unlist(strsplit(as.character(tab[,2]),";")))[c(FALSE,TRUE)])
        tab[,2]=as.character((unlist(strsplit(as.character(tab[,2]),";")))[c(TRUE,FALSE)])
        
        Seltab=tab[which(tab[,3]>=mm),]
        
        if(nrow(Seltab)< 6){
            output$ErrorMessage <- renderText({
                paste("Error: Not enough piRNAi fragments were found with the characterisitics described. Try to change to other parameters")
            })
            ErrorFlag=1
        }else{
            Seltab=Seltab[order(Seltab[,1]),]
            pos=quantile(Seltab[,1],c(0,.2,.4,.6,.8,1))
            idx=c()
            idx=append(idx,which.min(abs(Seltab[,1]-pos[1])))
            idx=append(idx,which.min(abs(Seltab[,1]-pos[2])))
            idx=append(idx,which.min(abs(Seltab[,1]-pos[3])))
            idx=append(idx,which.min(abs(Seltab[,1]-pos[4])))
            idx=append(idx,which.min(abs(Seltab[,1]-pos[5])))
            idx=append(idx,which.min(abs(Seltab[,1]-pos[6])))
            
            if(length(which(dist(Seltab[idx,1])<=20)) > 0){
                
                output$ErrorMessage <- renderText({
                    paste("Error: The sole piRNAi fragments found seem overlapping. Try to change its distribution and or the parameters")
                })
                ErrorFlag=1
                
                }
            }
        
        ##ProduceApE output
            if(ErrorFlag==0){
                output$ErrorMessage <- renderText({
                    paste("Sequences used:",Seltab[idx,2])
                })
                
                output$downloadseq <- renderUI({
                        write(paste(">Optimized_cDNA:Codon-",optsin,"\n",aaaads,SeqtoOpt,"\n",sep="",collapse=""),paste("WorkingSpace/users/",session_id,"/piRNAs.ape", sep=""))
                        
                    downloadButton('DownApeOut', 'Download Ape File')
                })
                
                
            }
        
    }, ignoreInit = T)
    
    
    ##Retrieve output ape
    output$DownApeOut <- downloadHandler(
        filename <- function() {
            paste("piRNAi", "ape", sep=".")
        },
        
        content <- function(file) {
            file.copy(paste("WorkingSpace/users/",session_id,"/piRNAs.ape", sep=""), file)
        },
    )
})  
