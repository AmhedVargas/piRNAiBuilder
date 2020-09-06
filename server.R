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
library(Biostrings)
library(DT)

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
    Intr=read.table("WorkingSpace/Introns_min.tsv",sep="\t",header=F)
    rownames(Intr)=as.character(Intr[,1])
    ##Extranames
    alias=read.table("WorkingSpace/Alias.txt",sep="\t", header=F, stringsAsFactors=F)
    rownames(alias)=as.character(alias[,1])
    
    ##Main search function
    observeEvent(input$actiongenesearch, {
        output$ErrorMessage <- renderText({})
        type="Not"
        wbid=""
        mygene =input$geneinput
        
        if(as.character(mygene) %in% rownames(alias)){
            mygene=alias[as.character(mygene),2]
        }
        
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
                
                selectInput("isoform", label = HTML("<b>Targeted isoform
                                                           [<a href=\"\" onclick=\"$('#explain_isoform').toggle(); return false;\">info</a>]
                                                           </b>"), 
                            unique(Genes[which(as.character(Genes[,4])==wbid),5])),
                
                HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_isoform\">
            Isoforms based on the WormBase 270 annotations</div></p>
                     "),
                selectInput("selectMM", label = HTML("<b>Synthetic piRNA off-target sites
                                                           [<a href=\"\" onclick=\"$('#explain_uniqueness').toggle(); return false;\">info</a>]
                                                           <br>(minimum # of mismatches)</b>"), 
                            choices = list(">= 5 mismatches to genome" = 1, ">= 4 mismatches to genome" = 2, ">= 3 mismatches to genome" = 3,
                                           ">= 5 mismatches to exome" = 4, ">= 4 mismatches to exome" = 5, ">= 3 mismatches to exome" = 6),
                            selected = 1),
                HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_uniqueness\">
            Select the specificity of piRNAi fragments. Each piRNAi fragment was mapped across the C. elegans genome and verified its uniqueness at least to n mismatches, i.e. even with n base pairs changed the sequence remains unique.
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
            checkboxInput("FlaControl", label = HTML("<b>Generate negative control
                                                               [<a href=\"\" onclick=\"$('#explain_control').toggle(); return false;\">info</a>]
                                                               </b>"), value = FALSE, width='100%'),
            HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_control\">
            Reverse complement sequences so they have the same orientation as the coding sequence
                                     </div></p>"),
            actionButton("actionPI", label = "Generate piRNAi cluster")
            
            )
        })
        }
    }, ignoreInit = T)
    
    
    ##Main search function
    observeEvent(input$actionPI, {
        ErrorFlag = 0
        output$ErrorMessage <- renderText({})
        
        matches=as.integer(input$selectMM)
        mm=c(5,4,3,5,4,3)[matches]
        isoform = input$isoform
        ControlEx = input$FlaControl
        
        wbid = as.character(unique(Genes[which(isoform==Genes[,5]),4]))
        loc = as.character(unique(Genes[which(isoform==Genes[,5]),6]))
        
        strand = as.character(unique(Genes[which(isoform==Genes[,5]),7]))
        genest = as.numeric(unique(Genes[which(isoform==Genes[,5]),2]))
        geneend = as.numeric(unique(Genes[which(isoform==Genes[,5]),3]))
        
        file=paste(c("DataBase/",as.character(wbid),"_",as.character(isoform),"_",as.character(loc),".txt"),sep="",collapse="")
        
        tab=read.table(file,sep="\t",header=F)
        tab[,5]=as.integer((unlist(strsplit(as.character(tab[,2]),";")))[c(FALSE,TRUE)])
        tab[,2]=as.character((unlist(strsplit(as.character(tab[,2]),";")))[c(TRUE,FALSE)])
        
        if(matches>3){
        tab[,3]=tab[,4]
        }
        
        Seltab=tab[which(tab[,3]>=mm),]
        
        ##Remove Inside Introns
        if(as.character(isoform) %in% rownames(Intr)){
            ItrS=as.numeric(unlist(strsplit(as.character(Intr[as.character(isoform),2]),",")))
            ItrE=as.numeric(unlist(strsplit(as.character(Intr[as.character(isoform),3]),",")))
            idx= IRanges(Seltab[,1], Seltab[,1] + 19) %over% IRanges(ItrS,ItrE)
            if(sum(idx) > 0){Seltab=Seltab[-which(idx),]}
            }
        
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
                    paste("Error: The program selected overlapping piRNAi sites. Try to change the parameters or use Advanced function")
                })
                ErrorFlag=1
                
                }
        }
        
        ##Error for at least 6 piRNAs
        if(length(idx) < 6){
            
            output$ErrorMessage <- renderText({
                paste("Error: Not enough piRNAi sites to create the cluster. Try to change the parameters or use advanced function")
            })
            ErrorFlag=1
            
        }
        
        #Remove if errors()
        if(ErrorFlag == 1){
            output$SelPiTabSummary <- renderUI({ HTML(paste0("<b>Try again!</b>",sep=""))})
            output$SelPiTab=renderTable({})
            output$downloadseq <- renderUI({})
            output$SimpleFragment <- renderText({})
            
            }
        ##Produce outputs
            if(ErrorFlag == 0){
                
                ##Table results
                Pitab=Seltab[idx,c(1,2,5)]
                
                #Render piRNA coordinates 1-based
                Pitab[,1] = Pitab[,1] + 1
                
                ##Get coords
                if(strand =="-"){
                    Pistrt=Pitab[,1]+19
                    Piedt=Pitab[,1]
                }else{
                    Pistrt=Pitab[,1]
                    Piedt=Pitab[,1]+19
                }
                
                
                ##Check if introns data
                if(as.character(isoform) %in% rownames(Intr)){
                    ItrS=as.numeric(unlist(strsplit(as.character(Intr[as.character(isoform),2]),",")))
                    ItrE=as.numeric(unlist(strsplit(as.character(Intr[as.character(isoform),3]),",")))
                }else{
                    ItrS=c(1)
                    ItrE=c(0)
                    }
                Ncoors=ConvCooTr2cDNA(Pistrt - genest +1 ,Piedt - genest +1 ,ItrS - genest +1 ,ItrE - genest+1)
                lengcdna=(geneend-genest+1 - sum(ItrE-ItrS + 1))
                if(strand =="-"){
                    Ncoors[,2]= lengcdna - Ncoors[,2] + 1
                    Ncoors[,1]= lengcdna - Ncoors[,1] + 1
                    }
                Pitab=cbind(paste(as.integer(Ncoors[,1]),"to", as.integer(Ncoors[,2])),Pitab[,c(2,3)])
                
                colnames(Pitab)=c("Location","Sequence (antisense to target)","%GC")
                colnames(Pitab)[1]=paste("cDNA location ","(",lengcdna,"bp long)",sep="")
                Pitab=Pitab[order(Ncoors[,1]),]
                
                output$SelPiTabSummary <- renderUI({ HTML(paste0("<b>Synthetic piRNAs</b>",sep=""))})
                output$SelPiTab=renderTable(Pitab)
                
                ##Ape output
                output$downloadseq <- renderUI({
                    ##If control experiment, invert sequences
                    if(ControlEx){
                        Seltab[idx[1],2]=as.character(reverseComplement(DNAString(as.character(Seltab[idx[1],2]))))
                        Seltab[idx[2],2]=as.character(reverseComplement(DNAString(as.character(Seltab[idx[2],2]))))
                        Seltab[idx[3],2]=as.character(reverseComplement(DNAString(as.character(Seltab[idx[3],2]))))
                        Seltab[idx[4],2]=as.character(reverseComplement(DNAString(as.character(Seltab[idx[4],2]))))
                        Seltab[idx[5],2]=as.character(reverseComplement(DNAString(as.character(Seltab[idx[5],2]))))
                        Seltab[idx[6],2]=as.character(reverseComplement(DNAString(as.character(Seltab[idx[6],2]))))
                        }
                    
                    uno="cgcgcttgacgcgctagtcaactaacataaaaaaggtgaaacattgcgaggatacatagaaaaaacaatacttcgaattcatttttcaattacaaatcctgaaatgtttcactgtgttcctataagaaaacattgaaacaaaatattAagT"
                    uno=tolower(uno)
                    seq1=as.character(Seltab[idx[1],2])
                    dos="ctaattttgattttgattttgaaatcgaatttgcaaatccaattaaaaatcattttctgataattagacagttccttatcgttaattttattatatctatcgagttagaaattgcaacgaagataatgtcttccaaatactgaaaatttgaaaatatgtt"
                    dos=tolower(dos)
                    seq2=as.character(reverseComplement(DNAString(as.character(Seltab[idx[2],2]))))
                    tres="AttGccagaactcaaaatatgaaatttttatagttttgttgaaacagtaagaaaatcttgtaattactgtaaactgtttgctttttttaaagtcaacctacttcaaatctacttcaaaaattataatgtttcaaattacataactgtgt"
                    tres=tolower(tres)
                    seq3= as.character(reverseComplement(DNAString(as.character(Seltab[idx[3],2]))))
                    cuatro="ActgtagagcttcaatgttgataagatttattaacacagtgaaacaggtaatagttgtttgttgcaaaatcggaaatctctacatttcatatggtttttaattacaggtttgttttataaaataattgtgtgatggatattattttcagacctcatactaatctgcaaaccttcaaacaatatgtgaagtctactctgtttcactcaaccattcatttcaatttggaaaaaaatcaaagaaatgttgaaaaattttcctgtttcaacattatgacaaaaatgttatgattttaataaaaaCaaT"
                    cuatro=tolower(cuatro)
                    seq4=as.character(Seltab[idx[4],2])
                    cinco="ttctgtttttcttagaagtgttttccggaaacgcgtaattggttttatcacaaatcgaaaacaaacaaaaatttttttaattatttctttgctagttttgtagttgaaaattcactataatcatgaataagtgagctgcccaagtaaacaaagaaaatttggcagcggccgacaactaccgggttgcccgatttatcagtggagga"
                    cinco=tolower(cinco)
                    seq5= as.character(reverseComplement(DNAString( as.character(Seltab[idx[5],2]))))
                    seis="AtcTaatgtgatgtacacggttttcatttaaaaacaaattgaaacagaaatgactacattttcaaattgtctatttttgctgtgtttattttgccaccaacaaT"
                    seis=tolower(seis)
                    seq6=as.character(Seltab[idx[6],2])
                    siete="tcaatctagtaaactcacttaatgcaattcctccagccacatatgtaaacgttgtatacatgcagaaaacggttttttggttttaatgggaacttttgacaaattgttcgaaaatcttaagctgtcccatttcagttgggtgatcgattt"
                    siete=tolower(siete)

                    xtracom=paste("Recoded 21ur-1224 locus. Parameters: Gene = ",wbid,"; Isoform = ",isoform,"; At least ",mm," mismatches")
                    binrev=c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
                    if(ControlEx){
                       xtracom = append(xtracom,"Control experiment: piRNAi Sequences have been reverse complemented. THIS FRAGMENT WONT SILENCE THE SELECTED GENE")
                       binrev=!binrev
                    }

                    Compseq=paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete),sep="",collapse="")
                    
                    toadd="WorkingSpace/Piconst2.txt"
                    
                    pats= c(seq1,seq2,seq3,seq4,seq5,seq6)
                    fwdc= c(rep("#00ff00",6))
                    revc= c(rep("#ff0000",6))
                    tooltis= paste("piRNA",1:6)

                    writeLines(PasteApe(paste(wbid,"_21ur_1224_",sep="",collapse=""),Compseq,pats,fwdc,revc,tooltis,xtracom,toadd,binrev,"C.elegans"),paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""))
                    
                    output$SimpleFragment <- renderText({
                        paste("Recoded 21ur-1224 piRNA cluster\n",paste(Compseq,sep="",collapse=""),sep="",collapse="")
                    })
                    
                    downloadButton('DownApeOut', 'Download annotated genbank file')
                })
                
                
            }
        
    }, ignoreInit = T)
    
    ###Advanced construct######
    ####Dirty coded addition of extra clusters
    observeEvent(input$actionconstruct, {
        AdvancedErrorFlag=0
        output$AdvancedErrorMessage <- renderText({})
        clust=as.character(input$clustercon)

        if(clust == "1"){

        pipi1=as.character(input$piRNAseq1_1)
        pipi2=as.character(input$piRNAseq2_1)
        pipi3=as.character(input$piRNAseq3_1)
        pipi4=as.character(input$piRNAseq4_1)
        pipi5=as.character(input$piRNAseq5_1)
        pipi6=as.character(input$piRNAseq6_1)
        
        ##Check for size
        if((nchar(pipi1) == 0) | (nchar(pipi2) == 0) | (nchar(pipi3) == 0) | (nchar(pipi4) == 0) | (nchar(pipi5) == 0) | (nchar(pipi6) == 0)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Please pick at least 6 synthetic piRNAs")
            })
        }
        
        if(AdvancedErrorFlag == 0){
        if((nchar(pipi1) != 20) | (nchar(pipi2) != 20) | (nchar(pipi3) != 20) | (nchar(pipi4) != 20) | (nchar(pipi5) != 20) | (nchar(pipi6) != 20)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Error: piRNAi sequences should be 20bp long")
            })
            }}
        
        #Check for input characters
        toto=paste(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,sep="",collapse="")
        if((AdvancedErrorFlag == 0) & (nchar(gsub("A|T|C|G","",toupper(toto))) != 0)){ ##Check for strange non ATCG characters
            output$AdvancedErrorMessage <- renderText({
                paste("Error: Unrecognized characters in piRNAi sequences")
            })
            AdvancedErrorFlag=1
        }
        
        #Checkfor GC
        if(AdvancedErrorFlag == 0){
        Gcvals=sapply(c(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6),CalculateGC)
        if((sum(Gcvals <.3)+sum(Gcvals >.7))>0){
            output$AdvancedErrorMessage <- renderText({
                paste("Warning: some sequences have a GC content lower to 30% or higher to 70%")
            })
        }
        }
        
        ##Main Routine
        if(AdvancedErrorFlag == 0){
            output$downloadconstruct <- renderUI({
            
                uno="cgcgcttgacgcgctagtcaactaacataaaaaaggtgaaacattgcgaggatacatagaaaaaacaatacttcgaattcatttttcaattacaaatcctgaaatgtttcactgtgttcctataagaaaacattgaaacaaaatattAagT"
                uno=tolower(uno)
                seq1=as.character(pipi1)
                dos="ctaattttgattttgattttgaaatcgaatttgcaaatccaattaaaaatcattttctgataattagacagttccttatcgttaattttattatatctatcgagttagaaattgcaacgaagataatgtcttccaaatactgaaaatttgaaaatatgtt"
                dos=tolower(dos)
                seq2=as.character(reverseComplement(DNAString(as.character(pipi2))))
                tres="AttGccagaactcaaaatatgaaatttttatagttttgttgaaacagtaagaaaatcttgtaattactgtaaactgtttgctttttttaaagtcaacctacttcaaatctacttcaaaaattataatgtttcaaattacataactgtgt"
                tres=tolower(tres)
                seq3= as.character(reverseComplement(DNAString(as.character(pipi3))))
                cuatro="ActgtagagcttcaatgttgataagatttattaacacagtgaaacaggtaatagttgtttgttgcaaaatcggaaatctctacatttcatatggtttttaattacaggtttgttttataaaataattgtgtgatggatattattttcagacctcatactaatctgcaaaccttcaaacaatatgtgaagtctactctgtttcactcaaccattcatttcaatttggaaaaaaatcaaagaaatgttgaaaaattttcctgtttcaacattatgacaaaaatgttatgattttaataaaaaCaaT"
                cuatro=tolower(cuatro)
                seq4=as.character(pipi4)
                cinco="ttctgtttttcttagaagtgttttccggaaacgcgtaattggttttatcacaaatcgaaaacaaacaaaaatttttttaattatttctttgctagttttgtagttgaaaattcactataatcatgaataagtgagctgcccaagtaaacaaagaaaatttggcagcggccgacaactaccgggttgcccgatttatcagtggagga"
                cinco=tolower(cinco)
                seq5= as.character(reverseComplement(DNAString( as.character(pipi5))))
                seis="AtcTaatgtgatgtacacggttttcatttaaaaacaaattgaaacagaaatgactacattttcaaattgtctatttttgctgtgtttattttgccaccaacaaT"
                seis=tolower(seis)
                seq6=as.character(pipi6)
                siete="tcaatctagtaaactcacttaatgcaattcctccagccacatatgtaaacgttgtatacatgcagaaaacggttttttggttttaatgggaacttttgacaaattgttcgaaaatcttaagctgtcccatttcagttgggtgatcgattt"
                siete=tolower(siete)
                
                xtracom=paste("Recoded 21ur-1224 locus via Advanced Search.")
                binrev=c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
                
                Compseq=paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete),sep="",collapse="")
                
                toadd="WorkingSpace/Piconst2.txt"
                
                pats= c(seq1,seq2,seq3,seq4,seq5,seq6)
                fwdc= c(rep("#00ff00",6))
                revc= c(rep("#ff0000",6))
                tooltis= paste("piRNA",1:6)
                
                writeLines(PasteApe("21ur_1224_",Compseq,pats,fwdc,revc,tooltis,xtracom,toadd,binrev,"C.elegans"),paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""))
                
                output$AdvancedFragment <- renderText({
                    paste("Recoded 21ur-1224 piRNA cluster\n",paste(Compseq,sep="",collapse=""),sep="",collapse="")
                })
                
            downloadButton('DownConOut', 'Download annotated genbank file')
            
                })
            
            }
        }

        ##Second cluster
        
        if(clust == "2"){

        pipi1=as.character(input$piRNAseq1_2)
        pipi2=as.character(input$piRNAseq2_2)
        pipi3=as.character(input$piRNAseq3_2)
        pipi4=as.character(input$piRNAseq4_2)
        pipi5=as.character(input$piRNAseq5_2)
        pipi6=as.character(input$piRNAseq6_2)
        
        ##Check for size
        if((nchar(pipi1) == 0) | (nchar(pipi2) == 0) | (nchar(pipi3) == 0) | (nchar(pipi4) == 0) | (nchar(pipi5) == 0) | (nchar(pipi6) == 0)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Please pick at least 6 synthetic piRNAs")
            })
        }
        
        if(AdvancedErrorFlag == 0){
        if((nchar(pipi1) != 20) | (nchar(pipi2) != 20) | (nchar(pipi3) != 20) | (nchar(pipi4) != 20) | (nchar(pipi5) != 20) | (nchar(pipi6) != 20)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Error: piRNAi sequences should be 20bp long")
            })
            }}
        
        #Check for input characters
        toto=paste(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,sep="",collapse="")
        if((AdvancedErrorFlag == 0) & (nchar(gsub("A|T|C|G","",toupper(toto))) != 0)){ ##Check for strange non ATCG characters
            output$AdvancedErrorMessage <- renderText({
                paste("Error: Unrecognized characters in piRNAi sequences")
            })
            AdvancedErrorFlag=1
        }
        
        #Checkfor GC
        if(AdvancedErrorFlag == 0){
        Gcvals=sapply(c(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6),CalculateGC)
        if((sum(Gcvals <.3)+sum(Gcvals >.7))>0){
            output$AdvancedErrorMessage <- renderText({
                paste("Warning: some sequences have a GC content lower to 30% or higher to 70%")
            })
        }
        }
        
        ##Main Routine
        if(AdvancedErrorFlag == 0){
            output$downloadconstruct <- renderUI({
            
                uno="ttcgtggtgcacttatctttctccttcaaattgaaaactcagtttttaattatagtcaaatctcttttgctgacaggtccaaagtactttattatttcatattatataaaattcattctcgaatttatttataaattttcgctgagtcaa"
                uno=tolower(uno)
                seq1=as.character(reverseComplement(DNAString(as.character(pipi1))))
                dos="ActgactaacaaaaacccctgtcaatttacttgtaatgtgaaactgtatcggtttcatattatctatgattcgagtacattgtttcaaattcaaT"
                dos=tolower(dos)
                seq2=as.character(pipi2)
                tres="attttacgctggtttgaaaatttgaaatattccaaaataaattatttagttttcgttttttgtacattgtcataaaacattttggttttttttaaca"
                tres=tolower(tres)
                seq3= as.character(reverseComplement(DNAString(as.character(pipi3))))
                cuatro="Acttaaaatcaaaaattgttacactttataacagttcattgaaactgaaaattattttcttttcccaaataataataccatcaaatgtcgtggtgtactcatcttttccttttcttcttttttttcaatttctccttcaaatctctacacactcttcactgccaatctttttttctttccttatccaaaagcacacttttgtgcagagtaaataatgcactttgtgaaaaaaaaactatttttaaaactgtatttttttaagtttggcaatttttgagaaaatttcaacaaaatctgatatagattggaatttaaatggttcaaatttg"
                cuatro=tolower(cuatro)
                seq4=as.character(reverseComplement(DNAString(as.character(pipi4))))
                cinco="AtctattcaaagttttattcgaagtttttaacagacacttgaaacagtgtaataattttctgacaaaaattaaaacaaatgttactactttgcttttcttactttatccgttttttatcacccttatttttcagtcaaccctagcaacgttaccgacggaatcggtaggactacaccgactgcatcaaatttgggaagaagccgtgagaatttgagtttcaatcaacaccgcccagaccatccttcgtcatattttgatagtttggagcatggtgagcattttatattaaaacagttgttttggtgttcatattactaatgtctgaatactaacttgcattaaaattggaaattaaaaaaattactgtttctcaaaagtattttcaatacctatatttttttgctacagT"
                cinco=tolower(cinco)
                seq5= as.character(pipi5)
                seis="caatattttcaaatattttataccagatttttcgaaaaagttgaattttcaattaacaataacgcatttatgcatttttcactcttttttgagatttaatgctgaaaaaatagttctgaaaatgacaaaagttatgttttcaatattttttatcaaactaaatttatttaatttgttaactgttgcttttttgtttttcttcaagt"
                seis=tolower(seis)
                seq6=as.character(reverseComplement(DNAString(as.character(pipi6))))
                siete="Atcttcgaagcaacttatttgatgttttataaacgacctgaaacatactggtgatgcccaataatgttttttttaaatttagtctcgtgaaaaaaataaaattaaaacagaaaattacatttgcgccgaagaaacttaagatctggaactt"
                siete=tolower(siete)
                
                xtracom=paste("Recoded 21ur-1692 locus via Advanced Search.")
                binrev=c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE)
                
                Compseq=paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete),sep="",collapse="")
                
                toadd=""
                
                pats= c(seq1,seq2,seq3,seq4,seq5,seq6)
                fwdc= c(rep("#00ff00",6))
                revc= c(rep("#ff0000",6))
                tooltis= paste("piRNA",1:6)
                
                writeLines(PasteApe("_21ur-1692_",Compseq,pats,fwdc,revc,tooltis,xtracom,toadd,binrev,"C.elegans"),paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""))
                
                output$AdvancedFragment <- renderText({
                    paste("Recoded 21ur-1692 piRNA cluster\n",paste(Compseq,sep="",collapse=""),sep="",collapse="")
                })
                
            downloadButton('DownConOut', 'Download annotated genbank file')
            
                })
            
            }
        }

        ###Third cluster

        if(clust == "3"){

        pipi1=as.character(input$piRNAseq1_3)
        pipi2=as.character(input$piRNAseq2_3)
        pipi3=as.character(input$piRNAseq3_3)
        pipi4=as.character(input$piRNAseq4_3)
        pipi5=as.character(input$piRNAseq5_3)
        pipi6=as.character(input$piRNAseq6_3)
        pipi7=as.character(input$piRNAseq7_3)
        
        ##Check for size
        if((nchar(pipi1) == 0) | (nchar(pipi2) == 0) | (nchar(pipi3) == 0) | (nchar(pipi4) == 0) | (nchar(pipi5) == 0) | (nchar(pipi6) == 0)| (nchar(pipi7) == 0)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Please pick at least 7 synthetic piRNAs")
            })
        }
        
        if(AdvancedErrorFlag == 0){
        if((nchar(pipi1) != 20) | (nchar(pipi2) != 20) | (nchar(pipi3) != 20) | (nchar(pipi4) != 20) | (nchar(pipi5) != 20) | (nchar(pipi6) != 20)| (nchar(pipi7) != 20)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Error: piRNAi sequences should be 20bp long")
            })
            }}
        
        #Check for input characters
        toto=paste(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,pipi7,sep="",collapse="")
        if((AdvancedErrorFlag == 0) & (nchar(gsub("A|T|C|G","",toupper(toto))) != 0)){ ##Check for strange non ATCG characters
            output$AdvancedErrorMessage <- renderText({
                paste("Error: Unrecognized characters in piRNAi sequences")
            })
            AdvancedErrorFlag=1
        }
        
        #Checkfor GC
        if(AdvancedErrorFlag == 0){
        Gcvals=sapply(c(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,pipi7),CalculateGC)
        if((sum(Gcvals <.3)+sum(Gcvals >.7))>0){
            output$AdvancedErrorMessage <- renderText({
                paste("Warning: some sequences have a GC content lower to 30% or higher to 70%")
            })
        }
        }
        
        ##Main Routine
        if(AdvancedErrorFlag == 0){
            output$downloadconstruct <- renderUI({
            
                uno="agggtggtcgcatagaagttgggcgcacttcatttacaaaaatatgttcaaattttgtgatttcatgttcaggcattttgttttgatgataaaacatagtgtgactgtttttatatgtttataaaatgtcttatgattaaacaagatcaaT"
                uno=tolower(uno)
                seq1=as.character(pipi1)
                dos="aaatgtcggtttatttcacactgataggaatttttcaaaaaaatatatcagcaaagtactgtattaaaatgtgaaaatctcataaaaagtttaagtttca"
                dos=tolower(dos)
                seq2=as.character(reverseComplement(DNAString(as.character(pipi2))))
                tres="ATtgtggcaaatagattttgacaatttttatcaaattcatgaaacagtagaattttttccaaaaaactcacaaaataaatacgaatttcaatttgcccactttatcaaataaatgtttacacaaaagtaggccgtgcaacgcgcctatcctagatgctacattccttggttttgagttgtgaaacgttggaataatgtacttcattttgtgacttacttttttgtaaggcaattgttttttttatttaataaaagtactttcctaaattcaaatatcaaatttgtcttcatttttgtgacaagtaaaa"
                tres=tolower(tres)
                seq3= as.character(reverseComplement(DNAString(as.character(pipi3))))
                cuatro="Agcttttttagaaaaaaaaagtcagtatttaagaacaatttgaaacagttcaatttttcagtgtaatttccaaccagaactttttgagtaaataattacagaaaacttatttagaaaataggactaaataatgcaaatattttccggactggcatttataaataagcaagtataag"
                cuatro=tolower(cuatro)
                seq4=as.character(reverseComplement(DNAString(as.character(pipi4))))
                cinco="AttgaggtaattttaaaaagcatacatataagcaagtcgtgaaacagtcgtttaaattttatttttcaaaaagttataacgcgacagcagtttcatctgtttcatattccctatttgttgaaatttgagacgtattttacgaT"
                cinco=tolower(cinco)
                seq5= as.character(pipi5)
                seis="tgccatccgaatcttgaactttgtatcaattgttcacatttttttccaaaaacgtattaactcactttca"
                seis=tolower(seis)
                seq6=as.character(reverseComplement(DNAString(as.character(pipi6))))
                siete="AtcttttgtttttaacaaagagatcatatatactcattgagaaacagtacattttttgaaagtacatttgccttggtcaaatatataaagttgacaaaagtttaaaaatgtttccaaaagttaattaaaaaatcaaatttattctagctcaaT"
                siete=tolower(siete)
                seq7=as.character(pipi7)
                ocho="cgtcaagtgatcaaaccatcatttttttcagattaagacctgatttgtcagtgaattgaaaaaaacgtgttcattgcgtgtttcgcattttttatatataaaaaagcaagtttcggcggcaataacgaagtattcccaacagatcaatag"
                ocho=tolower(ocho)
                
                xtracom=paste("Recoded 21ur-8831 locus via Advanced Search.")
                binrev=c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE)
                
                Compseq=paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete,seq7,ocho),sep="",collapse="")
                
                toadd=""
                
                pats= c(seq1,seq2,seq3,seq4,seq5,seq6,seq7)
                fwdc= c(rep("#00ff00",7))
                revc= c(rep("#ff0000",7))
                tooltis= paste("piRNA",1:7)
                
                writeLines(PasteApe("21ur-8831_",Compseq,pats,fwdc,revc,tooltis,xtracom,toadd,binrev,"C.elegans"),paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""))
                
                output$AdvancedFragment <- renderText({
                    paste("Recoded 21ur-8831 piRNA cluster\n",paste(Compseq,sep="",collapse=""),sep="",collapse="")
                })
                
            downloadButton('DownConOut', 'Download annotated genbank file')
            
                })
            
            }
        }

        ###Fourth
        
        if(clust == "4"){

        pipi1=as.character(input$piRNAseq1_4)
        pipi2=as.character(input$piRNAseq2_4)
        pipi3=as.character(input$piRNAseq3_4)
        pipi4=as.character(input$piRNAseq4_4)
        pipi5=as.character(input$piRNAseq5_4)
        pipi6=as.character(input$piRNAseq6_4)
        pipi7=as.character(input$piRNAseq7_4)
        pipi8=as.character(input$piRNAseq8_4)

        ##Check for size
        if((nchar(pipi1) == 0) | (nchar(pipi2) == 0) | (nchar(pipi3) == 0) | (nchar(pipi4) == 0) | (nchar(pipi5) == 0) | (nchar(pipi6) == 0)| (nchar(pipi7) == 0)| (nchar(pipi8) == 0)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Please pick at least 8 synthetic piRNAs")
            })
        }
        
        if(AdvancedErrorFlag == 0){
        if((nchar(pipi1) != 20) | (nchar(pipi2) != 20) | (nchar(pipi3) != 20) | (nchar(pipi4) != 20) | (nchar(pipi5) != 20) | (nchar(pipi6) != 20)| (nchar(pipi7) != 20)| (nchar(pipi8) != 20)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Error: piRNAi sequences should be 20bp long")
            })
            }}
        
        #Check for input characters
        toto=paste(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,pipi7,pipi8,sep="",collapse="")
        if((AdvancedErrorFlag == 0) & (nchar(gsub("A|T|C|G","",toupper(toto))) != 0)){ ##Check for strange non ATCG characters
            output$AdvancedErrorMessage <- renderText({
                paste("Error: Unrecognized characters in piRNAi sequences")
            })
            AdvancedErrorFlag=1
        }
        
        #Checkfor GC
        if(AdvancedErrorFlag == 0){
        Gcvals=sapply(c(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,pipi7,pipi8),CalculateGC)
        if((sum(Gcvals <.3)+sum(Gcvals >.7))>0){
            output$AdvancedErrorMessage <- renderText({
                paste("Warning: some sequences have a GC content lower to 30% or higher to 70%")
            })
        }
        }
        
        ##Main Routine
        if(AdvancedErrorFlag == 0){
            output$downloadconstruct <- renderUI({
            
                uno="caaaaaacaatacgtcccttatcttctggaatcagctcattgtgctcatcggagctatccgcaccgtcaactatactcgctagatcttccgtgttctgatcttgagtgtatagtggaggggggtcaacctgaaatttcagatttttgttg"
                uno=tolower(uno)
                seq1=as.character(reverseComplement(DNAString(as.character(pipi1))))
                dos="ActgttttagaagtgatgagtcttattataataacttgttgaaactgtggatttatattttttaaaaattaccggcgaaattgattcataatctcttattaccatagttaaagtctctagaataagcacaaaactactaaagtttgtaaaataattgaatatgccacaactgataagagactttttcctcttatcagcataaagtccaaagcgataaaattcaaaagagacaagtacaaatgtatattaatctgctttgttggaaaaaaattaaactttttatctaaacctgtcattgatccaaaagattaagtttcctgcaaaattgtttcgaaatattattgtgattgaaacttttgactttttcaacttatcaataagtcattggcttaagataaagtaatcaaT"
                dos=tolower(dos)
                seq2=as.character(pipi2)
                tres="cgcgctcagcactcaatttctgcccaaatagtt"
                tres=tolower(tres)
                seq3= as.character(reverseComplement(DNAString(as.character(pipi3))))
                cuatro="Attgagtatctaaatgaaaacctaaaatatgaacagttagaaacaggaaatttttgaaaagttaaaaaacaacctatacaattaatttccaagaaaaatttaacaatcgattttcatttctgaaatcccaaaatcggtgaattcttgatgaaaatgcatttgaaaatacaattttgtttta"
                cuatro=tolower(cuatro)
                seq4=as.character(reverseComplement(DNAString(as.character(pipi4))))
                cinco="Atctaattagatatgcaagcctaatatttgtatcattcttgaaactgtaaataaaaaatgtttgcaaaaaaaatcaattttttagcgaatgttaacataaaaccttaaatttttctgggttttgaccgtttctcatatttcaaaa"
                cinco=tolower(cinco)
                seq5= as.character(reverseComplement(DNAString(as.character(pipi5))))
                seis="AtcgaT"
                seis=tolower(seis)
                seq6=as.character(pipi6)
                siete="aaaaaatatgctgaaacgtgattgctttttgtgcttttttatacaagtttgcaatcgcacaaatcatatgaaaaattattaagcacgcttaaactatgtgatctgaaatacgaaaactagtatacgttaaacaggaaaaaaaaatcaactgtttcaaaatttgtgtttaatcaaattaaagttgctattccgaT"
                siete=tolower(siete)
                seq7=as.character(pipi7)
                ocho="taatagttaattttcaaaatagaaagttttaaaacatcctgtttcggtatgctgatttttacagactccactttgtagttaacT"
                ocho=tolower(ocho)
                seq8=as.character(pipi8)
                nueve="ggacagaaatattgatattttgccagttaccaggaaaataaattattctttgcaacatctgactttaagaataaaaactcacaaattccttttccatttctgaaatattttagtgtcctcttcccgcaaccactccctgtaaatcgaaaa"
                nueve=tolower(nueve)
                
                xtracom=paste("Recoded 21ur-1610 locus via Advanced Search.")
                binrev=c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
                
                Compseq=paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete,seq7,ocho,seq8,nueve),sep="",collapse="")
                
                toadd=""
                
                pats= c(seq1,seq2,seq3,seq4,seq5,seq6,seq7,seq8)
                fwdc= c(rep("#00ff00",8))
                revc= c(rep("#ff0000",8))
                tooltis= paste("piRNA",1:8)
                
                writeLines(PasteApe("21ur-1610_",Compseq,pats,fwdc,revc,tooltis,xtracom,toadd,binrev,"C.elegans"),paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""))
                
                output$AdvancedFragment <- renderText({
                    paste("Recoded 21ur-1610 piRNA cluster\n",paste(Compseq,sep="",collapse=""),sep="",collapse="")
                })
                
            downloadButton('DownConOut', 'Download annotated genbank file')
            
                })
            
            }
        }

        }, ignoreInit = T)
    
    ###Advanced searchform
    ###Only upon retrieval of a gene in database new form will appear
    observeEvent(input$actionAdvsearch, {
        
        type="Not"
        wbid=""
        mygene =input$Advancedgeneinput
        
        if(as.character(mygene) %in% rownames(alias)){
            mygene=alias[as.character(mygene),2]
        }
        
        if(as.character(mygene) %in% as.character(Genes[,4])){type=4}
        if(as.character(mygene) %in% as.character(Genes[,5])){type=5}
        if(as.character(mygene) %in% locus){type=6}
        
        if(type == "Not"){
            output$AdvDesignControls <- renderUI({
                HTML("<b>Gene not found</b>")
            })
        }else{
            wbid=as.character(unique(Genes[which(mygene==Genes[,type]),4]))
            
            ##Control for table
            output$AdvDesignControls <- renderUI({
                fluidRow(
                    
                    column(width = 3,selectInput("AdvIsoform", label = HTML("<b>Targeted isoform
                                                           [<a href=\"\" onclick=\"$('#explain_isoform_advanced').toggle(); return false;\">info</a>]
                                                           </b>"), 
                                unique(Genes[which(as.character(Genes[,4])==wbid),5]))),
                    
                    column(width = 3, selectInput("AdvSelectMM", label = HTML("<b>Synthetic piRNA off-target sites
                    [<a href=\"\" onclick=\"$('#explain_uniqueness_advanced').toggle(); return false;\">info</a>]<br>(minimum # of mismatches)
                                                           </b>"),
                                                  choices = list(">= 5 mismatches to genome" = 1, ">= 4 mismatches to genome" = 2, ">= 3 mismatches to genome" = 3,
                                                                 ">= 5 mismatches to exome" = 4, ">= 4 mismatches to exome" = 5, ">= 3 mismatches to exome" = 6),
                                                                      selected = 1)),
                    column(width = 3, sliderInput("Posslider", label = HTML("<b>Relative position in Genebody (%)
                                                                            [<a href=\"\" onclick=\"$('#explain_Posgene').toggle(); return false;\">info</a>]
                                                                            </b>
                                                                            "),
                                                  0, 100, c(0, 100), step = 5)),
                    column(width = 3, sliderInput("Gcslider", label = HTML("<b>GC content (%)
                                                                            [<a href=\"\" onclick=\"$('#explain_GCcont').toggle(); return false;\">info</a>]
                                                                            </b>
                                                                           ")
                                                  ,0, 100, c(30, 70), step = 5)),
                    
                    HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_isoform_advanced\">
            Isoforms based on the WormBase 270 annotations</div></p>
                     "),

                    HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_uniqueness_advanced\">
            Select the specificity of piRNAi fragments. Each piRNAi fragment was mapped across the C. elegans genome and verified its uniqueness at least to n mismatches, i.e. even with n base pairs changed the sequence remains unique.
                                                 </div></p>
                     "),
                    
                    HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_Posgene\">
            Select the preferential place for piRNAi targeting. Ordering is from the 5' to the 3' of the selected isoform. 
                                                 </div></p>
                     "),
                    
                    HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_GCcont\">
            We recommend piRNAi fragments with GC content between 30 to 70%.
                                                 </div></p>
                     ")
                    )
            })
            
            }
        
        
        }, ignoreInit = T)
    
    ###Observe function for adding fragments based on piRNA table
    ##############################################################
    
    observeEvent(input$AdvIsoform,{
        ##Now design table; NOt sure if it will work as table should be most of the time be out of observe functions
        ADVisoform = input$AdvIsoform
        
        
        wbid = as.character(unique(Genes[which(ADVisoform==Genes[,5]),4]))
        loc = as.character(unique(Genes[which(ADVisoform==Genes[,5]),6]))
        
        strand = as.character(unique(Genes[which(ADVisoform==Genes[,5]),7]))
        genest = as.numeric(unique(Genes[which(ADVisoform==Genes[,5]),2]))
        geneend = as.numeric(unique(Genes[which(ADVisoform==Genes[,5]),3]))
        
        file=paste(c("DataBase/",as.character(wbid),"_",as.character(ADVisoform),"_",as.character(loc),".txt"),sep="",collapse="")
        
        tab=read.table(file,sep="\t",header=F)
        tab[,5]=as.integer((unlist(strsplit(as.character(tab[,2]),";")))[c(FALSE,TRUE)])
        tab[,2]=as.character((unlist(strsplit(as.character(tab[,2]),";")))[c(TRUE,FALSE)])
        
        #if(strand =="-"){tab[,1]= tab[,1]+18}
        #tab[,1]=c(tab[,1]+2-genest)/(geneend-genest+1)
        
        #if(strand =="-"){tab[,1]= 1 - as.numeric(tab[,1])}
        
        ##Remove pis inside introns
        if(as.character(ADVisoform) %in% rownames(Intr)){
            ItrS=as.numeric(unlist(strsplit(as.character(Intr[as.character(ADVisoform),2]),",")))
            ItrE=as.numeric(unlist(strsplit(as.character(Intr[as.character(ADVisoform),3]),",")))
            idx= IRanges(tab[,1], tab[,1] + 19) %over% IRanges(ItrS,ItrE)
            if(sum(idx) > 0){tab=tab[-which(idx),]}
        }
        
        
        #output$AllPiTab <- DT::renderDataTable(DT::datatable({
        output$AllPiTab <- DT::renderDataTable({
            datatab = tab
            matches=as.integer(input$AdvSelectMM)
            mm=c(5,4,3,5,4,3)[matches]
            minGC=input$Gcslider[1]
            maxGC=input$Gcslider[2]
            minPos=input$Posslider[1]/100
            maxPos=input$Posslider[2]/100
            if(matches>3){
                datatab[,3]=datatab[,4]
                }
            
            datatab = datatab[which(datatab[,3]>=mm),]
            
            datatab = datatab[which((datatab[,5]>=minGC)&(datatab[,5]<=maxGC)),]
            
            relpos=datatab[,1]
            if(strand =="-"){relpos= relpos+18}
            relpos=c(relpos+2-genest)/(geneend-genest+1)
            if(strand =="-"){relpos= 1 - as.numeric(relpos)}
            
            datatab = datatab[which((relpos>=minPos)&(relpos<=maxPos)),]
            
            
            ##Table results
            Pitab=datatab
            
            #Render piRNA coordinates 1-based
            Pitab[,1] = Pitab[,1] + 1
            
            ##Get coords
            if(strand =="-"){
                Pistrt=Pitab[,1]+19
                Piedt=Pitab[,1]
            }else{
                Pistrt=Pitab[,1]
                Piedt=Pitab[,1]+19
            }
            
            
            ##Check if introns data
            if(as.character(ADVisoform) %in% rownames(Intr)){
                ItrS=as.numeric(unlist(strsplit(as.character(Intr[as.character(ADVisoform),2]),",")))
                ItrE=as.numeric(unlist(strsplit(as.character(Intr[as.character(ADVisoform),3]),",")))
            }else{
                ItrS=c(1)
                ItrE=c(0)
            }
            
            Ncoors=ConvCooTr2cDNA(Pistrt - genest +1 ,Piedt - genest +1 ,ItrS - genest +1 ,ItrE - genest+1)
            lengcdna=(geneend-genest+1 - sum(ItrE-ItrS + 1))
            if(strand =="-"){
                Ncoors[,2]= lengcdna - Ncoors[,2] + 1
                Ncoors[,1]= lengcdna - Ncoors[,1] + 1
            }
            
            datatab[,1]=Ncoors[,1]
            
            datatab=datatab[order(Ncoors[,1]),]
            

            Pdata=data.frame(
                Location = paste(as.integer(datatab[,1]), "to", as.integer(datatab[,1])+19),
                Sequence = datatab[,2],
                GCcontent = datatab[,5],
                Select = shinyInput(actionButton, as.character(datatab[,2]), 'button_', label = "Add to contruct", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' ),
                stringsAsFactors = FALSE,
                row.names = 1:nrow(datatab)
                )
            
            colnames(Pdata)[1]=paste("cDNA location ","(",(lengcdna),"bp long)",sep="")
            colnames(Pdata)[2]="Sequence (antisense to target)"
            rownames(Pdata)=1:nrow(Pdata)
            Pdata
        #},server = FALSE, escape = FALSE, selection = 'none'))
        },server = FALSE, escape = FALSE, selection = 'none')
        }, ignoreInit = F)
    
    ##Handle shiny to add dynamic button
    shinyInput <- function(FUN, seqs, id, ...) {
        inputs <- character(length(seqs))
        for (i in 1:length(seqs)) {
            inputs[i] <- as.character(FUN(paste0(id, seqs[i]), ...))
        }
        inputs
    }
    
    #Handle id-seq of dynamic button
    observeEvent(input$select_button, {
        clust=as.character(input$clustercon)
        fill=0
        selectedSeq <- as.character(strsplit(input$select_button, "_")[[1]][2])
        
        if(clust == "1"){
        if(fill == 0){
        if(as.character(input$piRNAseq1_1)==""){
            updateTextAreaInput(session, "piRNAseq1_1", value = selectedSeq)
            fill = 1
        }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq2_1)==""){
                updateTextAreaInput(session, "piRNAseq2_1", value = selectedSeq)
                fill = 1
            }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq3_1)==""){
                updateTextAreaInput(session, "piRNAseq3_1", value = selectedSeq)
                fill = 1
            }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq4_1)==""){
                updateTextAreaInput(session, "piRNAseq4_1", value = selectedSeq)
                fill = 1
            }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq5_1)==""){
                updateTextAreaInput(session, "piRNAseq5_1", value = selectedSeq)
                fill = 1
            }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq6_1)==""){
                updateTextAreaInput(session, "piRNAseq6_1", value = selectedSeq)
                fill = 1
            }}
        
        ##Send custom message
        if(fill == 0){
            showNotification("Construct has already 6 sequences.")
        }
        }
        
        ##Second
        if(clust == "2"){
            if(fill == 0){
                if(as.character(input$piRNAseq1_2)==""){
                    updateTextAreaInput(session, "piRNAseq1_2", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq2_2)==""){
                    updateTextAreaInput(session, "piRNAseq2_2", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq3_2)==""){
                    updateTextAreaInput(session, "piRNAseq3_2", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq4_2)==""){
                    updateTextAreaInput(session, "piRNAseq4_2", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq5_2)==""){
                    updateTextAreaInput(session, "piRNAseq5_2", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq6_2)==""){
                    updateTextAreaInput(session, "piRNAseq6_2", value = selectedSeq)
                    fill = 1
                }}
            
            ##Send custom message
            if(fill == 0){
                showNotification("Construct has already 6 sequences.")
            }
        }
        
        ##Third
        if(clust == "3"){
            if(fill == 0){
                if(as.character(input$piRNAseq1_3)==""){
                    updateTextAreaInput(session, "piRNAseq1_3", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq2_3)==""){
                    updateTextAreaInput(session, "piRNAseq2_3", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq3_3)==""){
                    updateTextAreaInput(session, "piRNAseq3_3", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq4_3)==""){
                    updateTextAreaInput(session, "piRNAseq4_3", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq5_3)==""){
                    updateTextAreaInput(session, "piRNAseq5_3", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq6_3)==""){
                    updateTextAreaInput(session, "piRNAseq6_3", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq7_3)==""){
                    updateTextAreaInput(session, "piRNAseq7_3", value = selectedSeq)
                    fill = 1
                }}
            
            ##Send custom message
            if(fill == 0){
                showNotification("Construct has already 7 sequences.")
            }
        }
        
        
        ###Four
        if(clust == "4"){
            if(fill == 0){
                if(as.character(input$piRNAseq1_4)==""){
                    updateTextAreaInput(session, "piRNAseq1_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq2_4)==""){
                    updateTextAreaInput(session, "piRNAseq2_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq3_4)==""){
                    updateTextAreaInput(session, "piRNAseq3_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq4_4)==""){
                    updateTextAreaInput(session, "piRNAseq4_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq5_4)==""){
                    updateTextAreaInput(session, "piRNAseq5_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq6_4)==""){
                    updateTextAreaInput(session, "piRNAseq6_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq7_4)==""){
                    updateTextAreaInput(session, "piRNAseq7_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq8_4)==""){
                    updateTextAreaInput(session, "piRNAseq8_4", value = selectedSeq)
                    fill = 1
                }}
            
            ##Send custom message
            if(fill == 0){
                showNotification("Construct has already 8 sequences.")
            }
        }
        
        
    })
    
    ####Other functions###########
    ##Convert coordinates
    ConvCooTr2cDNA = function(posS, posE, ExonS, ExonE){
        if(length(posS) != length(posE)){return(c())}
        if((length(posS)<1) | (length(posE)<1)){return(c())}
        if((length(ExonS)<1) | (length(ExonE)<1)){return(cbind(posS,posE))}
        dists=ExonE - ExonS + 1
        if(sum(dists)==0){return(cbind(posS, posE))}
        resS=c()
        resE=c()
        for(i in 1:length(posS)){
            idx=which(ExonS < posS[i])
            if(length(idx) > 0){
                resS=append(resS,posS[i] - sum(dists[idx]))
            }else{
                resS=append(resS,posS[i])
            }
            idx=which(ExonS < posE[i])
            if(length(idx) > 0){
                resE=append(resE,posE[i] - sum(dists[idx]))
            }else{
                resE=append(resE,posE[i])
            }
        }
        return(cbind(resS, resE))
    }
    
    ##############################
    
    ###Create ApeFIle as pasteLines####
    PasteApe = function(locus_name,sequence,patterns,FWDcolors,REVcolors,tooltips,xtraComments,xtraLines, BinRevComp, organism){
        if (is.null(sequence)){return(NULL)}
        if(!is.character(sequence)){return(c())}
        if(length(patterns) < 1 ){return(c(paste(sequence)))}
        if(length(patterns) != length(FWDcolors)){return(c())}
        if(length(REVcolors) != length(FWDcolors)){return(c())}
        if(length(tooltips) != length(FWDcolors)){return(c())}
        if(length(BinRevComp)==0){BinRevComp=rep(FALSE,length(patterns))}
        
        ##Save Lines
        FileLines=c()
        FileLines=append(FileLines,paste("LOCUS",paste(locus_name,sep="",collapse=""),paste(nchar(sequence),"bp ds-DNA", sep=""),"linear",paste(c(unlist(strsplit(date()," ")))[c(3,2,5)],sep="",collapse="-"),sep="     "))
        FileLines=append(FileLines,paste("DEFINITION",".",sep="     "))
        FileLines=append(FileLines,paste("ACCESSION",".",sep="     "))
        FileLines=append(FileLines,paste("VERSION",".",sep="     "))
        FileLines=append(FileLines,paste("SOURCE",".",sep="     "))
        FileLines=append(FileLines,paste("ORGANISM",organism,sep="     "))
        posipat=c()
        ##Match sequences
        for(i in 1:length(patterns)){
            stpos=start(matchPattern(DNAString(as.character(patterns[i])),DNAString(sequence),fixed=T))
            edpos=end(matchPattern(DNAString(as.character(patterns[i])),DNAString(sequence),fixed=T))
            if(length(stpos)>0){
                posipat=rbind(posipat, cbind(stpos,edpos,rep(tooltips[i],length(stpos)),rep(FWDcolors[i],length(stpos)),rep(REVcolors[i],length(stpos))))
            }
        }
        
        if(!(is.null(posipat))){
            colnames(posipat)=c("start","end","label","fwdc","revc")
        }
        
        if(xtraComments[1] != ""){
            FileLines=append(FileLines,paste("COMMENT",xtraComments,sep="     "))
        }
        
        if(!(is.null(posipat))){
            for(i in 1:length(patterns)){
                tempat=as.character(patterns[i])
                if(BinRevComp[i]){tempat=as.character(reverseComplement(DNAString(tempat)))}
                FileLines=append(FileLines,paste("COMMENT",paste(as.character(tooltips[i]),tempat),sep="     "))
            }
        }
        
        FileLines=append(FileLines,paste("COMMENT","Generated using wormbuilder.dev/piRNAi/",sep="     "))
        FileLines=append(FileLines,paste("COMMENT","ApEinfo:methylated:1",sep="     "))
        
        if(!(is.null(posipat))){
            FileLines=append(FileLines,paste("FEATURES             Location/Qualifiers",sep=""))
            for(n in 1:nrow(posipat)){
                nixt= which(posipat[n,3] == tooltips)
                if(BinRevComp[nixt]){
                    xnoteA="complement("
                    xnoteB=")"
                    }else{
                        xnoteA=""
                        xnoteB=""
                    }
                FileLines=append(FileLines,paste("     primer_bind     ",xnoteA,c(posipat[n,1]),"..",c(posipat[n,2]),xnoteB,"",sep=""))
                FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"",c(posipat[n,3]),"\"",sep="",collapse=""),sep="     "))
                FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"",c(posipat[n,3]),"\"",sep="",collapse=""),sep="     "))
                FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"",c(posipat[n,4]),"\"",sep=""))
                FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"",c(posipat[n,5]),"\"",sep=""))
                FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
            }
            
        }
        
        if(xtraLines != ""){
            FileLines=append(FileLines,readLines(xtraLines))
        }
        
        FileLines=append(FileLines,paste("ORIGIN"))
        
        Compseq=unlist(strsplit(sequence,""))
        
        partseq=c()
        
        for(i in seq(1,length(Compseq),10)){
            endseq=i+9
            if(length(Compseq)-i < 9){endseq=length(Compseq)}
            partseq=append(partseq,paste(Compseq[i:endseq],collapse=""))
            
        }
        
        i=1
        for(num in seq(1,length(Compseq),60)){
            index=as.character(num)
            spaces=paste(rep(" ",6-nchar(index)),collapse="")
            endseq=i+5
            if((length(partseq)-i) < 5){endseq=length(partseq)}
            FileLines=append(FileLines , paste(spaces,index," ",paste(partseq[i:(endseq)],collapse=" "),sep=""))
            
            i=i+6
        }
        
        FileLines=append(FileLines,paste("//"))
        
        return(FileLines)
    }
    
    ##Retrieve output ape
    output$DownApeOut <- downloadHandler(
        filename <- function() {
            paste("piRNAi", "gb", sep=".")
        },
        
        content <- function(file) {
            file.copy(paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), file)
        },
    )
    
    ##Retrieve construct ape
    output$DownConOut <- downloadHandler(
        filename <- function() {
            paste("piRNAi_construct", "gb", sep=".")
        },
        
        content <- function(file) {
            file.copy(paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), file)
        },
    )
    
    ##Calculate GC
    CalculateGC = function(x){
        if(!is.character(x)){return(c())}
        x=toupper(x)
        vecseq=unlist(strsplit(x,""))
        return((countPattern("C",x)+countPattern("G",x))/length(vecseq))
    }
    
    ##Clean boxes
    observeEvent(input$actionclean, {
        clust=as.character(input$clustercon)
        if(clust == "1"){
        updateTextAreaInput(session, "piRNAseq1_1", value = "")
        updateTextAreaInput(session, "piRNAseq2_1", value = "")
        updateTextAreaInput(session, "piRNAseq3_1", value = "")
        updateTextAreaInput(session, "piRNAseq4_1", value = "")
        updateTextAreaInput(session, "piRNAseq5_1", value = "")
        updateTextAreaInput(session, "piRNAseq6_1", value = "")
        }
        if(clust == "2"){
            updateTextAreaInput(session, "piRNAseq1_2", value = "")
            updateTextAreaInput(session, "piRNAseq2_2", value = "")
            updateTextAreaInput(session, "piRNAseq3_2", value = "")
            updateTextAreaInput(session, "piRNAseq4_2", value = "")
            updateTextAreaInput(session, "piRNAseq5_2", value = "")
            updateTextAreaInput(session, "piRNAseq6_2", value = "")
        }
        if(clust == "3"){
            updateTextAreaInput(session, "piRNAseq1_3", value = "")
            updateTextAreaInput(session, "piRNAseq2_3", value = "")
            updateTextAreaInput(session, "piRNAseq3_3", value = "")
            updateTextAreaInput(session, "piRNAseq4_3", value = "")
            updateTextAreaInput(session, "piRNAseq5_3", value = "")
            updateTextAreaInput(session, "piRNAseq6_3", value = "")
            updateTextAreaInput(session, "piRNAseq7_3", value = "")
        }
        if(clust == "4"){
            updateTextAreaInput(session, "piRNAseq1_4", value = "")
            updateTextAreaInput(session, "piRNAseq2_4", value = "")
            updateTextAreaInput(session, "piRNAseq3_4", value = "")
            updateTextAreaInput(session, "piRNAseq4_4", value = "")
            updateTextAreaInput(session, "piRNAseq5_4", value = "")
            updateTextAreaInput(session, "piRNAseq6_4", value = "")
            updateTextAreaInput(session, "piRNAseq7_4", value = "")
            updateTextAreaInput(session, "piRNAseq8_4", value = "")
        }
        })
})  
