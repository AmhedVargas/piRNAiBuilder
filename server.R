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
            Select the specificity of piRNAi fragments. Each piRNAi fragment was mapped across the C. elegans genome and verified its uniqueness up to n Mismatches.
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
        ErrorFlag = 0
        
        matches=as.integer(input$selectMM)
        mm=c(5,4,3,2,1,0)[matches]
        isoform = input$isoform
        ControlEx = input$FlaControl
        
        wbid = as.character(unique(Genes[which(isoform==Genes[,5]),4]))
        loc = as.character(unique(Genes[which(isoform==Genes[,5]),6]))
        
        strand = as.character(unique(Genes[which(isoform==Genes[,5]),7]))
        genest = as.numeric(unique(Genes[which(isoform==Genes[,5]),2]))
        geneend = as.numeric(unique(Genes[which(isoform==Genes[,5]),3]))
        
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
        
        ##Produce outputs
            if(ErrorFlag == 0){
                
                ##Table results
                Pitab=Seltab[idx,c(1,2,4)]
                
                
                Pitab[,1]=c(Pitab[,1]-genest)/(geneend-genest)
                if(strand =="-"){Pitab[,1]= 1 - as.numeric(Pitab[,1])}
                position=c()
                for(val in Pitab[,1]){
                    posssi="Unkwon"
                    if((val >= 0)&(val <= 0.25)){posssi="Proximal to 5' end"}
                    if((val >= 0.25)&(val <= 0.75)){posssi="Centered around isoform"}
                    if((val >= 0.75)&(val <= 1)){posssi="Proximal to 3' end"}
                    position=append(position,posssi)
                    }
                
                Pitab=cbind(position,Pitab[,c(2,3)])
                
                colnames(Pitab)=c("Location","Sequence","%GC")
                output$SelPiTabSummary <- renderUI({ HTML(paste0("<b>Selected sequences: </b>",sep=""))})
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
                    seq1=as.character(Seltab[idx[1],2])
                    dos="ctaattttgattttgattttgaaatcgaatttgcaaatccaattaaaaatcattttctgataattagacagttccttatcgttaattttattatatctatcgagttagaaattgcaacgaagataatgtcttccaaatactgaaaatttgaaaatatgtt"
                    seq2=as.character(reverseComplement(DNAString(as.character(Seltab[idx[2],2]))))
                    tres="AttGccagaactcaaaatatgaaatttttatagttttgttgaaacagtaagaaaatcttgtaattactgtaaactgtttgctttttttaaagtcaacctacttcaaatctacttcaaaaattataatgtttcaaattacataactgtgt"
                    seq3= as.character(reverseComplement(DNAString(as.character(Seltab[idx[3],2]))))
                    cuatro="ActgtagagcttcaatgttgataagatttattaacacagtgaaacaggtaatagttgtttgttgcaaaatcggaaatctctacatttcatatggtttttaattacaggtttgttttataaaataattgtgtgatggatattattttcagacctcatactaatctgcaaaccttcaaacaatatgtgaagtctactctgtttcactcaaccattcatttcaatttggaaaaaaatcaaagaaatgttgaaaaattttcctgtttcaacattatgacaaaaatgttatgattttaataaaaaCaaT"
                    seq4=as.character(Seltab[idx[4],2])
                    cinco="ttctgtttttcttagaagtgttttccggaaacgcgtaattggttttatcacaaatcgaaaacaaacaaaaatttttttaattatttctttgctagttttgtagttgaaaattcactataatcatgaataagtgagctgcccaagtaaacaaagaaaatttggcagcggccgacaactaccgggttgcccgatttatcagtggagga"
                    seq5= as.character(reverseComplement(DNAString( as.character(Seltab[idx[5],2]))))
                    seis="AtcTaatgtgatgtacacggttttcatttaaaaacaaattgaaacagaaatgactacattttcaaattgtctatttttgctgtgtttattttgccaccaacaaT"
                    seq6=as.character(Seltab[idx[6],2])
                    siete="tcaatctagtaaactcacttaatgcaattcctccagccacatatgtaaacgttgtatacatgcagaaaacggttttttggttttaatgggaacttttgacaaattgttcgaaaatcttaagctgtcccatttcagttgggtgatcgattt"
                        #write(paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete),sep="",collapse=""),paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""))

                    write(paste("LOCUS",paste(wbid,"_21ur_1224_",sep="",collapse=""),"1344 bp ds-DNA","linear",paste(c(unlist(strsplit(date()," ")))[c(3,2,5)],sep="",collapse="-"),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""))
                    write(paste("DEFINITION",".",sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("ACCESSION",".",sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("VERSION",".",sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("SOURCE",".",sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("ORGANISM","C.elegans",sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    if(ControlEx){
                        write(paste("COMMENT     ",paste("Control experiment: piRNAi Sequences have been reverse complemented. THIS FRAGMENT WONT SILENCE THE SELECTED GENE"),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    }
                    write(paste("COMMENT     ",paste("Recoded 21ur-1224 locus. Parameters: Gene =",wbid,"; Isoform = ",isoform,"; Up to",mm,"mismatches"),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("COMMENT",paste("piRNA1",as.character(Seltab[idx[1],2])),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("COMMENT",paste("piRNA2",as.character(Seltab[idx[2],2])),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("COMMENT",paste("piRNA3",as.character(Seltab[idx[3],2])),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("COMMENT",paste("piRNA4",as.character(Seltab[idx[4],2])),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("COMMENT",paste("piRNA5",as.character(Seltab[idx[5],2])),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("COMMENT",paste("piRNA6",as.character(Seltab[idx[6],2])),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("COMMENT","Generated using wormbuilder.dev/piRNAi/",sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("COMMENT","",sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("COMMENT","ApEinfo:methylated:1",sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("FEATURES             Location/Qualifiers",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("     primer_bind     152..171",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ",paste("/locus_tag=","\"",isoform," piRNA1\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ",paste("/ApEinfo_label=","\"",isoform," piRNA1\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_fwdcolor=\"#00ff00\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_revcolor=\"#ff0000\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("     primer_bind     complement(332..351)",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ",paste("/locus_tag=","\"",isoform," piRNA2\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ",paste("/ApEinfo_label=","\"",isoform," piRNA2\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_fwdcolor=\"#00ff00\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_revcolor=\"#ff0000\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("     primer_bind     complement(501..520)",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ",paste("/locus_tag=","\"",isoform," piRNA3\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ",paste("/ApEinfo_label=","\"",isoform," piRNA3\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_fwdcolor=\"#00ff00\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_revcolor=\"#ff0000\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("     primer_bind     825..844",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ",paste("/locus_tag=","\"",isoform," piRNA4\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ",paste("/ApEinfo_label=","\"",isoform," piRNA4\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_fwdcolor=\"#00ff00\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_revcolor=\"#ff0000\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("     primer_bind     complement(1051..1070)",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ",paste("/locus_tag=","\"",isoform," piRNA5\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ",paste("/ApEinfo_label=","\"",isoform," piRNA5\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_fwdcolor=\"#00ff00\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_revcolor=\"#ff0000\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("     primer_bind     1175..1194",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ",paste("/locus_tag=","\"",isoform," piRNA6\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ",paste("/ApEinfo_label=","\"",isoform," piRNA6\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_fwdcolor=\"#00ff00\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_revcolor=\"#ff0000\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(readLines("WorkingSpace/Piconst.txt"),paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    
                    Compseq=paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete),sep="",collapse="")
                    Compseq=unlist(strsplit(Compseq,""))
                    
                    partseq=c()
                    
                    for(i in seq(1,length(Compseq),10)){
                        endseq=i+9
                        if(length(Compseq)-i < 9){endseq=length(Compseq)}
                        partseq=append(partseq,paste(Compseq[i:endseq],collapse=""))
                        
                    }
                    
                    write(paste("        1 ",paste(partseq[1:6],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("       61 ",paste(partseq[7:12],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      121 ",paste(partseq[13:18],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      181 ",paste(partseq[19:24],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      241 ",paste(partseq[25:30],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      301 ",paste(partseq[31:36],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      361 ",paste(partseq[37:42],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      421 ",paste(partseq[43:48],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      481 ",paste(partseq[49:54],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      541 ",paste(partseq[55:60],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      601 ",paste(partseq[61:66],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      661 ",paste(partseq[67:72],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      721 ",paste(partseq[73:78],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      781 ",paste(partseq[79:84],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      841 ",paste(partseq[85:90],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      901 ",paste(partseq[91:96],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("      961 ",paste(partseq[97:102],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("     1021 ",paste(partseq[103:108],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("     1081 ",paste(partseq[109:114],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("     1141 ",paste(partseq[115:120],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("     1201 ",paste(partseq[121:126],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("     1261 ",paste(partseq[127:132],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write(paste("     1321 ",paste(partseq[133:length(partseq)],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    write("//",paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), append=T)
                    
                    downloadButton('DownApeOut', 'Download piRNAi fragment')
                })
                
                
            }
        
    }, ignoreInit = T)
    
    ###Advanced construct######
    observeEvent(input$actionconstruct, {
        AdvancedErrorFlag=0
        output$AdvancedErrorMessage <- renderText({})
        
        pipi1=as.character(input$piRNAseq1)
        pipi2=as.character(input$piRNAseq2)
        pipi3=as.character(input$piRNAseq3)
        pipi4=as.character(input$piRNAseq4)
        pipi5=as.character(input$piRNAseq5)
        pipi6=as.character(input$piRNAseq6)
        
        ##Check for size
        if((nchar(pipi1) != 20) | (nchar(pipi2) != 20) | (nchar(pipi3) != 20) | (nchar(pipi4) != 20) | (nchar(pipi5) != 20) | (nchar(pipi6) != 20)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Error: piRNAi sequences should be 20bp long")
            })
            }
        
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
                seq1=as.character(pipi1)
                dos="ctaattttgattttgattttgaaatcgaatttgcaaatccaattaaaaatcattttctgataattagacagttccttatcgttaattttattatatctatcgagttagaaattgcaacgaagataatgtcttccaaatactgaaaatttgaaaatatgtt"
                seq2=as.character(reverseComplement(DNAString(as.character(pipi2))))
                tres="AttGccagaactcaaaatatgaaatttttatagttttgttgaaacagtaagaaaatcttgtaattactgtaaactgtttgctttttttaaagtcaacctacttcaaatctacttcaaaaattataatgtttcaaattacataactgtgt"
                seq3= as.character(reverseComplement(DNAString(as.character(pipi3))))
                cuatro="ActgtagagcttcaatgttgataagatttattaacacagtgaaacaggtaatagttgtttgttgcaaaatcggaaatctctacatttcatatggtttttaattacaggtttgttttataaaataattgtgtgatggatattattttcagacctcatactaatctgcaaaccttcaaacaatatgtgaagtctactctgtttcactcaaccattcatttcaatttggaaaaaaatcaaagaaatgttgaaaaattttcctgtttcaacattatgacaaaaatgttatgattttaataaaaaCaaT"
                seq4=as.character(pipi4)
                cinco="ttctgtttttcttagaagtgttttccggaaacgcgtaattggttttatcacaaatcgaaaacaaacaaaaatttttttaattatttctttgctagttttgtagttgaaaattcactataatcatgaataagtgagctgcccaagtaaacaaagaaaatttggcagcggccgacaactaccgggttgcccgatttatcagtggagga"
                seq5= as.character(reverseComplement(DNAString( as.character(pipi5))))
                seis="AtcTaatgtgatgtacacggttttcatttaaaaacaaattgaaacagaaatgactacattttcaaattgtctatttttgctgtgtttattttgccaccaacaaT"
                seq6=as.character(pipi6)
                siete="tcaatctagtaaactcacttaatgcaattcctccagccacatatgtaaacgttgtatacatgcagaaaacggttttttggttttaatgggaacttttgacaaattgttcgaaaatcttaagctgtcccatttcagttgggtgatcgattt"
                
                write(paste("LOCUS",paste("Undefined_21ur_1224_",sep="",collapse=""),"1344 bp ds-DNA","linear",paste(c(unlist(strsplit(date()," ")))[c(3,2,5)],sep="",collapse="-"),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""))
                write(paste("DEFINITION",".",sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("ACCESSION",".",sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("VERSION",".",sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("SOURCE",".",sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("ORGANISM",".",sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("COMMENT     ",paste("Recoded 21ur-1224 locus."),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("COMMENT",paste("piRNA1",as.character(pipi1)),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("COMMENT",paste("piRNA2",as.character(pipi2)),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("COMMENT",paste("piRNA3",as.character(pipi3)),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("COMMENT",paste("piRNA4",as.character(pipi4)),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("COMMENT",paste("piRNA5",as.character(pipi5)),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("COMMENT",paste("piRNA6",as.character(pipi6)),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("COMMENT","Generated using wormbuilder.dev/piRNAi/",sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("COMMENT","",sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("COMMENT","ApEinfo:methylated:1",sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("FEATURES             Location/Qualifiers",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("     primer_bind     152..171",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ",paste("/locus_tag=","\"piRNA1\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ",paste("/ApEinfo_label=","\"piRNA1\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_fwdcolor=\"#00ff00\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_revcolor=\"#ff0000\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("     primer_bind     complement(332..351)",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ",paste("/locus_tag=","\"piRNA2\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ",paste("/ApEinfo_label=","\"piRNA2\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_fwdcolor=\"#00ff00\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_revcolor=\"#ff0000\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("     primer_bind     complement(501..520)",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ",paste("/locus_tag=","\"piRNA3\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ",paste("/ApEinfo_label=","\"piRNA3\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_fwdcolor=\"#00ff00\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_revcolor=\"#ff0000\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("     primer_bind     825..844",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ",paste("/locus_tag=","\"piRNA4\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ",paste("/ApEinfo_label=","\"piRNA4\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_fwdcolor=\"#00ff00\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_revcolor=\"#ff0000\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("     primer_bind     complement(1051..1070)",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ",paste("/locus_tag=","\"piRNA5\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ",paste("/ApEinfo_label=","\"piRNA5\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_fwdcolor=\"#00ff00\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_revcolor=\"#ff0000\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("     primer_bind     1175..1194",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ",paste("/locus_tag=","\"piRNA6\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ",paste("/ApEinfo_label=","\"piRNA6\"",sep="",collapse=""),sep="     "), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_fwdcolor=\"#00ff00\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_revcolor=\"#ff0000\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(readLines("WorkingSpace/Piconst.txt"),paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                
                Compseq=paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete),sep="",collapse="")
                Compseq=unlist(strsplit(Compseq,""))
                
                partseq=c()
                
                for(i in seq(1,length(Compseq),10)){
                    endseq=i+9
                    if(length(Compseq)-i < 9){endseq=length(Compseq)}
                    partseq=append(partseq,paste(Compseq[i:endseq],collapse=""))
                    
                }
                
                write(paste("        1 ",paste(partseq[1:6],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("       61 ",paste(partseq[7:12],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      121 ",paste(partseq[13:18],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      181 ",paste(partseq[19:24],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      241 ",paste(partseq[25:30],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      301 ",paste(partseq[31:36],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      361 ",paste(partseq[37:42],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      421 ",paste(partseq[43:48],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      481 ",paste(partseq[49:54],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      541 ",paste(partseq[55:60],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      601 ",paste(partseq[61:66],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      661 ",paste(partseq[67:72],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      721 ",paste(partseq[73:78],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      781 ",paste(partseq[79:84],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      841 ",paste(partseq[85:90],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      901 ",paste(partseq[91:96],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("      961 ",paste(partseq[97:102],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("     1021 ",paste(partseq[103:108],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("     1081 ",paste(partseq[109:114],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("     1141 ",paste(partseq[115:120],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("     1201 ",paste(partseq[121:126],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("     1261 ",paste(partseq[127:132],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write(paste("     1321 ",paste(partseq[133:length(partseq)],collapse=" "),sep=""), paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                write("//",paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), append=T)
                
                
            downloadButton('DownConOut', 'Download piRNAi fragment')
            
                })
            
            }
        }, ignoreInit = T)
    
    ###Advanced search
    observeEvent(input$actionAdvsearch, {
        
        type="Not"
        wbid=""
        mygene =input$Advancedgeneinput
        
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
                    
                    column(width = 3,selectInput("AdvIsoform", label = HTML("<b>Select Isoform
                                                           [<a href=\"\" onclick=\"$('#explain_isoform_advanced').toggle(); return false;\">info</a>]
                                                           </b>"), 
                                unique(Genes[which(as.character(Genes[,4])==wbid),5]))),
                    
                    column(width = 3, selectInput("AdvSelectMM", label = HTML("<b>Uniqueness of sequence
                    [<a href=\"\" onclick=\"$('#explain_uniqueness_advanced').toggle(); return false;\">info</a>]
                                                           </b>"),
                                                  choices = list("Up to 5 Mismatches across genome" = 1, "Up to 4 Mismatches across genome" = 2, "Up to 3 Mismatches across genome" = 3,
                                                                 "Up to 2 Mismatches across genome" = 4,"Up to 1 Mismatches across genome" = 5,"Up to 0 Mismatches across genome" = 6),
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
            Select the specificity of piRNAi fragments. Each piRNAi fragment was mapped across the C. elegans genome and verified its uniqueness up to n Mismatches.
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
        tab[,4]=as.integer((unlist(strsplit(as.character(tab[,2]),";")))[c(FALSE,TRUE)])
        tab[,2]=as.character((unlist(strsplit(as.character(tab[,2]),";")))[c(TRUE,FALSE)])
        tab[,1]=c(tab[,1]-genest)/(geneend-genest)
        
        if(strand =="-"){tab[,1]= 1 - as.numeric(tab[,1])}
        
        #output$AllPiTab <- DT::renderDataTable(DT::datatable({
        output$AllPiTab <- DT::renderDataTable({
            datatab = tab
            matches=as.integer(input$AdvSelectMM)
            mm=c(5,4,3,2,1,0)[matches]
            minGC=input$Gcslider[1]
            maxGC=input$Gcslider[2]
            minPos=input$Posslider[1]/100
            maxPos=input$Posslider[2]/100
            
            datatab = datatab[which(datatab[,3]>=mm),]
            
            datatab = datatab[which((datatab[,4]>=minGC)&(datatab[,4]<=maxGC)),]
            
            datatab = datatab[which((datatab[,1]>=minPos)&(datatab[,1]<=maxPos)),]
            
            position=c()
            for(val in datatab[,1]){
                posssi="Unkwon"
                if((val >= 0)&(val <= 0.25)){posssi="Proximal to 5' end"}
                if((val >= 0.25)&(val <= 0.75)){posssi="Centered around isoform"}
                if((val >= 0.75)&(val <= 1)){posssi="Proximal to 3' end"}
                position=append(position,posssi)
            }
            
            Pdata=data.frame(
                Location = position,
                Sequence = datatab[,2],
                GCcontent = datatab[,4],
                Select = shinyInput(actionButton, as.character(datatab[,2]), 'button_', label = "Add to contruct", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' ),
                stringsAsFactors = FALSE,
                row.names = 1:nrow(datatab)
                )
            
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
        fill=0
        selectedSeq <- as.character(strsplit(input$select_button, "_")[[1]][2])
        
        if(fill == 0){
        if(as.character(input$piRNAseq1)==""){
            updateTextAreaInput(session, "piRNAseq1", value = selectedSeq)
            fill = 1
        }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq2)==""){
                updateTextAreaInput(session, "piRNAseq2", value = selectedSeq)
                fill = 1
            }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq3)==""){
                updateTextAreaInput(session, "piRNAseq3", value = selectedSeq)
                fill = 1
            }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq4)==""){
                updateTextAreaInput(session, "piRNAseq4", value = selectedSeq)
                fill = 1
            }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq5)==""){
                updateTextAreaInput(session, "piRNAseq5", value = selectedSeq)
                fill = 1
            }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq6)==""){
                updateTextAreaInput(session, "piRNAseq6", value = selectedSeq)
                fill = 1
            }}
        
        ##Send custom message
        if(fill == 0){
            showNotification("Construct has already 6 sequences.")
            }
        
    })
    
    ####Other functions
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
})  
