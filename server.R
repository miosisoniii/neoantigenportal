library(foreach)
library(doMC)
registerDoMC(20)
library(shiny)

source("functions.R")

shinyServer(function(input, output) {
  
  selected_gene <- reactive({
    gene_seq_df %>%
      filter(gene == input$selectgene)
  }) 
  
  output$table_react <- renderTable({
    selected_gene()
  })
  
  #create directory and search file
  searchfile <- eventReactive(input$create_searchfile, {
    system(paste("mkdir ", selected_gene()$gene[1], sep=""))
    
    for (j in 1:nrow(selected_gene())){
      a <- data.frame(matrix(ncol = 1))
      for (i in 1:(nchar(as.vector(selected_gene()$seq[j]))-8)){
        a<-rbind(a,substr(selected_gene()$seq[j], i, i+8))  
      }
      sink(paste(selected_gene()$gene[j], "netmhc.txt", sep=""))
      for (i in 2:nrow(a)){
        cat(paste(">",selected_gene()$gene[j], "_", i-1,sep=""))
        cat("\n")
        cat(a$matrix.ncol...1.[i])
        cat("\n")
      }
      sink()
    }
    print(paste("Creation of Searchfile for ", selected_gene()$gene[1], " complete.", sep = ""))
  })
  
  #run netMHC
  run_netMHC <- eventReactive(input$run_netMHC, {
    #store as variable outside of running parallel loop
    storedgene <- selected_gene()$gene[1]
    foreach(i=1:ncol(hla_col)) %dopar% {
      system(paste("~/netMHC -f ", 
                   paste(storedgene, "netmhc.txt", sep=""), 
                   " -a ", hlanames[i], " > ", storedgene, "/", hlanames[i], ".txt", sep=""))
    }
    print(paste("NetMHC search for ", storedgene, " complete.", sep = ""))
  })
    output$searchfile_complete <- renderText({
    searchfile()
  })
  output$netmhc_complete <- renderText({
    run_netMHC()
  })
  
  #create list of netMHC files
  #this code creates list of .txt files, but full.name includes file path/connection
  readingfiles1 <- eventReactive(input$read_output, {
    HLAfiles <- list.files(path = paste(selected_gene()$gene[1], sep = ""), pattern = "*.txt", 
                           full.names = "TRUE") #can remove to obtain only file names in list
    HLAfiles <- HLAfiles[HLAfiles != "NA.txt"]
    HLAfiles
  })
  output$read1complete <- renderText({
    #print(paste("Total number of HLA's: ", length(readingfiles1())))
    readingfiles1()[1]
  })
  
  removeblank1 <- eventReactive(input$removeblank1_go, {
    system(paste("mkdir ", selected_gene()$gene[1], "/peptides/", sep = ""))
    for(f in readingfiles1()){
      x <- readLines(f)
      y <- gsub( "<= ", "<=", x )
      cat(y, file=f, sep="\n")
    }
  })
  output$removeblank1_test <- renderText({
    removeblank1()
  })
  output$removeblank1complete <- renderText({
    removeblank1()
    print("Removal of blanks complete.")
  })
  
  
  #Filter Peptides
  filterpep <- eventReactive(input$filter_pep, {
    for (f in readingfiles1()){
      x <- readLines(f)
      y <- x[-(1:40)]
      z <- c(x[39],y[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)])
      z <- data.frame(z)
      f <-  f %>% str_replace(".*/", "")
      write.table(z, file = paste(selected_gene()$gene[1], "/peptides/", f, sep=""), quote = F, row.names = F, col.names = F, sep ="\t")
    }
     print("Peptide filtering complete.")
  })
  output$filter_complete <- renderText({
    filterpep()
  })
  
  #read in files from peptides directory
  readingfiles2 <- eventReactive(input$read_output2, {
    HLAfiles <- list.files(path = paste(selected_gene()$gene[1], "/peptides", sep = ""), pattern = "*.txt", 
                           full.names = "TRUE") #can remove to obtain only file names in list
    HLAfiles <- HLAfiles[HLAfiles != "NA.txt"]
    HLAfiles
  })
  output$read2complete <- renderText({
    head(readingfiles2())
  })
  
  ###code does not loop through entire list of peptides to create list of binders!
  #Combine Binders into matrix and populate with binding score
  writebinders <- eventReactive(input$binders_go, {
    #create binders directory
    system(paste("mkdir ", selected_gene()$gene[1], "/peptides/binders/", sep = ""))
    filestest <- readingfiles2()
    for (f in filestest){
      a <- read.delim(f ,header=T,stringsAsFactors=F, sep="")
      b <- subset(a, a$BindLevel=="<=SB")
      f <-  f %>% str_replace(".*/", "")
      write.table(b, file=paste(selected_gene()$gene[1], "/peptides/binders/", f, sep=""), quote = F, row.names = F, sep ="\t")
    #this returns the table for a and b subsetted!
    #return(a)
    }
    print("Binder creation complete.")
  })
  # output$binders_complete <- renderText({
  #   nrow(writebinders())
  # })
  output$binders_delim_table <- renderTable({
    head(writebinders())
  })
  
  #read in files from binders directory
  readingfiles3 <- eventReactive(input$read_output3, {
  HLAfiles1 <- list.files(path = paste(selected_gene()$gene[1], "/peptides/binders", sep = ""), pattern = "*.txt",
                    full.names = "TRUE") #can remove to obtain only file names in list
  #HLAfiles1 <- HLAfiles1[HLAfiles1 != "NA.txt"]
  HLAfiles1
  })
  output$read3complete <- renderText({
    head(readingfiles3())
  })
  
  #combine binders pt 1/2
  combine_test1 <- eventReactive(input$combine_go, {
    files <- readingfiles3()
    
    #TESTING ----> WORKS!
    #do.call(rbind, lapply(files, read.delim, header = T, stringsAsFactors = F, sep = "")) -> a
    
    
    #b <- subset(a, a$BindLevel == "<=SB")
    
    # a %>% 
    #   distinct(Identity, .keep_all = TRUE) -> a
    ###TESTING

    #this produces table but does not include the gene/protein names as rows.
    #need to create df from read.delim of file names
    #combined <- data.frame(matrix(nrow = nrow(a), ncol = ncol(hla[,-1])))
    
    #cannot have duplicate rownames in DF
    #a <- writebinders()
    a <- read.delim(readingfiles2()[1], header = T, stringsAsFactors = F, sep = "")
    combined <- data.frame(matrix(nrow = nrow(a), ncol = ncol(hla[,-1])))
    
    #atttemot to use matrix
    # combined <- matrix(nrow = nrow(a), ncol = ncol(hla[,-1]))
    # 
    # #cannot have duplicate rownames
    row.names(combined) <- a$Identity
    colnames(combined) <- colnames(hla[,-1])
    colnames(combined) <- gsub("\\.", "-", colnames(combined))
    # 
   

    binders <- readingfiles3()
    for (b in binders){
      binder <- read.delim(b, header = T, stringsAsFactors = F, sep="")
      muts <- binder

      if (nrow(muts) >0){
        for (i in 1:nrow(muts)){
          combined[muts$Identity[i], muts$HLA[i]] <- muts$X.Rank[i]
        }
      }
    }
    
    combined
    
    
  })
  # output$combo_complete1 <- renderTable({
  #   head(combine_test1())
  # },rownames = TRUE)
  # 
  # output$testcol <- renderText({
  #   combine_test1()
  # })

  #combine binders pt 2/2
  combine_test2 <- eventReactive(input$combine_go2, {
    
    combine_test1() -> test1

    #Set HLA frequency
    #hla has 85 columns, 1st column has incorrect HLAfrequency title
    #test1["HLA_frequency", 1:84] <- hla[1,2:85]
    test1["HLA_frequency", 1:84] <- hla[1,2:84]
    #this changes all NA's to 0's
    test1[is.na(test1)] <- as.numeric(0)
    write.csv(test1, "testcombined.csv")

    test1

  })
  output$combo_complete2 <- renderTable({
    combine_test2()[53,]
  }, rownames = TRUE)
  
  #Calculate HLA binders and frequency for all peptides
  calculateHLA <- eventReactive(input$calc_HLA_go, {
    combine_test2() -> combined
    for (i in 1:(nrow(combined)-1)){
      
      neofreqTCGA <- c()
      hlabinders <- c()
      
      #calculate HLA frequency for rows/columns in combined table
      for (j in 1:84){
        if (as.numeric(combined[i,j]) > 0){
          neofreqTCGA <- c(neofreqTCGA, combined["HLA_frequency",j])
          hlabinders <- c(hlabinders, substring(colnames(combined[j]), 5,9))
        }
        
        probTCGA <- 1
        if (length(neofreqTCGA) > 0){
          for (k in 1:length(neofreqTCGA)){
            probTCGA <- probTCGA*(1-as.numeric(neofreqTCGA[k]))
          }
        }
      }
      probTCGA <- 1-probTCGA
      
      combined$"HLA_frequency"[i] <- probTCGA
      combined$"Alleles_bound"[i] <- length(hlabinders)
      combined$"HLA_binders"[i] <- paste(unlist(hlabinders), collapse = ", ")
    }
    
    combined
  })
  output$calc_HLA_out <- renderTable({
    #show last 3() columns of HLA calculations
    #head(calculateHLA())
    calculateHLA()[53, "HLA_frequency"]
  }, rownames = TRUE)
  
  
  #set AA postions
  score1test <- eventReactive(input$set_AA, {
    combined <- calculateHLA()
    combined$position <- row.names(combined)
    combined$gene <- gsub( "_.*$", "", row.names(combined))
    for (i in 1:nrow(combined)){
      #column that produces amino acid position
      combined$position[i] <- sub('.*\\_', '', row.names(combined)[i])
    }
    
    combined
  })
  output$scoring_table <- renderTable({
    score1test()[53, c("HLA_frequency", "gene")]
  },rownames = TRUE)
  
  
  
  #set peptide sequences for 9aa
  score2test <- eventReactive(input$scoring_table, {
    combined <- score1test()
    gentab_test <- selected_gene()[1,]
    
    for (i in 1:nrow(combined)){
      r <- match(gsub(" ", "",combined$gene[i]), gentab_test$gene)
      combined$aa[i] <- substr(gentab_test$seq[r], combined$position[i], combined$position[i])
      combined$pep[i] <- substr(gentab_test$seq[r], combined$position[i], (as.numeric(combined$position[i])+8))
    }
    write.csv(combined, "maptest.csv")
    combined
  })
  
  output$scoring_out <- renderTable({
    score2test()[53, c("HLA_frequency", "gene") ]
  }, rownames = TRUE)
  
  
  plot_9aa <- eventReactive(input$plot9aa, {
    combined <- score2test()
    barplot(combined$HLA_frequency, ylim=c(0,1))
  })
  
  output$plot_9aa_out <- renderPlot({
    plot_9aa()
  })
  

})





