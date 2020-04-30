library(shiny)
source("functions.R")

shinyUI(
  fluidPage(
    titlePanel("ShinyNAP (Neo-Antigen Portal)"),
    sidebarLayout(
      sidebarPanel(
        selectInput(inputId = 'selectgene',
                    label = 'Select Gene to Analyze',
                    choices = gene_seq_df$gene,
                    selected = "MYCN"),
        actionButton("create_searchfile", "Create Searchfile"),
        br(),
        actionButton("run_netMHC", "Run netMHC"),
        br(),
        actionButton("read_output", "Read Output 1"),
        br(),
        br(),
        br(),
        actionButton("removeblank1_go", "Remove Blanks"),
        br(),
        actionButton("filter_pep", "Filter Peptides"),
        br(),
        br(),
        actionButton("read_output2", "Read Peptide Files 2"),
        br(),
        actionButton("binders_go", "Create Binders"),
        br(),
        actionButton("read_output3", "Read Binders Files 3"),
        br(),
        actionButton("combine_go", "Combine Binders"),
        br(),
        actionButton("combine_go2", "Combine Binders 2"),
        br(),
        br(),
        actionButton("calc_HLA_go", "Calculate HLA Frequency"),
        br(),
        br(),
        actionButton("set_AA", "Set Amino Acids"),
        br(),
        actionButton("scoring_table", "Analyze Scores"),
        br(),
        actionButton("plot9aa", "Plot Top Scoring Region for 9 Amino Acids")
      ),
      mainPanel(
        tableOutput("table_react"),
        textOutput("searchfile_complete"),
        textOutput("netmhc_complete"),
        textOutput("read1complete"),
        textOutput("removeblank1_test"),
        textOutput("removeblank1complete"),
        textOutput("filter_complete"),
        textOutput("read2complete"),
        textOutput("binders_complete"),
        tableOutput("binders_delim_table"),
        textOutput("read3complete"),
        tableOutput("combo_complete1"),
        textOutput("testcol"),
        tableOutput("combo_complete2"),
        tableOutput("calc_HLA_out"),
        tableOutput("scoring_table"),
        tableOutput("scoring_out"),
        plotOutput("plot_9aa_out"), width = 6
      )
    )
  )
)

