library(tidyverse)
library(shiny)
library(plotly)
library(DT)
library(shinyBS)

percent = c("0%", "25%", "50%", "75%", "100%")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Subset the original full_cleaned.tsv so it fits in the 1GB memory on shinyapps.io
# read_tsv(
#     "data/full_cleaned.tsv",
#     col_types = cols(
#         .default = col_character(),
#         Start_Position = col_double(), End_Position = col_double(),
#         t_depth = col_double(), t_ref_count = col_double(), t_alt_count = col_double(),
#         n_depth = col_double(), n_ref_count = col_double(), n_alt_count = col_double()
#     )
# ) %>% select(
#     mc3_exome_barcode,
#     Hugo_Symbol, `Hugo_Symbol:1`,
#     Chromosome, `Chromosome:1`,
#     t_depth, t_alt_count,
#     i_VAF, 
#     Variant_Classification, HGVSp_Short, `Variant_Classification:1`,
#     FILTER
# ) %>% write_tsv('data/full_cleaned.for_shiny.tsv.gz')

# Read the final mutation overlap table
data_tbl <- read_tsv(
        "data/full_cleaned.for_shiny.tsv.gz",
        col_types = cols(
            .default = col_character(),
            t_depth = col_double(), t_alt_count = col_double()
        )
    ) %>%
    # Introduce the indicator variables for the overlap status
    mutate(
        match = ifelse(!is.na(Chromosome) & !is.na(`Chromosome:1`), 1, 0),
        PCAWG_only = ifelse(is.na(Chromosome) & !is.na(`Chromosome:1`), 1, 0),
        MC3_only = ifelse(!is.na(Chromosome) & is.na(`Chromosome:1`), 1, 0)
    ) %>%
    # Calculate VAFs
    mutate(
        m_VAFfixed = case_when(
            is.na(FILTER) ~ NA_real_,
            TRUE ~ t_alt_count/t_depth
        ),
        i_VAFfixed = sapply(
            strsplit(i_VAF, "|", fixed = TRUE), 
            function(x) mean(as.numeric(x))
        )
    )

# Select all the samples
samples <- data_tbl %>% 
    distinct(mc3_exome_barcode) %>% 
    mutate(sample = substr(mc3_exome_barcode, 1, 12)) %>%
    pull(sample)

# User interface
ui <- fluidPage(
    img(src='MAFit_logo2.png', align = "right"),
    titlePanel("ICGC MC3 Overlap"),
    
    fluidRow(
        column(3, 
                 
                       sliderInput("ivaf_slider", "ICGC VAF Cut-off", min=0, max=1, value= c(0.0,1)),
               
                       sliderInput("mvaf_slider", "MC3 VAF Cut-off", min=0, max=1, value= c(0.0,1)),
                 
                 
                       checkboxGroupInput("filts_var", "MC3 Variant-Level Filters:",
                                  c("PASS"="PASS",
#                                    "oxog"="oxog", #Any samples with oxog where removed from this app
                                    "common_in_exac"="common_in_exac",
                                    "StrandBias"="StrandBias")),
                 
                       checkboxGroupInput("filts_sample", "MC3 Sample-Level Filters:",
                                  c("nonpreferredpair"="nonpreferredpair",
                                    "native_wga_mix"="native_wga_mix",
                                    "wga"="wga",
                                    "gapfiller"="gapfiller")),
                 
                       downloadButton("downloadData", "Download")
        ),
        column(9, plotlyOutput("pid", height = "100%"))
    ),
    
    fluidRow(
        column(4, offset = 2, selectInput("bar", "TCGA-Barcode", c("All", samples))),
        column(4, selectizeInput("gene", "HugoSymbol", 
                              c("Type gene symbol", 
                                sort(unique(c(data_tbl$Hugo_Symbol, data_tbl$"Hugo_Symbol:1")))),
                              options = list(maxOptions = 25000)
        ))
    ),
    
    fluidRow(column(12, DT::dataTableOutput("table"))),

    bsTooltip(id = "ivaf_slider", title = "Select a variant allele fraction window for PCAWG vars", trigger = "hover"),
    bsTooltip(id = "mvaf_slider", title = "Select a variant allele fraction windows for MC3 vars", trigger = "hover"),
    bsTooltip(id = "filts_var", title = "Select a variant with a given MC3 filter (PASS = default)", trigger = "hover"),
    bsTooltip(id = "filts_sample", title = "Select samples with a given MC3 filter", trigger = "hover"),
    bsTooltip(id = "downloadData", title = "Download mutation data based on your filter criteria", trigger = "hover")
) 

server <- function(input, output, session) {
    # when the slider changes update the text box
    observe({
        updateTextInput(session, "mytext", value=c(input$ivaf_slider, input$mvaf_slider))
    })

   
    dat0 <- reactive({
        matches_var <- paste(input$filts_var, collapse="|")
        matches_sample <- paste(input$filts_sample, collapse="|")
        # VAF filter
        SHINE2 <- data_tbl %>%
            filter(
                (i_VAFfixed >= input$ivaf_slider[1] & i_VAFfixed < input$ivaf_slider[2]) | is.na(i_VAFfixed),
                (m_VAFfixed >= input$mvaf_slider[1] & m_VAFfixed < input$mvaf_slider[2]) | is.na(m_VAFfixed)
            )
        # variant-level filter
        SHINE3 <- SHINE2[which(grepl(matches_var, SHINE2$FILTER) | is.na(SHINE2$FILTER)),]
        # sample-level filter
        if (matches_sample == "") { 
            SHINE5 <- SHINE3 
        } else {
            SHINE4 <- SHINE2[which(grepl(matches_sample, SHINE2$FILTER)),]
            SHINE5 <- SHINE2[which(
                grepl(matches_sample, SHINE2$FILTER) |
                    grepl(matches_var, SHINE2$FILTER) | 
                    (is.na(SHINE2$FILTER) & SHINE2$mc3_exome_barcode %in% SHINE4$mc3_exome_barcode)),]
        }
        SHINE5
    })
    # when the slider changes update the dataset
    dat1 <- reactive({
        SHINE5 <- dat0()
        
        newdat <- SHINE5 %>% 
            group_by(mc3_exome_barcode) %>% 
            summarize(id_match = sum(match), id_mc3 = sum(MC3_only) , id_pcawg = sum(PCAWG_only))
        toPlot <- newdat %>% 
            filter(id_match + id_mc3 != 0) %>%
            mutate(
                perc_PCAWG_in_MC3 = id_match / (id_match + id_mc3),
                perc_MC3_in_PCAWG = id_match / (id_match + id_pcawg),
                sample = substr(mc3_exome_barcode, 1, 12)
            )
        toPlot 
    })
    
    dat2 <- reactive({
        SHINE6 <- dat0() %>% select(
                mc3_exome_barcode,
                Hugo_Symbol, `Hugo_Symbol:1`,
                m_VAFfixed, i_VAFfixed, 
                Variant_Classification, HGVSp_Short, `Variant_Classification:1`
            )
        SHINE7 <- SHINE6 %>% mutate_if(is.numeric, round, 3)
        SHINE7
        
    })
    
    dat3 <- reactive({
        tab_data <- dat2()
        if (input$bar == "All" & input$gene == "Type gene symbol" & nrow(tab_data) >= 2000 )  {
            # Show a truncated table if the current table contains too many rows
            tab_data <- data.frame(
                rbind(
                    c('(Table too large)', rep('...', 7)) 
                )
            ) %>% as_tibble()
        } else {
            # Filter the table by sample and gene
            if (input$bar != "All") {
                tab_data <- tab_data %>% 
                    filter(startsWith(mc3_exome_barcode, input$bar))
            }
            if (input$gene != "Type gene symbol") {
                tab_data <- tab_data %>% 
                    filter(Hugo_Symbol == input$gene | `Hugo_Symbol:1` == input$gene)
            }
        }
        colnames(tab_data) <- c(
            "TCGA Barcode", "MC3 gene", "PCAWG gene",
            "MC3 VAF", "PCAWG VAF", "MC3 var. Class", "MC3 HGVS short", " PCAWG var. class"
        )
        tab_data
    })
    
    output$downloadData <- downloadHandler(
        filename = function(){ paste("yourdata", ".csv", sep = "")},
        content = function(file) {write.csv(dat1(), file, row.names = FALSE)}
    )

    output$pid <- renderPlotly({ ##GD
        toPlot <- dat1() ##GD
        
        p <- ggplot(toPlot, aes(x=perc_PCAWG_in_MC3, y=perc_MC3_in_PCAWG)) + 
            geom_point(aes(text = paste("Sample ID:", sample,
                                        "\nMatch:", id_match,
                                        "\nMC3:", id_mc3,
                                        "\nPCAWG:", id_pcawg)), 
                       col="black") +
            scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = percent, limits = c(0,1)) + 
            scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = percent, limits = c(0,1)) +
            scale_color_manual(values = cbPalette) + 
            guides(color = guide_legend(reverse=TRUE)) + 
            theme_bw() + 
            theme(
                panel.border = element_blank(),
                axis.line.x = element_line(color="black"), 
                axis.line.y = element_line(color="black"),
                legend.position = c(0, 1), 
                legend.justification = c(0, 1)
            ) + 
            labs(x = "Matched/(Matched+MC3) variants", y = "Matched/(Matched+ICGC) variants") + 
            geom_vline(xintercept = .8, colour="#ff0066", linetype = "longdash") + 
            geom_hline(yintercept = .8, colour="#ff0066", linetype = "longdash")
        
        empty <- ggplot() + geom_blank() + theme_bw() + theme(panel.border = element_blank())
        hist_top <- ggplot(toPlot, aes(perc_PCAWG_in_MC3)) + 
            geom_histogram(fill="#5A80A6", breaks = seq(0, 1, 0.01)) +
            theme_bw() + 
            # remove axis ticks and tick mark labels 
            # REF: https://stackoverflow.com/questions/35090883/remove-all-of-x-axis-labels-in-ggplot
            theme(
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks = element_blank(),
                panel.border = element_blank()
            ) 
            
        hist_right <- ggplot(toPlot, aes(perc_MC3_in_PCAWG)) +
            geom_histogram(fill="#BE312D", breaks = seq(0, 1, 0.01)) +
            coord_flip()+
            theme_bw() + 
            theme(
                axis.title.x.top = element_text("Count"),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                panel.border = element_blank()
            ) 
            
        s <- subplot(
            hist_top, empty, p, hist_right, 
            nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), 
            shareX = TRUE, shareY = TRUE, titleX = TRUE, titleY = TRUE
        )
        ggplotly(s)
    })
    
    # Since dat() is generated from a reactive that 
    # is triggered by input$myslider this table will update
    # any time that input$myslider updates
    output$table <- DT::renderDataTable({
        datatable(dat3(), filter="top", selection="multiple", escape=FALSE, 
                  options = list(sDom  = '<"top">lrt<"bottom">ip'))
    })

    addTooltip(session, id = "ivaf_slider", title = "Select a variant allele fraction window for PCAWG vars", trigger = "hover")
    addTooltip(session, id = "mvaf_slider", title = "Select a variant allele fraction windows for MC3 vars", trigger = "hover")
    addTooltip(session, id = "filts_var", title = "Select a variant with a given MC3 filter (PASS = default)", trigger = "hover")
    addTooltip(session, id = "filts_sample", title = "Select samples with a given MC3 filter", trigger = "hover")
    addTooltip(session, id = "downloadData", title = "Download mutation data based on your filter criteria", trigger = "hover")
    
}

shinyApp(ui = ui, server = server)