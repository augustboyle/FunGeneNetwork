library(tidyverse)
library(shiny)
library(igraph)
library(visNetwork)
library(RColorBrewer)
library(shinythemes)
library(ComplexHeatmap)
library(circlize)
# markdown library is required

# Paths to network data
all_edges_rds = "cofunctional_network.edges.20181211.rds"
all_nodes_tsv = "cofunctional_network.genes.20181211.tsv.gz"
mat_scaled_rds = "cofunctional_scaled_matrix.20181211.rds"
screen_data_dir = "screen_data"

ui = fluidPage(
      theme = shinytheme("simplex"),
      # App title
      titlePanel("Co-functional network settings"),
      sidebarLayout(
         sidebarPanel(
          
          # choose query gene
          uiOutput("gene_query"),
          
          # layout of genes in network
          selectInput("net_layout", h4("Select layout"),
            choices = c("layout_randomly", "layout_in_circle", "layout_with_fr", "layout_with_lgl", "layout_with_kk"), selected = "layout_in_circle"),
          
          # preset metrics (files can be added to this list)
          selectInput("metric", h4("Gene score"),
            choices = c(
              "Custom",
              "KO_mean", 
              "strength",
              "degree",
              "community_score",
              "community_size",
              "Growth - CRISPRa (Jost 2017)",
              "Growth - CRISPRi (Jost 2017)",
              "Rigosertib - CRISPRa (Jost 2017)",
              "Rigosertib - CRISPRi (Jost 2017)"
              ),
            selected = "KO_mean"),
          
          # option to add your own file with 'gene_symbol' and 'score' columns in your browser
          fileInput("custom_data", "Upload gene scores",
                    multiple = FALSE
                    ),
          
          # gray out nodes that are part of a different community from the query gene
          checkboxInput("community_mask", "Mask by community", FALSE),
          
          # set aesthetics for number of genes, transparency, layout seed and edge threshold
          splitLayout(cellWidths = c("22%", "22%", "22%", "34%"), 
            numericInput("n_genes", h5("# genes"), 20, min = 1, max = 100, step = 1),      
            numericInput("alpha", h5("alpha"), 0.8, min = 0.0, max = 1, step = 0.01),
            numericInput("seed", h5("seed"), 1, min = 1, max = 100, step = 1),
            numericInput("min_cor", h5("min |R|"), 0.5, min = 0.0, max = 1, step = 0.01)
          ),

          # layer STRING database on top
          splitLayout(cellWidths = c("50%", "50%"), 
            selectizeInput("string_anno", h5("STRING interaction"), 
              choices = c("neighborhood",            
                "neighborhood_transferred", "fusion",                  
                "cooccurence",              "homology",                
                "coexpression",             "coexpression_transferred",
                "experiments",              "experiments_transferred", 
                "database",                 "database_transferred",    
                "textmining",               "textmining_transferred",  
                "combined_score"), 
              selected = "combined_score"),
            numericInput("string_thresh", h5("STRING confidence"), 250, min = 0, max = 999, step = 1)
          )
         ),

        # set dimensions of interactive plot
        mainPanel(      
  tabsetPanel(
    tabPanel("README",
          includeMarkdown("README.md")),
    tabPanel("Network visualization",
          visNetworkOutput("networkPlot", width = 750, height = 550)
        ),
    tabPanel("Heatmap visualization", 
          plotOutput("heatmap", width = 750, height = 550))
      )
    )
  )
)

server = function(input, output) {

  all_edges = readRDS(all_edges_rds) 
  all_nodes = read_tsv(all_nodes_tsv)
  mat_scaled = readRDS(mat_scaled_rds)

  # rust, white, blue
  node_colors = c("#F98400", "#ffffff", "#5BBCD6")
  # beige, white, steel
  edge_colors = c("#83889e", "#ffffff", "#e2ca96")
  # blue, white, gold
  fitness_colors = c("#5BBCD6", "#ffffff", "#F2AD00")
  # green, white, rust
  correlation_colors = c("#00A08A", "#ffffff", "#F98400")
  # black, white, red
  score_colors = c("#000000", "#ffffff", "#FF0000")

  all_genes = c(all_edges$gene_symbol_1,all_edges$gene_symbol_2) %>% unique

  output$gene_query = renderUI({
     selectizeInput("query_gene", h3("Query gene"), choices = all_genes, selected = "TP53") 
  })


  # acquire nodes co-functional to query gene
  query_nodes = reactive( {
      query_gene = input$query_gene
      connected_genes = bind_rows(
        tibble(gene_symbol = query_gene, cor = 1),
        filter(all_edges, gene_symbol_1 == query_gene) %>% select(gene_symbol_2, cor) %>% rename(gene_symbol = gene_symbol_2), 
        filter(all_edges, gene_symbol_2 == query_gene) %>% select(gene_symbol_1, cor) %>% rename(gene_symbol = gene_symbol_1)
        ) %>% arrange(abs(cor) %>% desc)
      
      # add gene annotations
      query_nodes = left_join(connected_genes, all_nodes) %>% rename(id = gene_symbol) %>% mutate(q_gene = query_gene == id)
      query_nodes
    })

  # get data for network aesthetics
  query_screen = reactive( {
      metric = input$metric
      custom_pheno_file = input$custom_data$datapath
      nodes = query_nodes()
      
      # use a provided metric
      if(metric %in% c("KO_mean", "strength", "degree", "community_score", "community_size"))
        query_screen = (nodes %>% rename_("score" = metric)) else
      
      # upload your own file with 'gene_symbol' and 'score' columns
      if(metric == "Custom")
        query_screen = read_tsv(custom_pheno_file) %>% left_join(nodes, ., by=c("id" = "gene_symbol")) else
        query_screen = read_tsv(paste0(screen_data_dir,"/", metric, ".tsv.gz")) %>% left_join(nodes, ., by=c("id" = "gene_symbol"))
      query_screen
    })

  # use visNetwork to plot
  output$networkPlot = renderVisNetwork({
      n_genes = input$n_genes
      query_screen = query_screen()
      query_gene = input$query_gene
      nodes = query_nodes() %>% mutate(
        value = with(query_screen, ifelse(is.na(score), median(score %>% abs, na.rm=TRUE), abs(score))),
        color = ifelse(is.na(query_screen$score), node_colors[2], node_colors[2 + sign(query_screen$score)]),
        shape = ifelse(q_gene, "diamond", "dot")) %>% head(n_genes)
      
      # Gray out genes that are members of different communities
      if(input$community_mask){
        query_community = nodes %>% filter(q_gene) %>% pull(community)
        nodes$color = ifelse((nodes$community != query_community) | is.na(nodes$community), "#bdbdbd", nodes$color)
      }
      net_layout = input$net_layout
      string_anno = input$string_anno
      string_thresh = input$string_thresh    
      opacity = input$alpha
      seed = input$seed
      min_cor = input$min_cor
      font = list(size=25, color="black")
      edges = filter(all_edges, gene_symbol_1 %in% nodes$id, gene_symbol_2 %in% nodes$id, abs(cor) > min_cor) %>% rename(from=gene_symbol_1, to = gene_symbol_2)
      edges$value = edges$cor**2
      edges$color = ifelse(is.na(edges$cor), edge_colors[2], edge_colors[2 + sign(edges$cor)]) 
      edges$dashes = ifelse((pull(edges, string_anno) > string_thresh) %>% is.finite, "false", "[10,20]")
      visNetwork(nodes, edges, width = "100%") %>% visIgraphLayout(net_layout, randomSeed = seed) %>% visNodes(font=font) %>% 
        visEdges(color = list(opacity = opacity), scaling = list(min=3.5)) %>% 
        visOptions(highlightNearest = list(enabled =TRUE, degree = 1, hover = T)) %>% 
        visExport(type = "png", name = query_gene)
    })

  output$heatmap = renderPlot({
    query_gene = input$query_gene
    n_genes = input$n_genes
    nodes = query_nodes()
    query_screen = query_screen()

    # pull relevant metric out for gene score and set range
    metric_matrix = as.matrix(select(query_screen, id, score) %>% head(n_genes) %>% column_to_rownames("id"))
    metric_max = max(metric_matrix, 0, na.rm = TRUE)
    metric_min = min(metric_matrix, 0, na.rm = TRUE)
    
    # use ComplexHeatmap to generate heatmap visualization
    Heatmap(mat_scaled[nodes$id %>% head(n_genes),], name = "KO phenotype", col = colorRamp2(c(-2, 0, 2), fitness_colors), 
      show_column_dend = FALSE,
      # top_annotation = ha, top_annotation_height = unit(22 / 2, "mm"), 
      show_row_names = FALSE, show_column_names = FALSE,
      clustering_distance_rows = function(a,b) 1 - abs(cor(a,b)),
      clustering_distance_columns = function(a,b) 1 - cor(a,b),
      clustering_method_columns = "mcquitty", 
      clustering_method_rows = "ward.D2") +
    Heatmap(metric_matrix, name = "Gene score", show_row_names = FALSE, width = unit(5, "mm"), col = colorRamp2(c(metric_min, 0, metric_max), score_colors),
      column_names_gp = gpar(fontsize = 10)) +
    Heatmap(as.matrix(select(query_screen, id, cor) %>% head(n_genes) %>% column_to_rownames("id")), name = "Correlation", show_row_names = TRUE, width = unit(15 / 3, "mm"), col = colorRamp2(c(-1, 0, 1), correlation_colors),
      clustering_distance_columns = function(a,b) 1 - cor(a,b),
      clustering_method_columns = "mcquitty",
      column_names_gp = gpar(fontsize = 10),
      row_names_gp = gpar(fontsize = 10)) 

  })  
}


shinyApp(ui = ui, server = server)