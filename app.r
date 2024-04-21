install.packages("shiny")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("reticulate")
install.packages("bslib")
install.packages("thematic")
install.packages("shinythemes")
install.packages("reshape2")
install.packages("DT")
install.packages("Peptides")
library(shiny)
library(tidyverse)
library(ggplot2)
library(reticulate)
library(bslib)
library(thematic)
library(shinythemes)
library(reshape2)
library(DT)
library(Peptides)

virtualenv_create("r-reticulate")
virtualenv_install("r-reticulate", "modlamp")
virtualenv_install("r-reticulate", "pandas")
virtualenv_install("r-reticulate", "numpy")
use_virtualenv("r-reticulate")

modlamp <- import("modlamp")
pandas <-import("pandas")
numpy <-import("numpy")

df = read.csv("./data/master_file_duplicated_removed.csv")
df = data.frame(df)


amino_acid_sequence_validator <- function(sequence) {
  
  if (nchar(sequence) == 0) {
    stop("Empty Input") # or you can customize the message or behavior here
  }
  
  valid_amino_acids <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
                         "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  sequence <- toupper(sequence)
  sequence_chars <- strsplit(sequence, "")[[1]]
  all_valid <- all(sequence_chars %in% valid_amino_acids)
  
  return(all_valid)
}
descriptor_data_frame <- function(sequence) {
  if (!amino_acid_sequence_validator(sequence)) {
    stop("Invalid amino acid sequence.")
  }
  sequence <- toupper(sequence)
  result <- tryCatch({
    desc <- modlamp$descriptors$GlobalDescriptor(sequence)
    desc$calculate_all()
    descriptors <- as.list(desc$featurenames)
    values <- as.list(desc$descriptor)
    dff <- data.frame(Descriptors = unlist(descriptors), Values = unlist(values), stringsAsFactors = FALSE)
    
  }, error = function(e) {
    dff <- data.frame(Descriptors = "Error", Values = e$message, stringsAsFactors = FALSE)
  })
  return(result)
}
amino_acid_distribution_plot_function <- function(sequence){
  if(amino_acid_sequence_validator(sequence)){
    sequence <- toupper(sequence)
    desc = modlamp$descriptors$GlobalDescriptor(sequence)
    desc$count_aa()
    amino_acids = as.list((desc$featurenames))
    values = as.list((desc$descriptor))
    dff = data.frame(unlist(amino_acids),unlist(values))
    names(dff) = c("Amino_Acids","Values")
    amino_acid_distribution_plot = ggplot(dff,aes(x=Amino_Acids,y=Values)) + 
      geom_bar(stat = "identity", fill = "#0071e3") +
      theme_classic() +
      ylab("Frequency") +
      xlab("Amino Acids") +
      labs(title = "Amino Acid Distribution Plot")+
      theme(axis.title = element_text(face="bold", size = 12),
            plot.title = element_text(hjust=0.5,face="bold", size = 18))+
      theme(text=element_text(size=20))
  }
  else{
    amino_acid_distribution_plot = plot.new()
  }
  return(amino_acid_distribution_plot)
  
}
length_plot_function <- function(sequence){
  if(amino_acid_sequence_validator(sequence)){
    sequence <- toupper(sequence)
    desc = modlamp$descriptors$GlobalDescriptor(sequence)
    desc$length()
    length = desc$descriptor
    length = as.list((desc$descriptor))
    length = unlist(length)
    length_plot = ggplot(df)+
      geom_boxplot(aes(y=log(Length)), colour = "black", outlier.alpha = 0.3, outlier.colour = "orange") + 
      geom_point(aes(x=0,y=log(length)),colour = "blue") +
      theme_classic() +
      labs(title = "Length Box Plot") + 
      xlab(NULL)+
      ylab("log(Length)") +
      theme(axis.title = element_text(face="bold", size = 12),
            plot.title = element_text(hjust=0.5,face="bold", size = 18))+
      theme(text=element_text(size=20))
  }
  else{
    length_plot = plot.new()
  }
  return(length_plot)
}
charge_plot_function <- function(sequence){
  if(amino_acid_sequence_validator(sequence)){
    sequence <- toupper(sequence)
    desc = modlamp$descriptors$GlobalDescriptor(sequence)
    desc$calculate_charge(ph=7.4)
    charge = desc$descriptor
    charge = as.list((desc$descriptor))
    charge = unlist(charge)
    charge_plot = ggplot(df)+
      geom_boxplot(aes(y=Charge), colour = "black", outlier.alpha = 0.3, outlier.colour = "orange") + 
      geom_point(aes(y=charge,x=0),colour = "blue") +
      theme_classic() +
      labs(title = "Charge Box Plot") + 
      ylab("Charge")+
      xlab(NULL) +
      theme(axis.title = element_text(face="bold", size = 12),
            plot.title = element_text(hjust=0.5, face="bold", size = 18))+
      theme(text=element_text(size=20))
  }
  else{
    charge_plot = plot.new()
  }
  return(charge_plot)
}
instability_index_plot_function <- function(sequence){
  if(amino_acid_sequence_validator(sequence)){
    sequence <- toupper(sequence)  
    desc = modlamp$descriptors$GlobalDescriptor(sequence)
    desc$instability_index()
    instability_index = desc$descriptor
    instability_index = as.list((desc$descriptor))
    instability_index = unlist(instability_index)
    instability_index_plot = ggplot(df)+
      geom_boxplot(aes(y=InstabilityInd), colour = "black", outlier.alpha = 0.3, outlier.colour = "orange") + 
      geom_point(aes(y=instability_index,x=0),colour = "blue") +
      theme_classic() +
      labs(title = "Instability Index Box Plot") + 
      ylab("Instability Index")+
      xlab(NULL) +
      theme(axis.title = element_text(face="bold", size = 12),
            plot.title = element_text(hjust=0.5, face="bold", size = 18))+
      theme(text=element_text(size=20))
  }
  else{
    instability_index_plot = plot.new()
  }
  return(instability_index_plot)
}
boman_index_plot_function <- function(sequence){
  if(amino_acid_sequence_validator(sequence)){
    sequence <- toupper(sequence)
    desc = modlamp$descriptors$GlobalDescriptor(sequence)
    desc$boman_index()
    boman_index = desc$descriptor
    boman_index = as.list((desc$descriptor))
    boman_index = unlist(boman_index)
    boman_index_plot = ggplot(df)+
      geom_boxplot(aes(y=BomanInd), colour = "black", outlier.alpha = 0.3, outlier.colour = "orange") + 
      geom_point(aes(y=boman_index,x=0),colour = "blue") +
      theme_classic() +
      labs(title = "Boman Index Box Plot") + 
      ylab("Boman Index")+
      xlab(NULL) +
      theme(axis.title = element_text(face="bold", size = 12),
            plot.title = element_text(hjust=0.5,face="bold", size = 18))+
      theme(text=element_text(size=20))
  }
  else{
    boman_index_plot = plot.new()
  }
  return(boman_index_plot)
}
aliphatic_plot_function <- function(sequence){
  if(amino_acid_sequence_validator(sequence)){
    sequence <- toupper(sequence)  
    desc = modlamp$descriptors$GlobalDescriptor(sequence)
    desc$aliphatic_index()
    aliphatic = desc$descriptor
    aliphatic = as.list((desc$descriptor))
    aliphatic = unlist(aliphatic)
    aliphatic_plot = ggplot(df)+
      geom_boxplot(aes(y=AliphaticInd), colour = "black", outlier.alpha = 0.3, outlier.colour = "orange") + 
      geom_point(aes(y=aliphatic,x=0),colour = "blue") +
      theme_classic() +
      labs(title = "Aliphatic Box Plot") + 
      ylab("Aliphatic")+
      xlab(NULL) +
      theme(axis.title = element_text(face="bold", size = 12),
            plot.title = element_text(hjust=0.5,face="bold", size = 18))+
      theme(text=element_text(size=20))
  }
  else{
    aliphatic_plot = plot.new()
  }
  return(aliphatic_plot)
}
aromaticity_plot_function <- function(sequence){
  if(amino_acid_sequence_validator(sequence)){
    sequence <- toupper(sequence)
    desc = modlamp$descriptors$GlobalDescriptor(sequence)
    desc$aromaticity()
    aromaticity = desc$descriptor
    aromaticity = as.list((desc$descriptor))
    aromaticity = unlist(aromaticity)
    aromaticity_plot = ggplot(df)+
      geom_boxplot(aes(y=Aromaticity), colour = "black", outlier.alpha = 0.3, outlier.colour = "orange") + 
      geom_point(aes(y=aromaticity,x=0),colour = "blue") +
      theme_classic() +
      labs(title = "Aromaticity Box Plot") + 
      ylab("Aromaticity")+
      xlab(NULL) +
      theme(axis.title = element_text(face="bold", size = 12),
            plot.title = element_text(hjust=0.5, face="bold", size = 18))+
      theme(text=element_text(size=20))
  }
  else{
    aromaticity_plot = plot.new()
  }
  return(aromaticity_plot)
}
hydrophobic_ratio_plot_function <- function(sequence){
  if(amino_acid_sequence_validator(sequence)){
    sequence <- toupper(sequence)
    desc = modlamp$descriptors$GlobalDescriptor(sequence)
    desc$hydrophobic_ratio()
    hydrophobic_ratio = desc$descriptor
    hydrophobic_ratio = as.list((desc$descriptor))
    hydrophobic_ratio = unlist(hydrophobic_ratio)
    hydrophobic_ratio_plot = ggplot(df)+
      geom_boxplot(aes(y=HydrophRatio), colour = "black", outlier.alpha = 0.3, outlier.colour = "orange") + 
      geom_point(aes(y=hydrophobic_ratio,x=0),colour = "blue") +
      theme_classic() +
      labs(title = "Hydrophobic-Ratio  Box Plot") + 
      ylab("HydrophRatio")+
      xlab(NULL) +
      theme(axis.title = element_text(face="bold", size = 12),
            plot.title = element_text(hjust=0.5, face="bold", size = 18))+theme(text=element_text(size=20))
  }
  else{
    hydrophobic_ratio_plot = plot.new()
  }
  return(hydrophobic_ratio_plot)
}
random_amino_acid_mutations <- function(sequence, num_mutations = 1) {
  if (!amino_acid_sequence_validator(sequence)) {
    stop("Invalid amino acid sequence.")
  }
  sequence <- toupper(sequence)
  amino_acids <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  mutations_df <- data.frame(Position = integer(num_mutations),
                             Previous_AA = character(num_mutations),
                             New_AA = character(num_mutations),
                             Mutated_Sequence = character(num_mutations), stringsAsFactors = FALSE)
  for (i in 1:num_mutations) {
    position <- sample.int(nchar(sequence), 1)
    previous_aa <- substr(sequence, position, position)
    new_aa <- sample(amino_acids[amino_acids != previous_aa], 1)
    mutated_sequence <- paste0(substr(sequence, 1, position - 1),
                               new_aa,
                               substr(sequence, position + 1, nchar(sequence)))
    mutations_df[i, ] <- c(position, previous_aa, new_aa, mutated_sequence)
  }
  return(mutations_df)
}
mutate_peptide_by_group <- function(sequence, group) {
  if (!amino_acid_sequence_validator(sequence)) {
    stop("Invalid amino acid sequence.")
  }
  sequence <- toupper(sequence)
  amino_acid_groups <- list(
    polar = c("S", "T", "Y", "N", "Q"),
    nonpolar = c("A", "V", "L", "I", "P", "F", "M", "W"),
    positively_charged = c("K", "R", "H"),
    negatively_charged = c("D", "E"),
    aromatic = c("F", "Y", "W", "H"),
    small = c("A", "G", "C", "S", "T"),
    hydrophobic = c("A", "V", "L", "I", "P", "F", "M", "W", "G", "C"),
    charged = c("K", "R", "H", "D", "E")
  )
  if (!(group %in% names(amino_acid_groups))) {
    stop("Invalid amino acid group.")
  }
  group_amino_acids <- amino_acid_groups[[group]]
  mutated_aa <- sample(group_amino_acids, 1)
  position <- sample.int(nchar(sequence), 1)
  mutated_sequence <- paste0(
    substr(sequence, 1, position - 1),
    mutated_aa,
    substr(sequence, position + 1, nchar(sequence))
  )
  mutations_df = data.frame(mutated_sequence,c("Mutated Sequence"))
  names(mutations_df) = c("Mutated_Sequence", "Description")
  return(mutations_df)
}
muated_peptide_group_plot <- function(sequence,group){
  df = mutate_peptide_by_group(sequence,group)
  sequence = df$Mutated_Sequence[1]
  san = as.data.frame(aaComp(sequence))
  plot = barplot(san$Number,names.arg=rownames(san))
  plot9 <- barplot(san$Number,names.arg=rownames(san), 
                   xlab="Physiochemical Properties",
                   ylab="Frequency",col="#0071e3", 
                   main="Amino Acid Property Chart",
                   border="black",
  )
  return(plot9)
}
generate_mutations <- function(amino_acid_sequence) {
  if(amino_acid_sequence_validator(amino_acid_sequence)){
    amino_acids <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", 
                     "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
    mutations_df <- data.frame(
      Position = integer(),
      Original = character(),
      Mutant = character(),
      Mutant_Sequence = character(),
      stringsAsFactors = FALSE
    )
    for (pos in 1:nchar(amino_acid_sequence)) {
      original_aa <- substr(amino_acid_sequence, pos, pos)
      for (mutant_aa in amino_acids) {
        if (mutant_aa != original_aa) { 
          mutant_sequence <- substring(amino_acid_sequence, 1, pos-1)
          mutant_sequence <- paste0(mutant_sequence, mutant_aa, 
                                    substring(amino_acid_sequence, pos+1))
          mutations_df <- rbind(mutations_df, data.frame(
            Position = pos,
            Original = original_aa,
            Mutant = mutant_aa,
            Mutant_Sequence = mutant_sequence,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  else{
    mutations_df <- data.frame(Wrong = double(), 
                               Input = integer())
  }
  return(mutations_df)
}
mutation_df_error_handler <- function(sequence) {
  if (!amino_acid_sequence_validator(sequence)) {
    stop("Invalid amino acid sequence.")
  }
  else{
    return(NULL)
  }
}
data_frame_analyser <- function(sequence){
  if (nchar(sequence) == 0) {
    return(data.frame(Empty = double(), 
                      Input = integer(), 
                      stringsAsFactors = FALSE) )
  }
  else if (!amino_acid_sequence_validator(sequence)) {
    return(data.frame(Wrong = double(), 
                      Input = integer(), 
                      stringsAsFactors = FALSE) )
  }
  sequence <- toupper(sequence)
  df = generate_mutations(sequence)
  desc = modlamp$descriptors$GlobalDescriptor(as.list(df$Mutant_Sequence))
  desc$calculate_all()
  dff = pandas$DataFrame(data=desc$descriptor, columns=desc$featurenames)
  dff["Sequence"] = df$Mutant_Sequence
  dff = data.frame(dff)
  dff = subset(dff, select = -c(Length,MW,ChargeDensity,pI) )
  dff = dff %>% relocate(Sequence)
  return(dff)
}
filter_var <- function(x, val) {
  if (is.numeric(x)) {
    !is.na(x) & x >= val[1] & x <= val[2]
  } else if (is.factor(x)) {
    x %in% val
  } else {
    # No control, so don't filter
    TRUE
  }
}
eda_plot <- function(df){
  g = modlamp$analysis$GlobalAnalysis(list(df))
  plot = g$plot_summary()
  
  return(plot)
}
boxplots <- function(df){
  df = subset(df, select = -c(MW,ChargeDensity,pI) )
  data_long <- melt(df)
  plot <- ggplot(data_long, aes(x = variable, y = log(value), color = variable)) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.title = element_text(face="bold", size = 12),
          plot.title = element_text(hjust=0.5,face="bold", size = 18))+
    theme(text=element_text(size=20),legend.key.size = unit(1.5, "cm"),
          legend.key.width = unit(0.5,"cm") )
  return(plot)
}

ui <- fluidPage(
  theme = bslib::bs_theme(
    bootswatch = "lux"
  ),
  tags$div(
    navbarPage(
      includeCSS("www/style.css"),
      title = span(tags$u(tags$b("AVP-DE")), style = "padding: 3px; color:magenta"),
      tabPanel(
        "Home",
        id = "Home_Page",
        tags$div(
          tags$h1("The Antiviral Peptides Designer and Explorer"),
          class = "text-center"
        ),
        bslib::layout_columns(
          bslib::card(
            tags$div(
              tags$h3("Descriptors Tab"),
              tags$h4("Antiviral Peptide Validator:"),
              tags$p("The inclusion of box plots and descriptive statistics in the Descriptors Tab provides users with a visual understanding of the distribution and characteristics of various peptide parameters. This feature is crucial for rational-based peptide design, allowing researchers to identify patterns that may contribute to peptide stability and activity.", style = "text-align:justify"),
              style = "margin: 2%"
            ),
            style = "margin: 2%"
          ),
          bslib::card(
            tags$div(
              tags$h3("Mutator Tab"),
              tags$h4("Antiviral Peptide Mutator:"),
              tags$p("The Mutator Tab offers a dynamic approach to peptide mutation, enabling users to generate random mutations, and group-based mutations, and explore all possible mutations. The interactive amino acid bar chart further enhances the understanding of how mutations affect peptide properties, aiding in the design of novel peptides with desired characteristics.", style = "text-align:justify"),
              style = "margin: 2%"
            ),
            style = "margin: 2%"
          ),
          bslib::card(
            tags$div(
              tags$h3("Analyzer Tab"),
              tags$h4("Antiviral Peptide Analyzer:"),
              tags$p("The Analyzer Tab's ability to provide a data frame for each possible mutation and render it in an interactive table format is a significant asset. Researchers can utilize dynamic filtering based on descriptors to narrow down sequences of interest, facilitating in-depth exploration and analysis of peptide variants.", style = "text-align:justify"),
              style = "margin: 2%"
            ),
            style = "margin: 2%"
          ),
          bslib::card(
            tags$div(
              tags$h3("EDA Tab"),
              tags$h4("Peptide Exploratory Data Analysis (EDA):"),
              tags$p("The EDA Tab plays a pivotal role in understanding the overall trends and characteristics of peptide sequences. By generating box plots and EDA plots, researchers gain valuable insights into the distribution of descriptor values and can identify outliers or uncommon patterns that may be of interest for further investigation.", style = "text-align:justify"),
              style = "margins: 2%"
            ),
            style = "margin: 2%"
          ),
          style = "margin: 2%"
        ),
      ),
      tabPanel(
        "Descriptors",
        tags$div(
          titlePanel("AntiViral Peptide Validator"),
          tags$hr(),
          class = "text-center"
        ),
        tags$div(
          sidebarLayout(
            sidebarPanel(
              fluidRow(
                textInput(
                  "sequence",
                  "Enter Your Sequence:",
                  "HGVSGHGQHGVHG",
                  width = "100%"
                )
              ),
              fluidRow(
                textOutput("Descriptor")
              ),
              fluidRow(
                tableOutput("Descriptor_table")
              ),
              fluidRow(
                textOutput("Descriptor_download")
              ),
              fluidRow(
                downloadButton("Download_descriptor_table", "Download Data")
              ),
              width = 3
            ),
            mainPanel(
              fluidRow(
                wellPanel(
                  tabsetPanel(
                    id = "tabset",
                    type = "pills",
                    tabPanel("Distribution",
                             fluidRow(
                               tags$div(
                                 plotOutput("plot1", height = "450px"),
                                 style = "width=98%"
                               )
                             ),
                             tags$div(
                               fluidRow(
                                 textOutput("plot1_download")
                               ),
                               class = "mt-3 mb-3 text-center",
                             ),
                             fluidRow(
                               downloadButton("Download_Plot1", "Download Plot")
                             ),
                    ),
                    tabPanel("Length",
                             plotOutput("plot2", height = "450px"),
                             tags$div(
                               fluidRow(
                                 textOutput("plot2_download")
                               ),
                               class = "mt-3 mb-3 text-center",
                             ),
                             fluidRow(
                               downloadButton("Download_Plot2", "Download Plot")
                             )
                    ),
                    tabPanel("Charge",
                             plotOutput("plot3", height = "450px"),
                             tags$div(
                               fluidRow(
                                 textOutput("plot3_download")
                               ),
                               class = "mt-3 mb-3 text-center",
                             ),
                             fluidRow(
                               downloadButton("Download_Plot3", "Download Plot")
                             ),
                    ),
                    tabPanel("Instability",
                             plotOutput("plot4", height = "450px"),
                             tags$div(
                               fluidRow(
                                 textOutput("plot4_download")
                               ),
                               class = "mt-3 mb-3 text-center",
                             ),
                             fluidRow(
                               downloadButton("Download_Plot4", "Download Plot")
                             ),
                    ),
                    tabPanel("BomanIndex",
                             plotOutput("plot5", height = "450px"),
                             tags$div(
                               fluidRow(
                                 textOutput("plot5_download")
                               ),
                               class = "mt-3 mb-3 text-center",
                             ),
                             fluidRow(
                               downloadButton("Download_Plot5", "Download Plot")
                             ),
                    ),
                    tabPanel("Aliphaticity",
                             plotOutput("plot6", height = "450px"),
                             tags$div(
                               fluidRow(
                                 textOutput("plot6_download")
                               ),
                               class = "mt-3 mb-3 text-center",
                             ),
                             fluidRow(
                               downloadButton("Download_Plot6", "Download Plot")
                             ),
                    ),
                    tabPanel("Aromaticity",
                             plotOutput("plot7", height = "450px"),
                             tags$div(
                               fluidRow(
                                 textOutput("plot7_download")
                               ),
                               class = "mt-3 mb-3 text-center",
                             ),
                             fluidRow(
                               downloadButton("Download_Plot7", height = "450px", "Download Plot")
                             ),
                    ),
                    tabPanel("Hydrophobicity",
                             plotOutput("plot8", height = "450px"),
                             tags$div(
                               fluidRow(
                                 textOutput("plot8_download")
                               ),
                               class = "mt-3 mb-3 text-center",
                             ),
                             fluidRow(
                               downloadButton("Download_Plot8", "Download Plot")
                             ),
                    ),
                  ),
                  width = "9"
                ),
              ),
              width = "9"
            )
          )
        ),
        tags$br(),
        tags$hr()
      ),
      tabPanel(
        "Mutator",
        tags$div(
          titlePanel("AntiViral Peptide Mutator"),
          class = "text-center",
          tags$hr()
        ),
        sidebarLayout(
          sidebarPanel(
            tags$div(
              textInput(
                "Mutator",
                "Enter Your Sequence:",
                width = "100%",
                value = "HGVSGHGQHGVHG",
              ),
            ),
            tags$div(
              numericInput(
                "Number_of_mutation",
                "Number of Random Muations:",
                value = 1,
                min = 1,
                max = 10,
                width = "100%",
              )
            ),
            tags$div(
              tableOutput("Mutated_table")
            ),
            width = 4
          ),
          mainPanel(
            fluidRow(
              tags$div(
                fluidRow(
                  column(7,
                         wellPanel(
                           fluidRow(
                             textInput(
                               "Mutator2",
                               "Enter Your Sequence:",
                               width = "100%",
                               value = "HGVSGHGQHGVHG",
                             ),
                           ),
                           fluidRow(
                             selectInput(
                               "Selecet_group",
                               label = "Select Group To Add As Mutation:",
                               choices = c("polar", "nonpolar", "positively_charged",
                                           "negatively_charged", "aromatic", "small",
                                           "hydrophobic", "charged"),
                               width = "100%"
                             )
                           ),
                           fluidRow(
                             tableOutput(
                               "Mutated_sequence"
                             ),
                           ),
                           fluidRow(
                             plotOutput(
                               "plot9",
                               width = "97%",
                             )
                           )
                         )
                  ),
                  column(5,
                         wellPanel(
                           fluidRow(
                             textInput(
                               "Mutator3",
                               "Enter Your Sequence:",
                               width = "100%",
                               value = "HGVSGHGQHGVHG",
                             ),
                           ),
                           fluidRow(
                             textOutput(
                               "Mutation_Download_df"
                             ),
                             class = "mb-3"
                           ),
                           fluidRow(
                             downloadButton(
                               "Download_Mutation_df",
                               "Download Data"
                             ),
                             class = "mb-3"
                           ),
                           fluidRow(
                             tableOutput("Mutation_df_error")
                           )
                         )
                  )
                ),
              ),
            ),
          ),
        ),
      ),
      tabPanel(
        "Analyzer",
        tags$div(
          titlePanel("AntiViral Peptide Analyzer"),
          class = "text-center",
          tags$hr()
        ),
        sidebarLayout(
          sidebarPanel(
            textInput("Sequence_analyzer", "Enter Your Sequence:", value = "HGVSGHGQHGVHG", width = "100%"),
            sliderInput("Charge", label = "Charge: ", min = -8, max = 9, value = c(-8,9)),
            sliderInput("InstabilityInd", label = "Instability Index: ", min = -56, max = 130, value = c(-56,130)),
            sliderInput("Aromaticity", label = "Aromaticity: ", min = 0, max = 0.3, value = c(0.,0.3)),
            sliderInput("AliphaticInd", label = "Aliphatic Index: ", min = -23, max =210, value = c(-23,210)),
            sliderInput("BomanInd", label = "Boman Index: ", min = -3, max = 6, value = c(-3,6)),
            sliderInput("HydrophRatio", label = "Hydrophobic Ratio: ", min = 0., max =0.8, value = c(0,0.8))
          ),
          mainPanel(
            wellPanel(
              DT::dataTableOutput("data") 
            )
          )
        )
      ),
      tabPanel(
        "EDA",
        tags$div(
          titlePanel("Peptide Exploratory Data Analysis"),
          class = "text-center",
          tags$hr()
        ),
        fluidRow(
          column(3,
                 wellPanel(
                   tags$div(h3("EDA Panel", style = "color: #f5f5f7;"), class = "text-center", style = "background-color:  rgba(22, 23, 23, .7); margin-right: 75px; margin-left: 75px;border-radius:10px 15px;"),
                   fileInput("upload", "Choose CSV File", accept = ".csv", buttonLabel = "Upload..."),
                 )
          ),
          column(9,
                 tags$div(
                   wellPanel(
                     tabsetPanel(
                       id = "tabsetEDA",
                       type = "pills",
                       tabPanel(
                         "Descriptors",
                         tags$div(
                           fluidRow(
                             DT::dataTableOutput("descriptor_data_head")
                           ),
                           fluidRow(
                             textOutput(
                               "Descriptor_download_EDA",
                             )
                           ),
                           class = "mt-3 mb-3 text-center",
                         ),
                         fluidRow(
                           downloadButton(
                             "Download_descriptor_table_EDA",
                             "Download Data"
                           ),
                         ),
                       ),
                       tabPanel(
                         "Summary",
                         tags$div(
                           fluidRow(
                             DT::dataTableOutput("descriptor_data_summary")
                           ),
                           fluidRow(
                             textOutput(
                               "Descriptor_download_EDA_summary",
                             )
                           ),
                           class = "mt-3 mb-3 text-center",
                         ),
                         fluidRow(
                           downloadButton(
                             "Download_descriptor_table_EDA_summary",
                             "Download Data"
                           ),
                         ),
                       ),
                       tabPanel(
                         "EDAPlot",
                         tags$div(
                           fluidRow(
                             plotOutput("eda_plot", width = "93%", height = "700") 
                           ),
                           class = "ms-4 ps-3"
                         ),
                         tags$div(
                           fluidRow(
                             textOutput(
                               "plot10_download",
                             )
                           ),
                           class = "mt-3 mb-3 text-center",
                         ),
                         fluidRow(
                           downloadButton(
                             "Download_Plot10",
                             "Download Plot"
                           )
                         ),
                       ),
                       tabPanel(
                         "BoxPlots",
                         tags$div(
                           fluidRow(
                             plotOutput("box_plots", width = "98%") 
                           ),
                         ),
                       )
                     )
                   )
                 )
          ),
        )
      )
    ),
    class = "font-monospace"
  )
)

server <- function(input, output, session){
  thematic::thematic_shiny()
  output$Descriptor = renderText("Descriptor Values:")
  output$Descriptor_table = renderTable(
    descriptor_data_frame(input$sequence),
    options = list(pageLength=10),
    width = "100%"
  )
  output$Descriptor_download = renderText("Download Your Descriptors:")
  output$Download_descriptor_table = downloadHandler(
    filename = function() {
      paste("Descriptors_data_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(descriptor_data_frame(input$sequence), row.names = FALSE,file)
    }
  )
  output$plot1 = renderPlot(amino_acid_distribution_plot_function(input$sequence))
  output$plot1_download = renderText("Download Your Amino Acid Distribution Plot:")
  output$Download_Plot1 <- downloadHandler(
    filename = function() {
      paste("Amino_Acid_Distribution_Plot", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      ggsave(file,amino_acid_distribution_plot_function(input$sequence), dpi = 300)
    }
  )
  output$plot2 = renderPlot(length_plot_function(input$sequence))
  output$plot2_download =  renderText("Download Your Length Plot:")
  output$Download_Plot2 <- downloadHandler(
    filename = function() {
      paste("Length_Plot_", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      ggsave(file,length_plot_function(input$sequence), dpi = 300)
    }
  )
  output$plot3 = renderPlot(charge_plot_function(input$sequence))
  output$plot3_download = renderText("Download Your Charge Plot:")
  output$Download_Plot3 <- downloadHandler(
    filename = function() {
      paste("Charge_Plot_", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      ggsave(file,charge_plot_function(input$sequence), dpi = 300)
    }
  )
  output$plot4 = renderPlot(instability_index_plot_function(input$sequence))
  output$plot4_download =renderText("Download Your InstabilityIndex Plot:")
  output$Download_Plot4 <- downloadHandler(
    filename = function() {
      paste("Instability_Index_Plot_", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      ggsave(file,instability_index_plot_function(input$sequence), dpi = 300)
    }
  )
  output$plot5 = renderPlot(boman_index_plot_function(input$sequence))
  output$plot5_download = renderText("Download Your BomanIndex Plot:")
  output$Download_Plot5 <- downloadHandler(
    filename = function() {
      paste("Boman_Index_Plot_", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      ggsave(file,boman_index_plot_function(input$sequence), dpi = 300)
    }
  )
  output$plot6 = renderPlot(aliphatic_plot_function(input$sequence))
  output$plot6_download =renderText("Download Your Aliphatic Plot:")
  output$Download_Plot6 <- downloadHandler(
    filename = function() {
      paste("Aliphatic_Plot_", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      ggsave(file,aliphatic_plot_function(input$sequence), dpi = 300)
    }
  )
  output$plot7 = renderPlot(aromaticity_plot_function(input$sequence))
  output$plot7_download = renderText("Download Your Aromaticity Plot:")
  output$Download_Plot7 <- downloadHandler(
    filename = function() {
      paste("Aromaticity_Plot_", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      ggsave(file,aromaticity_plot_function(input$sequence), dpi = 300)
    }
  )
  output$plot8 = renderPlot(hydrophobic_ratio_plot_function(input$sequence))
  output$plot8_download = renderText("Download Your Hydrophobicity Plot:")
  output$Download_Plot8 <- downloadHandler(
    filename = function() {
      paste("Hydrophobic_Ratio_Plot_", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      ggsave(file,hydrophobic_ratio_plot_function(input$sequence), dpi = 300)
    })
  output$Mutated_table =renderTable(
    random_amino_acid_mutations(
      input$Mutator, input$Number_of_mutation),
    width = "100%")
  
  output$Mutated_sequence =renderTable(
    mutate_peptide_by_group(input$Mutator2, input$Selecet_group),
    width = "100%"
  )
  output$plot9 = renderPlot(muated_peptide_group_plot(input$Mutator2,input$Selecet_group))
  output$Mutation_Download_df = renderText("Download Your All Possible Mutations:")
  output$Mutation_df_error = renderTable(mutation_df_error_handler(input$Mutator3))
  output$Download_Mutation_df = downloadHandler(
    filename = function() {
      paste("All_Possible_Mutation_data_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(generate_mutations(input$Mutator3), row.names = FALSE,file)
    }
  )
  filtering_data <- reactive({
    data_frame_analyser(input$Sequence_analyzer)
  })
  selected <- reactive({
    data <- filtering_data()
    filter_condition <- filter_var(data$Charge, input$Charge) &
      filter_var(data$InstabilityInd, input$InstabilityInd) &
      filter_var(data$Aromaticity, input$Aromaticity) &
      filter_var(data$AliphaticInd, input$AliphaticInd) &
      filter_var(data$BomanInd, input$BomanInd) &
      filter_var(data$HydrophRatio, input$HydrophRatio)
    data[filter_condition, ]
  })
  output$data <- DT::renderDataTable({
    DT::datatable(selected(), options = list(scrollX = TRUE))
  })
  output$Smiles = renderText(Peptides::aaSMILES(input$Sequence_EDA))
  output$Bonding = renderText(data.frame(value = (unlist(Peptides::crucianiProperties(input$Sequence_EDA))))$value[3])
  file_1 <- reactive({
    req(input$upload)
    df <- unlist(as.list(read.csv(input$upload$datapath,header=T)))
    rownames(df) <- NULL
    colnames(df) <- NULL
    desc <- modlamp$descriptors$GlobalDescriptor(df)
    desc$calculate_all()
    descriptors <- as.list(desc$featurenames)
    values <- desc$descriptor
    df = data.frame(values)
    names(df) = (c(descriptors))
    return(df)
  })
  file_plot_1 <- reactive({
    req(input$upload)
    df <- read.csv(input$upload$datapath, header = TRUE)
    df  = as.list(df$SEQUENCE)
    return(df)
  })
  output$descriptor_data_head <- DT::renderDataTable({
    DT::datatable(file_1(), options = list(scrollX = TRUE), class = 'table-xl')
  })
  output$Descriptor_download_EDA = renderText("Download Your Descriptors EDA:")
  output$Download_descriptor_table_EDA = downloadHandler(
    filename = function() {
      paste("Descriptors_data_EDA", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(file_1(), row.names = FALSE,file)
    }
  )
  output$descriptor_data_summary <- DT::renderDataTable({
    description <-  summary_df <- as.data.frame(summary(file_1()))
    transposed_summary <- t(summary_df)
    transposed_df <- as.data.frame(transposed_summary, stringsAsFactors = FALSE)
    DT::datatable(transposed_df, options = list(scrollX = TRUE), rownames = FALSE, colnames = NULL, width = "100%", class = 'cell-border stripe table-xl')
  })
  output$Descriptor_download_EDA_summary = renderText("Download Your Descriptors EDA Summry:")
  output$Download_descriptor_table_EDA_summary = downloadHandler(
    filename = function() {
      paste("Descriptors_data_EDA_summary", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(summary(file_1()), row.names = FALSE,file)
    }
  )
  output$eda_plot = renderPlot(eda_plot(file_plot_1()))
  output$plot10_download = renderText("Download Your EDA Plot")
  output$Download_Plot10 <- downloadHandler(
    filename = function() {
      paste("EDA_Plot", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      png(filename = file, width = 900, height = 700)
      plot <- eda_plot(file_plot_1())
      dev.off()  # Close the device after plotting
    }
  )
  output$box_plots = renderPlot(boxplots(file_1()))
}

shinyApp(ui,server);
