# plot mokapot output
library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(tictoc)
library(ggsci)
library(purrr)
library(scales)



#======================= INPUT PARAMETERS ==========================

args <- commandArgs()
print(args)

input_file = args[6]
file_root = args[7]

#------------ Define dirs ------------
project_dir <- "/mmfs1/gscratch/stergachislab/mwperez/ctcf-footprinting"

#ARGS
data_folder <- sprintf("%s/candidate_footprints", project_dir)
mokapot_dir <- sprintf("%s/mokapot_res", project_dir)
output_folder <- sprintf("%s/figures", mokapot_dir)
cat(data_folder, "\n")
cat(mokapot_dir, "\n")
cat(output_folder, "\n")

#---------- HELPER FUNCTIONS ----------

# Initiate PDF
openPDF <- function (output_file_fh, pdf_dims) {
  cat("Saving plots to: ", output_file_fh, "\n")
  # Open PDF
  pdf(output_file_fh, width = pdf_dims[1], height = pdf_dims[2],
      useDingbats=FALSE, family = "ArialMT")
}

# Close PDF
closePDF <- function () {
  while (!is.null(dev.list()))  dev.off()
  cat("File saved.\n")
}

#---------- PLOT & COLOR ARGS ----------

font_size <- 12
font_info <- element_text(size=font_size, family="ArialMT", color="black")

# define plot feature colors
# light grey, purple, darker grey, blue
vline_col <- "#CECDCC"
m6a_col <- "#7F3F98"
fiber_col <- "#D1D2D4"
motif_annot_col <- "#006738"

# colors from Mitchell's presentation
# orange, red, purple, light grey, navy (gene track)
mitch_cols <- c("#FF8C00", "#FF0000", "#9370DB", "#E6E6E6", "#0C0C78")

#---------- FORMAT DATA ----------

# load data file
data_file <- sprintf("%s/%s", mokapot_dir, input_file)
if (!file.exists(data_file)) {stop("File does not exist: %s", data_file)}

# read table
df <- fread(data_file)
cat(sprintf("m6a rows: %s", format(dim(df)[1], big.mark=",", scientific=FALSE)), "\n")

df <- df %>% mutate(Label = ifelse(Label == 1, "Positive", "Negative"))

# unique CTCF motifs, queries, & motif_queries
cat(sprintf("unique motifs: %s", format(length(unique(df$motif_name)), big.mark=",", scientific=FALSE)), "\n")
cat(sprintf("unique queries: %s", format(length(unique(df$query_name)), big.mark=",", scientific=FALSE)), "\n")
cat(sprintf("unique motif-query groups: %s", format(length(unique(df$motif_query)), big.mark=",", scientific=FALSE)), "\n")


#---------- AGG M6A DENSITY ----------

# output file
output_file_name <- sprintf("CTCF_%s.mokapot.m6a_fiberseq_density.pdf", file_root)
output_file <- sprintf("%s/%s", output_folder, output_file_name)
cat("Saving m6a density plots to: ", output_file, "\n")
pdf_dims <- c(10, 5)

# plot m6a density grouped by FDR
plot_m6a_density <- function(df, FDR) {

    df$FDR_group <- ifelse(df$FDR < FDR, "low_FDR", "high_FDR")
    df$FDR_group <- factor(df$FDR_group, levels=c("low_FDR", "high_FDR"))
    FDR_str <- strsplit(as.character(FDR), ".", fixed=TRUE)[[1]][2]
    
    # set limits
    y_limits <- c(0, NA)
    x_limits <- c(-100, 100)
    plot_title <- sprintf("Aggregate m6a density (FDR < %s) - %s", FDR, "L Pos & 5% Neg")
    plot_subtitle <- sprintf("FDR low: %s | FDR high: %s", 
                        format(table(df$FDR_group)["low_FDR"], big.mark=",", scientific=FALSE), 
                        format(table(df$FDR_group)["high_FDR"], big.mark=",", scientific=FALSE))
    x_axis_name <- "distance from motif start"
    y_axis_name <- "aggregate m6a methylation density"
    legend_title <- "FDR"
    legend_labels <- c(sprintf("FDR < %s", FDR), sprintf("FDR > %s", FDR))
    
    # color values - low FDR (tf w/o m6a's) red | high FDR (accessible w/ m6a's) blue
    df$FDR_colors <- ifelse(df$FDR_group == "low_FDR", "#FF0000", "#0C0C78")
    df$FDR_colors <- factor(df$FDR_colors, levels=c("#FF0000", "#0C0C78"))
    motif_annot_col <- "#006738"
    vline_col <- "#CECDCC"
    font_size <- 12
    font_info <- element_text(size=font_size, family="ArialMT", color="black")
    
    # plot
    p <- ggplot(df, aes(x=centered_start, group=FDR_group)) +
            # add box over motif location
            annotate("rect", xmin=0, xmax=35, ymin=-Inf, ymax=+Inf, alpha=0.1, fill=motif_annot_col) +
            # add vertical line at x=0
            geom_vline(xintercept=0, color=vline_col, show.legend=FALSE, linetype="dashed", linewidth=1) +
            # density plot
            geom_density(aes(col=FDR_colors), linewidth=1.2, adjust=0.5) +
            scale_color_identity(guide="legend", legend_title, 
                                 labels=legend_labels, breaks=levels(df$FDR_colors)) +
            #scale_color_manual(name="FDR", values=df$FDR_colors, labels=legend_labels) +
            scale_y_continuous(limits=y_limits, expand=expansion(mult=c(0,0.05)), name=y_axis_name,
                               labels = function(x) format(x, scientific=FALSE), 
                               guide=guide_axis(check.overlap=TRUE)) +
            ggtitle(label = plot_title, subtitle = plot_subtitle) +
            scale_x_continuous(name=x_axis_name) +
            theme_classic() +
            theme(text = font_info,
                  axis.ticks = element_line(color="#000000", lineend="square"),
                  axis.text = font_info,
                  panel.background = element_rect(fill="transparent", color="#000000", linewidth=1),
                  plot.title = element_text(family="ArialMT", size=font_size+2, hjust=0.5),
                  plot.subtitle = element_text(family="ArialMT", size=font_size+2, hjust=0.5))
    
    return(p)
    
}

# FDR values to plot
FDR_vals <- c(0.05, 0.01, 0.001)

openPDF(output_file, pdf_dims)

for (FDR in FDR_vals) {
    cat("Making agg. m6a density plot for FDR: ", FDR, "\n")
    p <- plot_m6a_density(subset(df, Label="Positive"), FDR=FDR)
    print(p)
    cat("\n")
}

closePDF()

print("Done with aggregate m6A density!")


