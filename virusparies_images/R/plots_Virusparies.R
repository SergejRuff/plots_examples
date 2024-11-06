rm(list=ls())

library(Virusparies)

# Import  -----------------

path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
file <- ImportVirusTable(path)

# import gatherer files
path2 <- system.file("extdata", "virusgatherer.tsv", package = "Virusparies")
vg_file <- ImportVirusTable(path2)



# VhgBoxplot ---------------

# plot 1 for E-values
plot1 <- VhgBoxplot(file, x_column = "best_query", y_column = "ViralRefSeq_E")
plot1

# plot 2 for identity
plot2 <- VhgBoxplot(file, x_column = "best_query", y_column = "ViralRefSeq_ident")
plot2

# plot 3 custom arguments used
plot3 <- VhgBoxplot(file,
                    x_column = "best_query",
                    y_column = "ViralRefSeq_E",
                    theme_choice = "grey",
                    subtitle = "Custom subtitle: Identity for custom query",
                    xlabel = "Custom x-axis label: Custom query",
                    ylabel = "Custom y-axis label: Viral Reference Evalue in -log10 scale",
                    legend_position = "right")
plot3



# plot 4: Virusgatherer plot for ViralRefSeq_taxonomy agains contig length
plot5 <- VhgBoxplot(vg_file,x_column = "ViralRefSeq_taxonomy",y_column = "contig_len")
plot5



# VhgIdentitiyfacScatterplot -------------------



scatefac_plot <- VhgIdenFacetedScatterPlot(file,cutoff = 1e-5)
scatefac_plot


# VhGIdentityScatterplot ------------------------

# Basic plot
scate_plot <- VhgIdentityScatterPlot(file,cutoff = 1e-5,
                                     highlight_groups = c("Gemini_Rep","Genomo_Rep"))

plot(scate_plot$plot)

# VhgRunsBarplot -----------------------------

Runsplot <- VhgRunsBarplot(file,cut = 1e-5)
Runsplot



# VhgSumHitsBarplot -----------------------


# plot 1: plot boxplot for reads
readplot <- VhgSumHitsBarplot(file,cut = 1e-5)
readplot

# plot 2: plot boxplot for micro_reads
plot_micro <- VhgSumHitsBarplot(file,cut = 1e-5,
                                y_column = "ViralRefSeq_contigs")
plot_micro




# plot 3: contigs in Gatherer
contig_plot <- VhgSumHitsBarplot(vg_file,groupby = "ViralRefSeq_taxonomy",
                                 y_column = "contig")
contig_plot



# VgConLenViolin -----------------------------

# create a violinplot.
violinplot <- VgConLenViolin(vg_file=vg_file,cut = 1e-5,log10_scale = TRUE)

violinplot


# RunsBarplot --------------------------------


table <- VhgRunsTable(file,cut = 1e-5)
table


# VhgTabularRasa -----------------------------

### Plot boxplot for "identity"

identity <- VhgBoxplot(file,y_column = "ViralRefSeq_ident")

# Generate table

tabtable <- VhgTabularRasa(identity$summary_stats)
tabtable


# Export ------------------------------

ExportVirusPlot(plot=plot1$boxp,file_name = "boxplot1.png",units = "in",
                width = 9,height = 9 )
ExportVirusPlot(plot=plot2$boxp,file_name = "boxplot2.png",units = "in",
                width = 10,height = 9 )
ExportVirusPlot(plot=plot3$boxp,file_name = "boxplot3.png",units = "in",
                width = 9,height = 9 )
ExportVirusPlot(plot=plot5$boxp,file_name = "boxplot4.png",units = "in",
                width = 10,height = 9 )

ExportVirusPlot(plot=scatefac_plot$plot,file_name = "scatefac_plot.png",units = "in",
                width = 10,height = 9 )
ExportVirusPlot(plot=scate_plot$plot,file_name = "scate_plot.png",units = "in",
                width = 10,height = 9 )
ExportVirusPlot(plot=Runsplot$plot,file_name = "Runsplot.png",units = "in",
                width = 12,height = 9 )
ExportVirusPlot(plot=readplot$plot,file_name = "readplot.png",units = "in",
                width = 14,height = 10 )
ExportVirusPlot(plot=plot_micro$plot,file_name = "plot_micro.png",units = "in",
                width = 12,height = 9 )
ExportVirusPlot(plot=contig_plot$plot,file_name = "contig_plot.png",units = "in",
                width = 12,height = 9 )
ExportVirusPlot(plot=violinplot$plot,file_name = "violinplot.png",units = "in",
                width = 10,height = 9 )


# Export tables ---------------------------

ExportVirusGt(gtable = table,filename = "table1.png")
ExportVirusGt(gtable = tabtable,filename = "table2.png")
