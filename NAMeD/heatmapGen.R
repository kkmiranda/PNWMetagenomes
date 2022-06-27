library(ggplot2)
library(plotly)
library(gplots)
library(plotrix)
library(tidyr)
library(dplyr)
library(cowplot)
library(egg)
library(stringr)

get_env = Vectorize(function(x) {
  s = stringr::str_split_fixed(x,"_",2)
  paste(s)
})

ur_labels = c("","","Urease","","")
nfix_labels = c("","","","Nitrogen Fixation","","","","")
ntrf_labels = c("","","","Nitrification","","")
dnr_labels = c("","","","","Dissimilatory nitrate reduction","","","","")
dntrf_labels = c("","","","","Denitrification","","","","","")
anr_labels = c("","","Assimilatory nitrate reduction","","","")
anmx_labels = c("","","","Annamox","","","")

MAG_order = c("NLE_BLD_Bin_00001","NLE_BLD_Bin_00002","NLE_BLD_Bin_00003","NLE_BLD_Bin_00004","NLE_BLD_Bin_00005","NLE_BLD_Bin_00006","NLE_BLD_Bin_00007","NLE_BLD_Bin_00008","NLE_BLD_Bin_00009","NLE_BLD_Bin_00011","NLE_BLD_Bin_00013","NLE_BLD_Bin_00014","NLE_BLD_Bin_00015","NLE_BLD_Bin_00016","NLE_BLD_Bin_00017","NLE_BLD_Bin_00019","NLE_BLD_Bin_00020","NLE_BLD_Bin_00022","NLE_BLD_Bin_00024","NLE_BLD_Bin_00025","NLE_BLD_Bin_00027","NLE_BLD_Bin_00028","NLE_BLD_Bin_00029","NLE_BLD_Bin_00030","NLE_BLD_Bin_00031","NLE_BLD_Bin_00032","NLE_BLD_Bin_00033","NLE_BLD_Bin_00034","NLE_BLD_Bin_00035","NLE_BLD_Bin_00036","NLE_BLD_Bin_00037","NLE_BLD_Bin_00038",
              "LSE_BLD_Bin_00001","LSE_BLD_Bin_00002","LSE_BLD_Bin_00003","LSE_BLD_Bin_00004","LSE_BLD_Bin_00005","LSE_BLD_Bin_00006","LSE_BLD_Bin_00007","LSE_BLD_Bin_00008","LSE_BLD_Bin_00009","LSE_BLD_Bin_00010","LSE_BLD_Bin_00011","LSE_BLD_Bin_00012","LSE_BLD_Bin_00013","LSE_BLD_Bin_00014","LSE_BLD_Bin_00017","LSE_BLD_Bin_00018",
              "PSC_BLD_Bin_00001","PSC_BLD_Bin_00002","PSC_BLD_Bin_00003","PSC_BLD_Bin_00004","PSC_BLD_Bin_00005","PSC_BLD_Bin_00006","PSC_BLD_Bin_00007","PSC_BLD_Bin_00008","PSC_BLD_Bin_00009","PSC_BLD_Bin_00010","PSC_BLD_Bin_00011","PSC_BLD_Bin_00012","PSC_BLD_Bin_00013","PSC_BLD_Bin_00014","PSC_BLD_Bin_00017","PSC_BLD_Bin_00018",
              "PSC_RHZ_Bin_00001","PSC_RHZ_Bin_00002","PSC_RHZ_Bin_00003","PSC_SED_Bin_00001","PSC_SED_Bin_00002","PSC_SED_Bin_00003","PSC_SED_Bin_00004","PSC_SED_Bin_00005",
              "PSE_RHZ_Bin_00001","PSE_RHZ_Bin_00002","PSE_RHZ_Bin_00003","PSE_RHZ_Bin_00004","PSE_RHZ_Bin_00005","PSE_RHZ_Bin_00006","PSE_RHZ_Bin_00007","PSE_RHZ_Bin_00010","PSE_RHZ_Bin_00011","PSE_RHZ_Bin_00012","PSE_RHZ_Bin_00015","PSE_RHZ_Bin_00016","PSE_RHZ_Bin_00017",
              "ZMA_RHZ_Bin_00001","ZMA_RHZ_Bin_00002","ZMA_RHZ_Bin_00003","ZMA_RHZ_Bin_00004","ZMA_RHZ_Bin_00005","ZMA_RHZ_Bin_00006","ZMA_RHZ_Bin_00007","ZMA_RHZ_Bin_00008","ZMA_RHZ_Bin_00009","ZMA_RHZ_Bin_00010","ZMA_RHZ_Bin_00011","ZMA_RHZ_Bin_00013",
              "ZMA_SED_Bin_00001","ZMA_SED_Bin_00002","ZMA_SED_Bin_00003","ZMA_SED_Bin_00004","ZMA_SED_Bin_00005","ZMA_SED_Bin_00006","ZMA_SED_Bin_00009",
              "NLU_BLB_Bin_00001")

genHeatMap = function(fileName, fmt) {
  hmData = read.csv(fileName,sep="\t",header=TRUE)
  hmData$MAG_names = as.character(hmData$MAG_names)
  hmData$MAG_names = factor(hmData$MAG_names, levels = MAG_order)
  # output = paste(fle, "heatmap.png", sep = "_")
  
  # For Wide-format
  if (fmt == "w") {
    hmData = pivot_longer(hmData, !MAG_names, names_to = "Gene", values_to = "Presence")
    hmData$dup = hmData$MAG_names
    hmData = as.data.frame(hmData %>% 
                             separate(dup, c("env", "bin"), sep="_Bin") %>% 
                             select(-bin))
    if (fileName=="heatmapFiles/finalNitrogen.txt") {
      return(hmData)
    }
    
    baseGraph = ggplot(hmData, aes(x = MAG_names, y=Gene, fill=env, alpha=Presence))
  } else {
    baseGraph = ggplot(hmData, aes(x = MAG_names, y=1, fill=env, alpha=Presence))
  }
  
  
  # Factoring out the names
  outputGraph =  baseGraph + 
    geom_tile() + 
    scale_fill_manual(values = c("#7B904B", "#E5F2C9", "#81DCD4", "#FFC857", "#E9724C", "#C5283D", "#481D24", "#E0BAD7", "#065858"), 
                                    breaks=c("LSE_BLD", "NLU_BLB", "PSC_BLD", "PSC_RHZ", "PSC_SED", "PSE_RHZ", "ZMA_RHZ", "ZMA_SED", "NLE_BLD")) + 
    theme(
      axis.text.x = element_text(angle=90), 
      legend.position = "none",
      panel.background = element_rect(fill = "transparent",colour = NA))
  
  return(outputGraph)
}

finalThiamin = genHeatMap("heatmapFiles/finalB1.txt", fmt="w")
pathwayThiamin = genHeatMap("heatmapFiles/pathwayThiamin.txt", fmt="l")
ggarrange(finalThiamin + theme(axis.text.x = element_blank()), 
          pathwayThiamin + theme(axis.text.y = element_blank()))

finalRiboflavin = genHeatMap("heatmapFiles/finalB2.txt", fmt="w")
pathwayRiboflavin = genHeatMap("heatmapFiles/pathwayRiboflavin.txt", fmt="l")
ggarrange(finalRiboflavin + theme(axis.text.x = element_blank()), 
          pathwayRiboflavin + theme(axis.text.y = element_blank()))

finalBiotin = genHeatMap("heatmapFiles/finalB7.txt", fmt="w")
pathwayBiotin = genHeatMap("heatmapFiles/pathwayBiotin.txt", fmt="l")
ggarrange(finalBiotin + theme(axis.text.x = element_blank()),
          pathwayBiotin + theme(axis.text.x = element_blank()))

finalCobalamin = genHeatMap("heatmapFiles/finalB12.txt", fmt="w")
pathwayCobalamin = genHeatMap("heatmapFiles/pathwayCobalamin.txt",fmt="w")
ggarrange(finalCobalamin+ theme(axis.text.x = element_blank()), 
          pathwayCobalamin)

#Vitamin Comparison
ggarrange(finalThiamin + theme(axis.text.x = element_blank()) + ylab("B1 P/A"),
          finalRiboflavin + theme(axis.text.x = element_blank()) + ylab("B2 P/A"),
          
          pathwayThiamin + theme(axis.text.y = element_blank(), axis.text.x = element_blank()) + ylab("B1 pathway"),
          pathwayRiboflavin + theme(axis.text.y = element_blank(), axis.text.x = element_blank()) + ylab("B2 pathway"),
          
          finalBiotin + theme(axis.text.x = element_blank()) + ylab("B7 P/A"),
          finalCobalamin + theme(axis.text.x = element_blank()) + ylab("B12 P/A"), 
          
          pathwayBiotin + theme(axis.text.x = element_blank()) + ylab("B7 pathway"),
          pathwayCobalamin + theme(axis.text.x = element_blank()) + ylab("B12 pathway"), 
          nrow = 4, ncol = 2
          )

genQuantMap = function(fileName, fmt) {
  hmData = read.csv(fileName,sep="\t",header=TRUE)
  hmData$MAG_names = as.character(hmData$MAG_names)
  hmData$MAG_names = factor(hmData$MAG_names, levels = MAG_order)
  # output = paste(fle, "heatmap.png", sep = "_")
  
  # For Wide-format
  if (fmt == "w") {
    hmData = pivot_longer(hmData, !MAG_names, names_to = "Gene", values_to = "Presence")
    hmData$dup = hmData$MAG_names
    hmData = as.data.frame(hmData %>% 
                             separate(dup, c("env", "bin"), sep="_Bin") %>% 
                             select(-bin))
    if (fileName=="heatmapFiles/finalNitrogen.txt") {
      return(hmData)
    }
    
    baseGraph = ggplot(hmData, aes(x = MAG_names, y=Gene, fill=Presence))
  } else {
    baseGraph = ggplot(hmData, aes(x = MAG_names, y=1, fill=Presence))
  }
  
  
  # Factoring out the names
  outputGraph =  baseGraph + 
    geom_tile() + scale_fill_gradientn(limits=c(0,100), colors = c('white','black')) +
    theme(
      axis.text.x = element_text(angle=90), 
      legend.position = 'none',
      panel.background = element_rect(fill = "transparent",colour = NA))
  
  return(outputGraph)
}

allHydrolases = genQuantMap("heatmapFiles/finalAmmonHydrolases.txt", fmt="l")
finalCarbon = genQuantMap("heatmapFiles/finalCarbon.txt", fmt="w") + 
  scale_y_discrete(position="right",
                   breaks=c("Carbohydrates_general","Carbohydrates_pentoses","Carboxylic_acids","Compatible_solutes"),
                   labels=c("Carbohydrates general","Carbohydrates pentoses","Carboxylic acids","Compatible solutes"))

#nitrogen
nitrogen = genHeatMap("heatmapFiles/finalNitrogen.txt", fmt="w")
nitrogen$dup = nitrogen$Gene
nitrogen = as.data.frame(nitrogen %>% 
                         separate(dup, c("pathway", "gene"), sep="_") %>% 
                         select(-gene))
nitrogen$pathway = as.character(nitrogen$pathway)
nitrogen$pathway = factor(nitrogen$pathway, levels=unique(nitrogen$pathway))

dnr = ggplot(filter(nitrogen, pathway=="Dissimilatory.nitrate.reduction"), aes(x = MAG_names, y=Gene, fill=env, alpha=Presence)) + geom_tile()
anr = ggplot(filter(nitrogen, pathway=="Assimilatory.nitrate.reduction"), aes(x = MAG_names, y=Gene, fill=env, alpha=Presence)) + geom_tile()
ur = ggplot(filter(nitrogen, pathway=="Urease"), aes(x = MAG_names, y=Gene, fill=env, alpha=Presence)) + geom_tile()
ntrf = ggplot(filter(nitrogen, pathway=="Nitrification"), aes(x = MAG_names, y=Gene, fill=env, alpha=Presence)) + geom_tile()
dntrf = ggplot(filter(nitrogen, pathway=="Denitrification"), aes(x = MAG_names, y=Gene, fill=env, alpha=Presence)) + geom_tile()
anmx = ggplot(filter(nitrogen, pathway=="Annamox"), aes(x = MAG_names, y=Gene, fill=env, alpha=Presence)) + geom_tile()
nfix = ggplot(filter(nitrogen, pathway=="Nitrogen.fixation"), aes(x = MAG_names, y=Gene, fill=env, alpha=Presence)) + geom_tile()

dnr = dnr + geom_tile() + 
  scale_fill_manual(values = c("#7B904B", "#E5F2C9", "#81DCD4", "#FFC857", "#E9724C", "#C5283D", "#481D24", "#E0BAD7", "#065858"), 
                    breaks=c("LSE_BLD", "NLU_BLB", "PSC_BLD", "PSC_RHZ", "PSC_SED", "PSE_RHZ", "ZMA_RHZ", "ZMA_SED", "NLE_BLD")) + 
  scale_y_discrete(labels=dnr_labels) + 
  theme(
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=12),
    plot.margin = margin(0,0,0,10,"pt"),
    legend.position = "none",
    panel.background = element_rect(fill = "transparent",colour = NA))

anr = anr + geom_tile() + 
  scale_fill_manual(values = c("#7B904B", "#E5F2C9", "#81DCD4", "#FFC857", "#E9724C", "#C5283D", "#481D24", "#E0BAD7", "#065858"), 
                    breaks=c("LSE_BLD", "NLU_BLB", "PSC_BLD", "PSC_RHZ", "PSC_SED", "PSE_RHZ", "ZMA_RHZ", "ZMA_SED", "NLE_BLD")) + 
  scale_y_discrete(labels=anr_labels) + 
  theme(
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size=12),
    plot.margin = margin(10,0,0,10,"pt"),
    legend.position = "none",
    panel.background = element_rect(fill = "transparent",colour = NA))

ur = ur + geom_tile() + 
  scale_fill_manual(values = c("#7B904B", "#E5F2C9", "#81DCD4", "#FFC857", "#E9724C", "#C5283D", "#481D24", "#E0BAD7", "#065858"), 
                    breaks=c("LSE_BLD", "NLU_BLB", "PSC_BLD", "PSC_RHZ", "PSC_SED", "PSE_RHZ", "ZMA_RHZ", "ZMA_SED", "NLE_BLD")) + 
  scale_y_discrete(labels=ur_labels) + 
  theme(
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size=12),
    plot.margin = margin(0,0,0,10,"pt"),
    legend.position = "none",
    panel.background = element_rect(fill = "transparent",colour = NA))

ntrf = ntrf + geom_tile() + 
  scale_fill_manual(values = c("#7B904B", "#E5F2C9", "#81DCD4", "#FFC857", "#E9724C", "#C5283D", "#481D24", "#E0BAD7", "#065858"), 
                    breaks=c("LSE_BLD", "NLU_BLB", "PSC_BLD", "PSC_RHZ", "PSC_SED", "PSE_RHZ", "ZMA_RHZ", "ZMA_SED", "NLE_BLD")) + 
  scale_y_discrete(labels=ntrf_labels) + 
  theme(
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size=12),
    plot.margin = margin(0,0,0,10,"pt"),
    legend.position = "none",
    panel.background = element_rect(fill = "transparent",colour = NA))

dntrf = dntrf + geom_tile() + 
  scale_fill_manual(values = c("#7B904B", "#E5F2C9", "#81DCD4", "#FFC857", "#E9724C", "#C5283D", "#481D24", "#E0BAD7", "#065858"), 
                    breaks=c("LSE_BLD", "NLU_BLB", "PSC_BLD", "PSC_RHZ", "PSC_SED", "PSE_RHZ", "ZMA_RHZ", "ZMA_SED", "NLE_BLD")) + 
  ylab("DNTRF") +
  scale_y_discrete(labels=dntrf_labels) + 
  theme(
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size=12),
    plot.margin = margin(0,0,0,10,"pt"),
    legend.position = "none",
    panel.background = element_rect(fill = "transparent",colour = NA))

anmx = anmx + geom_tile() + 
  scale_fill_manual(values = c("#7B904B", "#E5F2C9", "#81DCD4", "#FFC857", "#E9724C", "#C5283D", "#481D24", "#E0BAD7", "#065858"), 
                    breaks=c("LSE_BLD", "NLU_BLB", "PSC_BLD", "PSC_RHZ", "PSC_SED", "PSE_RHZ", "ZMA_RHZ", "ZMA_SED", "NLE_BLD")) + 
  ylab("ANMX") +
  scale_y_discrete(labels=anmx_labels) + 
  theme(
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size=12),
    plot.margin = margin(0,0,10,10,"pt"),
    legend.position = "none",
    panel.background = element_rect(fill = "transparent",colour = NA))

nfix = nfix + geom_tile() + 
  scale_fill_manual(values = c("#7B904B", "#E5F2C9", "#81DCD4", "#FFC857", "#E9724C", "#C5283D", "#481D24", "#E0BAD7", "#065858"), 
                    breaks=c("LSE_BLD", "NLU_BLB", "PSC_BLD", "PSC_RHZ", "PSC_SED", "PSE_RHZ", "ZMA_RHZ", "ZMA_SED", "NLE_BLD")) + 
  ylab("NFIX") +
  scale_y_discrete(labels=nfix_labels) + 
  theme(
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=12),
    plot.margin = margin(0,0,0,10,"pt"),
    legend.position = "none",
    panel.background = element_rect(fill = "transparent",colour = NA))

outputPathwayPlot = ggarrange(finalCarbon + 
                              theme(axis.text.x = element_blank(), 
                                    axis.title.x = element_blank(),
                                    axis.ticks.x = element_blank(),
                                    axis.text.y = element_text(size=12),
                                    axis.title.y = element_text(size=16)) + 
                              ylab("DOC Uptake"),
                            anr,
                            allHydrolases + 
                              theme(
                                axis.text = element_blank(),
                                axis.text.x = element_blank(),
                                axis.ticks = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_text(angle=0, vjust=0.5, size = 12, colour = "#3c3c3c"), 
                                plot.margin = margin(0,0,0,10,"pt")
                              ) + 
                              ylab("Ammonification Hydrolases"),
                            ur,dnr,ntrf,dntrf,nfix,anmx,
                            pathwayThiamin + 
                              theme(
                                axis.text = element_blank(),
                                axis.text.x = element_blank(),
                                axis.ticks = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_text(angle=0, vjust=0.5, size = 12, colour = "#3c3c3c"), 
                                plot.margin = margin(0,0,0,10,"pt")
                              ) + 
                              ylab("Vit. B1"),
                            pathwayRiboflavin + 
                              theme(
                                axis.text = element_blank(),
                                axis.text.x = element_blank(),
                                axis.ticks = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_text(angle=0, vjust=0.5, size = 12, colour = "#3c3c3c"), 
                                plot.margin = margin(0,0,0,10,"pt")
                              ) + 
                              ylab("Vit. B2"),
                            pathwayBiotin + 
                              theme(
                                axis.text = element_blank(),
                                axis.text.x = element_blank(),
                                axis.ticks = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_text(angle=0, vjust=0.5, size = 12, colour = "#3c3c3c"), 
                                plot.margin = margin(0,0,0,10,"pt")
                              ) + 
                              ylab("Vit. B7"),
                            pathwayCobalamin + 
                              theme(
                                axis.title.y = element_blank(),
                                axis.text.x = element_text(size = 8), 
                                axis.text.y = element_text(size = 12),
                                plot.margin = margin(0,0,0,10,"pt")
                              ),
                            ncol=1, heights = c(3,1,1,1,1,1,1,1,1,1,1,1,2))



outputFinalPlot = ggarrange(finalCarbon + 
            theme(axis.text.x = element_blank(), 
                  axis.title.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.text.y = element_text(size=12),
                  axis.title.y = element_text(size=16)) + 
            ylab("DOC Uptake"),
          anr,
          allHydrolases + 
            theme(
              axis.text = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_text(angle=0, vjust=0.5, size = 12, colour = "#3c3c3c"), 
              plot.margin = margin(0,0,0,10,"pt")
            ) + 
            ylab("Ammonification Hydrolases"),
          ur,dnr,ntrf,dntrf,nfix,anmx,
          finalThiamin + 
            theme(
              axis.text = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_text(angle=0, vjust=0.5, size = 12, colour = "#3c3c3c"), 
              plot.margin = margin(0,0,0,10,"pt")
            ) + 
            ylab("Vit. B1"),
          finalRiboflavin + 
            theme(
              axis.text = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_text(angle=0, vjust=0.5, size = 12, colour = "#3c3c3c"), 
              plot.margin = margin(0,0,0,10,"pt")
            ) + 
            ylab("Vit. B2"),
          finalBiotin + 
            theme(
              axis.text = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_text(angle=0, vjust=0.5, size = 12, colour = "#3c3c3c"), 
              plot.margin = margin(0,0,0,10,"pt")
            ) + 
            ylab("Vit. B7"),
          finalCobalamin + 
            theme(
              axis.text = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_text(angle=0, vjust=0.5, size = 12, colour = "#3c3c3c"), 
              plot.margin = margin(0,0,0,10,"pt")
              ) +
            ylab("Vit. B12"),
          ncol=1, heights = c(3,1,1,1,1,1,1,1,1,1,1,1,1))

ggsave(plot = outputFinalPlot, "output_heatmap.png", width = 12.2, height = 8.4, device="png", dpi=700)
ggsave(plot = outputPathwayPlot, 'output_pathway.png', width = 12.2, height = 8.4, device="png", dpi=700)
dev.off()
