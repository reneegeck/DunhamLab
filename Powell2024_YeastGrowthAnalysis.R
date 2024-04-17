#This is analysis for Powell et al. 2024, medRxiv doi 10.1101/2024.04.12.24305393
#Written by Renee Geck, January 2024
#____________________________

library(ggplot2)
library(dplyr)
library(gcplyr)
setwd("~/Library/CloudStorage/Dropbox/Dunham Lab/G6PD/14 Powell-Skaar Collab/230926_Dagua")

#______________________________
##PART 1: Format the input data
#This file is Supplemental Table S11 of 10.1101/2024.04.12.24305393

plate_input <- read.delim("AoU_yeast_data_final.txt", sep="\t", header=FALSE)
plate_headers <- plate_input %>% slice(1:6) 
plate_data <- plate_input %>% slice(-(1:6)) 
plate_time <- plate_data %>% pull(1) 
plate_time <- as.numeric(plate_time) 
plate_length <- length(plate_time)
all_data=data.frame()
col_iter <- c(2:ncol(plate_data))
for (iter in col_iter) {
  col_num <- iter
  well_strain <- plate_headers[1,col_num] 
  strain_ls <- rep(well_strain, plate_length) 
  well_name <- plate_headers[2,col_num] 
  name_ls <- rep(well_name, plate_length)
  well_var <- plate_headers[3,col_num] 
  var_ls <- rep(well_var, plate_length)
  well_cond <- plate_headers[4,col_num] 
  cond_ls <- rep(well_cond, plate_length)
  well_repb <- plate_headers[5,col_num] 
  repb_ls <- rep(well_repb, plate_length)
  well_rept <- plate_headers[6,col_num] 
  rept_ls <- rep(well_rept, plate_length)
  od_ls <- plate_data[,col_num]
  od_ls <- as.numeric(od_ls)
  well_data <- data.frame(plate_time, od_ls, strain_ls, name_ls, var_ls, cond_ls, repb_ls, rept_ls) 
  all_data <- rbind(all_data, well_data)
  iter <- iter+1
}
all_data <- all_data %>% rename(time=plate_time, dens=od_ls, strain=strain_ls, name=name_ls, variant=var_ls, cond=cond_ls, biolrep=repb_ls, techrep=rept_ls)
all_data$variant <- factor(all_data$variant, levels=c("reference", "c.376A>G", "c.337G>A", "c.[202G>A;376A>G]", "c.563C>T", "c.592C>T", "none",
                                              "c.430C>G", "c.433A>T", "c.582C>T", "c.595A>G","c.697G>A", "c.751G>A"))

#________________________________
##PART 2: Graph the growth curves

control_vars=c("reference", "c.376A>G", "c.337G>A", "c.[202G>A;376A>G]", "c.563C>T", "c.592C>T", "none")
vus_vars=c("c.430C>G", "c.433A>T", "c.582C>T", "c.595A>G","c.697G>A", "c.751G>A")

#plot control variants of known effect for Supplemental Figure S7A
control_curves <- ggplot(data=subset(all_data, variant %in% control_vars),
       aes(x=time, y=dens)) + #plot od/time by variant
  geom_point(aes(fill=biolrep, color=biolrep), size=0.2, alpha=0.5, show.legend=F) +
  scale_y_continuous(trans = "log2", limits=c(0.01,2), n.breaks = 15,
                     labels = function(x) ifelse(x == 0, "0", x)) + #make log2 yaxis
  scale_x_continuous(breaks=c(0, 24, 48)) +
  labs(x="time (hours)", y="optical density") +
  facet_grid(cond~variant+name) + #make a different plot for each media
  scale_color_manual(values=c("darkgray", "black"))+
  scale_fill_manual(values=c("darkgray", "black"))+
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("AoUvars_yeast_control_curves.pdf", plot=control_curves, device=pdf, width=10, height=3.5)

#plot variants of uncertain effect for Supplemental Figure S7A
vus_curves <- ggplot(data=subset(all_data, variant %in% vus_vars),
                         aes(x=time, y=dens)) + #plot od/time by variant
  geom_point(aes(fill=biolrep, color=biolrep), size=0.2, alpha=0.5, show.legend=F) +
  scale_y_continuous(trans = "log2", limits=c(0.01,2), n.breaks = 15,
                     labels = function(x) ifelse(x == 0, "0", x)) + #make log2 yaxis
  scale_x_continuous(breaks=c(0, 24, 48)) +
  labs(x="time (hours)", y="optical density") +
  facet_grid(cond~variant+name) + #make a different plot for each media
  scale_color_manual(values=c("darkgray", "black"))+
  scale_fill_manual(values=c("darkgray", "black"))+
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("AoUvars_yeast_vus_curves.pdf", plot=vus_curves, device=pdf, width=8.71, height=3.5)

#_______________________________
##PART 3: Area Under the Curve
#This is computed using package gcplyr, see doi 10.1101/2023.04.30.538883
data_auc <- all_data %>% group_by(strain, name, variant, cond, biolrep, techrep) %>%
  summarize(auc = auc(x = time, y = dens))
data_auc$variant <- factor(data_auc$variant, levels=c("none", "c.751G>A", "c.697G>A", "c.595A>G", "c.582C>T",
                                                      "c.433A>T", "c.430C>G", "c.592C>T", "c.563C>T",
                                                      "c.[202G>A;376A>G]", "c.337G>A", "c.376A>G", "reference"))

#plot AUC in selection media (c-ura-met+1mM H2O2) for Figure 4
selection_auc <- ggplot(data=subset(data_auc, cond=="selection"), mapping=aes(x=variant, y=auc)) +
  geom_boxplot(show.legend=F, size=0.5,  color="black", outlier.size=0) +
  geom_jitter(aes(shape=biolrep), show.legend=F, width=0.1, height=0) +
  labs(x="", y="AUC") +
  scale_fill_viridis_d(option="mako", begin=0.2, end=0.7)+
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = "white", linewidth=0.5),
        panel.background = element_rect(fill = '#e9e9e9')) +
  ylim(0, 50)+
  coord_flip()
ggsave("AoUvars_yeast_selection_auc.pdf", plot=selection_auc, device=pdf, width=8, height=4)

print(TukeyHSD(aov(auc ~ variant, data = subset(data_auc, cond=="selection"))))

#plot AUC in control (c-ura) media for Supplemental Figure S7B
control_auc <- ggplot(data=subset(data_auc, cond=="control"), mapping=aes(x=variant, y=auc)) +
  geom_boxplot(show.legend=F, size=0.5,  color="black", outlier.size=0) +
  geom_jitter(aes(shape=biolrep), show.legend=F, width=0.1, height=0) +
  labs(x="", y="AUC") +
  scale_fill_viridis_d(option="mako", begin=0.2, end=0.7)+
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = "white", linewidth=0.5),
        panel.background = element_rect(fill = '#e9e9e9')) +
  ylim(0,55)+
  coord_flip()
ggsave("AoUvars_yeast_control_auc.pdf", plot=control_auc, device=pdf, width=8, height=4)

print(TukeyHSD(aov(auc ~ variant, data = subset(data_auc, cond=="control"))))