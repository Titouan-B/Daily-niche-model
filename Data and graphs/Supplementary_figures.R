library(ggplot2)
library(dplyr)
library(readxl)
library(ggpattern)
library(truncnorm)
library(patchwork)
library(tidyverse)
# Figure S 1: Comparisons of simulation outcomes, after 500 vs. 2000 days.


Duration <- read_excel("Supplementary_data.xlsx", sheet = '500vs2000')

Duration$Emergence <- as.factor(Duration$Emergence)
Duration$Legend <- as.factor(Duration$Legend)

Duration_Bimodal <- Duration %>%  filter(Type=="Bimodal")

p1 <- ggplot(Duration_Bimodal, aes(y = Percentage*100, x = Emergence, group=Legend))+
  geom_point(aes(shape = Legend), size = 4)+
  geom_errorbar(aes(ymin=Percentage*100-SD, ymax=Percentage*100+SD), width=.2,
  ) +
  theme_classic() +
  labs(x = expression(paste("Emergence time (", italic("e"), ")")), y = "Percentage of simulations")+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")

Duration_Early <- Duration %>%  filter(Type=="Delayed-dawn")

p2 <- ggplot(Duration_Early, aes(y = Percentage*100, x = Emergence, group=Legend))+
  geom_point(aes(shape = Legend), size = 4)+
  geom_errorbar(aes(ymin=Percentage*100-SD, ymax=Percentage*100+SD), width=.2,
  ) +
  theme_classic() +
  labs(x = expression(paste("Emergence time (", italic("e"), ")")), y = "")+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")

Duration_Late <- Duration %>%  filter(Type=="Immediate")
p3 <- ggplot(Duration_Late, aes(y = Percentage*100, x = Emergence, group=Legend))+
  geom_point(aes(shape = Legend), size = 4)+
  geom_errorbar(aes(ymin=Percentage*100-SD, ymax=Percentage*100+SD), width=.2,
  ) +
  theme_classic() +
  labs(x = expression(paste("Emergence time (", italic("e"), ")")), y = "")+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")

fplot <- p1 + p2 + p3 + plot_annotation(tag_levels = 'A')
fplot


# Figure S 2: Temporal niches in the sub-populations, depending on the timing of adult emergence


Phenology <- read_excel("Supplementary_data.xlsx", sheet = 'Uni_vs_Bimodal')

yexpression <- (expression(atop("Average reproductive", paste("activity time ( ",bar("h"["a"])," )"))))

p1 <- ggplot(data=Phenology)+
  geom_abline(intercept = 0, linewidth=1.5, linetype = "dashed", colour="gray")+
  geom_point(data=Phenology, size=4,shape=16, colour="darkolivegreen3",aes(x=E,y=Uni))+
  geom_line(data=Phenology, linewidth=1, colour="darkolivegreen3",aes(x=E,y=Uni))+
  geom_errorbar(aes(E,ymin=Uni_min, ymax=Uni_max), width = 0.02) +
  xlab(expression(paste("Emergence time (",italic(e), ")"))) + ylab(yexpression) +
  theme_classic()+
  theme(text = element_text(size = 17))


p2 <- ggplot(data=Phenology)+
  geom_abline(intercept = 0, linewidth=1.5, linetype = "dashed", colour="gray")+
  geom_point(data=Phenology, size=4,shape=15, colour="#00bfc4",aes(x=E,y=Bi_Late))+
  geom_point(data=Phenology, size=4,shape=17, colour="#F8766D",aes(x=E,y=Bi_Early))+
  geom_line(data=Phenology, linewidth=1, colour="#00bfc4",aes(x=E,y=Bi_Late))+
  geom_line(data=Phenology, linewidth=1, colour="#F8766D",aes(x=E,y=Bi_Early))+
  geom_errorbar(aes(E,ymin=Bi_Late_min, ymax=Bi_Late_max), width = 0.02) +
  geom_errorbar(aes(E,ymin=Bi_Early_min, ymax=Bi_Early_max), width = 0.02) +
  xlab(expression(paste("Emergence time (",italic(e), ")"))) + ylab(" ")+
  theme_classic()+
  theme(text = element_text(size = 17))


final_plot <- p1 + p2 + plot_annotation(tag_levels = 'A')
final_plot


# Figure S 3: Effect of the proportion of loci required to generate incompatibility (G) between individuals on the evolution of sub-populations with different temporal niches.

G <- read_excel("Supplementary_data.xlsx", sheet = 'G_histogram')

G$G <- as.factor(G$G)

ggplot(data=G, aes(G))+
  geom_bar(aes(weight = BM_G))+
  xlab("G") + ylab("Percentage of bimodal distributions\nafter 500 days")+
  geom_errorbar(aes(ymin=Gymin, ymax=Gymax), width=.2)+
  theme_classic()+
  theme(text = element_text(size = 17))


# Figure S 4: Effect of the emergence time e on the average timing of reproductive activities in a seasonal model (i.e. assuming non-overlapping generations)

Seasonal_anc <- read_excel("Supplementary_data.xlsx", sheet = 'Seasonal_ancestral')

ggplot(data=Seasonal_anc, aes(x=Emerg, y=Activity))+
  geom_line(linewidth=1)+
  geom_point(size=3)+
  geom_abline(intercept=0,slope=1,linewidth=1,linetype="dotted")+
  scale_fill_grey() +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.05)+
  ylab("Average sexual activity time")+
  xlab("Emergence time")+
  theme_classic()+
  theme(text = element_text(size = 17))

# Figure S 5: Effect of low levels of male-male competition on the emergence of bimodal distribution in the activity time, in a seasonal model

Season <- read_excel("Supplementary_data.xlsx", sheet = 'Season_low_competition')

Season$dltc <- as.factor(Season$dltc)

ggplot(data=Season, aes(dltc))+
  geom_point(size=3, aes(x=dltc,y=bimodality))+
  ylab("Percentage of bimodal distribitions\nafter 500 generations")+
  xlab(expression("Strength of male-male competition (" * italic(delta["c"]) * ")"))+
  ylim(0,100)+
  geom_errorbar(aes(ymin=bimodality-sd*100, ymax=bimodality+sd*100), width=.2) +
  theme_classic()+
  theme(text = element_text(size = 17))


# Figure S 6: Effect of male-male competition on the evolution of differentiated sub-populations in a seasonal model


Season_high_dltc  <- read_excel("Supplementary_data.xlsx", sheet = 'Season_high_competition')

p1 <- ggplot(Season_high_dltc, aes(x=dlt_c, y=BM*100)) +
  geom_point() +
  geom_errorbar(aes(ymin=ymin*100, ymax=ymax*100), width=.05,
                position=position_dodge(.2)) +
  ylab("Percentage of bimodal distribitions\nafter 500 generations")+
  xlab(expression(delta*"c"))+
  theme_classic()+
  theme(text = element_text(size = 17))

p2 <- ggplot(Season_high_dltc, aes(x=dlt_c, y=Latepos)) +
  geom_point() +
  geom_errorbar(aes(ymin=ymn, ymax=ymx), width=.05,
                position=position_dodge(.2)) +
  ylab("Average activity time of the late\npopulation after 500 generations")+
  xlab(expression(delta*"c"))+
  theme_classic()+
  theme(text = element_text(size = 17))

fplot <- p1 + p2 + plot_annotation(tag_levels = "A")
fplot

# Figure S 7: Effect of the proportion of loci required to generate incompatibility on population differentiation (G) depending on the incompatibility mechanism
G <- read_excel("Supplementary_data.xlsx", sheet = 'Comparison_G')

G$G <- as.factor(G$G)

p1 <- ggplot(data=G, aes(G, fill=Type))+
  geom_bar(aes(weight = Percentage), position="dodge")+
  xlab("G") + ylab("Percentage of bimodal distributions\nafter 500 days")+
  geom_errorbar(aes(ymin=Percentage-std_percent*100, ymax=Percentage+std_percent*100), position = position_dodge(0.9), width=.4)+
  theme_classic()+
  theme(text = element_text(size = 17))+
  scale_fill_grey()+
  theme(legend.position = "none")


p2 <- ggplot(data=G, aes(G, fill=Type))+
  geom_bar(aes(weight = Fst), position="dodge")+
  xlab("G") + ylab("Average Fst values")+
  geom_errorbar(aes(ymin=Fst-sd_fst, ymax=Fst+sd_fst), position = position_dodge(0.9), width=.4)+
  theme_classic()+
  theme(text = element_text(size = 17))+
  scale_fill_grey()


fplot <- p1 + p2 + plot_annotation(tag_levels = "A")
fplot

# Figure S 8: Peaks caracterisation via distance and width

G <- read_excel("Supplementary_data.xlsx", sheet = 'Base_G')

G$G <- as.factor(G$G)

p1 <- ggplot(data=G, aes(x = G, y = distance_peaks))+
  geom_point()+
  xlab("G") + ylab("Average distance between peaks")+
  geom_errorbar(aes(ymin=distance_peaks-sd_dist, ymax=distance_peaks+sd_dist), width=.2)+
  theme_classic()+
  theme(text = element_text(size = 17))+
  scale_fill_grey()


p2 <- ggplot(data=G, aes(x = G, y = width_peaks))+
  geom_point()+
  xlab("G") + ylab("Average width at half height of peaks")+
  geom_errorbar(aes(ymin=width_peaks-sd_width, ymax=width_peaks+sd_width), width=.2)+
  theme_classic()+
  theme(text = element_text(size = 17))+
  scale_fill_grey()

fplot <- p1 + p2 + plot_annotation(tag_levels = "A")
fplot

# Figure S 9: Generational time 

G <- read_excel("Supplementary_data.xlsx", sheet = 'Base_G')

G$G <- as.factor(G$G)

ggplot(data=G, aes(x = G, y = generation_time))+
  geom_point()+
  xlab("G") + ylab("Average generation lenght in days")+
  geom_errorbar(aes(ymin=generation_time-sd_generation, ymax=generation_time+sd_generation), width=.2)+
  theme_classic()+
  theme(text = element_text(size = 17))+
  scale_fill_grey()

# Figure S 10: Link between male-male competition and Dawn-shifted population sizes

dltc <- read_excel("Supplementary_data.xlsx",sheet = 'Dltc_popsize')

p1 <- ggplot(dltc, aes(x = dltc, y = ratio)) +
  geom_point(size = 3) +
  xlab(expression("Strength of male-male competition (" * italic(delta["c"]) * ")")) +
  ylab("Ratio of the immediate population size\nby the delayed dawn population size") +
  theme_classic() +
  theme(text = element_text(size = 17)) 


p2 <- ggplot(dltc, aes(x = dltc, y = popsize_E)) +
  geom_point(size = 3) +
  xlab(expression("Strength of male-male competition (" * italic(delta["c"]) * ")")) +
  ylab("Dawn-shifted subpopulation size") +
  geom_errorbar(aes(ymin = popsize_E   - sd_E, ymax = popsize_E   + sd_E), width = 0.1) +
  theme_classic() +
  theme(text = element_text(size = 17)) 


fplot <- p2 + p1 + plot_annotation(tag_levels = "A")
fplot


# Figure S 11: Effect of the coevolution of sexual activity timings and emergence timing on the evolution of differentiated sub-populations in a daily model. 

G <- read_excel("Supplementary_data.xlsx",sheet = 'Co_evo_G')
G$G <- as.factor(G$G)
G$peak_type<- as.factor(G$peak_type)
p1 <- ggplot(data=G, aes(G))+
  geom_bar(aes(weight = Percentage), position="dodge")+
  xlab("G") + ylab("Percentage of bimodal\ndistributions after 500 days")+
  geom_errorbar(aes(ymin=Percentage-std*100, ymax=Percentage+std*100), position = position_dodge(0.9), width=.02)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()+
  theme(legend.position = "none")

p2 <- ggplot(G, aes(G ,peaks, shape = peak_type)) +
  geom_point(size = 3) +
  xlab("G") + ylab("Position of bimodal\npeaks after 500 days")+
  geom_errorbar(aes(ymin = peaks - std_peaks, ymax = peaks + std_peaks), width = 0.1) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = "none") +
  scale_shape_manual(values = c(16, 1))
  
p3 <- ggplot(data=G, aes(x = G, y = peaks_dist))+
  geom_point()+
  xlab("G") + ylab("Average distance\nbetween peaks")+
  geom_errorbar(aes(ymin=peaks_dist-p_std, ymax=peaks_dist+p_std), width=.2)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()


p4 <- ggplot(data=G, aes(x = G, y = peak_width))+
  geom_point()+
  xlab("G") + ylab("Average width at\nhalf height of peaks")+
  geom_errorbar(aes(ymin=peak_width-w_std, ymax=peak_width+w_std), width=.2)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()

 


fplot <- p1 + p2 + p3 + p4 + plot_annotation(tag_levels = "A")
fplot


# Figure S 12: Effect of differential evolution of the reproductive activity timing of males and females

G <- read_excel("Supplementary_data.xlsx", sheet = 'Sex_G')

G$G <- as.factor(G$G)

p1 <- ggplot(data=G, aes(G))+
  geom_bar(aes(weight = Percentage), position="dodge")+
  xlab("G") + ylab("Percentage of bimodal distributions\nafter 500 days")+
  geom_errorbar(aes(ymin=Percentage-std_percent*100, ymax=Percentage+std_percent*100), position = position_dodge(0.9), width=.4)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()+
  theme(legend.position = "none")

p2 <- ggplot(G, aes(G ,peaks, shape = peak_type)) +
  geom_point(size = 3) +
  xlab("G") + ylab("Position of bimodal\npeaks after 500 days")+
  geom_errorbar(aes(ymin = peaks - std_peaks, ymax = peaks + std_peaks), width = 0.1) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = "none") +
  scale_shape_manual(values = c(16, 1))

p3 <- ggplot(data=G, aes(x = G, y = peaks_dist))+
  geom_point()+
  xlab("G") + ylab("Average distance\nbetween peaks")+
  geom_errorbar(aes(ymin=peaks_dist-p_std, ymax=peaks_dist+p_std), width=.2)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()


p4 <- ggplot(data=G, aes(x = G, y = width, shape = width_type))+
  geom_point(size = 3,position = position_dodge(0.5)) +
  xlab("G") + ylab("Average width at half height of peaks")+
  geom_errorbar(aes(ymin=width-std_w, ymax=width+std_w), width=.2, position = position_dodge(0.5))+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()+
  scale_shape_manual(values = c(16, 1)) 


fplot <- p1 + p2 + p3 +p4 + plot_annotation(tag_levels = "A")
fplot

# Figure S 13: Effect of the mutational standard deviation mu of ha on population differentiation
G <- read_excel("Supplementary_data.xlsx", sheet = 'Mu')

G$G <- as.factor(G$mu)

p1 <- ggplot(data=G, aes(G))+
  geom_bar(aes(weight = Percentage), position="dodge")+
  xlab("Mu") + ylab("Percentage of bimodal distributions\nafter 500 days")+
  geom_errorbar(aes(ymin=Percentage-std*100, ymax=Percentage+std*100), width=.4)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()+
  theme(legend.position = "none")

p2 <- ggplot(G, aes(x = G, y = Fst)) +
  geom_point(size = 3) +
  xlab("G") + ylab("Average Fst values")+
  geom_errorbar(aes(ymin = Fst - fst_sd, ymax = Fst + fst_sd), width = 0.3) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = "none")

p3 <- ggplot(data=G, aes(x = G, y = width))+
  geom_point(size = 3,position = position_dodge(0.5)) +
  xlab("G") + ylab("Average width at half height of peaks")+
  geom_errorbar(aes(ymin=width-sd_w, ymax=width+sd_w), width=.2)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()

p4 <- ggplot(data=G, aes(x = G, y = distance))+
  geom_point(size = 3) +
  xlab("G") + ylab("Average distance\nbetween peaks")+
  geom_errorbar(aes(ymin=distance-sd_d, ymax=distance+sd_d), width=.2)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()



fplot <- p1 + p2 + p3 +p4 + plot_annotation(tag_levels = "A")
fplot

# Figure S 14: Caracterisation of linkage disequilibrium between activity time and the neutral genetic components

G <- read_excel("Supplementary_data.xlsx", sheet = 'Base_G')

G$G <- as.factor(G$G)

p1 <- ggplot(data=G, aes(G))+
  geom_bar(aes(weight = linkage_percent), position="dodge")+
  xlab("G") + ylab("Percentage of simulations\ndisplaying linkage disequilibrium")+
  geom_errorbar(aes(ymin=linkage_percent-sd_link_p*100, ymax=linkage_percent+sd_link_p*100), width=.4)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()+
  theme(legend.position = "none")

p2 <- ggplot(data=G, aes(G))+
  geom_bar(aes(weight = linkage_number), position="dodge")+
  xlab("G") + ylab("Average number of loci\nexhibiting linkage disequilibrium")+
  geom_errorbar(aes(ymin=linkage_number-sd_link_n, ymax=linkage_number+sd_link_n), width=.4)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()+
  theme(legend.position = "none")


fplot <- p1 + p2 + plot_annotation(tag_levels = "A")
fplot

