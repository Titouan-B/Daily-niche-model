library(ggplot2)
library(dplyr)
library(readxl)
library(ggpattern)
library(truncnorm)
library(patchwork)
library(tidyverse)
setwd("C:/Users/adm.violaine/Desktop/Data and graphs")



# Figure S 1: Effect of the proportion of loci required to generate incompatibility on population differentiation (G) depending on the incompatibility mechanism
G <- read_excel("Supplementary_data.xlsx", sheet = 'Comparison_G')

G$G <- as.factor(G$G)

p1 <- ggplot(data=G, aes(G, fill=Hypothesis))+
  geom_bar(aes(weight = Percentage), position="dodge")+
  xlab("G") + ylab("Percentage of bimodal distributions\nafter 500 days")+
  geom_errorbar(aes(ymin=Percentage-std_percent*100, ymax=Percentage+std_percent*100), position = position_dodge(0.9), width=.4)+
  theme_classic()+
  theme(text = element_text(size = 17))+
  scale_fill_grey()+
  theme(legend.position = "none")


p2 <- ggplot(data=G, aes(G, fill=Hypothesis))+
  geom_bar(aes(weight = Fst), position="dodge")+
  xlab("G") + ylab("Average Fst values")+
  geom_errorbar(aes(ymin=Fst-sd_fst, ymax=Fst+sd_fst), position = position_dodge(0.9), width=.4)+
  theme_classic()+
  theme(text = element_text(size = 17))+
  scale_fill_grey()


fplot <- p1 + p2 + plot_annotation(tag_levels = "A")
fplot



# Figure S 2: Comparisons of simulation outcomes, after 500 vs. 2000 days.


Duration <- read_excel("Supplementary_data.xlsx", sheet = '500vs2000')

Duration$Emergence <- as.factor(Duration$Emergence)
Duration$Legend <- as.factor(Duration$Legend)

Duration_Bimodal <- Duration %>%  filter(Type=="Bimodal")

p1 <- ggplot(Duration_Bimodal, aes(y = Percentage*100, x = Emergence, group=Legend))+
  geom_point(aes(shape = Legend), size = 4)+
  geom_errorbar(aes(ymin=Percentage*100-SD, ymax=Percentage*100+SD), width=.2,
  ) +
  theme_classic() +
  labs(x = expression(paste("Emergence time (", italic("e"), ")")), y = "Proportion of simulations ending with\ntwo differentiated sub-populations")+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")

Duration_Early <- Duration %>%  filter(Type=="Delayed-dawn")

p2 <- ggplot(Duration_Early, aes(y = Percentage*100, x = Emergence, group=Legend))+
  geom_point(aes(shape = Legend), size = 4)+
  geom_errorbar(aes(ymin=Percentage*100-SD, ymax=Percentage*100+SD), width=.2,
  ) +
  theme_classic() +
  labs(x = expression(paste("Emergence time (", italic("e"), ")")), y = "Proportion of simulations ending with\na single dawn-shifted population")+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")

Duration_Late <- Duration %>%  filter(Type=="Immediate")
p3 <- ggplot(Duration_Late, aes(y = Percentage*100, x = Emergence, group=Legend))+
  geom_point(aes(shape = Legend), size = 4)+
  geom_errorbar(aes(ymin=Percentage*100-SD, ymax=Percentage*100+SD), width=.2,
  ) +
  theme_classic() +
  labs(x = expression(paste("Emergence time (", italic("e"), ")")), y = "Proportion of simulations ending with\na single immediate population")+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")

fplot <- p1 + p2 + p3 + plot_annotation(tag_levels = 'A')
fplot


# Figure S 3: Link between male-male competition and Dawn-shifted population sizes

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



# Figure S 4: Effect of the mutational standard deviation mu of ha on population differentiation
G <- read_excel("Supplementary_data.xlsx", sheet = 'Mu')

G$G <- as.factor(G$mu)

p1 <- ggplot(data=G, aes(G))+
  geom_bar(aes(weight = Percentage), position="dodge")+
  xlab("µ") + ylab("Percentage of bimodal\ndistributions after 500 days")+
  geom_errorbar(aes(ymin=Percentage-std*100, ymax=Percentage+std*100), width=.4)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()+
  theme(legend.position = "none")

p2 <- ggplot(G, aes(x = G, y = Fst)) +
  geom_point(size = 3) +
  xlab("µ") + ylab("Average Fst values")+
  geom_errorbar(aes(ymin = Fst - fst_sd, ymax = Fst + fst_sd), width = 0.3) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = "none")

p3 <- ggplot(data=G, aes(x = G, y = width))+
  geom_point(size = 3,position = position_dodge(0.5)) +
  xlab("µ") + ylab("Average width at half height of peaks")+
  geom_errorbar(aes(ymin=width-sd_w, ymax=width+sd_w), width=.2)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()

p4 <- ggplot(data=G, aes(x = G, y = distance))+
  geom_point(size = 3) +
  xlab("µ") + ylab("Average distance\nbetween peaks")+
  geom_errorbar(aes(ymin=distance-sd_d, ymax=distance+sd_d), width=.2)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_fill_grey()



fplot <- p1 + p2 + p3 +p4 + plot_annotation(tag_levels = "A")
fplot

# Figure S 5: Generational time 

G <- read_excel("Supplementary_data.xlsx", sheet = 'Base_G')

G$G <- as.factor(G$G)

ggplot(data=G, aes(x = G, y = generation_time))+
  geom_point()+
  xlab("G") + ylab("Average generation lenght in days")+
  geom_errorbar(aes(ymin=generation_time-sd_generation, ymax=generation_time+sd_generation), width=.2)+
  theme_classic()+
  theme(text = element_text(size = 17))+
  scale_fill_grey()


# Figure S 6: Factors shaping the duration of coexistence between early and late sub-populations

K_c <- read_excel("Supplementary_data.xlsx", sheet = 'K_Coexistence')

K_c$K <- as.factor(K_c$K)

p1 <- ggplot(K_c, aes(x=K, y=Elen, group = 1)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin=Elen-sd, ymax=Elen+sd), width=.2) +
  theme_classic()+
  theme(text = element_text(size = 15))+
  ylim(0,200)+
  ylab("Average coexistence duration (in days)")+
  xlab(expression("Carrying capacity of the environment (" * italic(K) * ")"))


G_c <- read_excel("Supplementary_data.xlsx", sheet = 'G_Coexistence')

G_c$G <- as.factor(G_c$G)

p2 <- ggplot(G_c, aes(x=G, y=Elength, group = 1)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin=Elength-sd_L, ymax=Elength+sd_L), width=.2) +
  theme_classic()+
  theme(text = element_text(size = 15))+
  ylab(" ")+
  xlab(expression("Genetic threshold of reproductive isolation (" * italic(G) * ")"))


dltc_c <- read_excel("Supplementary_data.xlsx", sheet = 'Competition_Coexistence')

dltc_c$dltc <- dltc_c$dltc*100

p3 <- ggplot(dltc_c, aes(x=dltc, y=Elength)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin=yminEL, ymax=ymaxEL), width=.05) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))+
  theme_classic()+
  theme(text = element_text(size = 15))+
  ylab(" ")+
  xlab(expression("Strength of male-male competition (" * italic(delta["c"]) * ")"))


final_plot <- p1 + p2 + p3 + plot_annotation(tag_levels = 'A')
final_plot

# Figure S 7: Temporal niches in the sub-populations, depending on the timing of adult emergence


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


p2 <- ggplot(G, aes(G ,width_peaks, shape = peak_type)) +
  geom_point(size = 3, position = position_dodge(0.9)) +
  xlab("G") + ylab("Average width at half height of peaks")+
  geom_errorbar(aes(ymin = width_peaks - sd_width, ymax = width_peaks + sd_width), position = position_dodge(0.9), width = 0.3) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  scale_shape_manual(values = c(16, 1))

fplot <- p1 + p2 + plot_annotation(tag_levels = "A")
fplot

# Figure S 9: Effect of variance in the timing of adult emergence on the evolution of sub-populations with differentiated temporal niches

NG_ve <- read_excel("Data.xlsx", sheet = 'Var_Emerg_NG')

ggplot(NG_ve, aes(fill = Type, y = Percentage, x = ve, 
                  pattern = Type)) +
  geom_bar_pattern(aes(pattern = Type, pattern_angle = Type),
                   position = "stack", stat = "identity", 
                   pattern_fill = "black", fill = "white", 
                   colour = "black", pattern_spacing = 0.02,
  ) +
  scale_pattern_manual(values = c('Delayed-dawn' = 'stripe',
                                  'Immediate' = 'circle',
                                  'Bimodal' = 'none')) +
  scale_pattern_angle_manual(name="Type", values = c(0, 45,45))+
  theme_classic() +
  labs(x = expression("Variance in the emergence time (" * italic(V[e]) * ")"), y = "Percentage of each\ntype of sub-population")+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")

# Figure S 10: Effect of the proportion of loci required to generate incompatibility (G) between individuals on the evolution of sub-populations with different temporal niches.

G <- read_excel("Supplementary_data.xlsx", sheet = 'G_histogram')

G$G <- as.factor(G$G)

ggplot(data=G, aes(G))+
  geom_bar(aes(weight = BM_G))+
  xlab("G") + ylab("Percentage of bimodal distributions\nafter 500 days")+
  geom_errorbar(aes(ymin=Gymin, ymax=Gymax), width=.2)+
  theme_classic()+
  theme(text = element_text(size = 17))


# Figure S 11: Effect of the emergence time e on the average timing of reproductive activities in a seasonal model (i.e. assuming non-overlapping generations)

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

# Figure S 12: Effect of low levels of male-male competition on the emergence of bimodal distribution in the activity time, in a seasonal model

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


# Figure S 13: Effect of male-male competition on the evolution of differentiated sub-populations in a seasonal model


Season_high_dltc  <- read_excel("Supplementary_data.xlsx", sheet = 'Season_high_competition')

p1 <- ggplot(Season_high_dltc, aes(x=dlt_c, y=BM)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=BM - sd*100, ymax=BM + sd*100), width=.15,
                position=position_dodge(.2)) +
  ylab("Percentage of bimodal distribitions\nafter 500 generations")+
  xlab(expression(delta*"c"))+
  theme_classic()+
  theme(text = element_text(size = 17))


p3 <- ggplot(Season_high_dltc, aes(x=dlt_c, y=width, shape = peak_type)) +
  geom_point(size = 3, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin=width-w_sd, ymax=width+w_sd), width=.2, position = position_dodge(0.9)) +
  ylab("Average width at half height of peaks")+
  xlab(expression(delta*"c"))+
  theme_classic()+
  theme(text = element_text(size = 17), legend.position = "none") +
  scale_shape_manual(values = c(16, 1))

p4 <- ggplot(Season_high_dltc, aes(dlt_c ,peak_pos, shape = peak_pos_cat)) +
  geom_point(size = 3) +
  xlab(expression(delta*"c")) + ylab("Position of bimodal\npeaks after 500 days")+
  geom_errorbar(aes(ymin = peak_pos - sd_peak, ymax = peak_pos + sd_peak), width = 0.1) +
  theme_classic() +
  theme(text = element_text(size = 17), legend.position = "none") +
  scale_shape_manual(values = c(16, 1))


fplot <- p1 + p3 + p4 + plot_annotation(tag_levels = "A")
fplot



# Figure S 14: Effect of independent genetic control of the reproductive activity timing of males and females

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
  geom_point(size=3)+
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


# Figure S 15: Effect of the coevolution of sexual activity timings and emergence timing on the evolution of differentiated sub-populations in a daily model. 

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

