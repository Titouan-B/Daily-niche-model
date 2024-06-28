library(ggplot2)
library(dplyr)
library(readxl)
library(ggpattern)
library(truncnorm)
library(patchwork)

# Figure S 1: Comparisons of simulation outcomes, after 500 vs. 2000 generations.


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

Duration_Early <- Duration %>%  filter(Type=="Early")

p2 <- ggplot(Duration_Early, aes(y = Percentage*100, x = Emergence, group=Legend))+
  geom_point(aes(shape = Legend), size = 4)+
  geom_errorbar(aes(ymin=Percentage*100-SD, ymax=Percentage*100+SD), width=.2,
  ) +
  theme_classic() +
  labs(x = expression(paste("Emergence time (", italic("e"), ")")), y = "")+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")

Duration_Late <- Duration %>%  filter(Type=="Late")
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
  xlab("G") + ylab("Percentage of bimodal distributions\nafter 500 generations")+
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

