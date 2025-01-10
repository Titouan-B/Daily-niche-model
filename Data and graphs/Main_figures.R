library(ggplot2)
library(dplyr)
library(readxl)
library(ggpattern)
library(truncnorm)
library(patchwork)


# Figure 2: Emergence of sub-populations with different temporal niches depending on the timing of emergence of adults.

NG_e <- read_excel("Data.xlsx", sheet = 'Emergence_NG')

NG_e$Emergence <- as.factor(NG_e$Emergence)
ggplot(NG_e, aes(fill = Type, y = Percentage, x = Emergence, 
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
  labs(x = "Emergence (e)", y = "Percentage of each\ntype of sub-population")+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")


# Figure 3: Effects of the male-male competition and male-female asymmetry on the evolution of sub-populations with different temporal niches.

dltc <- read_excel("Data.xlsx", sheet = 'Competition_Bimodal')

dltc$dltc <- dltc$dltc*100

# Values are * 100 to make the log scale in R works, otherwise it breaks

p1 <- ggplot(dltc, aes(x=dltc, y=Percentage)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Percentage-sd100, ymax=Percentage+sd100), width=.1,
                position=position_dodge(.2)) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))+
  xlab(expression("Strength of male-male competition (" * italic(delta["c"]) * ")")) +ylab("Percentage of bimodal distribution\nafter 500 days")+
  ylim(0,100)+
  theme_classic()+
  theme(text = element_text(size = 17))


Femelles <- read_excel("Data.xlsx", sheet = 'Competition_Females')

Femelles$p <- as.factor(Femelles$p)
p2 <- ggplot(data=Femelles, aes(x=p,y=bimodality))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=bimodality-sd100, ymax=bimodality+sd100), width=.2,
                position=position_dodge(.2)) +
  ylab("")+
  xlab(expression("Maximum number of mates per female (" * italic(p) * ")"))+
  ylim(0,100)+
  theme_classic()+
  theme(text = element_text(size = 17))

final_plot <- p1 + p2 + plot_annotation(tag_levels = 'A')
final_plot


# Figure 4: Effect of genetic incompatibility and male-male competition on the genetic differentiation between sub-populations with different temporal niche

G <- read_excel("Data.xlsx", sheet = 'G_Fst')

G$G <- as.factor(G$G) 

p1 <- ggplot(G%>% na.omit(), aes(x=G, y=Fst_G)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin=Fstg_min, ymax=Fstg_max), width=.4,
                position=position_dodge(0)) +
  xlab(expression("Genetic threshold of reproductive isolation (" * italic("G") * ")")) + ylab("Average Fst values")+
  theme_classic()+
  theme(text = element_text(size = 17))


dltc <- read_excel("Data.xlsx", sheet = 'Competition_Fst')

p2 <- ggplot(dltc, aes(x=dltc, y=Fst_dltc)) +
  geom_point(size =3) +
  geom_errorbar(aes(ymin=Fstd_min, ymax=Fstd_max), width=.03) +
  xlab(expression("Strength of male−male competition (" * italic(delta*"c") * ")" )) + ylab("  ")+
  theme_classic()+
  theme(text = element_text(size = 17))


final_plot <- p1 + p2 + plot_annotation(tag_levels = 'A')
final_plot


# Figure 5: Caracterisation of linkage disequilibrium between activity time and the neutral genetic components

G <- read_excel("Data.xlsx", sheet = 'G_Fst')

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

