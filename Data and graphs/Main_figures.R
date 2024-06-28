library(ggplot2)
library(dplyr)
library(readxl)
library(ggpattern)
library(truncnorm)
library(patchwork)


# Figure 1: Emergence of sub-populations with different temporal niches depending on the timing of emergence of adults.

NG_e <- read_excel("Data.xlsx", sheet = 'Emergence_NG')

NG_e$Emergence <- as.factor(NG_e$Emergence)
ggplot(NG_e, aes(fill = Type, y = Percentage, x = Emergence, 
                pattern = Type)) +
  geom_bar_pattern(aes(pattern = Type, pattern_angle = Type),
                   position = "stack", stat = "identity", 
                   pattern_fill = "black", fill = "white", 
                   colour = "black", pattern_spacing = 0.02,
  ) +
  scale_pattern_manual(values = c('Early' = 'stripe',
                                  'Late' = 'circle',
                                  'Bimodal' = 'none')) +
  scale_pattern_angle_manual(name="Type", values = c(0, 45,45))+
  theme_classic() +
  labs(x = "Emergence (e)", y = "Percentage of each\ntype of sub-population")+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")



# Create a sequence of x values
x_valsE <- seq(0, 1, length.out = 100)
# Compute the corresponding y values using dtruncnorm
y_valsE <- dtruncnorm(x_valsE, mean = 0.15, sd = 0.05, a = 0, b = 1)

# Create the plot
Early <- ggplot(data.frame(x = x_valsE, y = y_valsE), aes(x = x, y = y)) + 
  geom_vline(xintercept = 0.3, linetype="dashed")+
  geom_area_pattern(aes(pattern = "stripe"), 
                    fill = "white", 
                    pattern_fill = "black", 
                    pattern_density = 0.1,
                    pattern_spacing = 0.02,
                    pattern_size = 0.2,
                    pattern = "stripe") +
  stat_function(fun=function(x) dtruncnorm(x, mean = 0.15, sd = 0.05, a = 0, b = 1))+
  theme_classic() +
  labs(x = "Day time", y = "Density of active individuals")+
  theme(text = element_text(size = 22))


x_valsL <- seq(0, 1, length.out = 100)
y_valsL <- dtruncnorm(x_valsL, mean = 0.5, sd = 0.05, a = 0, b = 1)

Late <- ggplot(data.frame(x = x_valsL, y = y_valsL), aes(x = x, y = y)) + 
  geom_vline(xintercept = 0.3, linetype="dashed")+
  geom_area_pattern(aes(pattern = "circle"), 
                    fill = "white", 
                    pattern_fill = "black", 
                    pattern_density = 0.1,
                    pattern_spacing = 0.02,
                    pattern_size = 0.2,
                    pattern = "circle") +
  stat_function(fun=function(x) dtruncnorm(x, mean = 0.5, sd = 0.05, a = 0, b = 1))+
  theme_classic() +
  labs(x = "Day time", y = "")+
  theme(text = element_text(size = 22))

Bimodal <- ggplot(data.frame(x = c(x_valsL,x_valsE), y = c(y_valsL,y_valsE)), aes(x = x, y = y)) + 
  geom_vline(xintercept = 0.3, linetype="dashed")+
  stat_function(fun=function(x) dtruncnorm(x, mean = 0.15, sd = 0.05, a = 0, b = 1) 
                + dtruncnorm(x, mean = 0.5, sd = 0.05, a = 0, b = 1))+
  theme_classic() +
  labs(x = "Day time", y = "")+
  theme(text = element_text(size = 22))


final_plot <- Early + Late + Bimodal+ plot_annotation(tag_levels = (tag_levels = list(c('B', 'C','D'), '1')))
final_plot


# Figure 2: Effect of variance in the timing of adult emergence on the evolution of sub-populations with differentiated temporal niches

NG_ve <- read_excel("Data.xlsx", sheet = 'Var_Emerg_NG')

ggplot(NG_ve, aes(fill = Type, y = Percentage, x = ve, 
               pattern = Type)) +
  geom_bar_pattern(aes(pattern = Type, pattern_angle = Type),
                   position = "stack", stat = "identity", 
                   pattern_fill = "black", fill = "white", 
                   colour = "black", pattern_spacing = 0.02,
  ) +
  scale_pattern_manual(values = c('Early' = 'stripe',
                                  'Late' = 'circle',
                                  'Bimodal' = 'none')) +
  scale_pattern_angle_manual(name="Type", values = c(0, 45,45))+
  theme_classic() +
  labs(x = expression("Variance in the emergence time (" * italic(V[e]) * ")"), y = "Percentage of each\ntype of sub-population")+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")

# Figure 3: Effects of the male-male competition and male-female asymmetry on the evolution of sub-populations with different temporal niches.

dltc <- read_excel("Data.xlsx", sheet = 'Competition_Bimodal')

dltc$dltc <- dltc$dltc*100
p1 <- ggplot(dltc, aes(x=dltc, y=Percentage)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Percentage-sd100, ymax=Percentage+sd100), width=.1,
                position=position_dodge(.2)) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))+
  xlab(expression("Strength of male-male competition (" * italic(delta["c"]) * ")")) +ylab("Percentage of bimodal distribution\nafter 500 generations")+
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

p1 <- ggplot(G, aes(x=G, y=Fst_G)) +
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


# Figure 5: Factors shaping the duration of coexistence between early and late sub-populations

K_c <- read_excel("Data.xlsx", sheet = 'K_Coexistence')

K_c$K <- as.factor(K_c$K)

p1 <- ggplot(K_c, aes(x=K, y=Elen, group = 1)) +
  geom_point() +
  geom_errorbar(aes(ymin=Elen-sd, ymax=Elen+sd), width=.2) +
  theme_classic()+
  theme(text = element_text(size = 15))+
  ylab("Average coexistence duration (in days)")+
  xlab(expression("Carrying capacity of the environment (" * italic(K) * ")"))


G_c <- read_excel("Data.xlsx", sheet = 'G_Coexistence')

G_c$G <- as.factor(G_c$G)

p2 <- ggplot(G_c, aes(x=G, y=Elength, group = 1)) +
  geom_point() +
  geom_errorbar(aes(ymin=Elength-sd_L, ymax=Elength+sd_L), width=.2) +
  theme_classic()+
  theme(text = element_text(size = 15))+
  ylab(" ")+
  xlab(expression("Genetic threshold of reproductive isolation (" * italic(G) * ")"))


dltc_c <- read_excel("Data.xlsx", sheet = 'Competition_Coexistence')

p3 <- ggplot(dltc_c, aes(x=dltc, y=Elength)) +
  geom_point() +
  geom_errorbar(aes(ymin=yminEL, ymax=ymaxEL), width=.05) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))+
  ylab(" ")+
  xlab(expression("Strength of male-male competition (" * italic(delta["c"]) * ")"))+
  theme_classic()
  
  
final_plot <- p1 + p2 + p3 + plot_annotation(tag_levels = 'A')
final_plot

