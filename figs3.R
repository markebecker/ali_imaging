# forthepape II

# libraries
library(ggplot2)
library(ggbeeswarm)
library(data.table)
library(dplyr)
library(MetBrewer)
library(lme4)
library(ggsignif)

####################################################################################
####################################################################################
####################################################################################
##################                    Figure 1                    ##################
####################################################################################
####################################################################################
####################################################################################

# qPCR over time, figure 1
qpcr_p3old <- read.csv("C:/Users/Mark/Desktop/figure_images/1/qpcr_p3old.csv")
keys <- colnames(qpcr_p3old)[!grepl('copies.mm.2',colnames(qpcr_p3old))]
X <- as.data.table(qpcr_p3old)
qpcrcurve = X[,list(mm= mean(copies.mm.2), sd = sd(copies.mm.2), hpi = hpi),'Sample.Name']
trimmedq = distinct(qpcrcurve)
#PVAL = 0.027 for 1:72 hpi
t.test(subset(trimmedq, hpi == 1)$mm, subset(trimmedq, hpi == 72)$mm)
shapiro.test(subset(trimmedq, hpi == 120)$mm)

#mean of the final timepoint
mean(subset(trimmedq, hpi == 120)$mm)

#the plot
ggplot(trimmedq, aes(x=hpi, y=mm)) + 
  labs(x = 'HPI', y = bquote('N copies per '~mm^2))+
  scale_y_log10(limits = c(-1, 8e8)) +
  geom_jitter(alpha = 0.5, color = 'deeppink2', size = 2) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               size = 0.25,
               linewidth = 0.5,
               geom = "pointrange") +
  theme_bw() + 
  #geom_hline(yintercept = 118, color = 'red') + 
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8), plot.margin = margin(5, 10, 0, 5))

ggsave("test.pdf", width = 2.3, height = 1.4, units = "in")

####################################################################################
####################################################################################
####################################################################################
##################                    Figure 2                    ##################
####################################################################################
####################################################################################
####################################################################################

# Figure 2c: fxn gfp+ culture area over time
relevant = subset(norspots, (condition == 'vanilla' | condition == 'static' | condition == 'shaken') & ((donor != 7 & donor != 8) | (time < 120)) & (virus == 'egfp'))
ggplot(subset(relevant), aes(time, as.numeric(norfreq), colour=factor(donor))) +
  labs(x = 'Hours post infection', y ="# of GFP+ spots", color = "Donor:")+
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8), 
  legend.key.size = unit(1, 'point'),
  legend.box.margin = margin(0, 0, 0, -10),
  legend.text = element_text(size=6), legend.title = element_text(size=8),
  plot.margin = margin(1, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_point(alpha = 0.6, size = 0.2) 
  #geom_point(alpha = 0.3) +
  #geom_smooth(aes(color = NULL))
ggsave("C:/Users/Mark/Desktop/figure_images/2/fig2c.pdf", width = 2.2, height = 1.8, units = "in")
  

# Figure 2d: peak # gfp+ spots vs. n rna copies per mm2
ggplot(subset(good.spots.summary.j, virus == "egfp" & (condition == 'vanilla' | condition == 'shaken' | condition == 'static')), (aes(y=as.numeric(n.mm.2mean), x= good.peak.spots,  color = factor(donor)))) + 
  labs(y = bquote('N copies per '~mm^2~''), x ="# GFP+ Spots", color = "Donor:")+
  #bquote('N copies per '~mm^2' culture')
  geom_point(show.legend = FALSE) +
  scale_y_log10()+
  theme_bw() + 
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/2/fig2d.pdf", width = 2, height = 1.8, units = "in")


# interesting to look at colored by rotary or by donor.
n.mm.2mean
trimrnafrac$ncap <- as.numeric(trimrnafrac$mm)
cor.test(log10(as.numeric(good.spots.summary.j$n.mm.2mean)), good.spots.summary.j$good.peak.spots, use = "complete.obs", method = 'pearson')
cor.test(as.numeric(good.spots.summary.j$n.mm.2mean), good.spots.summary.j$good.peak.spots, use = "complete.obs", method = 'pearson')

#Figure 2e: number of initial foci
level_order <- c('mock', 'egfp')
ggplot(subset(good.spots.summary.j, condition == "vanilla" | condition == 'static' | condition == 'shaken'), (aes(y=X20hpi.spots, x=factor(virus, level = level_order), color = factor(donor)))) + 
  labs(x = 'Condition', y ="# GFP+ spots at 20 HPI", color = "Donor:")+
  geom_hline(yintercept = 0, color = 'black') +
  geom_quasirandom(alpha = 0.8, show.legend = FALSE) +
  stat_summary(fun.y = median,
               fun.ymin = function(z) { quantile(z,0.25) }, 
               fun.ymax = function(z) { quantile(z,0.75) }, 
               geom = "pointrange",
               show.legend = FALSE,
               aes(group = factor(virus, level = level_order))) +
  scale_x_discrete(labels = c("Mock", "Infected"))+
  theme_bw() + 
  #scale_y_log10() +
  theme(axis.text = element_text(size=5), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/2/fig2e.pdf", width = 1, height = 1.8, units = "in")

relevant <- subset(good.spots.summary.j, virus == 'egfp' & (condition == "vanilla" | condition == 'static' | condition == 'shaken'))
median(relevant$X20hpi.spots, na.rm = 'TRUE')

#Figure 2f: time to peak egfp+ area
level_order <- c('mock', 'egfp') 
ggplot(subset(good.spots.summary.j, (condition == "vanilla" | condition == 'static' | condition == 'shaken') & (2*video.end.frame+1 > good.peak.spots.time) & !(expt.name == 'agarose_d3' & sample.no == 3)), (aes(y=good.peak.spots.time, x=factor(virus, level = level_order), color = good.peak.spots))) + 
  labs(x = 'Condition', y ="Time of Peak GFP+ Spots (HPI)", color = "# Spots")+
  geom_quasirandom(size = 0.3, alpha = 0.8) +
  stat_summary(fun.y = median,
               fun.ymin = function(z) { quantile(z,0.25) }, 
               fun.ymax = function(z) { quantile(z,0.75) }, 
               color = 'black',
               linewidth = 0.2,
               size = 0.3,
               geom = "crossbar") +
  scale_color_continuous(high = '#05e895', low = '#402375', limits = c(NA, NA)) +
  scale_x_discrete(labels = c("Mock", "Infected"))+
  theme_bw() +
  theme(axis.text = element_text(size=5), axis.title = element_text(size=8),
        plot.margin = margin(0, -5, 0, 0),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.1, 'in'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=4), legend.title = element_text(size=5))
ggsave("C:/Users/Mark/Desktop/figure_images/2/fig2f.pdf", width = 1.4, height = 2, units = "in")

relevant <- subset(good.spots.summary.j, condition == "vanilla" & (2*video.end.frame+1 > good.peak.spots.time) & !(expt.name == 'agarose_d3' & sample.no == 3) & virus == 'egfp')
2*median(relevant$peak.fxn.1100.frame, na.rm = 'TRUE')+1



#Figure 2h: correlation of focus type w/ mucus mvmt & peak egfp
level_order <- c('mock', 'egfp') 
level_order <- c('plaque', 'comet', 'diffuse', 'crypt', 'mock') 

relevant <- subset(good.spots.summary.j, (condition == 'vanilla' | condition == 'shaken' | condition == 'static') & virus == 'egfp')
level_order <- c('plaque', 'comet', 'diffuse') 
ggplot(relevant, (aes(y=as.numeric(good.peak.spots), x=factor(focus.simple, levels = level_order), color = (mucus.simple == 'rotary')))) + 
  labs(x = 'Focus Type', y ="Peak GFP+ Spots", color = "MCC Type:")+
  scale_color_manual(labels = c("Disorganized", "Rotary"), values = c('#4A4CA0', 'deeppink2'))+
  theme_bw() + 
  geom_quasirandom(alpha = 0.7, width = 0.2) +
  ylim(0, 40000) +
  scale_x_discrete(limits = level_order, labels = c('Plaque', 'Comet', 'Diffuse'))+
  stat_summary(fun = median,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "crossbar",
               width = 0.3,
               aes(group = FALSE, color = NULL)) +
  guides(shape = guide_legend(override.aes = list(size = 0.01))) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, -5, 0, 0),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.1, 'in'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=5), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/2/fig2i.pdf", width = 2, height = 1.35, units = "in")

 wilcox.test(subset(relevant, focus.simple == 'comet')$good.peak.spots, subset(relevant, focus.simple == 'plaque')$good.peak.spots)
wilcox.test(as.numeric(subset(relevant, focus.simple == 'plaque')$n.mm.2mean), as.numeric(subset(relevant, focus.simple == 'diffuse')$n.mm.2mean))

####################################################################################
####################################################################################
####################################################################################
##################                    Figure 3                    ##################
####################################################################################
####################################################################################
####################################################################################
relevant <- subset(norspots, expt.name == 'agarose' | expt.name == 'agarose_d3' | expt.name == 'd2d4ag' | expt.name == 'd3ag2')

# fig 3d. plot points over time from the agarose experiments
ggplot(subset(relevant, (condition == "vanilla" | condition == "agarose") & (donor == 2 | donor == 3 | donor == 4)), aes(time, as.numeric(norfreq), colour=factor(condition))) +
  labs(x = 'Hours post infection', y ="# of GFP+ spots", color = "Condition")+
  theme_bw() +
  scale_color_discrete(labels = c("Agarose", "Plain"), name = 'Condition:')+
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8), 
        legend.key.size = unit(1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_point(alpha = 0.6, size = 0.2) 
ggsave("C:/Users/Mark/Desktop/figure_images/3/fig3d.pdf", width = 2, height = 1.35, units = "in")

relevantt <- subset(good.spots.summary.j, (expt.name == 'agarose' | expt.name == 'agarose_d3' | expt.name == 'd2d4ag' | expt.name == 'd3ag2') & virus == 'egfp')
ggplot(subset(relevantt, (condition == "vanilla" | condition == 'agarose') & virus == 'egfp'), (aes(y=good.peak.spots, x=factor(condition)))) + 
  labs(x = 'Condition', y ="Peak # GFP+ spots", color = "Donor:")+
  geom_quasirandom(alpha = 0.8, show.legend = TRUE) +
  stat_summary(fun.y = median,
               fun.ymin = function(z) { quantile(z,0.25) }, 
               fun.ymax = function(z) { quantile(z,0.75) }, 
               geom = "pointrange",
               show.legend = FALSE,
               aes(group = factor(condition))) +
  scale_x_discrete(labels = c("Agarose", "Plain"))+
  theme_bw() + 
  ylim(0, 25000) +
  theme(axis.text = element_text(size=5), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        #legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/3/fig3e_NOoo.pdf", width = 1, height = 1.35, units = "in")
wilcox.test(subset(relevantt, condition == 'vanilla')$good.peak.spots, subset(relevantt, condition == 'agarose')$good.peak.spots)


relevantt <- subset(good.spots.summary.j, (expt.name == 'agarose' | expt.name == 'agarose_d3' | expt.name == 'd2d4ag' | expt.name == 'd3ag2') & virus == 'egfp')
ggplot(subset(relevant, (condition == "vanilla" | condition == 'agarose') & virus == 'egfp' & time == 95), (aes(y=norfreq, x=factor(condition)))) + 
  labs(x = 'Condition', y ="GFP+ spots at 4 DPI", color = "Donor:")+
  geom_quasirandom(alpha = 0.8, show.legend = TRUE) +
  stat_summary(fun.y = median,
               fun.ymin = function(z) { quantile(z,0.25) }, 
               fun.ymax = function(z) { quantile(z,0.75) }, 
               geom = "pointrange",
               show.legend = FALSE,
               aes(group = factor(condition))) +
  scale_x_discrete(labels = c("Agarose", "Plain"))+
  theme_bw() + 
  ylim(0, 25000) +
  theme(axis.text = element_text(size=5), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        #legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/3/fig3e.pdf", width = 1, height = 1.35, units = "in")
wilcox.test(subset(relevant, condition == 'vanilla' & time == 95 & virus == 'egfp')$norfreq, subset(relevant, condition == 'agarose' & time == 95 & virus == 'egfp')$norfreq)

####################################################################################
####################################################################################
####################################################################################
##################                    Figure 4                    ##################
####################################################################################
####################################################################################
####################################################################################

# fig 4d. the bee swarm.
ggplot(subset(cbfsall, (cbfsall$power > 40)), aes(as.factor(V2), cbf)) +
  xlab("KO") +
  ylim(0, 20) +
  ylab("Dominant pixel frequency (Hz)") +
  scale_x_discrete(labels = c("CYPA", "DNAH5", "DNAI1")) +
  scale_color_gradient(limits = c(0, 20)) +
  theme_bw() +
  geom_quasirandom(alpha = 0.006, size = 0.2) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange", position = position_dodge(width = 0.9), color = '#05e895', size = 0.25) +
  theme(axis.text = element_text(size=6, face = 'italic'), axis.title = element_text(size=8),
      plot.margin = margin(1, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, -10),
      legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/4/fig4d.pdf", width = 1.2, height = 2.42, units = "in")

sum(cbfsall$power > 40 & cbfsall$V2 == 'pt1')
sum(cbfsall$power > 40 & cbfsall$V2 == 'pt3')
sum(cbfsall$power > 40 & cbfsall$V2 == 'pt5')

fracbeatpt1 <- sum(cbfsall$power > 40 & cbfsall$V2 == 'pt1') / sum(cbfsall$V2 == 'pt1')
fracbeatpt3 <-sum(cbfsall$power > 40 & cbfsall$V2 == 'pt3') / sum(cbfsall$V2 == 'pt3')
fracbeatpt5 <-sum(cbfsall$power > 40 & cbfsall$V2 == 'pt5') / sum(cbfsall$V2 == 'pt5')
g1 <- subset(cbfsall, (cbfsall$power > 40 & V2 == 'pt5'))$cbf
g2 <- subset(cbfsall, (cbfsall$power > 40 & V2 == 'pt1'))$cbf
wilcox.test(g1, g2)
# fig. 4g. the spots over time.
mk <- subset(norspots, expt.name == 'ciliako2')
ggplot(mk, aes(time, as.numeric(norfreq), colour=factor(condition))) +
  #scale_color_met_d("Hiroshige", labels = c("1% agarose", "2% agarose", "Mock", "Plain"))+
  labs(x = 'Hours post infection', y ="# of GFP+ spots", color = "KO")+
  geom_point() +
  scale_color_discrete(labels = c("CYPA", "DNAH5", "DNAI1"))+
  theme_bw() + 
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_point(alpha = 0.6, size = 0.2) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6, face = 'italic'), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/4/fig4g.pdf", width = 2, height = 2.2, units = "in")


# fig 4h
ggplot(fig5_n, aes(culture, as.numeric(n))) +
  #scale_color_met_d("Hiroshige", labels = c("1% agarose", "2% agarose", "Mock", "Plain"))+
  labs(x = 'KO', y = bquote('N copies per'~mm^2))+
  geom_point() +
  scale_x_discrete(labels = c("CYPA", "DNAH5", "DNAI1"))+
  #scale_shape_discrete(labels = c("Donor 2", "Donor 3"))+
  #theme_linedraw() + 
  scale_y_log10() +
  theme_bw() +
  theme(axis.text = element_text(size=5, face = 'italic'), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6, face = 'italic'), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/4/fig4h.pdf", width = 1.5, height = 1, units = "in")



#fig 4i: correlation of good peak spots w/ fraction beating pixels
mkk <- subset(good.spots.summary.j, expt.name == 'ciliako2')

fracbeat <- c(fracbeatpt1, fracbeatpt3, fracbeatpt5)
mkk$fracbeat <- fracbeat
fracbeatpt1 <- sum(cbfsall$power > 40 & cbfsall$V2 == 'pt1') / sum(cbfsall$V2 == 'pt1')
fracbeatpt3 <-sum(cbfsall$power > 40 & cbfsall$V2 == 'pt3') / sum(cbfsall$V2 == 'pt3')
fracbeatpt5 <-sum(cbfsall$power > 40 & cbfsall$V2 == 'pt5') / sum(cbfsall$V2 == 'pt5')
medbeatpt1 <- median(subset(cbfsall, cbfsall$power > 40 & cbfsall$V2 == 'pt1')$cbf)
medbeatpt3 <- median(subset(cbfsall, cbfsall$power > 40 & cbfsall$V2 == 'pt3')$cbf)
medbeatpt5 <- median(subset(cbfsall, cbfsall$power > 40 & cbfsall$V2 == 'pt5')$cbf)
medbeat <- c(medbeatpt1, medbeatpt3, medbeatpt5)
mkk$medbeat <- medbeat
ggplot(mkk, aes(fracbeat, as.numeric(good.peak.spots))) +
  #scale_color_met_d("Hiroshige", labels = c("1% agarose", "2% agarose", "Mock", "Plain"))+
  labs(x = 'Fraction Beating Area', y ="# GFP+ spots at 5 DPI")+
  geom_point() +
  #scale_color_discrete(labels = c("1% agarose", "2% agarose", "Mock", "Plain"))+
  #scale_shape_discrete(labels = c("Donor 2", "Donor 3"))+
  #theme_linedraw() + 
  theme_bw() +
  stat_summary()+
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6, face = 'italic'), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/4/fig4i.pdf", width = 1.5, height = 1.2, units = "in")

cor.test(as.numeric(mkk$good.peak.spots), mkk$fracbeat, method = 'pearson')


####################################################################################
####################################################################################
####################################################################################
##################                    Figure 5                    ##################
####################################################################################
####################################################################################
####################################################################################

# 5a: mucus disc spinning cessation time
level_order <- c('mock', 'egfp') 
relevant <- subset(good.spots.summary.j, (condition == 'vanilla' | condition == 'static' | condition == 'shaken'))
ggplot(relevant, (aes(x = virus, y = 2*as.numeric(mucus.stop.frame)+1, color = factor(donor)))) + 
  labs(x = 'Condition', y ="Time to cessation of mucus movement \n (hours post infection)", color = "Donor:")+
  geom_quasirandom(size = 0.2) +
  scale_x_discrete(labels = c('Infected', "Mock"))+
  theme_bw() + 
  ylim(0, 150) +
  stat_summary(fun.y = median,
               fun.ymin = function(z) { quantile(z,0.25) }, 
               fun.ymax = function(z) { quantile(z,0.75) }, 
               geom = "pointrange",
               width = 0.1,
               size = .2,
               aes(group = FALSE, color = NULL)) +
  geom_hline(yintercept = 118, color = 'red') + 
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/5/fig5a.pdf", width = 1.6, height = 2, units = "in")

ok <- table(relevant$virus, (2*as.numeric(relevant$mucus.stop.frame)+1 <= 118))
chisq.test(relevant$virus, (2*as.numeric(relevant$mucus.stop.frame)+1 <= 118))

wilcox.test(subset(relevant, virus == 'mock')$mucus.stop.frame, subset(relevant, virus == 'egfp')$mucus.stop.frame)
# CBF Wrangling Time
resultbyfov_bulk <- allovem %>%
  group_by(hpi, vidpoint, donor, condition, experiment) %>%
  summarise(total = n(), powerful = sum(power > 40), ratio = powerful/total, slowpix = sum(slow), fastpix = sum(fast), domfreqsum = median(domfreq), gmean = exp(mean(log(domfreq))))

###testing something out but it was all for naught
powerful_resultbyfov <- subset(allovem, power > 40) %>%
  group_by(hpi, vidpoint, donor, condition, experiment) %>%
  summarise(total = n(), powerful = sum(power > 40), ratio = powerful/total, slowpix = sum(slow), fastpix = sum(fast), domfreqsum = median(domfreq), gmean = exp(mean(log(domfreq))))
ggplot(subset(powerful_resultbyfov, vidpoint != 34 & vidpoint != 68 & vidpoint != 65), aes(y=gmean, x = factor(hpi), color = condition, group = condition)) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16),
        legend.justification = c("right", "top"),
        #legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))+
  xlab("HPI") +
  theme_bw() +
  #ylim(0, 0.2) +
  ylab("Fraction of beating pixels per FOV") +
  geom_quasirandom(dodge.width = 0.8, alpha = 0.8) +
  stat_summary(fun.y = median,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange", position = position_dodge(width = 0.8), color = 'black')
wilcox.test(subset(resultbyfov_bulk, vidpoint != 34 & vidpoint != 68 & hpi == 72 & condition == 'infected')$powerful.frac, subset(resultbyfov_bulk, vidpoint != 34 & vidpoint != 68 & hpi == 72 & condition == 'mock')$powerful.frac)
wilcox.test(subset(resultbyfov_bulk, vidpoint != 34 & vidpoint != 68 & hpi == 48 & condition == 'infected')$powerful.frac, subset(resultbyfov_bulk, vidpoint != 34 & vidpoint != 68 & hpi == 48 & condition == 'mock')$powerful.frac)
wilcox.test(subset(resultbyfov_bulk, vidpoint != 34 & vidpoint != 68 & hpi == 24 & condition == 'infected')$powerful.frac, subset(resultbyfov_bulk, vidpoint != 34 & vidpoint != 68 & hpi == 24 & condition == 'mock')$powerful.frac)



resultbyfov_bulk$fastpix.frac <- resultbyfov_bulk$fastpix/65536
resultbyfov_bulk$slowpix.frac <- resultbyfov_bulk$slowpix/65536
resultbyfov_bulk$powerful.frac <- resultbyfov_bulk$powerful/65536

#Fig 5c: The beating fraction of the culture area.
# Slowly beating pixels don't really vary much
resultbyfov_bulk <- allovem %>%
  group_by(hpi, vidpoint, donor, condition, experiment) %>%
  summarise(total = n(), powerful = sum(power > 40), ratio = powerful/total, slowpix = sum(slow), fastpix = sum(fast), domfreqsum = median(domfreq))


ggplot(subset(resultbyfov_bulk, ratio < 1), aes(y=powerful.frac, x = factor(hpi), color = condition, group = condition)) +
  xlab("HPI") +
  theme_bw() +
  labs(x = 'Hours post infection', y ="Fraction beating pixels per FOV", color = "Condition")+
  scale_color_discrete(labels = c('Infected', 'Mock')) +
  ylim(0, 0.32) +
  ylab("Fraction of beating pixels") +
  geom_quasirandom(dodge.width = 0.8, alpha = 0.8, size = 0.2) +
  stat_summary(fun.y = median,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               size = 0.2,
               linewidth = 0.2,
               geom = "pointrange", position = position_dodge(width = 0.8), color = 'black') +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
      plot.margin = margin(1, 0, 0, 0),
      legend.key.size = unit(0.1, 'point'),
      legend.box.margin = margin(0, 0, 0, -10),
      legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/5/fig5c.pdf", width = 1.8, height = 1.35, units = "in")
wilcox.test(subset(resultbyfov_bulk, ratio < 1 & condition == 'infected' & hpi == 72)$powerful.frac, subset(resultbyfov_bulk, ratio < 1 & condition == 'mock' & hpi ==72)$powerful.frac)
wilcox.test(subset(resultbyfov_bulk, vidpoint != 34 & vidpoint != 68 & hpi == 48 & condition == 'infected')$powerful.frac, subset(resultbyfov_bulk, vidpoint != 34 & vidpoint != 68 & hpi == 48 & condition == 'mock')$powerful.frac)
wilcox.test(subset(resultbyfov_bulk, vidpoint != 34 & vidpoint != 68 & hpi == 24 & condition == 'infected')$powerful.frac, subset(resultbyfov_bulk, vidpoint != 34 & vidpoint != 68 & hpi == 24 & condition == 'mock')$powerful.frac)


# fig 5d: the cbf of all powerful pixels
# using 'beating' which i made for fig 5f down a ways
cleam = subset(allovem, power > 40 & vidpoint != 34 & vidpoint != 65 & vidpoint != 68)
ggplot(subset(beating, power > 40), aes(y=domfreq, x = factor(hpi), fill = condition)) +
  labs(x = 'Hours post infection', y ="Dominant pixel frequency (Hz)", color = "Condition")+
  ylim(0, 20) +
  scale_fill_discrete(labels = c("Infected", "Mock"), name = NULL)+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size = 0.2) +
  theme_bw() +
  # stat_summary(fun.y = mean,
  #              fun.ymin = function(x) mean(x) - sd(x), 
  #              fun.ymax = function(x) mean(x) + sd(x), 
  #              size = 0.1,
  #              linewidth = 0.1,
  #              geom = "pointrange", position = position_dodge(width = 0.9), color = 'black') +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
      plot.margin = margin(1, 0, 0, 0),
      legend.key.size = unit(4, 'point'),
      legend.box.margin = margin(0, 0, 0, -10),
      legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/5/fig5d.pdf", width = 1.8, height = 1.4, units = "in")
wilcox.test(subset(beating, hpi == 72 & condition == 'infected')$domfreq, subset(beating, hpi ==72 & condition == 'mock')$domfreq)

#Fig 5e: The pixels that are beating make up the same fxn of the
# Slowly beating pixels don't really vary much


gfpfov.bulk <- allovem %>%
  group_by(hpi, vidpoint, donor, condition, experiment, gfppos = ifelse(gfp > 900, 'yes', 'no')) %>%
  summarise(total = n(), powerful = sum(power > 40), ratio = powerful/total, slowpix = sum(slow), fastpix = sum(fast), domfreqsum = median(domfreq), gmean = exp(mean(log(domfreq))), gfp = median(gfp))
gfpfov.bulk$fastpix.frac <- gfpfov.bulk$fastpix/65536
gfpfov.bulk$slowpix.frac <- gfpfov.bulk$slowpix/65536
gfpfov.bulk$powerful.frac <- gfpfov.bulk$powerful/65536

ggplot(subset(gfpfov.4cat, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37)), aes(y=ratio, x = as.factor(hpi), color = status)) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16),
        legend.justification = c("right", "top"),
        #legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  xlab("HPI") +
  theme_bw() +
  #ylim(0, 0.3) +
  ylab("Powerful fxn") +
  #scale_color_manual(values=c('darkorchid4', 'chartreuse3'), labels = c("GFP-", "GFP+"), name = NULL)+
  #scale_color_discrete(labels = c("GFP-", "GFP+"), name = NULL)+
  geom_quasirandom(dodge.width = 0.8, alpha = 0.8) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange", position = position_dodge(width = 0.8), aes(group = status), color = 'black')

t.test(subset(gfpfov.onlypos, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'infected' & hpi == 24)$ratio, subset(gfpfov.onlypos, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'mock' & hpi == 24)$ratio)
t.test(subset(gfpfov.onlypos, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'infected' & hpi == 48)$ratio, subset(gfpfov.onlypos, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'mock' & hpi == 48)$ratio)
t.test(subset(gfpfov.onlypos, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'infected' & hpi == 72)$ratio, subset(gfpfov.onlypos, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'mock' & hpi == 72)$ratio)

wilcox.test(subset(gfpfov.bulk, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'infected' & hpi == 24)$ratio, subset(gfpfov.bulk, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'mock' & hpi == 24)$ratio)
wilcox.test(subset(gfpfov.bulk, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'infected' & hpi == 48)$ratio, subset(gfpfov.bulk, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'mock' & hpi == 48)$ratio)
wilcox.test(subset(gfpfov.bulk, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'infected' & hpi == 72)$ratio, subset(gfpfov.bulk, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'mock' & hpi == 72)$ratio)



#NEW 5F OR WHATEVER
gfpfov.onlypos <- allovem %>%
  group_by(hpi, vidpoint, donor, condition, experiment) %>%
  mutate(gfppos = ifelse(gfp > 900, 'yes', 'no')) %>%
  group_by(hpi, vidpoint, donor, condition, experiment, gfppos) %>%
  summarise(total = n(), powerful = sum(power > 40), ratio = powerful/total, slowpix = sum(slow), fastpix = sum(fast), domfreqsum = median(domfreq), gmean = exp(mean(log(domfreq))), gfp = median(gfp)) %>%
  filter(all(c('yes', 'no') %in% gfppos))

gfpfov.onlypos$fastpix.frac <- gfpfov.onlypos$fastpix/65536
gfpfov.onlypos$slowpix.frac <- gfpfov.onlypos$slowpix/65536
gfpfov.onlypos$powerful.frac <- gfpfov.onlypos$powerful/65536

relevant <- subset(gfpfov.onlypos, condition == 'infected')

ggplot(subset(relevant), aes(y=ratio, x = gfppos, group = vidpoint, color = gfppos)) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16),
        legend.justification = c("right", "top"),
        #legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  xlab("GFP+") +
  theme_bw() +
  ylim(0, 0.3) +
  ylab("Fraction of pixels above power threshold") +
  geom_line(color = 'black', alpha = 0.5) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values=c('darkorchid4', 'chartreuse3'), labels = c("GFP-", "GFP+"), name = NULL)+
  #scale_x_discrete(labels = c("CYPA", "DNAH5", "DNAI1")) +
  #scale_color_gradient(limits = c(0, 5000)) +
  facet_grid(rows = vars(hpi))
t.test(subset(gfpfov.onlypos, vidpoint != 34 & vidpoint != 68 & hpi == 72 & condition == 'infected' & gfppos == 'no')$ratio, subset(gfpfov.onlypos, vidpoint != 34 & vidpoint != 68 & hpi == 72 & condition == 'infected' & gfppos == 'yes')$ratio)  # Replace group1 and group2 with your actual data
wilcox.test(subset(gfpfov.onlypos, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'infected' & hpi == 72 & gfppos == 'no')$ratio, subset(gfpfov.onlypos, (vidpoint != 34) & (vidpoint != 68) & (vidpoint != 65) & (vidpoint != 37) & condition == 'infected' & hpi == 72 & gfppos == 'yes')$ratio)


deltas <- allovem %>%
  group_by(hpi, vidpoint, donor, condition, experiment) %>%
  mutate(gfppos = ifelse(gfp > 900, 'yes', 'no')) %>%
  group_by(hpi, vidpoint, donor, condition, experiment, gfppos) %>%
  summarise(total = n(), powerful = sum(power > 40), ratio = powerful/total, greenpowerful = sum(subset(poweslowpix = sum(slow), fastpix = sum(fast), domfreqsum = median(domfreq), gmean = exp(mean(log(domfreq))), gfp = median(gfp)) %>%
  filter(all(c('yes', 'no') %in% gfppos))


hmm <- allovem %>%
  group_by(hpi, vidpoint, donor, condition, experiment) %>%
  mutate(gfppos = ifelse(gfp > 900, 'yes', 'no')) %>%
  filter(all(c('yes', 'no') %in% gfppos)) %>%
  summarise(total = n(), powerful = sum(power > 40), ratio = powerful/total, slowpix = sum(slow), fastpix = sum(fast), domfreqsum = median(domfreq), gmean = exp(mean(log(domfreq))), gfppos, gfp = median(gfp)) %>%
  pivot_wider(names_from = gfppos, values_from = ratio, names_prefix = "ratio_") %>%
  mutate(delta_ratio = ratio_yes - ratio_no)

spread(gfppos, ratio) %>%
  mutate(delta_ratio = yes - no)

ggplot(subset(relevant), aes(y=ratio, x = hpi, group = vidpoint, color = gfppos)) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16),
        legend.justification = c("right", "top"),
        #legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  xlab("GFP+") +
  theme_bw() +
  ylim(0, 0.3) +
  ylab("Fraction of pixels above power threshold") +
  geom_line(color = 'black') +
  geom_point() +
  scale_color_manual(values=c('darkorchid4', 'chartreuse3'), labels = c("GFP-", "GFP+"), name = NULL)+
  #scale_x_discrete(labels = c("CYPA", "DNAH5", "DNAI1")) +
  #scale_color_gradient(limits = c(0, 5000)) +
  facet_grid(rows = vars(hpi))
wilcox.test(subset(gfpfov, vidpoint != 34 & vidpoint != 68 & hpi == 24 & condition == 'infected' & gfppos == 'no')$ratio, subset(gfpfov, vidpoint != 34 & vidpoint != 68 & hpi == 24 & condition == 'infected' & gfppos == 'yes')$ratio)  # Replace group1 and group2 with your actual data



ggplot(subset(cleam, power > 40 & condition == 'infected' & hpi != 24), aes(y=domfreq, x = factor(hpi), fill = (gfp > 900))) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16),
        legend.justification = c("right", "top"),
        #legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))+
  xlab("HPI") +
  ylim(0, 20) +
  ylab("Dominant pixel frequency (Hz)") +
  theme_bw() +
  scale_fill_manual(values=c('darkorchid4', 'chartreuse3'), labels = c("GFP-", "GFP+"), name = NULL)+
  #scale_x_discrete(labels = c("CYPA", "DNAH5", "DNAI1")) +
  #scale_color_gradient(limits = c(0, 5000)) +
  #facet_grid(rows = vars(experiment)) +
  #geom_boxplot(outlier.shape = NA) +
  #geom_point(position = position_jitterdodge()) +
  geom_violin() +
  #geom_quasirandom(alpha = 0.05, dodge_width = 0.8) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange", position = position_dodge(width = 0.9), color = 'black')


whatever <- sample_n(cleam, 600000)




####################################################################################
####################################################################################
####################################################################################
##################          Figure 5 NEW AND IMPROVED!            ##################
####################################################################################
####################################################################################
####################################################################################
gfpfov.onlypost <- allovem %>%
  group_by(video, hpi, vidpoint, donor, condition, experiment) %>%
  mutate(gfppos = ifelse(gfp > 900, 'yes', 'no')) %>%
  group_by(video, hpi, vidpoint, donor, condition, experiment, gfppos) %>%
  summarise(total = n(), powerful = sum(power > 40), ratio = powerful/total) %>%
  filter(all(c('yes', 'no') %in% gfppos))
gfpfov.onlypost <- gfpfov.onlypost %>% mutate(status =
                                              case_when(gfppos == 'yes' ~ 'infected', 
                                                        gfppos == 'no' ~ "near")
)

gfpfov.onlynegt <- allovem %>%
  group_by(video, hpi, vidpoint, donor, condition, experiment) %>%
  mutate(gfppos = ifelse(gfp > 900, 'yes', 'no')) %>%
  group_by(video, hpi, vidpoint, donor, condition, experiment, gfppos) %>%
  summarise(total = n(), powerful = sum(power > 40), ratio = powerful/total) %>%
  filter(all(!c('yes') %in% gfppos))
gfpfov.onlynegt$status <- 'distant'

gfpfov.all <- rbind(gfpfov.onlypost, gfpfov.onlynegt)
gfpfov.4cat <- gfpfov.all %>% mutate(status = case_when(condition == 'mock' ~ 'mock',
                                                        condition == 'infected' ~ status))
clean <- subset(gfpfov.4cat, ratio < 1 & total > 20)
clean$jointstatus <- paste(clean$hpi, clean$status) 
level_order <- c('infected', 'near', 'distant', 'mock') 
ggplot(clean, aes(y=ratio, x = as.factor(hpi), color = factor(status, levels = level_order))) +
  labs(x = 'Hours post infection', y ="Beating fraction", color = NULL)+
  theme_bw() +
  ylim(0, 0.4) +
  scale_color_manual(values=c('chartreuse3', 'darkorchid4', 'darkorchid2', 'deeppink1'), labels = c("GFP+", "Near GFP-", "Far GFP-", "Mock"), name = NULL)+
  geom_quasirandom(dodge.width = 0.8, alpha = 0.8, size =0.2, show.legend = FALSE) +
  #facet_grid(rows = vars(donor), cols = vars(experiment)) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               size = 0.2,
               linewidth = 0.2,
               geom = "pointrange", position = position_dodge(width = 0.8), aes(group = factor(status, levels = level_order)), color = 'black') +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
      plot.margin = margin(1, 1, 0, 0),
      legend.key.size = unit(4, 'point'),
      legend.box.margin = margin(0, 0, 0, -10),
      legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/5/fig5e.pdf", width = 1.8, height = 1.4, units = "in")




# Stats
#is it normal?
hist(log(clean$ratio), breaks = 50)
qqnorm(clean$ratio)
qqline(clean$ratio)
shapiro.test(clean$ratio)
#nope not at all
# is there any difference?
kruskal.test(ratio ~ jointstatus, data = clean)
# ya. but where?
library(dunn.test)
dunn.test(clean$ratio, clean$jointstatus, method = 'bh')
tclean = subset(clean, hpi == 24)
midclean = subset(clean, hpi == 48)
lateclean = subset(clean, hpi == 72)
dunn.test(tclean$ratio, tclean$status, method = 'bh')
dunn.test(midclean$ratio, midclean$status, method = 'bh')
dunn.test(lateclean$ratio, lateclean$status, method = 'bh')
# most places tbh


########################################################
##########        CBFS for each status        ##########
########################################################
allovem <- read.csv("D:/cbf/allovem_norux.csv")

labelled.pixels <- allovem %>%
  mutate(
    status = case_when(
      condition == 'mock' ~ 'mock',
      video %in% gfpfov.onlynegt$video ~ 'distant',
      video %in% gfpfov.onlypost$video & gfp > 900 ~ 'infected',
      video %in% gfpfov.onlypost$video & gfp <= 900 ~ 'near',
      TRUE ~ NA_character_  # To handle cases when no condition is met
    )
  )
bad.videos <- c('mb_220618_ruxcilia_t5cilia_65_R3D.dv', 'mb_220516_cilia_t7_10x_ciliavid_34_R3D.dv', 'mb_220516_cilia_t6_10x_ciliavid_68_R3D.dv')

good.pixels <- labelled.pixels %>%
  filter(!video %in% bad.videos)

beating = subset(good.pixels, power > 40)
level_order <- c('infected', 'near', 'distant', 'mock') 
ggplot(subset(beating), aes(y=domfreq, x = factor(hpi), fill = factor(status, levels = level_order))) +
  labs(x = 'Hours post infection', y ="Dominant pixel frequency (Hz)", color = NULL)+
  #ylim(0, 20) +
  theme_bw() +
  scale_fill_manual(values=c('chartreuse3', 'darkorchid4', 'darkorchid2', 'deeppink1'), labels = c("GFP+", "Near GFP-", "Far GFP-", "Mock"), name = NULL)+
  #facet_grid(rows = vars(experiment)) +
  #geom_boxplot(outlier.shape = NA) +
  #geom_point(position = position_jitterdodge()) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size = 0.2) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(4, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/5/fig5f.pdf", width = 2.4, height = 1.4, units = "in")






fbeat = subset(beating, hpi == 48)
kruskal.test(domfreq ~ status, data = fbeat)
dunn.test(fbeat$domfreq, fbeat$status, method = 'hochberg')
ks.test(subset(beating, hpi == '72' & status == 'infected')$domfreq, subset(beating, hpi == '72' & status == 'near')$domfreq)

####################################################################################
####################################################################################
####################################################################################
##################                    Figure 6                    ##################
####################################################################################
####################################################################################
####################################################################################
d7d8 <- norspots
# # gfp+ spots keeps rising after rinse. 6
rinsed = subset(d7d8, virus == 'egfp')
ggplot(rinsed, aes(time, as.numeric(norfreq), colour=factor(donor))) +
  labs(x = 'Hours post infection', y ="# of GFP+ spots", color = "Donor")+
  geom_point() +
  theme_bw() + 
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_point(alpha = 0.6, size = 0.01) +
  geom_vline(xintercept = 120, color = 'black', linetype = 'dashed') + 
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/6/fig6c.pdf", width = 2, height = 2.2, units = "in")

