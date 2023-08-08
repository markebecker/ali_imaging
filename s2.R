ggplot(subset(good.spots.summary.j, virus == "egfp" & (condition == 'vanilla' | condition == 'shaken' | condition == 'static')), (aes(x=donor, y= good.peak.spots,  color = factor(dyes)))) + 
  labs(y = bquote('N copies per '~mm^2~''), x ="# GFP+ Spots", color = "Donor:")+
  #bquote('N copies per '~mm^2' culture')
  geom_point(show.legend = TRUE) +
  scale_y_log10()+
  theme_bw() + 
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
rell <- subset(good.spots.summary.j, virus == 'egfp')
cor.test(as.numeric(rell$X20hpi.spots), as.numeric(rell$n.mm.2mean), use = "complete.obs", method = 'pearson')

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

#supfig 2. is cmo vs. spytub different?
# 2a: peak timing
ggplot(subset(good.spots.summary.j, (condition == "vanilla" | condition == 'static' | condition == 'shaken') & (donor == 2 | donor == 3 | donor == 4) & (dyes == 'spy + nuc' | dyes == 'cmo') & (2*video.end.frame-1 > good.peak.spots.time) & virus == 'egfp'), (aes(y=good.peak.spots.time, x=factor(dyes), color = factor(donor)))) + 
  labs(x = 'Dye', y ="Time of Peak GFP+ Spots (HPI)", color = "Donor")+
  stat_summary(fun.y = median,
               fun.ymin = function(z) { quantile(z,0.25) }, 
               fun.ymax = function(z) { quantile(z,0.75) }, 
               color = 'black',
               linewidth = 0.5,
               size = 0.3,
               geom = "crossbar") +
  geom_quasirandom() +
  scale_x_discrete(labels = c("CMO", "SPYtub + NV"))+
  theme_bw() + 
  ylim(0, 200) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(0, -5, 0, 0),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.1, 'in'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/supfigs/s2/s2a.pdf", width = 2, height = 2, units = "in")

rel <-subset(good.spots.summary.j, (condition == "vanilla" | condition == 'static' | condition == 'shaken') & (donor == 2 | donor == 3 | donor == 4) & (dyes == 'spy + nuc' | dyes == 'cmo') & (2*video.end.frame-1 > good.peak.spots.time) & virus == 'egfp')
wilcox.test(subset(rel, dyes == 'cmo')$good.peak.spots.time, subset(rel, dyes == 'spy + nuc')$good.peak.spots.time)


# 2b: # gfp spots at 20 hpi
ggplot(subset(good.spots.summary.j, (condition == "vanilla" | condition == 'static' | condition == 'shaken') & (donor == 2 | donor == 3 | donor == 4) & (dyes == 'spy + nuc' | dyes == 'cmo') & virus == 'egfp'), 
       (aes(y=X20hpi.spots, x=factor(dyes), color = factor(donor)))) + 
  labs(x = 'Dyes', y ="GFP+ Spots at 20 HPI", color = "Donor")+
  stat_summary(fun.y = median,
               fun.ymin = function(z) { quantile(z,0.25) }, 
               fun.ymax = function(z) { quantile(z,0.75) }, 
               color = 'black',
               linewidth = 0.5,
               size = 0.3,
               geom = "crossbar") +
  geom_quasirandom() +
  #scale_color_continuous(high = '#05e895', low = '#402375', limits = c(NA, 20)) +
  scale_x_discrete(labels = c("CMO", "SPYtub + NV"))+
  theme_bw() + 
  ylim(0, 350) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(0, -5, 0, 0),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.1, 'in'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/supfigs/s2/s2b.pdf", width = 2, height = 2, units = "in")

rel <- subset(good.spots.summary.j, (condition == "vanilla" | condition == 'static' | condition == 'shaken') & (donor == 2 | donor == 3 | donor == 4) & (dyes == 'spy + nuc' | dyes == 'cmo') & virus == 'egfp')
wilcox.test(subset(rel, dyes == 'cmo')$X20hpi.spots, subset(rel, dyes == 'spy + nuc')$X20hpi.spots)




# 2c: # gfp spots at peak
ggplot(subset(good.spots.summary.j, (2*(video.end.frame)-1 > good.peak.spots.time) & (condition == "vanilla" | condition == 'static' | condition == 'shaken') & (donor == 2 | donor == 3 | donor == 4) & (dyes == 'spy + nuc' | dyes == 'cmo') & virus == 'egfp'), 
       (aes(y=good.peak.spots, x=factor(dyes), color = factor(donor)))) + 
  labs(x = 'Dyes', y ="Peak # of GFP+ Spots", color = "Donor")+
  stat_summary(fun.y = median,
               fun.ymin = function(z) { quantile(z,0.25) }, 
               fun.ymax = function(z) { quantile(z,0.75) }, 
               color = 'black',
               linewidth = 0.5,
               size = 0.3,
               geom = "crossbar") +
  geom_quasirandom() +
  #scale_color_continuous(high = '#05e895', low = '#402375', limits = c(NA, 20)) +
  scale_x_discrete(labels = c("CMO", "SPYtub + NV"))+
  theme_bw() + 
  ylim(0, 38000) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(0, -5, 0, 0),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.1, 'in'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/supfigs/s2/s2c.pdf", width = 2, height = 2, units = "in")

rel <- subset(good.spots.summary.j,  (2*(video.end.frame)-1 > good.peak.spots.time) & (condition == "vanilla" | condition == 'static' | condition == 'shaken') & (donor == 2 | donor == 3 | donor == 4) & (dyes == 'spy + nuc' | dyes == 'cmo') & virus == 'egfp')
wilcox.test(subset(rel, dyes == 'cmo')$good.peak.spots, subset(rel, dyes == 'spy + nuc')$good.peak.spots)


# 2d: mcc cessation time
ggplot(subset(good.spots.summary.j, (video.end.frame > mucus.stop.frame) & (condition == "vanilla" | condition == 'static' | condition == 'shaken') & (donor == 2 | donor == 3 | donor == 4) & (dyes == 'spy + nuc' | dyes == 'cmo') & virus == 'egfp'), 
       (aes(y=2*(mucus.stop.frame-1)+1, x=factor(dyes), color = factor(donor)))) + 
  labs(x = 'Dyes', y ="MCC Cessation Time (HPI)", color = "Donor")+
  stat_summary(fun.y = median,
               fun.ymin = function(z) { quantile(z,0.25) }, 
               fun.ymax = function(z) { quantile(z,0.75) }, 
               color = 'black',
               linewidth = 0.5,
               size = 0.3,
               geom = "crossbar") +
  geom_quasirandom() +
  #scale_color_continuous(high = '#05e895', low = '#402375', limits = c(NA, 20)) +
  scale_x_discrete(labels = c("CMO", "SPYtub + NV"))+
  theme_bw() + 
  ylim(0, 140) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(0, -5, 0, 0),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.1, 'in'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/supfigs/s2/s2d.pdf", width = 2, height = 2, units = "in")

rel <- subset(good.spots.summary.j, (video.end.frame > mucus.stop.frame) & (condition == "vanilla" | condition == 'static' | condition == 'shaken') & (donor == 2 | donor == 3 | donor == 4) & (dyes == 'spy + nuc' | dyes == 'cmo') & virus == 'egfp')
wilcox.test(subset(rel, dyes == 'cmo')$mucus.stop.frame, subset(rel, dyes == 'spy + nuc')$mucus.stop.frame)


