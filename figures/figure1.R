# figure 1
# revised

#libraries
library(readr)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

######################
## 1d. % ciliation. ##
######################
goodspots_j <- read_csv("D:/000_rebuttal/240606_goodspots_summary.csv")
dataa = subset(goodspots_j, (condition == 'vanilla' | condition == 'shaken' | condition == 'static') 
               & (dyes != "cmo" & dyes != 'undyed') 
               & (usable == 1)
               # these experiments had spytub added like right before imaging instead of >48 hours preinfection
               & (expt.name != 'd2shake') & (expt.name != 'd1d9') & expt.name != 'multitest')

ggplot(dataa, aes(x=as.numeric(dataa$t0spy350)/as.numeric(dataa$t0spy0))) +
  geom_histogram(bins = 30) +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  labs(x = "Fraction Tubulin+\n Area",
       y = "Count") +
  theme(axis.text = element_text(size=6), 
        axis.title = element_text(size=6), 
        plot.margin = margin(5, 10, 0, 5))
ggsave("C:/Users/Mark/Desktop/figure_images/1/fig1dnew.pdf", width = 0.9, height = 1.1, units = "in")

median(as.numeric(dataa$t0spy350)/as.numeric(dataa$t0spy0))

#############################
## 1g. N copies over time. ##
#############################

qpcr_j <- read_csv("C:/Users/Mark/Desktop/figure_images/1/qpcr_p3oldandflippies.csv")
keys <- colnames(qpcr_j)[!grepl('copies.mm.2',colnames(qpcr_j))]
X <- as.data.table(qpcr_j)
X$copies.mm.2 <- as.numeric(X$copies.mm.2)
qpcrcurve = X[,list(mm= mean(copies.mm.2), sd = sd(copies.mm.2), hpi = hpi, type = type, donor = donor),'sample.name']
trimmedq = distinct(qpcrcurve)

#p 0.054 at 48hpi for grouped
t.test(subset(trimmedq, hpi == 1)$mm, subset(trimmedq, hpi == 48)$mm)
#p 0.0035 at 72hpi for grouped
t.test(subset(trimmedq, hpi == 168)$mm, subset(trimmedq, hpi == 72)$mm)
# no significant differences between timepoints with matched inverted and conventional points
t.test(subset(trimmedq, hpi == 24 & type == "inv")$mm, subset(trimmedq, hpi == 24 & type == "conv")$mm)
t.test(subset(trimmedq, hpi == 48 & type == "inv")$mm, subset(trimmedq, hpi == 48 & type == "conv")$mm)
t.test(subset(trimmedq, hpi == 72 & type == "inv")$mm, subset(trimmedq, hpi == 72 & type == "conv")$mm)
t.test(subset(trimmedq, hpi == 144 & type == "inv")$mm, subset(trimmedq, hpi == 144 & type == "conv")$mm)
t.test(subset(trimmedq, hpi == 168 & type == "inv")$mm, subset(trimmedq, hpi == 168 & type == "conv")$mm)

t.test(subset(trimmedq, hpi == 168)$mm, subset(trimmedq, hpi == image)$mm)

shapiro.test(subset(trimmedq, hpi == 120)$mm)

#mean of the final timepoint
mean(subset(trimmedq, hpi == 168)$mm)

#the plot
colors <- c("inv" = "deeppink2", "conv" = "chartreuse3")
ggplot(trimmedq, aes(x=hpi, y=mm, color = factor(type))) + 
  labs(x = 'Hours Post Infection', y = bquote('N copies per'~mm^2), color = "ALI type")+
  scale_y_log10(limits = c(-1, 3e9)) +
  geom_point(alpha = 0.4, size = 1, position = position_dodge(width = 1)) +
  scale_color_manual(values = colors, labels = c("conv" = "Conventional", "inv" = "Inverted")) +
  #geom_smooth(method = "loess", se = FALSE)+
  stat_summary(aes(group = factor(type)), 
               fun = mean, 
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               size = 0.25, 
               alpha = 0.8,
               position = position_dodge(width = 10),
               linewidth = 0.5, 
               geom = "pointrange") +
  theme_bw() + 
  #geom_hline(yintercept = 118, color = 'red') + 
  theme(axis.text = element_text(size=6), 
        axis.title = element_text(size=8), 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.5, 'lines'),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, 'lines'),
        legend.box.margin = margin(-8, -8, -8, -8),
        #legend.justification = c("right", "bottom"),
        legend.position = c(0.8, 0.15),
        plot.margin = margin(5, 10, 0, 5))
ggsave("C:/Users/Mark/Desktop/figure_images/1/fig1f.pdf", width = 3.1, height = 1.85, units = "in")


############################
## 1j & i. line profiles. ##
############################

cell <- read.csv("C:/Users/Mark/Desktop/figure_images/1/PUNCTA/revised/Values.csv")
cellslong <- cell %>%
  pivot_longer(cols = c("N", "S", "dsRNA"), 
               names_to = "channel", 
               values_to = "pixel_value")
colors <- c("N" = "deeppink1", "S" = "chartreuse3", "dsRNA" = "blue", "nuclei" = "grey")

ggplot(cellslong, aes(x = distance, y = pixel_value, color = channel)) +
  geom_line() +
  scale_color_manual(values = colors) +
  coord_flip() +
  scale_y_continuous(breaks = c(0, 50000)) +
  labs(x = "Distance (µm)",
       y = "Pixel Value",
       color = "Channel") +
  theme_bw() + 
  #geom_hline(yintercept = 118, color = 'red') + 
  theme(axis.text = element_text(size=6), 
        axis.title = element_text(size=8), 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.5, 'lines'),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, 'lines'),
        legend.box.margin = margin(-8, -8, -8, -8),
        plot.margin = margin(5, 10, 0, 5))
ggsave("C:/Users/Mark/Desktop/figure_images/1/fig1g.pdf", width = 2, height = 1.15, units = "in")

point <- read.csv("C:/Users/Mark/Desktop/figure_images/1/PUNCTA/revised/pointvalues.csv")
pointslong <- point %>%
  pivot_longer(cols = c("N", "S", "dsRNA"), 
               names_to = "channel", 
               values_to = "pixel_value")
colors <- c("N" = "deeppink1", "S" = "chartreuse3", "dsRNA" = "blue", "nuclei" = "grey")

ggplot(pointslong, aes(x = distance, y = pixel_value, color = channel)) +
  geom_line() +
  scale_color_manual(values = colors) +
  coord_flip() +
  xlim(0, 4) +
  labs(x = "Distance (µm)",
       y = "Pixel Value",
       color = "Channel") +
  theme_bw() + 
  scale_y_continuous(breaks = c(0, 3000)) +
  #geom_hline(yintercept = 118, color = 'red') + 
  theme(axis.text = element_text(size=6), 
        axis.title = element_text(size=8), 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.5, 'lines'),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, 'lines'),
        legend.box.margin = margin(-8, -8, -8, -8),
        plot.margin = margin(5, 10, 0, 5))
ggsave("C:/Users/Mark/Desktop/figure_images/1/fig1gb.pdf", width = 2, height = 1.15, units = "in")

#####################################
## 1 Text Furrow & cyst prevalence ##
#####################################
furrowcount <- table(dataa$crypts)
print(cryptcount)
cystcount <- table(dataa$cyst)
print(cystcount)

jointcryptcount <- table(dataa$crypts, dataa$cyst)
print(jointcryptcount)

########################################################
## S4. Unjamming prevalence; jamming x crypts x mucus ##
########################################################

exptlog_iii <- read_csv("D:/000_rebuttal/exptlog_ii.csv")
spots <- read_csv("D:/000_rebuttal/240702_goodspots.csv")

# take the spots over time csv and extract specific values of interest
daily.spots.summary <- spots %>%
  group_by(csv) %>%
  summarise(good.peak.spots = max(norfreq), 
            good.peak.spots.time = time[which.max(norfreq)], 
            onedayspots = norfreq[time == 23],
            twodayspots = norfreq[time == 47],
            threedayspots = norfreq[time == 71],
            fourdayspots = norfreq[time == 95],
            fivedayspots = norfreq[time == 119],
            aucspots = sum(norfreq)) 
# cross back with the info i want
daily.spots.summary.j <- inner_join(daily.spots.summary, exptlog_iii)


daily.spots.summary.j <- daily.spots.summary.j %>%
  mutate(jammed = case_when(
    jammed == 'jammed' ~ 'Jammed',
    jammed == 'edge' ~ 'Edge Migration',
    jammed == 'unjammed' ~ 'Widespread Migration',
    jammed == 'hypermobile' ~ 'Hypermobile',
    TRUE ~ jammed
  ))

# this is the correct dataset- including the KOs, I have inconsistent level names for crypts & cysts.
dataa = subset(daily.spots.summary.j, ((condition == 'vanilla') | (condition == 'shaken') | (condition == 'static') | (condition == 'rinsed') & (usable == 1)))
dataa$jammed <- factor(dataa$jammed, levels = c("Jammed", "Edge Migration", "Widespread Migration", "Hypermobile"))
dataa <- dataa %>%
  mutate(crypts = case_when(
    crypts == 'clean' ~ 'No Furrows',
    crypts == 'edgefurrow' ~ 'Edge Furrows',
    crypts == 'furrows' ~ 'Center Furrows',
    crypts == 'pointycrypts' ~ 'No Furrows'
  )) %>%
  mutate(mucus.simple = case_when(
    mucus.simple == 'rotary' ~ 'Rotary',
    mucus.simple == 'disorganized' ~ 'Disorganized'))
dataa$crypts <- factor(dataa$crypts, levels = c("No Furrows", "Edge Furrows", "Center Furrows"))

ggplot(dataa, aes(x = jammed, y = as.numeric(piv.median.speed), color = as.factor(crypts))) +
  geom_jitter() + 
  geom_hline(yintercept = 3, linetype = 'dashed', color = 'black') +
  labs(x = "Jamming Category (over entire imaging period)",
       y = "Median SPY650-tubulin speed\nat 4 HPI (um/hr)",
       color = "Crypts:") +
  theme_bw() + 
  #geom_hline(yintercept = 118, color = 'red') + 
  theme(axis.text = element_text(size=5), 
        axis.title = element_text(size=8), 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.5, 'lines'),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, 'lines'),
        legend.box.margin = margin(-8, -8, -8, -8),
        legend.position = c(0.2, 0.8),
        plot.margin = margin(5, 10, 0, 5))
ggsave("a.pdf", width = 3.2, height = 3, units = "in")
largs <- list(set_varnames = list(jammed = "Jammed Axis Title", mucus.simple = "Mucus Simple Axis Title"))
library(vcd)


table_data <- table(dataa$jammed, dataa$mucus.simple)
largs <- list(
  set_varnames = c(jammed = "Migration Extent", mucus.simple = "Mucus Pattern")
)
largs <- list(
  set_varnames = c(jammed = "Migration Extent", mucus.simple = "Mucus Pattern", crypts = 'Furrows'),
  rot_labels = c(0, 20, 0, 90),  # Adjust the rotation angles (top, right, bottom, left)
  label_args = list(
    cex = 0.2,  # Adjust the size of the labels
    font = 2    # Font style (1 = plain, 2 = bold, etc.)
  )
)

mosaic(~ jammed + mucus.simple, data = dataa,  
       labeling_args = largs, 
       shade = TRUE, 
       #gp = shading_max,
       legend = TRUE)
# save as pdf 7.5 x 7.5 in
mosaic(~ jammed + crypts, data = dataa,  
       labeling_args = largs, 
       shade = TRUE,
       #gp = shading_max,
       legend = TRUE)
mosaic(~ jammed + crypts, data = dataa, shade = TRUE, legend = TRUE, rot_labels = c(0, 20, 0, 90))
mosaic(~ jammed + crypts, data = dataa, shade = TRUE, legend = TRUE, rot_labels = c(0, 20, 0, 90))
mosaic(~ jammed + mucus.simple, data = dataa, shade = TRUE, legend = TRUE, rot_labels = c(0, 20, 90, 90))


##############################################################################
## How jammed are these puppies? How crypt-ridden? How do these impact MCC? ##
##############################################################################
library(readr)
spots <- read_csv("D:/000_rebuttal/240702_goodspots.csv")
# for info about each culture
exptlog_iii <- read_csv("D:/000_rebuttal/exptlog_ii.csv")
# take the spots over time csv and extract specific values of interest
daily.spots.summary <- spots %>%
  group_by(csv) %>%
  summarise(good.peak.spots = max(norfreq), 
            good.peak.spots.time = time[which.max(norfreq)],
            aucspots = sum(norfreq)) 
# cross back with the info i want
daily.spots.summary.j <- inner_join(daily.spots.summary, exptlog_iii)
daily.spots.summary.j <- daily.spots.summary.j %>%
  mutate(crypts = case_when(
    crypts == 'pointycrypts' ~ 'clean',
    TRUE ~ crypts
  ))
daily.spots.summary.j <- daily.spots.summary.j %>%
  mutate(mucus.simple = case_when(
    mucus.simple == 'static' ~ 'disorganized',
    TRUE ~ mucus.simple
  ))

##############################################
# how crypt ridden?
library(vcd)

dataa = subset(daily.spots.summary.j, ((condition == 'vanilla') | (condition == 'shaken') | (condition == 'static') | (condition == 'rinsed') & (usable == 1)))

# comparing crypts and jammedness, p = 2.2e-16
contingency_table <- table(dataa$crypts, dataa$jammed)
chi_squared_test <- chisq.test(contingency_table)
print(chi_squared_test)
print(contingency_table)

# comparing cysts and jammedness, p = 0.12
contingency_table <- table(dataa$cyst, dataa$jammed)
chi_squared_test <- chisq.test(contingency_table)
print(chi_squared_test)
print(contingency_table)

# comparing mucus & jammedness, p = 5e-6
contingency_table <- table(dataa$mucus.simple, dataa$jammed)
chi_squared_test <- chisq.test(contingency_table)
print(chi_squared_test)
print(contingency_table)

# comparing mucus & crypts, p = 6e-5
contingency_table <- table(dataa$mucus.simple, dataa$crypts)
chi_squared_test <- chisq.test(contingency_table)
print(chi_squared_test)
print(contingency_table)

# so it looks like jammedness, crypts, and mucus motion pattern are related somehow. how?
# log linear analysis for what is up with all this.
dataa$crypts <- as.factor(dataa$crypts)
dataa$jammed <- as.factor(dataa$jammed)
dataa$mucus.simple <- as.factor(dataa$mucus.simple)
contingency_table <- table(dataa$crypts, dataa$jammed, dataa$mucus.simple)

# Fit the log-linear model
log_linear_model <- loglm(~ 1 * 2 * 3, data = contingency_table)

# Print the summary of the log-linear model
summary(log_linear_model)

ff <- dataa %>%
  group_by(jammed, mucus.simple, crypts) %>%
  summarise(freq = n())
glmtest <- glm(freq ~ mucus.simple + crypts + jammed + mucus.simple:crypts + mucus.simple:jammed + crypts:jammed, data = ff, family=poisson)

mosaic(glmtest, residuals_type ='pearson', shade = TRUE, legend = TRUE)

contingency_table <- table(dataa$jammed, dataa$crypts, dataa$mucus.simple)

dimnames(contingency_table) <- list(
  'jammed',
  'crypts',
  'mucus.simple')
)


dataa$crypts <- as.factor(dataa$crypts)
dataa$jammed <- as.factor(dataa$jammed)
dataa$mucus.simple <- as.factor(dataa$mucus.simple)

# Create the contingency table
contingency_table <- table(dataa$crypts, dataa$jammed, dataa$mucus.simple)

# Set the dimension names explicitly
dimnames(contingency_table) <- list(
  crypts = levels(dataa$crypts),
  jammed = levels(dataa$jammed),
  mucus.simple = levels(dataa$mucus.simple)
)

log_linear_model <- loglm(~ crypts + jammed + mucus.simple, data = contingency_table)

mosaic(log_linear_model, residuals_type = 'pearson', shade = TRUE, legend = TRUE)


# Fit a log-linear model
log_linear_model <- loglm(~ crypts + jammed + mucus.simple + crypts:jammed + crypts:mucus.simple + jammed:mucus.simple, data = contingency_table)

# Print the summary of the log-linear model
summary(log_linear_model)

# Load the necessary library
install.packages('vcd')
library(vcd)

ggplot(dataa, aes(x = crypts, y = as.numeric(piv.median.speed), color = mucus.simple)) +
  geom_jitter() + 
  geom_hline(yintercept = 3)

 # Create a mosaic plot
dataa$crypts <- factor(dataa$crypts, levels = c("clean", "edgefurrow", "furrows"))
dataa$jammed <- factor(dataa$jammed, levels = c("jammed", "edge", "unjammed", 'hypermobile'))
mosaic(~ jammed + mucus.simple + crypts, data = dataa, shade = TRUE, legend = TRUE, rot_labels = c(0, 20, 0, 90))
mosaic(~ jammed + crypts, data = dataa, shade = TRUE, legend = TRUE, rot_labels = c(0, 20, 0, 90))
mosaic(~ jammed + mucus.simple, data = dataa, shade = TRUE, legend = TRUE, rot_labels = c(0, 20, 90, 90))
mosaic(~ crypts + mucus.simple, data = dataa, shade = TRUE, legend = TRUE, rot_labels = c(0, 20, 90, 90))

model <- glm(focus.type ~ crypts + cyst + jammed + mucus.simple, data = dataa, family = poisson)

residuals <- chisq.test(contingency_table)$residuals

# Calculate standardized residuals
std_residuals <- chisq.test(contingency_table)$stdres

# Calculate odds ratio
odds_ratio <- oddsratio(contingency_table)

# Print the results
print("Residuals:")
print(residuals)
print("Standardized Residuals:")
print(std_residuals)
print("Odds Ratio:")
print(odds_ratio)
print(contingency_table)
# Print the results
print(fisher_test)
dataa = subset(good.spots.summary.j, (usable == 1))

# these experiments had spytub added like right before imaging instead of >48 hours preinfection
xvar = as.numeric(dataa$t20spy400)/as.numeric(dataa$t0spy0)
fxn_crypt <- dataa %>%
  count(jammed) %>%
  mutate(fraction = n / sum(n))

ggplot(subset(dataa, virus == 'egfp'), aes(x=crypts, y = as.numeric(good.peak.spots))) +
  geom_violin()

contingency_table <- table(dataa$mucus.simple, dataa$jammed)
#contingency_table <- table(dataa$focus.simple, dataa$jammed)
print(contingency_table)
chi_squared_test <- chisq.test(contingency_table)
print(chi_squared_test)
###########################
# how becysted?
