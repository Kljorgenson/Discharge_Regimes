### Explore relationships between discharge regime metrics from anomaly calculation with catchment characteristics
library(tidyverse)
library(tidyhydat)
library(ggfortify)
library(factoextra)
library(ggforce)
library(cowplot)
library(easystats)
library(see)
library(corrplot)
library(ggExtra)
library(grid)
library(gridExtra)
library(ggpubr)
library(ggpattern)
library(vegan)
library(lme4)
library(car)
library(ggpubr)


##### Format all discharge regime metrics, climate variables, and catchment characteristics
## Metadata/catchment characteristics
meta_full <- read.csv("Output_data/Metadata_AK_CAN.csv")

### Read in metric data
metric.data <- read.csv("Output_data/DFFT_metrics_20yrs_ak_can.csv") %>%
  rename(name = site)

head(metric.data)

## Take average of annual regime metrics (HSAM, HSAM timing, flood duration)
# SAM annual -> average SAM over period
SAM <- read.csv("Output_data/SAM_20_ak_can.csv") %>% group_by(site) %>%
  summarise(SAM = mean(resid.sig), SAM_t = mean(timing), peak = unique(peak))
names(SAM) <- c("name", "HSAM", "HSAM_t",'peak')

# Flood event duration and number of events
# Flood is assigned to the year when the flood peaked
events <- read.csv("Output_data/events.csv") %>% select(site, event.duration,extreme.this.event,year,jday) %>%
  group_by(site,year) %>% summarise(event.duration = mean(event.duration)) 
events.d <- events %>% group_by(site) %>% summarise(event.duration.h = mean(event.duration)) %>% rename(name = site)

## Join with metadata
metric.d <- full_join(metric.data, SAM) %>% full_join(events.d) %>% left_join(meta_full) %>%
  filter(name !='SOLOMON') # Remove river with poor DFFT fit during summer

# Save compiled metric data
write.csv(metric.d, 'Output_data/Metrics_recent.csv')



##### Principal Component Aanalysis (PCA) on discharge regime metrics from DFFT
# Select data and run PCA
PCA_dat <- metric.d %>% dplyr::select(name, snr, rms.signal, sigma.hfb, HSAM, HSAM_t, event.duration.h, ecoregion1, glacial, permafrost_class)
PCA_1 <- prcomp(PCA_dat[,c(2:7)],scale. = TRUE) # Run PCA

# View PCA indices
ind <- get_pca_ind(PCA_1)
inds <- as.data.frame(ind$coord)
inds$name <- PCA_dat$name

# View variance explained by PCs
var.exp <- summary(PCA_1)$importance[2,]

# Export variable contribution
var <- get_pca_var(PCA_1)
var$contrib
write.csv(var$contrib, 'Output_data/PCA_contrib.csv', row.names = F)

# Variable coordinates
coords <- var$coord

# Join metric data and PCA data
metric.dat <- full_join(metric.d, inds[,c(1,2,3,7)]) 


#### PCA plots
# Explore autoplots
autoplot(PCA_1, data = PCA_dat, colour = 'glacial', loadings = TRUE, loadings.label = TRUE, loadings.colour = 1, size = 2) +
  theme_classic() + geom_mark_ellipse(aes(color = glacial), size = 0.1) 


### PCA plots with rivers grouped by ecoregions
# Ecoregion colors
cols = c("#4E79A7", "#D95F02", "#66A61E", "#A6761D", "#F5C710")
#colorblindcheck::palette_check(cols, plot = TRUE)

## Plot of PC1 vs PC2
p1 <- metric.dat %>% filter(is.na(ecoregion1) == F) %>% ggplot(aes(Dim.1, Dim.2, colour = ecoregion1)) + theme_bw() +
  stat_chull(geom = "polygon", aes(fill = after_scale(alpha(colour, 0.2))), size = 0, lty = "blank") +
  geom_point(aes(shape = glacial), size = 1.5) +
  scale_color_manual(name = NULL,values = cols) + ylab("PC2 (26%)") + xlab("PC1 (48%)") +
  scale_shape_manual(name = NULL, values = c(17,16), labels = c("Glacial", "Non-glacial")) +
  geom_segment(data = coords,
               aes(x = 0, y = 0, xend = Dim.1*3, yend = Dim.2*3),
               color = 1, arrow = arrow(angle = 25, length = unit(3, "mm"))) +
  geom_text(aes( x=0.4, y=-2.5, label="SNR"), col = 1, size=4) +
  geom_text(aes( x=3.2, y=-1.3, label=paste("A[RMS]")), parse = T, col = 1, size=4) +
  geom_text(aes( x=1.6, y=1.6, label="HSAM"), col = 1, size=4) +
  geom_text(aes( x=1.6, y=1.1, label="timing"), col = 1, size=4) +
  geom_text(aes( x=0.8, y=2.9, label="flood duration"), col = 1, size=4) +
  geom_text(aes( x=3.7, y=-0.4, label="HSAM"), col = 1, size=4) +
  geom_text(aes( x=3.4, y=.5, label=paste('sigma[hf]')), parse = TRUE, col = 1, size=4) +
  scale_fill_manual(name = NULL,values = cols) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                       axis.title = element_text(size = 14), axis.text = element_text(size = 12))

# Extract legend
legend <- cowplot::get_legend(p1)
pl <- grid.draw(legend)

# Add density curves to PC1 vs PC2 plot axis
p2 <- p1 + guides(col = 'none', shape = 'none', fill = 'none')
p3 <- ggMarginal(p2, groupColour = TRUE, groupFill = TRUE, alpha = 0)

ggdraw(plot_grid(plot_grid(p3, ncol=1, align='v'),
                 plot_grid(NULL),
                 plot_grid(legend),
                 plot_grid(NULL),
                 rel_widths=c(1, 0.05,0.5, 0.05), ncol = 4)) + theme(plot.background = element_rect(fill="white", color = NA))

## Plot PC2 vs PC3
p11 <- metric.dat %>% filter(is.na(ecoregion1) == F) %>% ggplot(aes(Dim.3, Dim.2, colour = ecoregion1)) + theme_bw() +
  stat_chull(geom = "polygon", aes(fill = after_scale(alpha(colour, 0.2))), size = 0, lty = "blank") +
  geom_point(aes(shape = glacial), size = 1.5) +
  scale_color_manual(name = NULL,values = cols) + ylab("PC2 (26%)") + xlab("PC3 (16%)") +
  scale_shape_manual(name = NULL, values = c(17,16), labels = c("Glacial", "Non-glacial")) +
  geom_segment(data = coords,
               aes(x = 0, y = 0, xend = Dim.3*3, yend = Dim.2*3),
               color = 1, arrow = arrow(angle = 25, length = unit(3, "mm"))) +
  geom_text(aes( x=-0.7, y=-2.4, label="SNR"), col = 1, size=4) +
  geom_text(aes( x=0.6, y=-1.2, label=paste("A[RMS]")), parse = T, col = 1, size=4) +
  geom_text(aes( x=-3, y=1.8, label="HSAM"), col = 1, size=4) +
  geom_text(aes( x=-3, y=1.3, label="timing"), col = 1, size=4) +
  geom_text(aes( x=0.7, y=2.9, label="flood duration"), col = 1, size=4) +
  geom_text(aes( x=1.1, y=0.2, label="HSAM"), col = 1, size=4) +
  geom_text(aes( x=0.6, y=.7, label=paste('sigma[hf]')), parse = TRUE, col = 1, size=4) + 
  scale_fill_manual(name = NULL,values = cols) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                       axis.title = element_text(size = 14), axis.text = element_text(size = 12))

# Extract legend
legend <- cowplot::get_legend(p11)
pl <- grid.draw(legend)

# Add density curves to PC1 vs PC2 plot axis
p12 <- p11 + guides(col = 'none', shape = 'none', fill = 'none')
p13 <- ggMarginal(p12, groupColour = TRUE, groupFill = TRUE, alpha = 0)

ggdraw(plot_grid(plot_grid(p13, ncol=1, align='v'),
                 plot_grid(NULL),
                 plot_grid(legend),
                 plot_grid(NULL),
                 rel_widths=c(1, 0.05,0.5, 0.05), ncol = 4)) + theme(plot.background = element_rect(fill="white", color = NA))


## Combine PC1 vs PC2 and PC2 vs PC3 plots
# Extract legend
p11.1 <- p11 + theme(legend.position = 'bottom', legend.margin=margin(t=-20), legend.text=element_text(size=13)) + guides(shape = 'none', col = guide_legend(nrow = 2)) # Reformat legend
legend = cowplot::get_plot_component(p11.1, 'guide-box-bottom', return_all = TRUE)

# Alter PC1 vs PC2 plot
p2.1 <- p1 + guides(col = 'none',fill = 'none') + theme(plot.margin = margin(b =0), legend.position = c(.77,.87), legend.box.background = element_rect(colour = "black"),legend.text=element_text(size=12))
p3.2 <- ggMarginal(p2.1, groupColour = TRUE, groupFill = TRUE, alpha = 0, margins = 'x')

# Alter PC2 vs PC3 plot
p12.2 <- p11 + guides(col = 'none', shape = 'none', fill = 'none') + 
  labs(y = NULL) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),plot.margin = margin(b = 0)) 
p13.2 <- ggMarginal(p12.2, groupColour = TRUE, groupFill = TRUE, alpha = 0)

# Figure 5
ggdraw(plot_grid(plot_grid(p3.2, p13.2, nrow=1, rel_widths = c(0.97,1)),
                 plot_grid(legend),
                 rel_heights=c(1, 0.2), nrow = 2)) + theme(plot.background = element_rect(fill="white", color = NA))

ggsave("Figures/Figure_5.png", width = 7.7, height = 5.1)



#### Plots of discharge regime metrics

# Boxplots of metrics by ecoregion 
dat <- metric.dat %>% dplyr::select(ecoregion1, sigma.hfb, snr, HSAM, HSAM_t, rms.signal, event.duration.h) %>% pivot_longer(cols = 2:7, names_to = "variable")
dat$variable = factor(dat$variable, labels = c('Flood~duration~(days)',"HSAM","HSAM~timing~(days)","A[RMS]", 
  "sigma[hf]","SNR")) 

dat %>% ggplot(aes(ecoregion1, value, col = ecoregion1)) + 
  geom_boxplot(outliers = F) + geom_jitter(size = 0.5) + facet_wrap(~variable, scales = "free", labeller = label_parsed) + 
  ylab(NULL) + theme_bw() + xlab(NULL) +
  theme(strip.background =element_rect(fill="white"), axis.ticks.x=element_blank(), text=element_text(size=13),legend.text=element_text(size=8.9),
        axis.text.y = element_text(size = 8),plot.margin = unit(c(0.5,0.7,0.2,0.5), "cm"),
        axis.text.x = element_blank()) + scale_colour_manual(name = NULL,values = cols) +
   theme(strip.background =element_rect(fill="white"), legend.margin = margin(l=0),legend.spacing.y = unit(100.0, 'cm'), legend.key.size = unit(0.5, 'cm'), legend.key.height = unit(1, 'cm'),
         axis.text.x = element_blank())+ scale_color_manual(name = "Ecoregion",values = cols, labels = function(x) str_wrap(x, width = 10))

ggsave("Figures/Figure_6.png", width = 6.8, height = 3.5)

# Boxplots of metrics by glacial vs non-glacial
dat2 <- metric.dat %>% dplyr::select(glacial, ecoregion1, sigma.hfb, snr, HSAM, HSAM_t, rms.signal, event.duration.h) %>% pivot_longer(cols = 3:8, names_to = "variable")
dat2$variable = factor(dat2$variable, labels = c('Flood~duration~(days)',"HSAM","HSAM~timing~(days)","A[RMS]", 
                                               "sigma[hf]","SNR")) 
dat2  %>% 
  ggplot(aes(glacial, value, col = glacial)) + 
  geom_boxplot(outliers = F) + geom_jitter(size = 0.5) + facet_wrap(~variable, scales = "free", labeller = label_parsed) + 
  ylab(NULL) + theme_bw() + xlab(NULL) +
  theme(strip.background =element_rect(fill="white"), legend.margin = margin(l=0),legend.spacing.y = unit(100.0, 'cm'), legend.key.size = unit(0.5, 'cm'), legend.key.height = unit(1, 'cm'),
        axis.text.x = element_blank(), axis.ticks.x=element_blank()) + labs(col = NULL) +
  scale_color_manual(values = c('grey55', 'black'), labels = c( 'Glacial (n = 28)','Non-glacial (n = 144)'))

ggsave("Figures/Figure_S3.png", width = 6.3, height = 3)

# Boxplots of metrics by permafrost categories
cols.p <- rev(c("lightgoldenrod3", "#7FCDBB", "#1D91C0", "#0C2C84"))
dat <- metric.dat %>% dplyr::select(permafrost_class, sigma.hfb, snr, HSAM, HSAM_t, rms.signal, event.duration.h) %>% pivot_longer(cols = 2:7, names_to = "variable")
dat$variable = factor(dat$variable, labels = c('Flood~duration~(days)',"HSAM","HSAM~timing~(days)","A[RMS]", 
                                               "sigma[hf]","SNR")) 

dat %>% ggplot(aes(permafrost_class, value, col = permafrost_class)) + 
  geom_boxplot(outliers = F) + geom_jitter(size = 0.5) + facet_wrap(~variable, scales = "free", labeller = label_parsed) + 
  ylab(NULL) + theme_bw() + xlab(NULL) +
  theme(strip.background =element_rect(fill="white"), axis.ticks.x=element_blank(), text=element_text(size=13),legend.text=element_text(size=8.9),
        axis.text.y = element_text(size = 8),plot.margin = unit(c(0.5,0.7,0.2,0.5), "cm"),
        axis.text.x = element_blank()) + scale_colour_manual(name = NULL,values = cols.p, labels = c("Continuous (n = 16)", "Discontinuous (n = 50)", "Sporadic (n = 73)", "Unfrozen (n = 33)")) +
  theme(strip.background =element_rect(fill="white"), legend.margin = margin(l=0),legend.spacing.y = unit(100.0, 'cm'), legend.key.size = unit(0.5, 'cm'), legend.key.height = unit(1, 'cm'),
        axis.text.x = element_blank())

ggsave("Figures/Figure S4.png", width = 7.2, height = 3.5)




#### Catchment attribute plot

dat <- metric.dat %>% filter(is.na(ecoregion1) == F) %>% dplyr::select(name, ecoregion1, latitude.c, longitude.c, area_km2, rain.20, snowfall.20,temp.20, snowmelt.t.20, SI.20) %>% pivot_longer(cols = 3:10, names_to = "variable")


themes <- theme_bw()  +
  theme(legend.margin = margin(l=0),legend.spacing.y = unit(100.0, 'cm'), legend.key.size = unit(0.5, 'cm'), legend.key.height = unit(1, 'cm'),
        axis.text.x = element_blank(), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10)) 

a <- dat %>% filter(variable == 'area_km2') %>% ggplot(aes(ecoregion1, value, col = ecoregion1)) + 
  geom_boxplot(outliers = F) + geom_jitter(size = 0.5) + 
  ylab(expression(Area~(km^{2}))) + xlab(NULL)+
  scale_color_manual(name = NULL,values = cols, labels = function(x) str_wrap(x, width = 20)) +
  scale_y_continuous(trans = "log10", labels = scales::label_log()) + themes +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.15), "cm"))

b<- dat %>% filter(variable == 'latitude.c') %>% ggplot(aes(ecoregion1, value, col = ecoregion1)) + 
  geom_boxplot(outliers = F) + geom_jitter(size = 0.5) + 
  ylab(expression(Latitude~('°'))) + xlab(NULL) + themes +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.3), "cm"))+ 
  scale_color_manual(name = "Ecoregion",values = cols, labels = function(x) str_wrap(x, width = 10))

c<- dat %>% filter(variable == 'longitude.c') %>% ggplot(aes(ecoregion1, value, col = ecoregion1)) + 
  geom_boxplot(outliers = F) + geom_jitter(size = 0.5) + 
  ylab(expression(Longitude~('°'))) + xlab(NULL) + themes +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.15), "cm")) +
  scale_color_manual(name = "Ecoregion",values = cols, labels = function(x) str_wrap(x, width = 10))

d<- dat %>% filter(variable == 'temp.20') %>% ggplot(aes(ecoregion1, value, col = ecoregion1)) + 
  geom_boxplot(outliers = F) + geom_jitter(size = 0.5) + 
  ylab(expression(Temperature~('°C'))) + xlab(NULL) + themes +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.15), "cm")) +
  scale_color_manual(name = "Ecoregion",values = cols, labels = function(x) str_wrap(x, width = 10))

e<- dat %>% filter(variable == 'rain.20') %>% ggplot(aes(ecoregion1, value, col = ecoregion1)) + 
  geom_boxplot(outliers = F) + geom_jitter(size = 0.5) + 
  ylab('Rain (mm)') + xlab(NULL) +themes +
  scale_color_manual(name = "Ecoregion",values = cols, labels = function(x) str_wrap(x, width = 10))

f<- dat %>% filter(variable == 'snowfall.20') %>% ggplot(aes(ecoregion1, value, col = ecoregion1)) + 
  geom_boxplot(outliers = F) + geom_jitter(size = 0.5) + 
  ylab('Snowfall (mm)') + xlab(NULL) + themes +
  scale_color_manual(name = "Ecoregion",values = cols, labels = function(x) str_wrap(x, width = 10)) +
  ylim(0,2621)

g<- dat %>% filter(variable == 'snowmelt.t.20') %>% ggplot(aes(ecoregion1, value, col = ecoregion1)) + 
  geom_boxplot(outliers = F) + geom_jitter(size = 0.5) + 
  ylab("Snowmelt timing \n(day of year)") + xlab(NULL) +themes +
  scale_color_manual(name = "Ecoregion",values = cols, labels = function(x) str_wrap(x, width = 10))

h <- dat %>% filter(variable == 'SI.20') %>% ggplot(aes(ecoregion1, value, col = ecoregion1)) + 
  geom_boxplot(outliers = F) + geom_jitter(size = 0.5) + 
  ylab("Rainfall Seasonality \n Index (Jun-Oct)") + xlab(NULL) +themes +
  scale_color_manual(name = "Ecoregion",values = cols, labels = function(x) str_wrap(x, width = 10))


ggarrange(a,b,c,d,e,f,g,h, common.legend = T, widths = c(1,1,1.05,1), legend = 'top', ncol = 4, nrow = 2) + bgcolor("white") + border('white')   
ggsave("Figures/Figure_6.png", width = 7, height = 3)



#### Plot of % of rivers in glacial and permafrost categories by ecoregions
metric.dat %>% group_by(ecoregion1) %>% summarise(n = length(unique(name)))
cols.p <- rev(c("lightgoldenrod3", "#7FCDBB", "#1D91C0", "#0C2C84"))
metric.dat %>% ggplot() + geom_bar(aes(y = ecoregion1, fill = permafrost_class), position="fill") + theme_bw() + ylab("") + xlab("") +
  scale_fill_manual(name = "Permafrost extent",values = cols.p, labels = c("Continuous", "Discontinuous", "Sporadic", "Unfrozen")) +
  geom_text(aes(x=.2,y=5, label = "10"), col = 'white', size = 3) +
  geom_text(aes(x=.2,y=4, label = "37"), size = 3)+
  geom_text(aes(x=.2,y=3, label = "69"), size = 3)+
  geom_text(aes(x=.2,y=2, label = "34"), size = 3)+
  geom_text(aes(x=.2,y=1, label = "22"), size = 3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_discrete(labels = ~str_wrap(.x, 15)) +
  scale_x_continuous(name = "",
                     breaks = c(0, 0.5, 1), 
                     labels = scales::percent(c(0, 0.5, 1))) +
  geom_bar(aes(y = ecoregion1, col = glacial), position="fill", fill = NA, size = 1) +
  scale_color_manual(name = NULL, labels = c('Glacial', 'Non-glacial'), values = c(1, 'grey55'))
  #geom_bar_pattern(aes(pattern = glacial), position="fill", pattern_spacing = 0.02, width = 0.5) + scale_pattern_manual(values = c('stripe', 'crosshatch'))

ggsave("Figures/Figure_2.png", width = 4.5, height = 2.5)


### Correlation plot of metrics and catchment characteristics
cordat <- metric.dat %>% dplyr::select(latitude.c, longitude.c, area_km2_log, temp.20, snowfall.20, rain.20, snowmelt.t.20, SI.20, snr, rms.signal, sigma.hfb, HSAM, HSAM_t, event.duration.h)
names(cordat) <- c("latitude", "longitude", "area", "temperature", "snowfall", "rainfall", 'snowmelt timing', 'rainfall seasonality index', "SNR", "$A[RMS]", "$sigma[hf]", "HSAM", "HSAM timing", "flood duration")
cor1 <-  cordat %>% 
  correlation()

png(filename="Figures/Figure_S2.png", width = 4000, height = 4000, res = 400)
corrplot(cor(cor1), type = 'lower', tl.col = 'black', tl.srt = 45, method = 'color', tl.cex = 1,
         addCoef.col = 'black', cl.pos = 'n', col = c("#C6413E", "#D96752", "#E88B6E", "#F5AD8C", "#FAC9B1", "#FDE2D2", "#FEF5F0", "#F2F8FB", "#DAEAF3", "#BDDAEA", "#9BCAE0", "#74B2D4", "#4B98C5",                                                                                                                                  "#3480B9"))
dev.off()






########################## Analysis
### PERMANOVA on PCA axes
# Select data
dat <- metric.dat %>% dplyr::select(permafrost_class, glacial, ecoregion1, Dim.1, Dim.2, Dim.3) %>% na.omit()

## PERMANOVA including all levels
adonis2(dat[,4:6]~dat[,1], method = "euc") # Permafrost class
adonis2(dat[,4:6]~dat[,2], method = "euc") # Glacial versus non-glacial
adonis2(dat[,4:6]~dat[,3], method = "euc") # Ecoregion


## Pairwise comparison of ecoregions
dat2 <- metric.dat %>% dplyr::select(ecoregion1, Dim.1, Dim.2) %>% na.omit()

# List of all pairwise comparisons
vars <- expand.grid(unique(dat2$ecoregion1), unique(dat2$ecoregion1)) %>% filter(Var1 != Var2)
vars <- vars %>% filter(duplicated(t(apply(vars, 1, sort))) == F)

# Loop through pairwise ecoregion comparisons
p_values <- c() 
psuedo.F <- c()
DFr <- c()
DFt <- c()

for(i in 1:nrow(vars)){
  var <- vars[i,]
  dat2 <- metric.dat %>% dplyr::select(ecoregion1, Dim.1, Dim.2, Dim.3) %>% na.omit() %>% filter(ecoregion1 %in% c(var$Var1, var$Var2))
  pairwise_anosim <- adonis2(dat2[,2:3]~dat2$ecoregion1, method = "euc")
  # Save the p-value
  p_values <- c(p_values, pairwise_anosim$`Pr(>F)`[1])
  psuedo.F <- c(psuedo.F, pairwise_anosim$`F`[1])
  DFr <- c(DFr, pairwise_anosim$Df[2])
  DFt <- c(DFt, pairwise_anosim$Df[3])
  }

# Adjust for multiple comparisons using the Holms method
adjusted_p_values <- p.adjust(p_values, method = "holm")

anosim.pairwise.ecoregion1 <- cbind( vars, DFr, DFt, psuedo.F, adjusted_p_values, p_values)
anosim.pairwise.ecoregion1


## Export all PERMANOVA results
glac <- data.frame(Var1 = 'Glacial', Var2 = 'Non-glacial', p_values = 0.001, psuedo.F = 6.2683, DFr=164, DFt = 168)
manova <- full_join(anosim.pairwise.ecoregion1, glac)
write.csv(manova, "Output_data/Table_S2.csv", row.names = F)




### Linear regression to explore influence of climate variables on discharge regime metrics
## Linear model of PC1
mix1 <- lm(data = metric.dat, Dim.1 ~ glacial + scale(rain.20)*glacial + scale(snowfall.20)*glacial + scale(temp.20)*glacial + scale(area_km2_log)*glacial)
mix_sum1 <- summary(mix1)

# Check variance inflation
vif(mix1, type = "predictor")

## Extract effect size and confidence intervals for plotting
coefs <- mix_sum1$coefficients[,1]
coefs

X <- model.matrix(mix1)
dof <- nrow(X)-ncol(X)
coefs_var <- vcov(mix1)

cnfnt <- confint(mix1)

# Confidence intervals for non-glacial
# rain
halfCI <- qt(0.975, dof) * sqrt(coefs_var[7,7]+coefs_var[3,3]+2*coefs_var[7,3])
rain <- as.vector(c(coefs[7]+coefs[3]-halfCI, coefs[7]+coefs[3]+halfCI))

# snowfall
halfCI <- qt(0.975, dof) * sqrt(coefs_var[8,8]+coefs_var[4,4]+2*coefs_var[8,4])
snow <- as.vector(c(coefs[8]+coefs[4]-halfCI, coefs[8]+coefs[4]+halfCI))

# temp
halfCI <- qt(0.975, dof) * sqrt(coefs_var[9,9]+coefs_var[5,5]+2*coefs_var[9,5])
temp <- as.vector(c(coefs[9]+coefs[5]-halfCI, coefs[9]+coefs[5]+halfCI))

# area
halfCI <- qt(0.975, dof) * sqrt(coefs_var[10,10]+coefs_var[6,6]+2*coefs_var[10,6])
area <- as.vector(c(coefs[10]+coefs[6]-halfCI, coefs[10]+coefs[6]+halfCI))

cofs1 <- data.frame(glacial = c(rep("glacial", 4), rep("non-glacial", 4)), 
                   term = c(rep(c("Rainfall", "Snowfall", "Temperature", "Area"), 2)),
                   estimate = c(coefs[3:6], c(coefs[3:6] + coefs[7:10])),
                   upperci = c(cnfnt[3:6,2], rain[2], snow[2], temp[2], area[2]),
                   lowerci = c(cnfnt[3:6,1], rain[1], snow[1], temp[1], area[1]))

cofs1
cofs1$PC <- 'PC1'

## Linear model for PC2
mix1 <- lm(data = metric.dat, Dim.2 ~ glacial + scale(rain.20)*glacial + scale(snowfall.20)*glacial + scale(temp.20)*glacial + scale(area_km2_log)*glacial)
mix_sum2 <- summary(mix1)

## Extract effect size and confidence intervals for plotting
coefs <- mix_sum2$coefficients[,1]
coefs

X <- model.matrix(mix1)
dof <- nrow(X)-ncol(X)
coefs_var <- vcov(mix1)

cnfnt <- confint(mix1)

# Confidence intervals for non-glacial
# rain
halfCI <- qt(0.975, dof) * sqrt(coefs_var[7,7]+coefs_var[3,3]+2*coefs_var[7,3])
rain <- as.vector(c(coefs[7]+coefs[3]-halfCI, coefs[7]+coefs[3]+halfCI))

# snowfall
halfCI <- qt(0.975, dof) * sqrt(coefs_var[8,8]+coefs_var[4,4]+2*coefs_var[8,4])
snow <- as.vector(c(coefs[8]+coefs[4]-halfCI, coefs[8]+coefs[4]+halfCI))

# temp
halfCI <- qt(0.975, dof) * sqrt(coefs_var[9,9]+coefs_var[5,5]+2*coefs_var[9,5])
temp <- as.vector(c(coefs[9]+coefs[5]-halfCI, coefs[9]+coefs[5]+halfCI))

# area
halfCI <- qt(0.975, dof) * sqrt(coefs_var[10,10]+coefs_var[6,6]+2*coefs_var[10,6])
area <- as.vector(c(coefs[10]+coefs[6]-halfCI, coefs[10]+coefs[6]+halfCI))

cofs2 <- data.frame(glacial = c(rep("glacial", 4), rep("non-glacial", 4)), 
                    term = c(rep(c("Rainfall", "Snowfall", "Temperature", "Area"), 2)),
                    estimate = c(coefs[3:6], c(coefs[3:6] + coefs[7:10])),
                    upperci = c(cnfnt[3:6,2], rain[2], snow[2], temp[2], area[2]),
                    lowerci = c(cnfnt[3:6,1], rain[1], snow[1], temp[1], area[1]))

cofs2
cofs2$PC <- 'PC2'


## Linear model for HSAM timing
mix1 <- lm(data = metric.dat, scale(HSAM_t) ~ glacial + scale(rain.20)*glacial + scale(snowfall.20)*glacial + scale(temp.20)*glacial + scale(area_km2_log)*glacial)
mix_sum3 <- summary(mix1)

## Emix1## Extract effect size and confidence intervals for plotting
coefs <- mix_sum3$coefficients[,1]
coefs

X <- model.matrix(mix1)
dof <- nrow(X)-ncol(X)
coefs_var <- vcov(mix1)

cnfnt <- confint(mix1)

# Confidence intervals for non-glacial
# rain
halfCI <- qt(0.975, dof) * sqrt(coefs_var[7,7]+coefs_var[3,3]+2*coefs_var[7,3])
rain <- as.vector(c(coefs[7]+coefs[3]-halfCI, coefs[7]+coefs[3]+halfCI))

# snowfall
halfCI <- qt(0.975, dof) * sqrt(coefs_var[8,8]+coefs_var[4,4]+2*coefs_var[8,4])
snow <- as.vector(c(coefs[8]+coefs[4]-halfCI, coefs[8]+coefs[4]+halfCI))

# temp
halfCI <- qt(0.975, dof) * sqrt(coefs_var[9,9]+coefs_var[5,5]+2*coefs_var[9,5])
temp <- as.vector(c(coefs[9]+coefs[5]-halfCI, coefs[9]+coefs[5]+halfCI))

# area
halfCI <- qt(0.975, dof) * sqrt(coefs_var[10,10]+coefs_var[6,6]+2*coefs_var[10,6])
area <- as.vector(c(coefs[10]+coefs[6]-halfCI, coefs[10]+coefs[6]+halfCI))

cofs3 <- data.frame(glacial = c(rep("glacial", 4), rep("non-glacial", 4)), 
                    term = c(rep(c("Rainfall", "Snowfall", "Temperature", "Area"), 2)),
                    estimate = c(coefs[3:6], c(coefs[3:6] + coefs[7:10])),
                    upperci = c(cnfnt[3:6,2], rain[2], snow[2], temp[2], area[2]),
                    lowerci = c(cnfnt[3:6,1], rain[1], snow[1], temp[1], area[1]))

cofs3
cofs3$PC <- 'HSAM_t'



## Combine results of 3 linear models
cofs <- rbind(cofs1,cofs2,cofs3)

sum1 <- as.data.frame(mix_sum1$coefficients)
sum1$response <- 'PC1'
sum2 <- as.data.frame(mix_sum2$coefficients)
sum2$response <- 'PC2'
sumt <- as.data.frame(mix_sum3$coefficients)
sumt$response <- 'HSAM.t'

sum <- rbind(sum1,sum2,sumt)
sum$term <- row.names(sum)
row.names(sum) <- NULL
names(sum) <- c('estimate','sd','t','p','response','term') 
sum <- sum %>%
  mutate(glacial = ifelse(grepl('glacial', term), 'non-glacial', 'glacial'))
head(sum)

cofs$PC <- factor(cofs$PC, levels = c("PC1","PC2", "HSAM_t"), labels = c("PC1","PC2", "HSAM~timing"))


### Plot effects of climate on regime metrics
## Plot effect sizes
p.efs <- cofs %>% ggplot(aes(estimate,term, shape = glacial, col = glacial))+ 
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_errorbar(aes(xmin = lowerci, xmax = upperci, y = term, width = 0.7), position = position_dodge(width = 0.8))+ 
  geom_point(size = 2.5, position = position_dodge(width = 0.8))  + theme_bw() +
  ylab("")  + theme(strip.background =element_rect(fill="white"), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), legend.margin=margin(t = 0, b = -.3, unit='cm'), 
                    plot.margin = unit(c(60,5.5,5.5,0), "pt"), text=element_text(size=15), legend.position = c(0.85,1.7),
                    axis.text.y=element_text(angle=-45, hjust = 1, vjust = 1)) + facet_wrap(~PC, labeller = label_parsed) +
  scale_color_manual(values = c(1, 'grey55'), labels = c("Glacial", "Non-glacial")) +
  scale_shape_manual(name = NULL, values = c(17,16), labels = c("Glacial", "Non-glacial")) +
  labs(col = NULL, shape = NULL) + xlim(-1.8,1.6) + xlab("Effect size") +scale_y_discrete(limits=rev)
p.efs

## PCA plots with rivers colored by scaled climate variable values
# Format data
dat <- metric.dat %>% filter(is.na(ecoregion1) == F) %>% mutate(Rain = scale(log(rain.20)), Snowfall = scale(log(snowfall.20)),
                                                                Temperature = scale(temp.20), Area = scale(area_km2_log)) %>%
  dplyr::select(name, ecoregion1, Dim.1, Dim.2, Area, Rain, Snowfall, Temperature, glacial) %>%
  pivot_longer(cols = 5:8, names_to = "variable") 
dat$variable = factor(dat$variable, labels = c(
  "Area~(km^{2})", 
  "Rainfall~(mm)", "Snowfall~(mm)", "Temperature~('°C')")) 

# Plot
plt <- dat %>% ggplot(aes(Dim.1, Dim.2, col = value, shape = glacial)) + geom_point(size = 1.2) + facet_wrap(~variable, labeller = label_parsed) +
  theme_bw() + labs(col = NULL) + scale_shape_manual(name = NULL, values = c(17,16), labels = c("Glacial", "Non-glacial")) +
  theme(strip.background =element_rect(fill="white"), legend.direction="horizontal",text=element_text(size=15),
        legend.position = c(.8,0.9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5,0.2,0,1.1), "cm")) + scale_color_gradientn(
          colors = c("darkgoldenrod","goldenrod1", "blue", "#00008B"),
          values = scales::rescale(c(-3.205624,-1,1, 4.120276), from = c(-3.205624, 4.120276)),
          limits = c(-3.205624, 4.120276)) + xlab("PC1") + ylab("PC2") +
  guides(shape = 'none')
plt

# Combine plots
fig <- ggarrange(plt,p.efs, nrow = 2, heights = c(1.2,0.8))
fig2 <- annotate_figure(fig, fig.lab = "\n      a)", fig.lab.pos = 'top.left', fig.lab.face = 'bold')
annotate_figure(fig2, fig.lab = "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\     b)", fig.lab.pos = 'top.left', fig.lab.face = 'bold')

ggsave("Figures/Figure_7.png", width = 4.8, height = 7.9)
