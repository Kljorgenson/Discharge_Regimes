### Explore relationships between flow regime metrics from anomaly claculation with catchment characteristics

library(tidyverse)
library(tidyhydat)
library(ggfortify)
library(factoextra)
library(ggforce)
library(cowplot)
library(ggExtra)
library(grid)
library(gridExtra)
library(vegan)
library(RVAideMemoire)
library(emmeans)
library(ggpubr)
library(stats)

### Read in metric data for recent period
metric.data.20 <- read.csv("Output_data/Metrics_recent.csv") %>% select(name,seasonal,rms.signal,rms.noise,snr,sigma.hfb,start.date,end.date,HSAM,HSAM_t,peak,event.duration.h)
metric.data.20$period <- "recent"
head(metric.data.20)

### Compile metric data from historical period (1970-1992)
metric.data.e <- read.csv("Output_data/DFFT_metrics_1970-1992_ak_can.csv")%>% rename(name = site)

# SAM annual -> average SAM over period
SAM.e <- read.csv("Output_data/SAM_1970-1992_ak_can.csv") %>% group_by(site) %>%
  summarise(SAM = mean(resid.sig), SAM_t = mean(timing), peak = unique(peak))
names(SAM.e) <- c("name", "HSAM", "HSAM_t",'peak')



# Event duration and number of events
events <- read.csv("Output_data/events_1970-1992.csv") %>% select(site, event.duration,extreme.this.event,year,jday) %>%
  group_by(site,year) %>% summarise(event.duration = mean(event.duration)) 
events.e <- events %>% group_by(site) %>% summarise(event.duration.h = mean(event.duration)) %>% rename(name = site)


# Combine all historical period metrics
metric.e <- full_join(metric.data.e, SAM.e) %>% full_join(events.e)
metric.e$period <- "historical"

## Join periods and add metadata
meta_full <- read.csv("Output_data/Metadata_AK_CAN.csv")
metric.d.20 <- metric.data.20 %>% filter(name %in% metric.e$name) %>%
  full_join(metric.e) %>% left_join(meta_full)



##### Principal Component Aanalysis (PCA) on discharge regime metrics separated by period
# Select data
PCA_dat.20 <- metric.d.20 %>% dplyr::select(name, snr, rms.signal, sigma.hfb, HSAM, HSAM_t, event.duration.h, ecoregion1, glacial, permafrost_class, period) %>% na.omit() %>% ungroup()

# Run PCA
PCA_1.20 <- prcomp(PCA_dat.20[,c(2:7)],scale. = TRUE)

# View PCA indices
ind.20 <- get_pca_ind(PCA_1.20)
ind.20
inds.20 <- as.data.frame(ind.20$coord)
inds.20$name <- PCA_dat.20$name
inds.20$period <- PCA_dat.20$period

# View variable contribution
var.20 <- get_pca_var(PCA_1.20)
var.20$contrib

# Variable coordinates
coords.20 <- as.data.frame(var.20$coord)
coords.20$term <- row.names(coords.20)

# View variance explained by PCs
var.exp.20 <- summary(PCA_1.20)$importance[2,]

# Join metric data, metadata, and PCA data
metric.dat.20 <- full_join(PCA_dat.20, inds.20[,c(1,2,3,7,8)])



#### PCA plots
# Explore autoplots
autoplot(PCA_1.20, data = PCA_dat.20, colour = 'period', loadings = TRUE, loadings.label = TRUE, loadings.colour = 1, size = 2) +
  theme_classic() + geom_mark_ellipse(aes(color = period), size = 0.1)
autoplot(PCA_1.20, data = PCA_dat.20, x = 3, y = 2, colour = 'period', loadings = TRUE, loadings.label = TRUE, loadings.colour = 1, size = 2) +
  theme_classic() + geom_mark_ellipse(aes(color = period), size = 0.1)

### PCA plots with rivers grouped by ecoregions and period
# Ecoregion colors
cols = c("#4E79A7", "#D95F02", "#66A61E", "#A6761D", "#F5C710")
#colorblindcheck::palette_check(cols, plot = TRUE)

# Separate historical and recent periods
old <- metric.dat.20 %>% filter(period == 'historical', ecoregion1 != "Tundra") 
new <- metric.dat.20 %>% filter(period == 'recent') %>% filter(name %in% old$name)

## Plot of PC1 vs PC2
p1 <- metric.dat.20 %>% filter(name %in% old$name | name %in% c("COPPERMINE RIVER AT OUTLET OF POINT LAKE", "KUPARUK")) %>% mutate(shape = paste(glacial,period,sep="")) %>%
  ggplot(aes(Dim.1, Dim.2, colour = ecoregion1)) + theme_bw() +
  stat_chull(data = new, geom = "polygon", aes(fill = after_scale(alpha(colour, 0.2))), size = 0, lty = "blank") + 
  stat_chull(data = old, geom = "polygon", aes(fill = after_scale(alpha(colour, 0))), size = 0.5, lty = 2) +
  geom_point(aes(fill = ecoregion1, shape = shape), size = 2) +
  scale_color_manual(name = NULL,values = cols) + ylab("PC2 (27%)") + xlab("PC1 (45%)") +
  scale_fill_manual(name = NULL,values = cols) + theme_bw()+
  scale_shape_manual(name = NULL, values = c(2,17,1,16), labels = c("1970-1992, glacial", "2000-2022, glacial", "1970-1992, non-glacial", "2000-2022, non-glacial")) +
   geom_segment(data = coords.20,
               aes(x = 0, y = 0, xend = Dim.1*2.5, yend = Dim.2*2.5),
               color = 1, arrow = arrow(angle = 25, length = unit(4, "mm"))) +
  geom_text(aes( x=0.5, y=-1.8, label="SNR"), col = 1, size=4) +
  geom_text(aes( x=2.7, y=-0.8, label = paste("A[RMS]")), parse = T, col = 1, size=4) +
  geom_text(aes( x=1.2, y=-0.6, label="HSAM"), col = 1, size=4) +
  geom_text(aes( x=1.3, y=-0.9, label="timing"), col = 1, size=4) +
  geom_text(aes( x=-0, y=2.7, label="flood duration"), col = 1, size=4) +
  geom_text(aes( x=3.2, y=0, label="HSAM"), col = 1, size=4) +
  geom_text(aes( x=2.6, y=.9, label=paste('sigma[hf]')), parse = T, col = 1, size=4) +
 #guides(color = guide_legend(override.aes = list(size = 0.5)),
  #      fill = guide_legend(override.aes = list(size = 0.5))) +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       axis.title = element_text(size = 14),
       axis.text = element_text(size = 12))

p1

## Plot of PC2 vs PC3
p2 <- metric.dat.20 %>% filter(name %in% old$name | name %in% c("COPPERMINE RIVER AT OUTLET OF POINT LAKE", "KUPARUK")) %>% mutate(shape = paste(glacial,period,sep="")) %>%
  ggplot(aes(Dim.3, Dim.2, colour = ecoregion1)) +
  stat_chull(data = new, geom = "polygon", aes(fill = after_scale(alpha(colour, 0.2))), size = 0, lty = "blank") + theme_bw() +
  stat_chull(data = old, geom = "polygon", aes(fill = after_scale(alpha(colour, 0))), size = 0.5, lty = 2) +
  geom_point(aes(fill = ecoregion1, shape = shape), size = 2) +
  scale_color_manual(name = NULL,values = cols) + ylab("PC2 (27%)") + xlab("PC3 (18%)") +
  scale_fill_manual(name = NULL,values = cols) + theme_bw()+
  scale_shape_manual(name = NULL, values = c(2,17,1,16), labels = c("1970-1992, glacial", "2000-2022, glacial", "1970-1992, non-glacial", "2000-2022, non-glacial")) +
  geom_segment(data = coords.20,
               aes(x = 0, y = 0, xend = Dim.3*2.5, yend = Dim.2*2.5),
               color = 1, arrow = arrow(angle = 25, length = unit(4, "mm"))) +
  geom_text(aes( x=0, y=-2, label="SNR"), col = 1, size=4) +
  geom_text(aes( x=1, y=-0.4, label = paste("A[RMS]")), parse = T, col = 1, size=4) +
  geom_text(aes( x=-3, y=-0.2, label="HSAM"), col = 1, size=4) +
  geom_text(aes( x=-3, y=-0.5, label="timing"), col = 1, size=4) +
  geom_text(aes( x=-0.5, y=2.4, label="flood duration"), col = 1, size=4) +
  geom_text(aes( x=0.9, y=0.2, label="HSAM"), col = 1, size=4) +
  geom_text(aes( x=0.5, y=0.8, label=paste('sigma[hf]')), parse = T, col = 1, size=4) +
  #guides(color = guide_legend(override.aes = list(size = 0.5)),
  #       fill = guide_legend(override.aes = list(size = 0.5))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
p2

# Combine plots
p2.2 <- p2 + theme(legend.position = 'bottom', axis.text.y=element_blank(), axis.ticks.y=element_blank()) + ylab("") + guides(shape = guide_legend(nrow=1), col = 'none',fill = 'none') 
p1.1 <- p1 + theme(legend.position = 'bottom', legend.spacing.x = unit(-0.2, 'cm')) + guides(shape = 'none')

# Extract the shape and color legends individually
legend_shape <- get_legend(p2.2+ theme(legend.text=element_text(size=12, margin = margin(l = 0, r =-4, t = 1, unit = "pt")), legend.key.size = unit(0.8,'cm'), legend.margin = margin(t = -.7, unit="cm")))
legend_color <- get_legend(p1.1+ theme(legend.text=element_text(size=12, margin = margin(l = 1, r = -1, unit = "pt")), legend.key.size = unit(0.4,'cm'), legend.margin = margin(b = -.3, t = -0.5, unit="cm")))

# Combine the legends vertically
combined_legends <- plot_grid(legend_color, legend_shape, ncol = 1)

# Combine the plots (without legends)
plots <- plot_grid(
  p1.1 + theme(legend.position = "none"),
  p2.2 + theme(legend.position = "none"),
  ncol = 2
)

# Add the stacked legends at the bottom
plot_grid(
  plots,
  combined_legends,
  ncol = 1,
  rel_heights = c(1, 0.2)
) + theme(plot.background = element_rect(fill="white", color = NA))

ggsave("Figures/Figure_8.png", width = 7.7, height = 4.7)






########################## Analysis
### PERMANOVA for differences between periods
# Select data
dat <- metric.dat.20 %>% dplyr::select(name, permafrost_class, glacial, ecoregion1, Dim.1, Dim.2, Dim.3, period) %>% na.omit()

## PERMANOVA for all rivers
ad1 <- adonis2(dat[,5:7]~dat$period, method = "euc", strata = dat$name)
ad1

## PERMANOVA by group
# Glacial and non-glacial
dat2 <- metric.dat.20 %>% dplyr::select(ecoregion1, glacial, permafrost_class, period, name, Dim.1, Dim.2, Dim.3) %>% na.omit() %>% filter(glacial == "glacial")
ad2 <- adonis2(dat2[,6:8]~dat2$period, method = "euc", strata = dat2$name)
dat3 <- metric.dat.20 %>% dplyr::select(ecoregion1, glacial, permafrost_class, period, name, Dim.1, Dim.2, Dim.3) %>% na.omit() %>% filter(glacial == "non-glacial")
ad3 <- adonis2(dat3[,6:8]~dat3$period, method = "euc", strata = dat3$name)

p.values <- c(ad2$`Pr(>F)`[1],ad3$`Pr(>F)`[1])
p.g <- data.frame(p = p.values, group = c('glacial', 'non-glacial'))
p.g$p.adj <- p.adjust(p.g$p, method = "holm") # Adjust p-values for multiple tests using the Holms method
p.g

# Ecoregions
dat4 <- metric.dat.20 %>% dplyr::select(ecoregion1, glacial, permafrost_class, period, name, Dim.1, Dim.2, Dim.3) %>% na.omit() %>% filter(ecoregion1 == "Northern Forests")
ad4 <- adonis2(dat4[,6:8]~dat4$period, method = "euc", strata = dat4$name)
dat5 <- metric.dat.20 %>% dplyr::select(ecoregion1, glacial, permafrost_class, period, name, Dim.1, Dim.2, Dim.3) %>% na.omit() %>% filter(ecoregion1 == "Northwestern Forested Mountains")
ad5 <- adonis2(dat5[,6:8]~dat5$period, method = "euc", strata = dat5$name)
dat7 <- metric.dat.20 %>% dplyr::select(ecoregion1, glacial, permafrost_class, period, name, Dim.1, Dim.2, Dim.3) %>% na.omit() %>% filter(ecoregion1 == "Taiga")
ad7 <- adonis2(dat7[,6:8]~dat7$period, method = "euc", strata = dat7$name)

p.values <- c(ad4$`Pr(>F)`[1],ad5$`Pr(>F)`[1],ad7$`Pr(>F)`[1])
p.e <- data.frame(p = p.values, group = c('Northern Forests', 'Northwestern Forested Mountains', 'Taiga'))
p.e$p.adj <- p.adjust(p.e$p, method = "holm")
p.e

# Permafrost categories
dat8 <- metric.dat.20 %>% dplyr::select(ecoregion1, glacial, permafrost_class, period, name, Dim.1, Dim.2, Dim.3) %>% na.omit() %>% filter(permafrost_class == "C")
ad8 <- adonis2(dat8[,6:8]~dat8$period, method = "euc", strata = dat8$name)
dat9 <- metric.dat.20 %>% dplyr::select(ecoregion1, glacial, permafrost_class, period, name, Dim.1, Dim.2, Dim.3) %>% na.omit() %>% filter(permafrost_class %in% c("D", "C"))
ad9 <- adonis2(dat9[,6:8]~dat9$period, method = "euc", strata = dat9$name)
dat10 <- metric.dat.20 %>% dplyr::select(ecoregion1, glacial, permafrost_class, period, name, Dim.1, Dim.2, Dim.3) %>% na.omit() %>% filter(permafrost_class == "S")
ad10 <- adonis2(dat10[,6:8]~dat10$period, method = "euc", strata = dat10$name)
dat11 <- metric.dat.20 %>% dplyr::select(ecoregion1, glacial, permafrost_class, period, name, Dim.1, Dim.2, Dim.3) %>% na.omit() %>% filter(permafrost_class == "Unfrozen")
ad11 <- adonis2(dat11[,6:8]~dat11$period, method = "euc", strata = dat11$name)

p.values <- c(ad9$`Pr(>F)`[1],ad10$`Pr(>F)`[1],ad11$`Pr(>F)`[1])
p.p <- data.frame(p = p.values, group = c('D', 'S', 'Unfrozen'))
p.p$p.adj <- p.adjust(p.p$p, method = "holm")
p.p


### Stat output dataframe
tot <- data.frame(p = ad1$`Pr(>F)`[1], group = 'All')
stats <- full_join(tot, p.e) %>% full_join(p.p) %>% full_join(p.g)
stats$F <- c(ad1$F[1], ad4$F[1],ad5$F[1],ad6$F[1],ad7$F[1],ad8$F[1],ad9$F[1],ad10$F[1],ad11$F[1], ad2$F[1],ad3$F[1])
stats$DFr <- c(ad1$Df[2], ad4$Df[2],ad5$Df[2],ad6$Df[2],ad7$Df[2],ad8$Df[2],ad9$Df[2],ad10$Df[2],ad11$Df[2], ad2$Df[2],ad3$Df[2])
stats$DFt <- c(ad1$Df[3], ad4$Df[3],ad5$Df[3],ad6$Df[3],ad7$Df[3],ad8$Df[3],ad9$Df[3],ad10$Df[3],ad11$Df[3], ad2$Df[3],ad3$Df[3])
write.csv(stats, "C:/Users/kljorgenson/Documents/Repos/AK_discharge/Data/DFFT metrics/PERMANOVA_periods.csv", row.names= F)





#### Test differences in climate variables between periods
# Reformat data
dat <- meta_full %>% select(station, name, ppt.20, ppt.60, rain.20, rain.60,snowfall.20,snowfall.60,temp.20,temp.60,SI.20,SI.60,snowmelt.t.20,snowmelt.t.60,ecoregion1) %>% pivot_longer(cols = 3:14, names_to = "variable") %>% 
  separate(col = variable, into = c('variable', 'period', 'period2')) %>% 
  mutate(period = ifelse(is.na(period2),period,period2)) %>%
    mutate(period = case_when(period == 20 ~ '2000-2023', period == 60 ~ '1970-1992')) %>%
  filter(is.na(value) == F) %>% group_by(name,variable) %>% filter(length(name) == 2) %>% 
  mutate(variable = case_when(variable == 'rain' ~ "Rain", variable == 'snowfall' ~ 'Snowfall',
                             variable == 'temp' ~ "Temperature", variable == 'ppt' ~ 'Precipitation',
                             variable == 'snowmelt' ~ 'Snowmelt timing', variable == 'SI' ~ 'Rainfall Seasonality Index')) %>% select(-c(period2))
dat

#3 PERMANOVA by ecoregion
# Rain
d1 <- dat %>% filter(variable == 'Rain', ecoregion1 == "Northern Forests")
ad1 <- adonis2(d1[,6]~d1$period, method = "euc", strata = d1$name)
d2 <- dat %>% filter(variable == 'Rain', ecoregion1 == "Northwestern Forested Mountains")
ad2 <- adonis2(d2[,6]~d2$period, method = "euc", strata = d2$name)
d3 <- dat %>% filter(variable == 'Rain', ecoregion1 == "Taiga")
ad3 <- adonis2(d3[,6]~d3$period, method = "euc", strata = d3$name)

p.values <- c(ad1$`Pr(>F)`[1],ad2$`Pr(>F)`[1],ad3$`Pr(>F)`[1])
p.r <- data.frame(p = p.values, group = c('Northern Forests', 'Northwestern Forested Mountains','Taiga'), variable = rep('Rain',3))
p.r$p.adj <- p.adjust(p.r$p, method = "holm")
p.r

# Snowfall
d1 <- dat %>% filter(variable == 'Snowfall', ecoregion1 == "Northern Forests")
ad4 <- adonis2(d1[,6]~d1$period, method = "euc", strata = d1$name)
d2 <- dat %>% filter(variable == 'Snowfall', ecoregion1 == "Northwestern Forested Mountains")
ad5 <- adonis2(d2[,6]~d2$period, method = "euc", strata = d2$name)
d3 <- dat %>% filter(variable == 'Snowfall', ecoregion1 == "Taiga")
ad6 <- adonis2(d3[,6]~d3$period, method = "euc", strata = d3$name)

p.values <- c(ad4$`Pr(>F)`[1],ad5$`Pr(>F)`[1],ad6$`Pr(>F)`[1])
p.s <- data.frame(p = p.values, group = c('Northern Forests', 'Northwestern Forested Mountains','Taiga'), variable = rep('Snowfall',3))
p.s$p.adj <- p.adjust(p.s$p, method = "holm")
p.s

# Temperature
d1 <- dat %>% filter(variable == 'Temperature', ecoregion1 == "Northern Forests")
ad7 <- adonis2(d1[,6]~d1$period, method = "euc", strata = d1$name)
d2 <- dat %>% filter(variable == 'Temperature', ecoregion1 == "Northwestern Forested Mountains")
ad8 <- adonis2(d2[,6]~d2$period, method = "euc", strata = d2$name)
d3 <- dat %>% filter(variable == 'Temperature', ecoregion1 == "Taiga")
ad9 <- adonis2(d3[,6]~d3$period, method = "euc", strata = d3$name)

p.values <- c(ad7$`Pr(>F)`[1],ad8$`Pr(>F)`[1],ad9$`Pr(>F)`[1])
p.t <- data.frame(p = p.values, group = c('Northern Forests', 'Northwestern Forested Mountains','Taiga'), variable = rep('Temperature',3))
p.t$p.adj <- p.adjust(p.t$p, method = "holm")
p.t

# Snowmelt timing
d1 <- dat %>% filter(variable == 'Snowmelt timing', ecoregion1 == "Northern Forests")
ad10 <- adonis2(d1[,6]~d1$period, method = "euc", strata = d1$name)
d2 <- dat %>% filter(variable == 'Snowmelt timing', ecoregion1 == "Northwestern Forested Mountains")
ad11 <- adonis2(d2[,6]~d2$period, method = "euc", strata = d2$name)
d3 <- dat %>% filter(variable == 'Snowmelt timing', ecoregion1 == "Taiga")
ad12 <- adonis2(d3[,6]~d3$period, method = "euc", strata = d3$name)

p.values <- c(ad10$`Pr(>F)`[1],ad11$`Pr(>F)`[1],ad12$`Pr(>F)`[1])
p.st <- data.frame(p = p.values, group = c('Northern Forests', 'Northwestern Forested Mountains','Taiga'), variable = rep('Snowmelt timing',3))
p.st$p.adj <- p.adjust(p.st$p, method = "holm")
p.st

# Snowmelt timing
d1 <- dat %>% filter(variable == 'Rainfall Seasonality Index', ecoregion1 == "Northern Forests")
ad13 <- adonis2(d1[,6]~d1$period, method = "euc", strata = d1$name)
d2 <- dat %>% filter(variable == 'Rainfall Seasonality Index', ecoregion1 == "Northwestern Forested Mountains")
ad14 <- adonis2(d2[,6]~d2$period, method = "euc", strata = d2$name)
d3 <- dat %>% filter(variable == 'Rainfall Seasonality Index', ecoregion1 == "Taiga")
ad15 <- adonis2(d3[,6]~d3$period, method = "euc", strata = d3$name)

p.values <- c(ad13$`Pr(>F)`[1],ad14$`Pr(>F)`[1],ad15$`Pr(>F)`[1])
p.si <- data.frame(p = p.values, group = c('Northern Forests', 'Northwestern Forested Mountains','Taiga'), variable = rep('Rainfall Seasonality Index',3))
p.si$p.adj <- p.adjust(p.si$p, method = "holm")
p.si


### Stat output dataframe
stats <- full_join(p.r, p.s) %>% full_join(p.t) %>% full_join(p.st) %>% full_join(p.si)
stats$F <- c(ad1$F[1], ad2$F[1],ad3$F[1],ad4$F[1],ad5$F[1],ad6$F[1],ad7$F[1],ad8$F[1],ad9$F[1],ad10$F[1],ad11$F[1],ad12$F[1],ad13$F[1],ad14$F[1],ad15$F[1])
stats$DFr <- c(ad1$Df[2], ad2$Df[2],ad3$Df[2],ad4$Df[2],ad5$Df[2],ad6$Df[2],ad7$Df[2],ad8$Df[2],ad9$Df[2],ad10$Df[2],ad11$Df[2],ad12$Df[2],ad13$Df[2],ad14$Df[2],ad15$Df[2])
stats$DFt <- c(ad1$Df[3], ad2$Df[3],ad3$Df[3],ad4$Df[3],ad5$Df[3],ad6$Df[3],ad7$Df[3],ad8$Df[3],ad9$Df[3],ad10$Df[3],ad11$Df[3],ad12$Df[3],ad13$Df[3],ad14$Df[3],ad15$Df[3])
write.csv(stats, "C:/Users/kljorgenson/Documents/Repos/AK_discharge/Data/DFFT metrics/PERMANOVA_periods_climate.csv", row.names= F)



### Plot of difference between time periods
# Reformat data
dat %>% mutate(value2 = ifelse(variable %in% c('rain', 'snowfall'), log(value),value)) %>% group_by(variable) %>%
  mutate(value3 = scale(value2)) %>% ggplot(aes(variable, value3, col = period)) + geom_boxplot() + theme_bw() +
  xlab("") + ylab("Scaled value") + labs(col = NULL) + scale_color_manual(values = c(1,'grey55')) +
  theme(legend.position = 'top', legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-8,-10))
dat2 <- dat %>% 
  pivot_wider(names_from = 'period', values_from = 'value') 
names(dat2) <- c('station', 'name', 'ecoregion1', 'variable', 'recent', 'historical')
dat3 <- dat2 %>% mutate(diff.p = (recent - historical)/historical*100, diff = recent - historical) %>%
  mutate(ecoregion1 = ifelse(ecoregion1 == 'Hudson Plain', 'Northern Forests',ecoregion1)) %>% filter(is.na(ecoregion1) == F)


# Ecoregion colors
cols = c("#4E79A7", "#D95F02", "#66A61E", "#A6761D", "#F5C710")

# Plot
a <- dat3 %>% filter(!variable %in% c('Temperature', 'Precipitation', 'Snowmelt timing')) %>% 
  mutate(variable2 = factor(variable, levels = c('Rain', 'Snowfall', 'Rainfall Seasonality Index'), labels = c(c('Rain', 'Snowfall', 'Rainfall Seasonality Index (Jun-Oct)'))),
         variable2 = fct_relabel(variable2, ~str_wrap(.x, width = 15))) %>%
  ggplot(aes(variable2, diff.p, col = ecoregion1))+
  geom_hline(yintercept = 0, linetype = 2, col = 'grey') + geom_boxplot(outliers = F) + theme_bw() +
  xlab("") + ylab("Difference (%)") + labs(col = NULL) + scale_color_manual(values = cols, labels = function(x) str_wrap(x, width = 20)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 7), plot.margin = unit(c(0.1,0.2,0,0.5), 'cm')) +
  geom_point(size = 0.7, position = position_jitterdodge(jitter.width = 0.1), shape = 1)  + theme(legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=8)) +
  geom_point(aes(x = 0.85,y=40), shape = 8, col = "#D95F02") +
  geom_point(aes(x = 1,y=40), shape = 8, col = "#66A61E") +
  geom_point(aes(x = 1.85,y=40), shape = 8, col = "#D95F02") +
  geom_point(aes(x = 2,y=40), shape = 8, col = "#66A61E") +
  geom_point(aes(x = 2.15,y=40), shape = 8, col = "#A6761D")+
geom_point(aes(x = 3,y=40), shape = 8, col = "#66A61E")
a

b <- dat3 %>% filter(variable == 'Temperature') %>% ggplot(aes(variable, diff, col = ecoregion1))+
  geom_hline(yintercept = 0, linetype = 2, col = 'grey') + geom_boxplot(outliers = F) + theme_bw() +
  xlab("") + ylab("Difference (°C)") + labs(col = NULL) + scale_color_manual(values = cols, labels = function(x) str_wrap(x, width = 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 7), plot.margin = unit(c(0.1,0.3,0.65,0), 'cm')) +
  geom_point(size = 0.7, position = position_jitterdodge(jitter.width = 0.1), shape = 1) + ylim(-0.3,0.3) + scale_y_continuous(position = "left") + 
  theme(legend.key.size = unit(0.7, 'cm')) +
  geom_point(aes(x = 0.85,y=0.2), shape = 8, col = "#D95F02") +
  geom_point(aes(x = 1,y=0.2), shape = 8, col = "#66A61E") +
  geom_point(aes(x = 1.15,y=0.2), shape = 8, col = "#A6761D")
b

c <- dat3 %>% filter(variable == 'Snowmelt timing') %>% ggplot(aes(str_wrap(variable, width = 10), diff, col = ecoregion1))+
  geom_hline(yintercept = 0, linetype = 2, col = 'grey') + geom_boxplot(outliers = F) + theme_bw() +
  xlab("") + ylab("Difference (days)") + labs(col = NULL) + scale_color_manual(values = cols, labels = function(x) str_wrap(x, width = 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 7), plot.margin = unit(c(0.1,0.3,0.3,0.1), 'cm')) +
  geom_point(size = 0.7, position = position_jitterdodge(jitter.width = 0.1), shape = 1) + ylim(-0.3,0.3) + scale_y_continuous(position = "left") + 
  theme(legend.key.size = unit(0.7, 'cm'))  +
  geom_point(aes(x = 1,y=-15), shape = 8, col = "#66A61E") 
c

ggarrange(a,c,b, widths = c(1,0.44,0.42), common.legend = T, legend = 'top', nrow = 1) 

ggsave("Figures/Figure_9.png", width = 6, height = 2.5, bg = "white")
