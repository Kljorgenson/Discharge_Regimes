

########################## Example plot of rivers contrasting in seasonality and flow variability
### High and low seasonality rivers
flow.dat <- read.csv("C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/Discharge_data/DFFT_20_ak_can.csv") %>% 
  mutate(site2 = case_when(site == "ALSEK RIVER ABOVE BATES RIVER" ~ "Alsek (glacial)",
                           site == "WATERFOUND RIVER BELOW THERIAU LAKE" ~ "Waterfound (non-glacial)",
                           site == "RUSSELL" ~ "Russell (non-glacial)",
                           site == "ANTLER" ~ "Antler (glacial)")) %>%
  filter(site2 %in% c("Alsek (glacial)", "Waterfound (non-glacial)", "Russell (non-glacial)", "Antler (glacial)"))

# Seasonal signal
sea <- flow.dat %>% group_by(site2,jday) %>% summarise(pred = mean(pred2, na.rm = T))

mets <- metric.d %>% rename(site = name) %>% mutate(site2 = case_when(site == "ALSEK RIVER ABOVE BATES RIVER" ~ "Alsek (glacial)",
                                                                      site == "WATERFOUND RIVER BELOW THERIAU LAKE" ~ "Waterfound (non-glacial)",
                                                                      site == "RUSSELL" ~ "Russell (non-glacial)",
                                                                      site == "ANTLER" ~ "Antler (glacial)")) %>%
  filter(site2 %in% c("Alsek (glacial)", "Waterfound (non-glacial)", "Russell (non-glacial)", "Antler (glacial)"))


flow.dat$site2 <- factor(flow.dat$site2, levels = c("Alsek (glacial)", "Antler (glacial)", "Waterfound (non-glacial)", "Russell (non-glacial)"))
sea$site2 <- factor(sea$site2, levels = c("Alsek (glacial)", "Antler (glacial)", "Waterfound (non-glacial)", "Russell (non-glacial)"))
mets$site2 <- factor(mets$site2, levels = c("Alsek (glacial)", "Antler (glacial)", "Waterfound (non-glacial)", "Russell (non-glacial)"))


SAM <- mets %>% select(HSAM, HSAM_t, site2) %>% mutate(jday = round(HSAM_t))
SAMs <- left_join(SAM,sea) %>% mutate(val = HSAM+pred)



SAM <- flow.dat %>%
  filter(site2 %in% c("Alsek (glacial)", "Antler (glacial)", "Waterfound (non-glacial)", "Russell (non-glacial)")) %>% 
  filter(is.na(resid.sig) == F, jday >= 100 & jday <= 304) %>% group_by(site2, year)  %>% 
  filter(length(resid.sig) >= 200, # Not missing big chunk of summer
         resid.sig == max(resid.sig, na.rm = T) | resid.sig == min(resid.sig, na.rm = T)) %>% 
  mutate(type = ifelse(resid.sig > 0, "HSAM", "LSAM")) %>% dplyr::select(date,year,jday,resid.sig,type, pred2)%>% mutate(val = resid.sig + pred2) %>% filter(type == "HSAM")

labs <- data.frame(site2 = c("Alsek (glacial)", "Antler (glacial)", "Waterfound (non-glacial)", "Russell (non-glacial)"), xlab = c("Low stochasticity of high flows", "High stochasticity of high flows", "", ""),
                   ylab = c("", " High signal:noise \nHigh seasonality", "", " Low signal:noise \nLow seasonality"))
labs$site2 <- factor(labs$site2, levels = c("Alsek (glacial)", "Antler (glacial)", "Waterfound (non-glacial)", "Russell (non-glacial)"))


flow.dat   %>% ggplot(aes()) + 
  geom_line(aes(jday,pred2+resid.sig, group = as.factor(year)), size = 0.3,col = "grey50") + geom_line(data = sea, aes(jday,pred), col = "red") + facet_wrap(~site2) +
  theme_bw() + theme(strip.background = element_rect(fill = 'white'),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("Normalized discharge") + xlab("Day of year") +
  geom_text(data = mets, aes(x = 57, y = 2, label = paste("=", round(sigma.hfb, digits = 2))), size = 3) +
  geom_text(data = mets, x = 8, y = 2, label = expression(sigma ["hf"]), size = 3) +
  geom_text(data = mets, aes(x = 50, y = 2.4, label = paste("SNR =", round(snr, digits = 2))), size = 3) +
  geom_text(data = mets, x = 18, y = 2.2, label = expression(A[RMS]), size = 3) +
  geom_text(data = mets, aes(x = 80, y = 2.2, label = paste("=", round(rms.signal, digits = 2))), size = 3) +
  geom_point(data = SAM, aes(jday,val), size = 1, col = "blue", shape = 1) +
  geom_text(data = labs, y = 1.3, x = 420,
            size = 4.1, aes(label = ylab), angle = 270) +
  geom_text(data = labs, y = 3.3, x = 180,
            size = 4.1, aes(label = xlab)) +
  coord_cartesian(xlim = c(0, 360), ylim = c(0.25, 2.5),# This focuses the x-axis on the range of interest
                  clip = 'off') +   # This keeps the labels from disappearing
  theme(plot.margin = unit(c(3,3.5,1,1), "lines"), strip.text.x = element_text(size = 10))

ggsave("Figures/Regimes paper/Example rivers.png", width = 6, height = 5)


