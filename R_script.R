### Load the packages
library(tidyverse)# For the packages dplyr and ggplot2
library(MuMIn)# For the function AICc
library(glmmTMB)# For the function glmmTMB
library(DHARMa)# For the function simulateResiduals
library(emmeans)# For the function emmeans
library(ggcorrplot)# For the function ggcorrplot

### Load the data
dta <- read.table("exposure_dry_phase.txt", h = T, sep = ";")

### Variables as factor
dta$flume <- factor(dta$flume)
dta$dry_timef <- factor(dta$dry_time)

######## 1. Hatching proportion (mean +- SD) #########################
# Dry phases: 0 - 5 h
dta %>%
  filter(dry_time %in% c(0:5))%>%
  group_by(dry_time) %>%
  summarise(moy=round(mean(prop),2),
            SD = round(sd(prop),2))

dta %>%
  filter(dry_time %in% c(0:5))%>%
  summarise(moy=round(mean(prop),2),
            SD = round(sd(prop),2))

# Dry phase: 5 h for the slate 29
dta %>%
  filter(dry_time == 5 & sub == "s29")%>%
  group_by(sub)%>%
  summarise(moy=round(mean(prop),2),
            SD = round(sd(prop),2))

# Dry phase: 6 h
dta %>%
  filter(dry_time %in% 6)%>%
  summarise(moy=round(mean(prop),2),
            SD = round(sd(prop),2))

# Dry phase: 6 h for the slates 48, 40, 22
dta %>%
  filter(sub %in% c("s48", "s40", "s22"))%>%
  group_by(sub) %>%
  summarise(moy=round(mean(prop),2),
            SD = round(sd(prop),2))

# Dry phases: 9 - 18 h
dta %>%
  filter(dry_time %in% c(9:18))%>%
  summarise(moy=round(mean(prop),2),
            SD = round(sd(prop),2))

# Dry phases: 9 - 18 h with the slates 17, 35, 50
dta %>%
  filter(sub %in% c("s17", "s35", "s50"))%>%
  group_by(sub) %>%
  summarise(moy=round(mean(prop),2),
            SD = round(sd(prop),2))


######## 2. Modelling #########################

### Add 'egg_mass’ as random effect nested within the ‘flume’ random effect
dta2 <- dta %>% 
  group_by(flume) %>% 
  mutate(egg_mass = sprintf("m%d", 1:n())) %>% 
  ungroup() %>% 
  mutate(egg_mass = factor(egg_mass))

  ### 2.1 Model with 'dry_time' as factor ###
s1d <- glmmTMB(cbind(hatch, unhatch) ~ 
                 dry_timef + (1|flume/egg_mass),
               data = dta2,
               family = binomial())


## Look at the residuals of the model
res1d <- simulateResiduals(s1d, refit = T, n = 20)
plot(res1d)

  
### 2.2 Model with 'dry_time' as integer ###
sinteger <- glmmTMB(cbind(hatch, unhatch) ~ 
                      dry_time + (1|flume/egg_mass),
                    data = dta2,
                    family = binomial())

### 2.3 Comparison of the AICc between models ###
AICc(sinteger, s1d)
              # AICc of the model with 'dry_time' as factor is INFERIOR to the AICc of the model with 'dry_time' as integer
              # Thus, the model with 'dry_time' as factor was used in the following analysis

### 2.4 Difference in hatching egg probability (on the logit scale) between dry-phase duration
p1 <- pairs(emmeans(s1d,  ~ dry_timef), adjust = "fdr")

### Create a summary table
p1 <- summary(p1) %>%
  as.data.frame(p1)
p1 <- do.call(rbind, strsplit(p1$contrast, " - ")) %>% 
  as.data.frame() %>%
  mutate_all(~gsub("\\D+", "", .)) %>%
  cbind(p1[,-1])

### 2.5 Table of the z.ratio for each difference in hatching egg probability between dry-phase duration
nc <- nlevels(dta2$dry_timef)
zratios <- matrix(numeric(), nc, nc)
diag(zratios) <- 0
colnames(zratios) <- rownames(zratios) <- 
  levels(dta2$dry_timef)
zratios[as.matrix(p1[1:2])] <- p1$z.ratio
zratios[as.matrix(p1[2:1])] <- p1$z.ratio
dimnames(zratios) <- lapply(dimnames(zratios),
                         function(x) sprintf("d%s", x))

### 2.6 Table of the p values for each difference in hatching egg probability between dry-phase duration
pvalues <- matrix(numeric(), nc, nc)
diag(pvalues) <- 0
colnames(pvalues) <- rownames(pvalues) <- 
  levels(dta2$dry_timef)
pvalues[as.matrix(p1[1:2])] <- p1$p.value
pvalues[as.matrix(p1[2:1])] <- p1$p.value
dimnames(pvalues) <- lapply(dimnames(pvalues),
                         function(x) sprintf("d%s", x))

### P values in detail
tab_pvalue <- p1[,c("V1","V2", "p.value")]
tab_pvalue$p.value <- round(tab_pvalue$p.value, 2)
tab_pvalue <- tab_pvalue %>% 
  rename(Dry_time1 = V1, Dry_time2 = V2)

  

# 0 - 5 h of dry-phases
tab_pvalue %>%
  filter(Dry_time1 %in% c(0:5) & Dry_time2 %in% c(0:5)) %>%
  summarise(min(p.value), max(p.value))

# 6 h and the other dry phase durations
tab_pvalue %>%
  filter(Dry_time1 %in% 6 | Dry_time2 %in% 6) %>%
  summarise(min(p.value), max(p.value))

# 9 - 18 h of dry phase durations
tab_pvalue %>%
  filter(Dry_time1 %in% 9:18 & Dry_time2 %in% 9:18) %>%
  summarise(min(p.value), max(p.value))

# 9 - 18 h and the other dry phase durations
tab_pvalue %>%
  filter(Dry_time1 %in% 0:6 & Dry_time2 %in% c(9,12,18)) %>%
  summarise(min(p.value), max(p.value))


### 2.7 Table of the difference in hatching egg probability
# between dry-phase duration
nc <- nlevels(dta2$dry_timef)
differences <- matrix(numeric(), nc, nc)
diag(differences) <- 0
colnames(differences) <- rownames(differences) <- 
  levels(dta2$dry_timef)
differences[as.matrix(p1[1:2])] <- p1$estimate
differences[as.matrix(p1[2:1])] <- p1$estimate
dimnames(differences) <- lapply(dimnames(differences),
                         function(x) sprintf("d%s", x))

### 2.8 Look at the z values
# The cross for significative p-values
ggcorrplot(zratios, hc.order = FALSE,
           type = "lower", p.mat = (pvalues < 0.05)*1,
           pch.col = "white") +
  scale_fill_viridis_c() + 
  labs(fill = "z value")

### 2.9 Look at the difference in hatching egg probability (on the logit scale) between dry-phase duration
# The cross for significative p-values
ggcorrplot(differences, hc.order = FALSE,
           type = "lower", p.mat = (pvalues<0.05)*1,
           pch.col="white") +
  scale_fill_viridis_c() + 
  labs(fill = "difference")

######## 3. Figures #########################

### 3.1 Figure 6
ggplot(dta2, aes(dry_time, prop, group = dry_timef))+
  geom_boxplot(outlier.colour =  "transparent")+
  geom_jitter(width = 0.2, alpha=0.15, size = 3)+
  xlab("Dry-phase duration (h)")+
  ylab("Proportion hatched")+
  annotate("text",
           x = c(0,1,2,3,4,5,6,9,12,18),
           y = c(0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07),
           label = paste0("N=",table(dta2$dry_time)),
           col = "black", 
           size=2.9)+
  theme(
    panel.grid.major = element_line(color="#CCCCCC"),
    panel.grid.minor = element_line(color="white"),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14))+
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 9, 12, 18),
                     labels = c(0, 1, 2, 3, 4, 5, 6, 9, 12, 18))

### 3.2 Figure S5
count_df <- dta2 %>%
  group_by(sub, dry_time) %>%
  summarize(count = n())

ggplot(data = dta2, aes(as.factor(sub), prop, fill = sub)) + 
  geom_boxplot(width = 0.2, outlier.colour =  "transparent", fill = NA) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha = 0.15, size=3) +
  ylab("Proportion hatched") +
  facet_wrap(~dry_time, ncol = 5, scales = "free_x") +
  theme(aspect.ratio = 1) +
  labs(x = "Slates")+
  geom_text(data = count_df, 
            aes(label = paste("N =", count), 
                x = as.factor(sub), 
                y = 0.09), 
            position = position_dodge(width = 0.8), vjust = -1, size = 3.5, color = "black")+
  ggtitle("Dry-phase duration (h)") +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, vjust = -1.5), 
    panel.grid.major = element_line(color="#CCCCCC"),
    panel.grid.minor = element_line(color="white"),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.position = "none")
