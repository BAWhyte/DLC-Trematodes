### - MAX PLANCK INSTITUTE OF ANIMAL BEHAVIOR PROJECT [WORM - ###
# AUTHOR: Brian Whyte, PhD, UC Berkeley
# EMAIL: ba.whyte@berkeley.edu
# START DATE: Nov. 18, 2021
# ---------------------------------------------------------------------------------------------------------------- 

# BLOCK 0: Set up + clean up -------------------------------------------------------------------------------------

## Housekeeping
cat("\014") # Clear console
rm(list=ls()) # Remove all variables
library(ggplot2) # for most graphs
library(tidyverse) # for multiple data transforming packages
library(gganimate) # for animated ggplots
library(lubridate) # for time conversion into hh:mm:ss
library(amt) # for spatial data analysis w/ dplyr commands
library(ggpubr) # for gg density

## Set working directory + load data
setwd("/Example/path")
data <- read.csv("YourDeepLabCutOutput.csv")
df <- tibble(data)

## Function for treating frame as seconds and converting to "D HH:MM:SS" format
dhms <- function(t){
  paste(t %/% (60*60*24) 
        ,paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0")
               ,formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0")
               ,formatC(t %% 60, width = 2, format = "d", flag = "0")
               ,sep = ":"
        )
  )
}

## Transform raw DLC .csv data frame into variables vs. observations
{ 
colnames(df) <- paste(df[1,],df[2,],df[3,],sep="-") # store character info from rows into colnames
colnames(df)[1] <- "frame" # rename first column to what it is (frame number)
df <- df[-c(1:3),] # remove character info rows
df <- df %>% # move character info into ID and BP columns and stack all data
  pivot_longer(cols = starts_with("HIMA"), # this stacks all data into columns ("long", not "wide")
               names_to = c("ID","bp",".value"), # find this pattern in the colnames
               names_pattern = "(.*)-(.*)-(.*)") %>%  # use all info inbetween "-" seperators
  arrange(ID, bp) %>% # rearrange row order by ID and bp
  mutate_at(c("frame","x","y","likelihood"), as.numeric) # make numeric columns numeric
df$frame <- df$frame %>% as.integer # make frame column have no decimals
df <- na.omit(df) # remove NA rows
df <- df[df$likelihood == 1,] # use only data points with high likelihood
df["dhms"] <- dhms(df$frame) # add conversions as their own column
df$dhms <- substring(df$dhms, 3) # remove 0 in front of HH:MM:SS time format
df$dhms <- as.POSIXct(strptime(df$dhms, "%H:%M:%S")) # assume date is the same as system 
} 

## Pixel dimensions according to file resolution
width <- 1024
height <- 768

## BLOCK 1: Density plots (heat maps) ----------------------------------------------------------------------------
# Source: https://www.r-graph-gallery.com/2d-density-plot-with-ggplot2.html
xr <- c(200,700) # arbitrary axis ranges
yr <- c(0,650)

## Heat map with contounrs, resolution limits
{
gd1 <- ggplot(df[df$ID == "HIMA1" & df$bp == "mouth",], aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_distiller(palette="Blues", direction=-1) +
  xlim(0,width) +
  ylim(0,height)
gd2 <- ggplot(df[df$ID == "HIMA2" & df$bp == "mouth",], aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_distiller(palette="Greens", direction=-1) +
  xlim(0,width) +
  ylim(0,height)
gd3 <- ggplot(df[df$ID == "HIMA3" & df$bp == "mouth",], aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_distiller(palette="PuRd", direction=-1) +
  xlim(0,width) +
  ylim(0,height)
}

## Heat map with white = hot intensity, and no contours
#gd <- ggplot(df[df$ID == "HIMA3" & df$bp == "mouth",], aes(x=x, y=y) ) +
gd1 <- ggplot(hbb1, aes(x=x_, y=y_) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette="Greys", direction=-1) +
  #scale_fill_manual(gray.colors(7)) +
  theme_linedraw() +
  xlim(xr) +
  ylim(yr)  
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0))
gd2 <- ggplot(hbb2, aes(x=x_, y=y_) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette="Greys", direction=-1) +
  theme_linedraw() +
  xlim(xr) +
  ylim(yr)  
#scale_x_continuous(expand = c(0, 0)) +
#scale_y_continuous(expand = c(0, 0))
gd3 <- ggplot(hbb3, aes(x=x_, y=y_) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette="Greys", direction=-1) +
  theme_linedraw() +
  xlim(xr) +
  ylim(yr)  
#scale_x_continuous(expand = c(0, 0)) +
#scale_y_continuous(expand = c(0, 0))

## BLOCK 2: Animated scatter plots ----------------------------------------------------------------------------

## One HIMA, only mouth
{
p <- ggplot(df[df$ID == "HIMA3" & df$bp == "mouth",], aes(x=x, y=y, color=ID) ) +
  geom_point(alpha = 0.7, show.legend = FALSE) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size = 20)) +
  scale_color_manual(values = c("#A550A0")) +
  # Here comes the gganimate specific bits
  labs(title = 'Frame: {frame_time}') +
  transition_time(frame) +
  shadow_wake(wake_length = 0.1, alpha = FALSE)
animate(p, fps = 5, nframes = 200)
}

## Re-sampling and displacement calculations for each ID (only mouth)
hbb <- read_csv("HABIBI.csv")
hbb1 <- hbb %>%
  filter(ID == "HIMA1", bp == "mouth") %>%
  make_track(x, y, frame, id = ID) %>% # needed to swap "dhms" with "frame"
  track_resample(., rate = seconds (1), tolerance = seconds(0.01)) %>% # use only 1 frame per 1 sec?
  mutate(step = step_lengths (.)) %>% # Displacement calculation?
  filter(step <= quantile(step, .95, na.rm = TRUE)) # removing anything outside .95 freq. quantile?
hbb2 <- hbb %>%
  filter(ID == "HIMA2", bp == "mouth") %>%
  make_track(x, y, frame, id = ID) %>% # needed to swap "dhms" with "frame"
  track_resample(., rate = seconds (1), tolerance = seconds(0.01)) %>% # use only 1 frame per 1 sec?
  mutate(step = step_lengths (.)) %>% # Displacement calculation?
  filter(step <= quantile(step, .95, na.rm = TRUE)) # removing anything outside .95 freq. quantile?
hbb3 <- hbb %>%
  filter(ID == "HIMA3", bp == "mouth") %>%
  make_track(x, y, frame, id = ID) %>% # needed to swap "dhms" with "frame"
  track_resample(., rate = seconds (1), tolerance = seconds(0.01)) %>% # use only 1 frame per 1 sec?
  mutate(step = step_lengths (.)) %>% # Displacement calculation?
  filter(step <= quantile(step, .95, na.rm = TRUE)) # removing anything outside .95 freq. quantile?
df <- rbind(hbb1,hbb2,hbb3) # combine all hbb's by stacking them
colnames(df) <- c("x","y","frame","ID","burst","step")
df$frame <- df$frame %>% as.integer # make frame column have no decimals

## All HIMA, only mouth (w/ filtered outlier displacements)
p <- ggplot(df, aes(x=x, y=y, color=ID) ) +
  geom_point(alpha = 0.7, show.legend = TRUE) +
  theme_linedraw() +
  #theme(axis.text.x = element_blank()) +
  #theme(axis.text.y = element_blank()) +
  xlab("") +
  ylab("") +
  theme(legend.position = "") +
  scale_color_manual(values = c("#1ACCEF","#6A8EDC","#D044DA")) +
  # Here comes the gganimate specific bits
  labs(title = 'Frame: {frame_time}') +
  transition_time(frame) +
  shadow_wake(wake_length = 0.1, alpha = FALSE)
#animate(p, fps = 8, nframes = 200, end_pause = 8)
animate(p, fps=7, height=4, width=4, units="in", res=100, end_pause=8)
anim_save("HIMAS_scatter.gif")

## One HIMA, all body parts
p <- ggplot(df[df$ID == "HIMA3",], aes(x=x, y=y, color=ID) ) +
  geom_point(alpha = 1, show.legend = TRUE) +
  theme_linedraw() +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  xlab("") +
  ylab("") +
  scale_color_manual(values = c("#A550A0")) +
  # Here comes the gganimate specific bits
  labs(title = 'Frame: {frame_time}') +
  transition_time(frame) +
  shadow_wake(wake_length = 0.05, alpha = 0.5)
#animate(plot=p,nframes=length(unique(df$frame)),detail = 8,fps = 8,end_pause = 8)
animate(plot=p, detail=8, fps=8, end_pause=8)
anim_save("HIMA3_WholeBody.gif") # saves current animation in viewer

## BLOCK 3: Traj object and displacements (Dist) -----------------------------------------------------------------

## Transform into a trajectory and plot these trajectories
coord <- cbind(df[df$ID == "HIMA1" & df$bp == "mouth",]$x, df[df$ID == "HIMA1" & df$bp == "mouth",]$y)
traj <- as.ltraj(xy = coord, date = df[df$ID == "HIMA1" & df$bp == "mouth",]$dhms, id = "HIMA1")
#plot(traj) # default plot reads unique ltraj object
#plotltr(traj, 'dist') # shows displacement per frame (is this just one column vs. time?)
tf <- tibble(traj[[1]]) # within ltraj is a usable data frame though
tf <- tf[-c(nrow(tf)),] # remove last row with NA dist
H1M <- tf[tf$dist < 5,] # only displacements per frame of 10 units or less

## Same as above but for HIMA3 mouth
coord <- cbind(df[df$ID == "HIMA3" & df$bp == "mouth",]$x, df[df$ID == "HIMA3" & df$bp == "mouth",]$y)
traj <- as.ltraj(xy = coord, date = df[df$ID == "HIMA3" & df$bp == "mouth",]$dhms, id = "HIMA3")
#plot(traj) # default plot reads unique ltraj object
#plotltr(traj, 'dist') # shows displacement per frame (is this just one column vs. time?)
tf <- tibble(traj[[1]]) # within ltraj is a usable data frame though
H3M <- tf[tf$dist < 5,] # only displacements per frame of 10 units or less

## Dist vs. time (filtered)
ggplot(H3M, aes(x=date, y=dist)) +
  geom_point()

## Animated scatter plot (HIMA 1 mouth using traj frame instead of df)
H1M$ID <- rep("A",nrow(H1M)) # put placeholder column to call for point colors
p <- ggplot(H1M, aes(x=x, y=y, color=ID) ) +
  geom_point(alpha = 0.7, show.legend = FALSE) +
  theme_linedraw() +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  xlab("") +
  ylab("") +
  scale_color_manual(values = c("#0FB1D2")) +
  # Here comes the gganimate specific bits
  labs(title = 'Date: {frame_time}') +
  transition_time(date) +
  shadow_wake(wake_length = 0.1, alpha = FALSE)
animate(p, fps = 5, nframes = 200)

## Histogram of dist values
ggdensity(hbb1, x="step", fill="ID") +
  scale_fill_manual(values=c("#0FB1D2")) +
  labs(title = "Displacement per frame for HIMA 2 mouth")

## BLOCK 3: Alternative method using "amt" -----------------------------------------------------------------------

## Read habibi and convert to simple df
hbb <- read_csv("HABIBI.csv")
hbb1 <- hbb %>%
  filter(ID == "HIMA1", bp == "mouth") %>%
  make_track(x, y, frame, id = ID) %>% # needed to swap "dhms" with "frame"
  track_resample(., rate = seconds (1), tolerance = seconds(0.01)) %>% # use only 1 frame per 1 sec?
  mutate(step = step_lengths (.)) %>% # displacement calculation?
  filter(step <= quantile(step, .95, na.rm = TRUE)) # removing anything outside .95 freq. quantile?
hbb1$sum <- cumsum(hbb1$step) # cumulative distance
hbb2 <- hbb %>%
  filter(ID == "HIMA2", bp == "mouth") %>%
  make_track(x, y, frame, id = ID) %>% # needed to swap "dhms" with "frame"
  track_resample(., rate = seconds (1), tolerance = seconds(0.01)) %>% # use only 1 frame per 1 sec?
  mutate(step = step_lengths (.)) %>% # displacement calculation?
  filter(step <= quantile(step, .95, na.rm = TRUE)) # removing anything outside .95 freq. quantile?
hbb2$sum <- cumsum(hbb2$step) # cumulative distance
hbb3 <- hbb %>%
  filter(ID == "HIMA3", bp == "mouth") %>%
  make_track(x, y, frame, id = ID) %>% # needed to swap "dhms" with "frame"
  track_resample(., rate = seconds (1), tolerance = seconds(0.01)) %>% # use only 1 frame per 1 sec?
  mutate(step = step_lengths (.)) %>% # displacement calculation?
  filter(step <= quantile(step, .95, na.rm = TRUE)) # removing anything outside .95 freq. quantile?
hbb3$sum <- cumsum(hbb3$step) # cumulative distance
lf <- rbind(hbb1,hbb2,hbb3) # combine all hbb's by stacking them
colnames(lf) <- c("x","y","frame","ID","burst","step","sum")

## Cumulative movement line graph
p <- ggplot(lf, aes(frame, sum, color=ID)) +
  geom_line() +
  theme_linedraw() +
  scale_color_manual(values = c("#1ACCEF","#6A8EDC","#D044DA")) +
  labs(x = "Frame", y = "Cumulative movement") +
  theme(legend.position = "top") +
  geom_point() + 
  transition_reveal(frame) 
animate(p, fps=7, height=4, width=4, units="in", res=100, end_pause=8)
anim_save("HIMAS_Line2.gif")

## Step scatter w/ linear trend
hbb1 %>%
  ggplot(aes(x=t_, y= step )) + 
  geom_point(alpha=0.1) + 
  geom_smooth(method = "lm")

## Relative angles (histogram/density plot)
hbb1 <- hbb %>%
  filter(ID == "HIMA1", bp == "mouth") %>%
  make_track(x, y, frame, id = ID) %>% # needed to swap "dhms" with "frame"
  track_resample(., rate = seconds (4), tolerance = seconds(0.01)) %>% # use only 1 frame per 1 sec?
  mutate(rel_ang = direction_rel(.)) # %>%
  #ggplot(aes(x=t_,y=rel_ang)) + geom_point(alpha=0.1)
hbb2 <- hbb %>%
  filter(ID == "HIMA2", bp == "mouth") %>%
  make_track(x, y, frame, id = ID) %>% # needed to swap "dhms" with "frame"
  track_resample(., rate = seconds (4), tolerance = seconds(0.01)) %>% # use only 1 frame per 1 sec?
  mutate(rel_ang = direction_rel(.)) # %>%
#ggplot(aes(x=t_,y=rel_ang)) + geom_point(alpha=0.1)
hbb3 <- hbb %>%
  filter(ID == "HIMA3", bp == "mouth") %>%
  make_track(x, y, frame, id = ID) %>% # needed to swap "dhms" with "frame"
  track_resample(., rate = seconds (4), tolerance = seconds(0.01)) %>% # use only 1 frame per 1 sec?
  mutate(rel_ang = direction_rel(.)) # %>%
#ggplot(aes(x=t_,y=rel_ang)) + geom_point(alpha=0.1)
df <- rbind(hbb1,hbb2,hbb3) # combine all hbb's by stacking them
colnames(df) <- c("x","y","frame","ID","burst","rel_ang")

ggdensity(df, x="rel_ang", fill="ID") +
  scale_fill_manual(values = c("#1ACCEF","#6A8EDC","#D044DA")) +
  labs(title = "Distribution of rel_ang per 4 sec (HIMA 1 mouth)")

## Calculating (angle?) correlations in every minute using moving window method
#hbb1 <- read_csv("HABIBI.csv")
hbb <- read_csv("HABIBI_WX.csv")
hbb1 <- hbb %>% 
  filter(ID == "HIMA1", bp == "mouth") %>%
  make_track(x, y, dhms, id = ID) %>%
  track_resample(., rate = seconds (1), tolerance = seconds(0.01))
#hbb1 <- na.omit(hbb1) # remove NA rows

cor = data.frame()
for (i in 1:(nrow(hbb1)-60)) {
  cor.i = tac(hbb1[i:(i+60),]) # for every 60 rows?
  cor = rbind(cor, data.frame(time = hbb1[i,3], cor = cor.i))
  #print(paste("i = ",i, "/ 85,096"))
}
colnames(cor) <- c("t_","ang_cor")
cor %>%
  ggplot(aes(x=t_, y= ang_cor)) + 
  geom_point(alpha=0.1) + 
  geom_smooth(method = "lm")


