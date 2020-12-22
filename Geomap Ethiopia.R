####
setwd("~/Dropbox/PGML_Projects")
weather <- read.csv("Weather data.csv", header = TRUE)
weather$Date <- as.Date(weather$Date, "%m/%d/%y")
####
library(ggplot2)
library(ggmap)
library(grid)

labs <- data.frame(env=c("Kobo","Meiso","Sheraro"),
                   lat=c(12.015, 9.233, 14.383),
                   long=c(39.633, 40.75, 37.767))
  
world <- map_data("world")

Ethiopia_bbox <- c(left = 31.5, bottom = 2, right = 71.5, top = 18)
Ethiopia_main_map <- get_stamenmap(Ethiopia_bbox, zoom = 5, maptype = "terrain-background")

tiff("Environmental Fig.tiff", width = 6.5, height = 3.25, units = "in", res = 300)

# Basic geographical map
p_main <- ggmap(Ethiopia_main_map) +
  theme_void() +
  theme(axis.title = element_blank(), 
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_hline(yintercept=15, linetype="dashed", color = "floralwhite") +
  geom_hline(yintercept=10, linetype="dashed", color = "floralwhite") +
  geom_hline(yintercept=5, linetype="dashed", color = "floralwhite") +
  geom_point(data = labs, aes(x = long, y = lat), color = "black", size = 1) +
  geom_text(data = labs, aes(x = long, y = lat, label = env), 
            size = 3.5, vjust = 0.5, hjust = -0.2) +
  geom_text(x=32, y=15, label="15° N", size=3.5, vjust = 0.5, hjust=0) +
  geom_text(x=32, y=10, label="10° N", size=3.5, vjust = 0.5, hjust=0) +
  geom_text(x=32, y=5, label="5° N", size=3.5, vjust = 0.5, hjust=0)

# Raninfall and temp data for three envs
p_rainfall <- ggplot(weather, aes(x = Date)) +
  geom_line(aes(y = MaxT), color = "red") + 
  geom_line(aes(y = MinT), color = "orange") + 
  geom_col(aes(Date, Rainfall/3), color="steelblue1")+
  geom_vline(xintercept = as.numeric(weather$Date[182]), linetype="dotted", color = "blue", size=1) +
  scale_y_continuous(
    # Features of the first axis
    name = "Temp (°C)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*3, name="Precipitation (mm)")) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b") +
  facet_grid(vars(Env), scales = "free") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=8, color = "black"),
        axis.text.y = element_text(size=8, color = "black"),
        axis.title.x = element_text(size=8, color = "black"),
        axis.title.y = element_text(size=8, color = "black"),
        strip.text.y = element_text(size=8, color="black"),
        strip.background = element_rect(colour="black", fill="#CCCCFF"))

p_main + 
  inset(ggplotGrob(p_rainfall), xmin = 48.5, xmax = 71.5, ymin = 2, ymax = 18)

dev.off()





