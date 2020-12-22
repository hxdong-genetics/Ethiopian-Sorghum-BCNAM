tiff("Marker-Trait Associations.tiff", units = "in", width = 11, height = 11, res = 300)

p1 <- ggplot(don1, aes(x=BPcum)) +
  facet_wrap(~Trait, ncol = 1, strip.position = "right", scales = "free_y") +
  geom_vline(stay.green, mapping = aes(xintercept = BPcum), lty=1, color="yellowgreen", alpha=0.3) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=1, color="black") +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  
  #Ma6
  geom_vline(data = filter(QTL.GWAS, Trait=="DF" | Trait=="DM" | Trait=="HE" | Trait=="NL" | Trait=="NT"), aes(xintercept = 356714466), 
             lty=1, lwd=1, alpha=1, color="orange") +
  #annotate("text", x=356714466, y = "Sheraro_JL", label="Ma6") +
  
  # Ma5
  #geom_vline(data = filter(QTL.GWAS, Trait=="DF" | Trait=="GYP" | Trait=="HE" | Trait=="NL" | Trait=="PH" | Trait=="TSW"), aes(xintercept = 6745069), 
  #           lty=1, lwd=1, alpha=1, color="orange") +
  #annotate("text", x=23145069, y =6.5, label="Ma5") +
  
  # NAT/CAX7
  geom_vline(data = filter(QTL.GWAS, Trait=="DF" | Trait=="GYP" | Trait=="HE" | Trait=="NL" | Trait=="PH" | Trait=="TSW"), aes(xintercept = 11158453), 
             lty=1, lwd=1, alpha=1, color="orange") +
  #annotate("text", x=23145069, y =6.5, label="Ma5") +
  
  # Ma3
  geom_vline(data = filter(QTL.GWAS, Trait=="DF" | Trait=="DM"| Trait=="PH" | Trait=="TSW"), aes(xintercept = 60852774), 
             lty=1, lwd=1, alpha=1, color="orange") +
  #annotate("text", x=23145069, y =6.5, label="Ma3") +
  
  # LHY
  geom_vline(data = filter(QTL.GWAS, Trait=="DF" | Trait=="DM"| Trait=="GYP"), aes(xintercept = 422952336), 
             lty=1, lwd=1, alpha=1, color="orange") + # LHY
  #annotate("text", x=23145069, y =6.5, label="LHY") +
  
  # P5CS2
  geom_vline(data = filter(QTL.GWAS, Trait=="DM"| Trait=="GYP" | Trait=="HE" | Trait=="LS" | Trait=="NL" | Trait=="NT" | Trait=="PH" |Trait=="TSW"), aes(xintercept = 218863379), 
             lty=1, lwd=1, alpha=1, color="orange") + # LHY
  #annotate("text", x=23145069, y =6.5, label="P5CS2") +
  
  # SbCN8
  geom_vline(data = filter(QTL.GWAS, Trait=="DF" | Trait=="LS" |Trait=="PH"), aes(xintercept = 592828477), 
             lty=1, lwd=1, alpha=1, color="orange") + # LHY
  #annotate("text", x=23145069, y =6.5, label="SbCN8") +
  
  # DREB1A
  geom_vline(data = filter(QTL.GWAS, Trait=="DF" | Trait=="GYP" | Trait=="NL" | Trait=="PH" | Trait=="TSW"), aes(xintercept = 478434065), 
             lty=1, lwd=1, alpha=1, color="orange") +
  
  # KCS12
  #geom_vline(data = filter(QTL.GWAS, Trait=="DF" | Trait=="DM" |Trait=="GYP" | Trait=="HE" |Trait=="NL" | Trait=="NT" |Trait=="PH" | Trait=="TSW"), aes(xintercept = 292412818), 
  #           lty=1, lwd=1, alpha=1, color="orange") +
  
  # LEA
  geom_vline(data = filter(QTL.GWAS, Trait=="HE" | Trait=="NL" | Trait=="NT" |Trait=="PH" | Trait=="TSW"), aes(xintercept = 295802230), 
             lty=1, lwd=1, alpha=1, color="orange") +
  
  
  #geom_vline(data = NULL, xintercept = 280595797, 
  #lty=1, lwd=1.5, alpha=0.4, color="orange") + # TOC1
  #annotate("text", x=23145069, y =6.5, label="TOC1") +
  
  #geom_vline(data = NULL, aes(xintercept = 154284883), 
  #lty=1, lwd=1.5, alpha=0.4, color="orange") + # SbGI
  #annotate("text", x=23145069, y =6.5, label="SbGI") +
  
  #geom_vline(data = NULL, aes(xintercept = 607890125), 
#lty=1, lwd=2, alpha=0.4, color="orange") + # SbCO
#annotate("text", x=23145069, y =6.5, label="SbCO") +

# Dw3
#geom_vline(data = filter(QTL.GWAS, Trait=="DF" | Trait=="GYP" | Trait=="NL" | Trait=="PH" | Trait=="TSW"), aes(xintercept = 476794611), 
#         lty=1, lwd=1, alpha=1, color="orange") +
#annotate("text", x=375515140, y =6.5, label="Dw3") +

# Dw2
#geom_vline(data = filter(QTL.GWAS, Trait=="DF" | Trait=="TSW"), aes(xintercept = 398790736), 
#           lty=1, lwd=1, alpha=1, color="orange") +
#annotate("text", x=375515140, y =6.5, label="Dw3") +


labs(x="", y= "") +
  geom_point(QTL.GWAS, mapping = aes(x=Peak_BPcum, y=Analysis), color="dodgerblue3", shape=15, size=1.2) + 
  theme_bw(base_size = 12)+ 
  theme(
    plot.margin = margin(0.5,0.1,0,0.1,"in"),
    #strip.background = element_rect(fill = "cadetblue1"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 9)) +
  guides(shape=guide_legend(override.aes = list(size=2)), color=guide_legend(override.aes = list(size=2)))

p2 <- ggplot(don1, aes(x=BPcum)) +
  facet_wrap(~Stat, ncol = 1, strip.position = "right", scales = "free_y") +
  geom_point(Fst, mapping = aes(x=BPcum, y=FST), shape=17, size=1.2) +
  geom_vline(stay.green, mapping = aes(xintercept = BPcum), lty=1, color="yellowgreen", alpha=0.3) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=1, color="black") +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  labs(x="Chromosome", y= "Fst") +
  theme_bw(base_size = 12)+ 
  theme(
    plot.margin = margin(0,0.1,0,0.1,"in"),
    #strip.background = element_rect(fill = "cadetblue1"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(margin = margin(0,-40,0,0)))

ggarrange(p1, p2, ncol = 1, heights = c(9,1), align = "v", labels = c("A","B"))
dev.off()









tiff("Marker-Trait Associations.tiff", units = "in", width = 12.945, height = 7.14, res = 300)

p1 <- ggplot(don1, aes(x=BPcum)) +
  facet_wrap(~Trait, ncol = 1, strip.position = "right", scales = "free_y") +
  geom_vline(stay.green, mapping = aes(xintercept = BPcum), lty=1, color="yellowgreen", alpha=0.3) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=1, color="black") +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  
  #Ma6
  geom_vline(data = QTL.GWAS, aes(xintercept = 356714466), 
             lty=1, lwd=1, alpha=1, color="orange") +
  #annotate("text", x=356714466, y = "Sheraro_JL", label="Ma6") +
  
  # Ma5
  #geom_vline(data = QTL.GWAS, aes(xintercept = 6745069), 
  #           lty=1, lwd=1, alpha=1, color="orange") +
  #annotate("text", x=23145069, y =6.5, label="Ma5") +
  
  # NAT/CAX7
  geom_vline(data = QTL.GWAS, aes(xintercept = 11158453), 
             lty=1, lwd=1, alpha=1, color="orange") +
  
  # Ma3
  geom_vline(data = QTL.GWAS, aes(xintercept = 60852774), 
             lty=1, lwd=1, alpha=1, color="orange") +
  #annotate("text", x=23145069, y =6.5, label="Ma3") +
  
  # LHY
  geom_vline(data = QTL.GWAS, aes(xintercept = 422952336), 
             lty=1, lwd=1, alpha=1, color="orange") + # LHY
  #annotate("text", x=23145069, y =6.5, label="LHY") +
  
  # P5CS2
  geom_vline(data = QTL.GWAS, aes(xintercept = 218863379), 
             lty=1, lwd=1, alpha=1, color="orange") + # LHY
  #annotate("text", x=23145069, y =6.5, label="P5CS2") +
  
  # SbCN8
  geom_vline(data = QTL.GWAS, aes(xintercept = 592828477), 
             lty=1, lwd=1, alpha=1, color="orange") + # LHY
  #annotate("text", x=23145069, y =6.5, label="SbCN8") +
  
  # DREB1A
  geom_vline(data = QTL.GWAS, aes(xintercept = 478434065), 
             lty=1, lwd=1, alpha=1, color="orange") +
  
  # KCS12
  #geom_vline(data = filter(QTL.GWAS, Trait=="DF" | Trait=="DM" |Trait=="GYP" | Trait=="HE" |Trait=="NL" | Trait=="NT" |Trait=="PH" | Trait=="TSW"), aes(xintercept = 292412818), 
  #           lty=1, lwd=1, alpha=1, color="orange") +
  
  # LEA
  geom_vline(data = QTL.GWAS, aes(xintercept = 295802230), 
             lty=1, lwd=1, alpha=1, color="orange") +
  
  
  
  #geom_vline(data = NULL, xintercept = 280595797, 
  #lty=1, lwd=1.5, alpha=0.4, color="orange") + # TOC1
  #annotate("text", x=23145069, y =6.5, label="TOC1") +
  
  #geom_vline(data = NULL, aes(xintercept = 154284883), 
  #lty=1, lwd=1.5, alpha=0.4, color="orange") + # SbGI
  #annotate("text", x=23145069, y =6.5, label="SbGI") +
  
#geom_vline(data = NULL, aes(xintercept = 607890125), 
#lty=1, lwd=2, alpha=0.4, color="orange") + # SbCO
#annotate("text", x=23145069, y =6.5, label="SbCO") +

# Dw3
#geom_vline(data = filter(QTL.GWAS, Trait=="DF" | Trait=="GYP" | Trait=="NL" | Trait=="PH" | Trait=="TSW"), aes(xintercept = 476794611), 
#         lty=1, lwd=1, alpha=1, color="orange") +
#annotate("text", x=375515140, y =6.5, label="Dw3") +

# Dw2
#geom_vline(data = filter(QTL.GWAS, Trait=="DF" | Trait=="TSW"), aes(xintercept = 398790736), 
#           lty=1, lwd=1, alpha=1, color="orange") +
#annotate("text", x=375515140, y =6.5, label="Dw3") +


labs(x="", y= "") +
  geom_point(QTL.GWAS, mapping = aes(x=Peak_BPcum, y=Analysis), color="dodgerblue3", shape=15, size=1.2) + 
  theme_bw(base_size = 12)+ 
  theme(
    plot.margin = margin(0.3,0.1,0,0.1,"in"),
    #strip.background = element_rect(fill = "cadetblue1"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 9)) +
  guides(shape=guide_legend(override.aes = list(size=2)), color=guide_legend(override.aes = list(size=2)))

p2 <- ggplot(don1, aes(x=BPcum)) +
  facet_wrap(~Stat, ncol = 1, strip.position = "right", scales = "free_y") +
  geom_vline(stay.green, mapping = aes(xintercept = BPcum), lty=1, color="yellowgreen", alpha=0.3) +
  
  # NAT/CAX7
  geom_vline(data = NULL, aes(xintercept = 11158453), lty=1, lwd=1, alpha=1, color="orange") +
  # P5CS2
  geom_vline(data = NULL, aes(xintercept = 218863379), lty=1, lwd=1, alpha=1, color="orange") + # LHY
  # DREB1A
  geom_vline(data = NULL, aes(xintercept = 478434065), lty=1, lwd=1, alpha=1, color="orange") +
  # LEA
  geom_vline(data = NULL, aes(xintercept = 295802230), lty=1, lwd=1, alpha=1, color="orange") +
  
  geom_point(Fst, mapping = aes(x=BPcum, y=FST), shape=17, size=1.2, color = "dodgerblue3") +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=1, color="black") +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  labs(x="Chromosome", y= "Fst") +
  theme_bw(base_size = 12)+ 
  theme(
    plot.margin = margin(0,0.1,0,0.1,"in"),
    #strip.background = element_rect(fill = "cadetblue1"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(margin = margin(0,-40,0,0)))

ggarrange(p1, p2, ncol = 1, heights = c(7,1), align = "v", labels = c("A","B"))
dev.off()








