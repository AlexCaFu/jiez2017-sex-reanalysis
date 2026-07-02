#!/usr/bin/env Rscript
###############################################################################
# Figure 2 - Source-resolution evidence map - HD 1200 dpi, stronger palette,
# ALL text bold. Arial (Liberation Sans, Arial-metric). Outputs PNG (then JPG).
###############################################################################
suppressMessages({library(ggplot2); library(systemfonts); library(ragg)})
LIB <- "/usr/share/fonts/truetype/liberation"
try(register_font(name="Arial",
  plain=file.path(LIB,"LiberationSans-Regular.ttf"),
  bold =file.path(LIB,"LiberationSans-Bold.ttf"),
  italic=file.path(LIB,"LiberationSans-Italic.ttf"),
  bolditalic=file.path(LIB,"LiberationSans-BoldItalic.ttf")), silent=TRUE)
FAM <- "Arial"

nodes <- c("Substrate","Microbial\nProduction","Host\nConversion",
           "Transport/\nUptake","Clearance","Vascular\nResponse")
paths <- c("TMAO","PAGln","Bile acids","Cholesterol\nConversion")
tier <- rbind(TMAO=c(1,3,3,0,1,0), PAGln=c(0,1,0,0,0,1),
              `Bile acids`=c(0,2,1,2,0,0), Cholesterol=c(0,1,0,0,0,0))
df <- expand.grid(node=factor(nodes,levels=nodes), path=factor(paths,levels=rev(paths)))
df$tier<-as.vector(t(tier)); df$glyph<-ifelse(df$tier==0,"\u00b7",as.character(df$tier))
df$tf<-factor(df$tier,levels=c(0,1,2,3))

# STRONGER palette (deeper/more saturated); black bold text stays readable
FILL <- c("0"="#E59696","1"="#A6C8EF","2"="#5A97D6","3"="#2E6EBE")

p <- ggplot(df, aes(node,path)) +
  geom_tile(aes(fill=tf), colour=NA, width=.9, height=.9) +
  geom_text(aes(label=glyph), colour="black", family=FAM, fontface="bold", size=6.8) +
  scale_fill_manual(values=FILL, drop=FALSE, name=NULL,
    labels=c("\u00b7  no sex-resolved evidence identified",
             "1  mechanistic, not sex-resolved",
             "2  sex-resolved (animal)","3  sex-resolved (human)")) +
  scale_x_discrete(position="top") + coord_equal() +
  guides(fill=guide_legend(nrow=2, byrow=TRUE, override.aes=list(colour=NA))) +
  theme_minimal(base_size=12, base_family=FAM) +
  theme(
    axis.title=element_blank(),
    axis.text.x.top=element_text(colour="black", family=FAM, size=13, face="bold",
                                 lineheight=.9, vjust=0),
    axis.text.y=element_text(colour="black", family=FAM, size=15, face="bold",
                             hjust=1, lineheight=.9),
    panel.grid=element_blank(),
    legend.position="bottom",
    legend.text=element_text(colour="black", family=FAM, size=11, face="bold"),
    legend.key.size=unit(15,"pt"), legend.spacing.x=unit(6,"pt"),
    legend.margin=margin(t=10),
    plot.background=element_rect(fill="white",colour=NA),
    panel.background=element_rect(fill="white",colour=NA),
    plot.margin=margin(10,14,10,10))

ggsave("Figure2_evidence_map_HD1200.png", p, width=8.2, height=5.9, dpi=1200,
       device=ragg::agg_png, bg="white")
cat("Saved 1200-dpi PNG master\n")
