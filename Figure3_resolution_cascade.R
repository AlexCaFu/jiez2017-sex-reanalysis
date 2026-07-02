#!/usr/bin/env Rscript
###############################################################################
# Figure 3  -  Resolution cascade: a compositional sex signal that does not
# localize (JieZ_2017 atherosclerosis cohort).  Publication vector, v2.
# Palette ties to Figure 2: BLUE = sex effect present; RED = absent.
# Left axis carries the coarse->fine direction (no in-figure sentence).
# Run: Rscript Figure3_resolution_cascade.R
###############################################################################
suppressMessages(library(ggplot2))

BLUE <- "#378ADD"; BLUE_D <- "#0C447C"; RED <- "#F7C1C1"
INK_W <- "#FFFFFF"; INK_R <- "#791F1F"
B_BLUE <- "#185FA5"; B_RED <- "#CE7C7C"; AXIS <- "#8A897F"; ARR <- "#6B6A64"

box <- function(x, ymid, w, h, fill, border) data.frame(
  xmin=x-w/2, xmax=x+w/2, ymin=ymid-h/2, ymax=ymid+h/2, fill=fill, border=border)
B <- rbind(
  box(0, 5.10, 9.4, 0.78, BLUE,   B_BLUE),
  box(0, 3.85, 6.6, 0.78, RED,    B_RED),
  box(0, 2.60, 5.0, 0.78, RED,    B_RED),
  box(0, 1.30, 9.9, 0.92, BLUE_D, BLUE_D))
arrows <- data.frame(x=0, xend=0,
  y=c(5.10-0.39, 3.85-0.39, 2.60-0.39), yend=c(3.85+0.39, 2.60+0.39, 1.30+0.46))

p <- ggplot() +
  # left resolution axis (replaces the old top sentence)
  geom_segment(aes(x=-5.6, y=5.55, xend=-5.6, yend=2.25), colour=AXIS, linewidth=0.6,
               arrow=arrow(length=unit(0.16,"cm"), type="closed")) +
  annotate("text", x=-5.6, y=5.75, label="coarse", size=2.5, colour=AXIS) +
  annotate("text", x=-5.6, y=2.05, label="fine",   size=2.5, colour=AXIS) +
  annotate("text", x=-6.15, y=3.85, label="mechanistic resolution", angle=90,
           size=2.9, colour="#5F5E5A") +
  # tier boxes
  geom_rect(data=B, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill=B$fill, colour=B$border, linewidth=0.5) +
  geom_segment(data=arrows, aes(x=x, y=y, xend=xend, yend=yend),
               colour=ARR, linewidth=0.6, arrow=arrow(length=unit(0.14,"cm"), type="closed")) +
  # tier 1
  annotate("text", x=0, y=5.24, label="Community composition", fontface=2, size=3.6, colour=INK_W) +
  annotate("text", x=0, y=4.98, label="sex detectable  \u00b7  p = 0.034  \u00b7  0.4% of variance", size=2.8, colour=INK_W) +
  # tier 2
  annotate("text", x=0, y=3.99, label="Individual taxa", fontface=2, size=3.6, colour=INK_R) +
  annotate("text", x=0, y=3.73, label="0 species survive FDR (lowest 0.35)", size=2.8, colour=INK_R) +
  # tier 3
  annotate("text", x=0, y=2.74, label="Metabolic pathways", fontface=2, size=3.6, colour=INK_R) +
  annotate("text", x=0, y=2.48, label="no TMA / bile / cholesterol pathway", size=2.8, colour=INK_R) +
  # conclusion
  annotate("text", x=0, y=1.44, label="A compositional sex difference that is real", fontface=2, size=3.5, colour=INK_W) +
  annotate("text", x=0, y=1.16, label="but does not localize to a co-metabolic node", fontface=2, size=3.5, colour=INK_W) +
  scale_x_continuous(limits=c(-6.5,5.3), expand=c(0,0)) +
  scale_y_continuous(limits=c(0.75,6.0), expand=c(0,0)) +
  theme_void()

ggsave("Figure3_resolution_cascade.pdf", p, width=7.6, height=4.7, device="pdf")
ptype <- if (isTRUE(capabilities("aqua"))) "quartz" else if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib"
try(ggsave("Figure3_resolution_cascade.png", p, width=7.6, height=4.7, dpi=600, device="png", type=ptype), silent=TRUE)
cat("Saved Figure 3 v2 (PDF + PNG).\n")
