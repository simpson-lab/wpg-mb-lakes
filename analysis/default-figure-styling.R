library('ggplot2')   # for plotting
library('cowplot')   # for plotting
library('extrafont') # for plot fonts
loadfonts(device = 'win', quiet = TRUE)
theme_set(theme_bw(base_family = 'serif') +
            theme(panel.grid = element_blank()))

# using Limnology and Oceanography requirements:
LABELS <- paste0(letters, '.') # top left, outside the figure
DPI <- 1200
WIDTH <- 3.5 # inches
HEIGHT <- WIDTH / 4 * 3 # 4:3 format

p2pdf <- function(filename, p, width = WIDTH, height = HEIGHT, scale = 1, x.plots = 1,
                  y.plots = 1) {
  if(width > 5 | height > 6)
    warning('Figure sizes should be no more than 5" in width and 6" in height')
  if((x.plots + y.plots != 2) & (width != WIDTH | height != HEIGHT))
    warning('Non-standard scaling! (Figure sizes & number of plots have been modified.)')
  
  ggsave(filename = filename,
         plot = p,
         device = cairo_pdf,
         path = 'figures/',
         scale = scale,
         width = width * x.plots,
         height = height * y.plots,
         units = 'in',
         dpi = DPI,
         limitsize = TRUE) # prevent saving images > 50x50 inches
}

# Most illustrations, except some maps and very wide graphs, should be 1-column size
# (3.5 inches) and a resolution of 300 dpi. The font size on the x and y axes should not
# be larger than that of the title, and the same font (Arial or Times New Roman is
# preferred) should be used throughout. Numbers on the x and y axes should be smaller
# than the descriptive title, which should be 12-point font. Fonts smaller than 12 points
# are generally not legible when reduced to 1 column size. Use boldface type with care;
# if illustrations are to be reduced, the letters with open spaces will disappear.
# Use sentence case (capitalize the first word ONLY) for axis titles, labels, and legends.

# Symbols and Lines: Avoid very small symbols (no smaller than 2 mm) on line graphs;
# print all elements of the graph with the same degree of intensity. In addition to the
# above guidelines, color figures must be submitted in the CMYK colorspace.

