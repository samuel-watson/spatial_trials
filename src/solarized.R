# SOLARISED COLOR PALETTE
library(extrafont)
#font_import() # only needs to be run once
loadfonts(device = "win")

solar_color <- c(
    base03 =   "#002b36",
    base02=    "#073642",
    base01=    "#586e75",
    base00=    "#657b83",
    base0=     "#839496",
    base1=     "#93a1a1",
    base2=     "#eee8d5",
    base3=     "#fdf6e3",
    yellow=    "#b58900",
    orange=    "#cb4b16",
    red=       "#dc322f",
    magenta=   "#d33682",
    violet=    "#6c71c4",
    blue=      "#268bd2",
    cyan=      "#2aa198",
    green=     "#859900"
)

theme_solar <- function(base_size = 11, base_family = "", base_line_size = base_size/22, 
                        base_rect_size = base_size/22){
  font <- "Tw Cen MT"
  theme_bw(base_size = base_size, base_family = base_family, 
             base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,colour = "#073642"), 
          panel.grid = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          strip.background = element_rect(fill = "white", colour = "#073642"), 
          legend.key = element_rect(fill = "white", colour = NA),
          plot.title = element_text(             #title
            family = font,            #set font family
            size = 12,                #set font size
            face = 'bold',            #bold typeface
            hjust = 0,                #left align
            vjust = 2),
          
          plot.subtitle = element_text(          #subtitle
            family = font,            #font family
            size = 14),               #font size
          
          plot.caption = element_text(           #caption
            family = font,            #font family
            size = 9,                 #font size
            hjust = 1),               #right align
          
          axis.title = element_text(             #axis titles
            family = font,            #font family
            size = 10),               #font size
          
          axis.text = element_text(              #axis text
            family = font,            #axis famuly
            size = 10),                #font size
          
          axis.text.x = element_text(            #margin for axis text
            margin=margin(5, b = 10)),
          
          complete = TRUE)
}
