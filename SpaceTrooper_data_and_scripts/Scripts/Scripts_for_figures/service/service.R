library(ggplot2)

### Define plot theme
my.theme <- theme(aspect.ratio = 1,
                  panel.border=element_blank(),
                  panel.background=element_rect(fill="transparent", color=NA),
                  plot.background=element_rect(fill="transparent", color=NA),
                  plot.title = element_blank(),
                  axis.line=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.text.y=element_blank(),
                  axis.text.x=element_blank(),
                  legend.background=element_rect(fill="transparent", color=NA),
                  legend.title= element_text(color = "black", size = 10),
                  legend.text=element_text(color="black", size = 6),
                  legend.position = "right",
                  legend.justification = c(0.95, 0.5),
                  legend.box.just = "right",
                  legend.margin = margin(t = -15)  ,
                  panel.grid=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.grid.major=element_blank())

# SpaceTrooper plotScaleBar

plotScaleBar <- function(p, spe, micronConvFact=0.12) {
  stopifnot(is(spe, "SpatialExperiment"))
  if (!inherits(p, "ggplot")) {
    stop("`p` must be a ggplot object (e.g. created by ggplot()).")
  }
  if (!is.numeric(micronConvFact) || length(micronConvFact) != 1 || micronConvFact <= 0) {
    stop("`micronConvFact` must be a single positive numeric value.")
  }
  x_vals <- spatialCoords(spe)[,1]
  y_vals <- spatialCoords(spe)[,2]
  x_range <- range(x_vals, na.rm = TRUE)
  x_length <- diff(x_range)
  if(metadata(spe)$technology %in% c("Nanostring_CosMx",
                                     "Nanostring_CosMx_Protein")){
    target <- x_length / 4 * micronConvFact
  } else {
    target <- x_length / 4
  }
  order <- 10^floor(log10(target))
  # guard rails: if target is NA or non-positive, fallback to 1*10^0
  if (!is.finite(order) || order <= 0) order <- 1
  mantissa <- target / order
  nice_mantissa <- if (mantissa < 1.5) {
    1
  } else if (mantissa < 3.5) {
    2
  } else if (mantissa < 7.5) {
    5
  } else {
    10
  }
  scale_length <- nice_mantissa * order
  scale_label <- if (scale_length >= 1000) {
    c(as.character(scale_length/2000), paste0(scale_length / 2000, " mm"))
  } else {
    micro <- "\u00B5"
    c(as.character(scale_length/2), paste0(scale_length/2, " ", micro, "m"))
  }
  if(metadata(spe)$technology %in% c("Nanostring_CosMx",
                                     "Nanostring_CosMx_Protein")){
    scale_length <- scale_length / micronConvFact
  }
  # Set coordinates for the scale bar
  x_start <- x_range[2] - scale_length * 1.1
  x_end <- x_range[2] - scale_length * 0.1
  y_pos <- min(y_vals, na.rm = TRUE) + scale_length * 0.1
  box_data <- data.frame(
    xmin = x_start - scale_length * 0.05,
    xmax = x_end + scale_length * 0.05,
    ymin = y_pos - scale_length * 0.05,
    ymax = y_pos + scale_length * 0.3
  )
  
  scale_bar_data <- data.frame(
    xmin = c(x_start, x_start + (scale_length/2)),
    xmax = c(x_start + (scale_length/2), x_start + (scale_length)),
    ymin = c(y_pos, y_pos),
    ymax = c(y_pos + (scale_length * 0.05), y_pos + (scale_length * 0.05))
  )
  
  label_data <- data.frame(
    x = c(x_start + (scale_length/4), x_start + 0.7*scale_length),
    y = c(y_pos + (scale_length * 0.2), y_pos + (scale_length * 0.2)),
    label = scale_label)
  
  p <- p+geom_rect(data = box_data, aes(xmin = xmin, xmax = xmax,
                                        ymin = ymin, ymax = ymax),
                   fill = "white", color = "grey", alpha = 0.8,
                   inherit.aes = FALSE) +
    geom_rect(data = scale_bar_data, aes(xmin = xmin, xmax = xmax,
                                         ymin = ymin, ymax = ymax),
              fill = c("black", "white"), color = "grey",
              inherit.aes = FALSE) +
    geom_text(data = label_data,
              aes(x = x, y = y, label = label), color = "grey40",
              inherit.aes = FALSE, size = 3, hjust = 0.2)
  return(p)
}

