suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(showtext)
})

font_add("TeX Gyre Pagella",
         regular    = "/home/miguel/.fonts/texgyrepagella-regular.otf",
         bold       = "/home/miguel/.fonts/texgyrepagella-bold.otf",
         italic     = "/home/miguel/.fonts/texgyrepagella-italic.otf",
         bolditalic = "/home/miguel/.fonts/texgyrepagella-bolditalic.otf")
showtext_auto()
showtext_opts(dpi = 300)
FONT <- "TeX Gyre Pagella"

# RUTAS
base    <- "/home/miguel/Abaum_Resistome_Network"
tables  <- file.path(base, "results/tables")
out_dir <- file.path(base, "results/figures/main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

trends  <- fread(file.path(tables, "temporal_trends.csv"))
prev_yr <- fread(file.path(tables, "temporal_prevalence_by_year.csv"))

# FILTRO TEMPORAL
prev_yr <- prev_yr[year >= 2006 & year <= 2023]

# GENES SIGNIFICATIVOS p < 0.05
sig     <- trends[pvalue < 0.05]
sig_inc <- sig[trend == "increasing"][order(-slope), gene]
sig_dec <- sig[trend == "decreasing"][order(slope),  gene]

gene_order <- c(sig_dec, sig_inc)   # dec=bottom, inc=top
gcols      <- intersect(gene_order, colnames(prev_yr))
gcols_inc  <- intersect(sig_inc, gcols)
gcols_dec  <- intersect(sig_dec, gcols)
n_inc <- length(gcols_inc)
n_dec <- length(gcols_dec)
n_tot <- length(gcols)

# PREVALENCIA MEDIA POR GEN
prev_mean <- prev_yr[, lapply(.SD, mean, na.rm = TRUE),
                     .SDcols = gcols]
prev_mean_long <- data.table(
  gene     = gcols,
  mean_prev = as.numeric(prev_mean[1, ])
)

# GENES ROBUSTOS: R² ≥ 0.60 AND prevalencia media ≥ 10%
robust_genes <- trends[pvalue < 0.05 & r2 >= 0.60, gene]
robust_genes <- intersect(robust_genes, gcols)
robust_in    <- prev_mean_long[gene %in% robust_genes &
                                 mean_prev >= 0.10, gene]

# CLASE DE RESISTENCIA
res_class <- c(
  "msr(e)"        = "Macrolide",
  "arma"          = "Aminoglycoside",
  "mph(e)"        = "Macrolide",
  "ftsi"          = "Beta-lactam",
  "blaoxa_23like" = "Carbapenem",
  "tet(b)"        = "Tetracycline",
  "blandm"        = "Carbapenem",
  "blaoxa_24like" = "Carbapenem",
  "arr-2"         = "Rifampicin",
  "pmrb"          = "Colistin",
  "aac(3)"        = "Aminoglycoside",
  "ant(2'')"      = "Aminoglycoside",
  "qacedelta"     = "Disinfectant",
  "sul"           = "Sulfonamide",
  "aph(3')"       = "Aminoglycoside",
  "aada"          = "Aminoglycoside",
  "adec"          = "Efflux",
  "blatem"        = "Beta-lactam",
  "mert"          = "Metal",
  "merr"          = "Metal",
  "blaoxa_58like" = "Carbapenem"
)

# FORMATO LARGO
prev_long <- melt(
  prev_yr[, c("year", gcols), with = FALSE],
  id.vars       = "year",
  variable.name = "gene",
  value.name    = "prev"
)
prev_long[, gene    := factor(gene, levels = gene_order)]
prev_long[, pct     := round(prev * 100, 0)]
prev_long[, lab     := ifelse(prev > 0, paste0(pct, "%"), "")]
prev_long[, txt_col := ifelse(prev >= 0.52, "white", "#1a1a1a")]

# Etiquetas eje Y: ★ solo en genes con R² ≥ 0.60 y prevalencia ≥ 10%
gene_labels <- setNames(
  ifelse(gene_order %in% robust_in,
         paste0(gene_order, "  *"),
         gene_order),
  gene_order
)

# POSICIONES
y_dec_center <- n_dec / 2 + 0.5
y_inc_center <- n_dec + n_inc / 2 + 0.5
years_available <- sort(unique(prev_long$year))
x_min <- min(years_available)
x_max <- max(years_available)

# PALETA
pal_colors <- c(
  "#FFFEF5", "#FFF7BC", "#FEC44F",
  "#FE9929", "#EC7014", "#CC4C02", "#8C2D04"
)
pal_values <- rescale(c(0, 0.10, 0.30, 0.50, 0.65, 0.80, 1.0))

# HEATMAP
p_heat <- ggplot(prev_long, aes(x = year, y = gene, fill = prev)) +

  # Franja Increasing
  annotate("rect",
           xmin = x_min - 1.62, xmax = x_min - 0.92,
           ymin = n_dec + 0.5, ymax = n_tot + 0.5,
           fill = "#CC4C02", color = NA, alpha = 0.90) +
  annotate("text",
           x = x_min - 1.27, y = y_inc_center,
           label = "Increasing \u2191",
           angle = 90, hjust = 0.5, size = 6.2,
           fontface = "bold", family = FONT, color = "white") +

  # Franja Decreasing
  annotate("rect",
           xmin = x_min - 1.62, xmax = x_min - 0.92,
           ymin = 0.5, ymax = n_dec + 0.5,
           fill = "#1B4F72", color = NA, alpha = 0.90) +
  annotate("text",
           x = x_min - 1.27, y = y_dec_center,
           label = "Decreasing \u2193",
           angle = 90, hjust = 0.5, size = 6.2,
           fontface = "bold", family = FONT, color = "white") +

  # Separador horizontal (solo en el área del heatmap, no sobre la barra lateral)
  annotate("rect",
           xmin = x_min - 0.5, xmax = x_max + 0.52,
           ymin = n_dec + 0.485, ymax = n_dec + 0.515,
           fill = "#222222", color = NA) +

  # Celdas
  geom_tile(color = "white", linewidth = 0.28) +

  # Porcentajes
  geom_text(aes(label = lab, color = txt_col),
            size = 6.2, family = FONT, show.legend = FALSE) +
  scale_color_identity() +

  scale_fill_gradientn(
    colors = pal_colors, values = pal_values,
    limits = c(0, 1), name = "Prevalence",
    labels = percent_format(accuracy = 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.0),
    guide  = guide_colorbar(
      barheight = unit(5.0, "cm"), barwidth = unit(0.42, "cm"),
      ticks.colour = "grey50", frame.colour = "grey60",
      title.position = "top", title.hjust = 0.5
    )
  ) +

  scale_x_continuous(
    breaks = years_available,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_discrete(
    labels = gene_labels,
    expand = expansion(add = c(0.4, 0.4))
  ) +

  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 22, base_family = FONT) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1,
                                    size = 18.0, color = "#333333", family = FONT),
    axis.text.y      = element_text(size = 18.0, color = "#111111",
                                    face = "italic", family = FONT),
    axis.title.x     = element_text(size = 17.5, family = FONT, margin = margin(t = 6)),
    axis.title.y     = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 19.0, face = "bold", family = FONT),
    legend.text      = element_text(size = 17.5, family = FONT),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(t = 10, r = 4, b = 10, l = 55)
  ) +
  labs(x = "Collection year", y = NULL)

# EXPORTAR
# EXPORTAR PNG + TIFF
png_path  <- file.path(out_dir, "Fig2_temporal_trend_heatmap.png")
tiff_path <- file.path(out_dir, "Fig2_temporal_trend_heatmap.tiff")
ggsave(png_path,  p_heat, width = 14, height = 10, dpi = 300, bg = "white")
ggsave(tiff_path, p_heat, width = 14, height = 10, dpi = 300, bg = "white",
       device = "tiff", compression = "lzw")
cat("\u2713 PNG :", png_path,  "\n")
cat("\u2713 TIFF:", tiff_path, "\n")
