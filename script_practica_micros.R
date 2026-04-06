
# Script : Práctica de filogeografía con microsatélites 
  # Author: Sofía Zorrilla
  # Date: 2026-04-06
  # Description: Análisis de diversidad y estructura genética con microsatélites
  # Arguments:
  #   - Input:
  #       - micros_hetaerina.csv
  #       - coordinates_micros.csv


# Preparación de datos ----------------------------------------------------


  ## Instalar paquetes ------------------------------------------------------
  
  # Lista de paquetes necesarios
  paquetes <- c(
    "adegenet",
    "pegas",
    "poppr",
    "hierfstat",
    "dplyr",
    "tibble",
    "reshape2",
    "broom",
    "ggplot2",
    "vegan",
    "PopGenReport",
    "tidyverse"
  )
  
  # Función para instalar paquetes faltantes
  instalar_si_falta <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
  
  # Instalar paquetes que no estén disponibles
  invisible(lapply(paquetes, instalar_si_falta))
  
  # Cargar todos los paquetes
  invisible(lapply(paquetes, library, character.only = TRUE))
  
  
  ## Cargar librerías --------------------------------------------------------
  
  # Análisis genético
  library(adegenet) # Estructura genética y análisis multivariado
  library(pegas) # Estadísticos poblacionales y haplotipos
  library(poppr) # Diversidad genética y clonación
  library(hierfstat) # Estadísticos F y diferenciación genética
  
  # Manipulación de datos
  library(dplyr) # Transformación de datos
  library(tibble) # Manejo moderno de data frames
  library(reshape2) # Reorganización de datos
  library(broom) # Conversión de resultados a formatos tabulares
  
  # Visualización y análisis multivariado
  library(ggplot2) # Gráficos
  library(vegan) # Análisis multivariado (ordenación, diversidad)
  library(PopGenReport) # Reportes automáticos de genética de poblaciones
  
  # Colección integrada
  library(tidyverse) # Conjunto de paquetes para ciencia de datos (dplyr, ggplot2, etc.)
  library(forcats)
  
  ## Organizacion de carpetas ------------------------------------------------
  
  # Crear estructura de carpetas del proyecto
  dirs <- c(
    "hetaerina_practica_micros/00_raw_data",
    "hetaerina_practica_micros/01_metadata",
    "hetaerina_practica_micros/03_outputs"
  )
  
  for (d in dirs) {
    if (!dir.exists(d))
      dir.create(d, recursive = TRUE)
  }
  
  fs::dir_tree("hetaerina_practica_micros")
  
  
  
  ## Descargar datos ---------------------------------------------------------
  
  micros <- read.csv(
    "https://raw.githubusercontent.com/sofiazorrilla/Ecologia-teorica-II-Filogeografia/refs/heads/main/hetaerina_practica_micros/00_raw_data/micros_hetaerina.csv"
  )
  
  write.csv(micros,
            "hetaerina_practica_micros/00_raw_data/micros_hetaerina.csv",
            row.names = F)
  
  coords <- read.csv(
    "https://raw.githubusercontent.com/sofiazorrilla/Ecologia-teorica-II-Filogeografia/refs/heads/main/hetaerina_practica_micros/01_metadata/coordinates_micros.csv"
  )
  
  write.csv(
    coords,
    "hetaerina_practica_micros/01_metadata/coordinates_micros.csv",
    row.names = F
  )
  
  ## Leer datos de microsatélites --------------------------------------------
  
  micros <- read.csv("hetaerina_practica_micros/00_raw_data/micros_hetaerina.csv") |>
    mutate(across(c(Poblacion, morfo, pops, Individuo), factor))
  
  rownames(micros) <- micros$Individuo
  
  as_tibble(micros)
  
  ## Convertir datos a genind ------------------------------------------------
  
  micros_genind <- adegenet::df2genind(
    X = micros[, -c(1:4)],
    sep = ":",
    NA.char = "NA",
    ind.names = micros$Individuo,
    ploidy = 2,
    type = "codom",
  )
  
  # Asignar metadatos de poblaciones
  strata(micros_genind) <- micros[, 1:4]
  setPop(micros_genind) <- ~ pops
  
  micros_genind
  
  # Resumen de los datos
  summary(micros_genind)
  
  # Espacios de memoria que tiene el objeto genind
  str(micros_genind, 2)
  
  #Explorar la tabla de genotipos
  head(micros_genind@tab)
  

#  Análisis ---------------------------------------------------------------

  
  ## Equilibrio H-W ----------------------------------------------------------
  
  round(pegas::hw.test(micros_genind, B = 1000), digits = 3)
  
  HWE.test <- data.frame(sapply(seppop(micros_genind), function(ls)
    pegas::hw.test(ls, B = 1000)[, 3]))
  
  HWE.test.chisq <- t(data.matrix(HWE.test))
  
  {
    cat("Chi-squared test (p-values):", "\n")
    round(HWE.test.chisq, 3)
  }
  
  ## ¿Qué proporción de loci están fuera de equilibrio
  
  # Define el nivel de significancia. Un locus se considera fuera de HWE si p < α.
  alpha = 0.05
  
  # Aqui vamos a utilizar un truco:
  
  # Para cada fila revisamos cuántas celdas cumplen la condición HWE.test.chisq <
  # α. Esto genera un vector de valores lógicos (TRUE o FALSE). Si convertimos ese
  # vector a valores numéricos, TRUE = 1 y FALSE = 0.
  
  # Para obtener la proporción de loci que cumplen la condición, podríamos contar
  # los que sí la cumplen y dividirlos entre el número total de loci. Sin embargo,
  # esta operación es equivalente a sumar los valores TRUE (que ahora valen 1) y
  # calcular su promedio.
  
  Prop.pops.out.of.HWE <- data.frame(Chisq = apply(HWE.test.chisq < alpha, 1, mean))
  
  Prop.pops.out.of.HWE
  
  
  ## Riqueza y rarefacción ---------------------------------------------------
  
  summary <- summary(micros_genind)
  
  # Riqueza absoluta (sin rarefaccionar)
  summary$pop.n.all
  
  summary_rarefied <- PopGenReport::allel.rich(micros_genind, min.alleles = NULL)
  
  knitr::kable(data.frame(
    N = summary$n.by.pop,
    Na = summary$pop.n.all,
    Ar = round(summary_rarefied$sum.richness, digits = 2)
  ))
  
  
  ## Heterocigocidad esperada y observada ------------------------------------
  
  heterocigocity <- data.frame(
    Ho = round(hierfstat::Ho(micros_genind), 3),
    He = round(adegenet::Hs(micros_genind), 3)
  ) |>
    mutate(
      diferencia = ifelse(
        Ho > He,
        "Exceso",
        ifelse(Ho < He, "Déficit", "Iguales")
      )
    )
  
  heterocigocity |> kable()
  
  
  ## Coeficiente de endogamia ------------------------------------------------
  
  basic_stats_micros <- hierfstat::basic.stats(micros_genind, diploid = TRUE)
  
  Fis_micros <- apply(
    basic_stats_micros$Fis,
    MARGIN = 2,
    FUN = mean,
    na.rm = TRUE
  ) %>%
    round(digits = 3)
  
  Fis_micros
  
  
  ## Tabla de medidas de diversidad ------------------------------------------
  
  micros_diversity <- data.frame(
    N = summary$n.by.pop,
    Na = summary$pop.n.all,
    Ar = round(summary_rarefied$sum.richness, digits = 2),
    Ho = heterocigocity$Ho,
    He = heterocigocity$He,
    Fis = Fis_micros
  ) |>
    rownames_to_column(var = "Pop")
  
  micros_diversity |> kable()
  

# Estructura genética -----------------------------------------------------


  ## Estimar número de grupos -----------------------------------------------
  
  kclusters <- find.clusters(
    micros_genind,
    max.n.clust = 26,
    n.pca = 20,
    # n.clust = 4,
    # n.iter = 10
  )
  
  mat <- tab(micros_genind, NA.method = "mean")
  grp <- pop(micros_genind)
  xval <- xvalDapc(
    mat,
    grp,
    n.pca.max = 100,
    training.set = 0.9,
    result = "groupMean",
    center = TRUE,
    scale = FALSE,
    n.pca = NULL,
    n.rep = 90,
    xval.plot = TRUE
  )
  
  xval[2:6]
  
  
  ## DAPC -------------------------------------------------------------------
  
  set.seed(12454)
  
  kclusters <- find.clusters(
    micros_genind,
    max.n.clust = 26,
    n.pca = 15,
    n.clust = 3,
    n.iter = 10
  )
  
  dapc1 <- dapc(micros_genind, n.pca = 15, n.da = 2, kclusters$grp)
  
  spp.label = c(
    expression(italic("H. nov. sp")),
    expression(italic("H. calverti")),
    expression(italic("H. americana"))
  )
  
  colors <- c("red", "blue", "green")
  
  scatter(
    dapc1,
    scree.da = FALSE,
    bg = "white",
    col = colors,
    pch = 20,
    cstar = 0,
    cex = 4,
    clabel = 0,
    leg = TRUE,
    posi.leg = "topleft",
    txt.leg = spp.label,
    scree.pca = TRUE,
    posi.pca = "bottomleft"
  )
  
  # Barplot
  
  compoplot.dapc(
    dapc1,
    posi = "bottomright",
    leg = FALSE,
    txt.leg = spp.label,
    lab = "a",
    xlab = "Individuos",
    col = colors
  )
  
  
  ## Mapa --------------------------------------------------------------------
  
  info <- data.frame(
    ind = names(dapc1$grp),
    grp = as.factor(dapc1$grp),
    pop = as.factor(micros_genind$pop)
  )
  
  # ¿cuantos individuos de cada poblacion pertenecen a cada grupo?
  tabla_pertenencia <- table(info$grp, info$pop) |> t()
  
  colnames(tabla_pertenencia) <- make.names(c(
    "H.nov.sp",
    "H.calverti",
    "H.americana"
  ))
  
  # Leer archivo de coordenadas
  coords <- read.csv("hetaerina_lab/01_metadata/coordinates_micros.csv")
  
  frecuencias_por_pop <- tabla_pertenencia |>
    as.data.frame() |>
    rename(pop = "Var1", sp = "Var2") |>
    group_by(pop) |>
    mutate(Freq = Freq / sum(Freq)) |>
    ungroup() |>
    pivot_wider(id_cols = pop, names_from = sp, values_from = Freq) |>
    left_join(coords, by = c("pop" = "Code"))
  
  # Generar poligono del mundo
  mapa_mundo <- ggplot2::map_data("world")
  
  # Mapa
  ggplot() +
    geom_polygon(
      data = mapa_mundo,
      aes(x = long, y = lat, group = group),
      fill = "grey85",
      color = "white",
      linewidth = 0.2
    ) +
    geom_point(data = frecuencias_por_pop, aes(x = Lon, y = Lat)) +
    scatterpie::geom_scatterpie(
      data = frecuencias_por_pop,
      aes(x = Lon, y = Lat, r = 0.7),
      alpha = 0.7,
      cols = colnames(frecuencias_por_pop)[startsWith(
        colnames(frecuencias_por_pop),
        "H"
      )],
      long_format = F
    ) +
    scale_fill_manual(values = colors) +
    coord_fixed(xlim = c(-120, -80), ylim = c(10, 42)) +
    labs(
      title = "Distribución de grupos por localidad",
      x = "Longitud",
      y = "Latitud",
      fill = "Grupo"
    ) +
    theme_classic(base_size = 13)
  
  
  ## Diferenciación Fst ------------------------------------------------------
  
  fst_micros <- genet.dist(micros_genind, method = "WC84") %>% round(digits = 3)
  fst_micros
  
  
  ## Heatmap
  
  #Ordenamos los datos de la matriz resultante
  lab_order <- c(
    "Bo",
    "CC",
    "Cm",
    "Co",
    "Cp",
    "Ct",
    "EF",
    "LC",
    "Mn",
    "Mo",
    "NR",
    "Ny",
    "RCh",
    "SA",
    "SG",
    "SM",
    "SN",
    "Ve",
    "Vi",
    "Za"
  )
  
  df_fst <- tidy(as.dist(fst_micros))
  names(df_fst) <- c("pop1", "pop2", "fst")
  
  df_fst <- df_fst |>
    mutate(
      pop1 = factor(pop1, levels = lab_order),
      pop2 = factor(pop2, levels = lab_order),
      fst = ifelse(fst < 0, 0, fst)
    )
  
  
  #Etiqueta en itálica
  fst.label = expression(italic("F")[ST])
  
  #Extraemos valores de Fst medios para hacer nuestro gradiente
  mid = max(df_fst$fst) / 2
  
  #Plotteamos el Heatmap
  fst_micros_graf <- ggplot(
    data = df_fst,
    aes(x = pop1, y = fct_rev(pop2), fill = fst)
  ) +
    geom_tile(colour = "black") +
    geom_text(aes(label = fst), color = "black", size = 3) +
    scale_fill_gradient2(
      low = "blue",
      mid = "pink",
      high = "red",
      midpoint = mid,
      name = fst.label,
      limits = c(0, max(df_fst$fst)),
      breaks = c(0, 0.25, 0.5, 0.75, 1)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), position = "right") +
    theme(
      axis.text = element_text(colour = "black", size = 10, face = "bold"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  fst_micros_graf
  