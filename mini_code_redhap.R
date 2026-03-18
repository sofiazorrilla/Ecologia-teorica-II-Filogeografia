####
# Script: Comandos para ordenar la red de haplotipos en una consola de R
# Author: Sofía Zorrilla
# Fecha: 11-03-2026
####

# Este script asume que en la práctica de filogeogrfia_libelulas.qmd ya estás en la sección de Red de haplotipos y ya la pudiste correr. Ahora lo que vamos a hacer es organizar las coordenadas de los haplotipos para que se vean mejor las relaciones.

library(pegas)
library(dplyr)
library(tidyr)

coiAlign_trim <- read.dna("hetaerina_lab/02_alignment/alignment_coi_trimmed.fasta", format = "fasta")

coiHaplo <- haplotype(coiAlign_trim)

metadata_COI <- read.csv("hetaerina_lab/01_metadata/metadata_COI.csv")

red <- haploNet(coiHaplo)

freq_hap <- summary(coiHaplo)

# ¿Qué haplotipo tiene cada individuo?
tmp <- stack(setNames(attr(coiHaplo, "index"), rownames(coiHaplo))) # recuperar qué secuencia pertenece a qué haplotipo
ind_hap <- tapply(as.integer(tmp$ind), tmp$values, identity) # Agrupar los haplotipos usando como índice los indiviuos, es decir, generar un objeto en el que podamos saber qué haplotipo tiene cada individuo.

# Asignar colores a cada pais
countries <- metadata_COI$Country # Vamos a colorear los haplotipos por pais
countries <- factor(countries, levels = unique(countries))
n_countries <- length(unique(countries)) # Numero de localidades
paleta <- RColorBrewer::brewer.pal(min(n_countries, 8), "Set1") # Asignar un color a cada pais
colores_countries <- setNames(
  paleta[as.integer(as.factor(unique(n_countries)))],
  unique(n_countries)
)

# Matriz de proporción por haplotipo
prop_matriz <- ind_hap |>
  lapply(function(idx) table(countries[idx]) / length(idx)) |>
  bind_rows() |>
  as.matrix()

coords_manual <- readRDS("hetaerina_lab/03_outputs/coords_red_haplotipos.rds")

# Graficar
plot(
  red,
  size = freq_hap * 2, # tamaño de los haplotipos proporcional a su frecuencia
  pie = prop_matriz, # Asignar los colores por pais (usando graficas de pie)
  bg = paleta,
  xy = coords_manual,
  legend = c(100, -10),
  show.mutation = 1, # mostrar número de mutaciones en las ramas
  main = "Red de haplotipos COI – Libélula",
  threshold = 0
)

# Capturar coordenadas actuales
coords_red <- replot()

coords_red <- saveRDS(coords_red, "hetaerina_lab/03_outputs/coords_red_haplotipos_modificadas.rds")