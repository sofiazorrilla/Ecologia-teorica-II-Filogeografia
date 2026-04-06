######
# Script : Práctica de filogeografía con microsatélites
# Author: Sofía Zorrilla
# Date: 2026-04-06
# Description: Análisis de diversidad y estructura genética con microsatélites
# Arguments:
#   - Input: 
#       - 
#   - Output: 
#######


# Instalar paquetes -------------------------------------------------------

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


# Cargar librerías --------------------------------------------------------

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

# Organizacion de carpetas ------------------------------------------------

# Crear estructura de carpetas del proyecto
dirs <- c(
  "hetaerina_practica_micros/00_raw_data",
  "hetaerina_practica_micros/01_metadata",
  "hetaerina_practica_micros/03_outputs"
)

for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

fs::dir_tree("hetaerina_practica_micros")

