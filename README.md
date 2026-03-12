# Ecología Teórica II
Repositorio que contiene las prácticas de filogeografía para las clases de Ecología Teórica II de la licenciatura de Ecología en la Escuela Nacional de Estudios Superiores Unidad Morelia

## Instrucciones previas a la práctica

1. Instalar TODOS los paquetes mencionados a continuación 

**Nota:** Algunos tardan bastante tiempo en instalar, háganlo con suficiente anticipación. Si tienen problemas pueden venir a verme en el Laboratorio de Biología Neotropical para buscar una solución

2. Generar una [cuenta en GenBank](https://www.ncbi.nlm.nih.gov/myncbi/) 

## Instalación de paquetes

```r
# Instalar paquetes desde CRAN
install.packages(c(
  "bold",        
  "rentrez",     
  "ape",         
  "seqinr",      
  "stringr",     
  "dplyr",       
  "tidyr",       
  "ggplot2",     
  "knitr",       
  "kableExtra",  
  "here",        
  "fs",          
  "pegas",       
  "adegenet",    
  "poppr",       
  "RColorBrewer",
  "maps",        
  "scatterpie",  
  "gtools"       
))

# Instalar paquetes desde Bioconductor
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "Biostrings",  
  "msa"          
))
```
