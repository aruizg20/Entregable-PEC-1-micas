# Se cargan las librerías necesarias
library(SummarizedExperiment)
library(S4Vectors)  # Para la clase DataFrame

# Se cargan los datos del estudio desde los archivos, usando la función de R File/Import Dataset.


# Se verifican las dimensiones de los datos cargados
cat("Dimensiones de features:", dim(features), "\n")  # Se esepra 1541 x 45
cat("Dimensiones de metadata:", dim(metadata), "\n")  # Se esepra 45 x 2
cat("Dimensiones de metaboliteNames:", dim(metaboliteNames), "\n")  # Se esepra 1541 x 3


# Verifica las dimensiones antes de crear el objeto
if (ncol(features) != nrow(metadata)) {
  stop("El número de columnas en el archivo features y el número de filas en el archivo metadata no coincide.")
}

## Se comparan los nombres de las columnas para detectar porsibles errores

# Verifica los nombres de las muestras en features
cat("Nombres de muestras en features:", colnames(features), "\n")

# Verifica los nombres de las muestras en metadata
cat("Nombres de muestras en metadata:", metadata$ID, "\n")  


#Se modifica el nombre de las columnas de Metadata
rownames(metadata)<-metadata$ID

library(SummarizedExperiment)

# Se crea el objeto SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = as.matrix(features)),
  colData = DataFrame(metadata)
)


# Se verifica que el objeto se haya creado correctamente
print(se)


##Exploración

# Verificación del número de variables (metabolitos) y muestras (sujetos)
dim(se)

# Revisión de las primeras filas del conteo de metabolitos
head(assay(se, "counts"))

# Se exploran los metadatos de las muestras
colData(se)

# Resumen estadístico de los valores de metabolitos
summary(assay(se, "counts"))

# Boxplot para cada muestra
boxplot(assay(se, "counts"), main = "Distribución de Concentraciones por Muestra", xlab = "Muestras", ylab = " ", las = 2, outline = FALSE)

  
install.packages("pheatmap")

# Se carga la biblioteca
library(pheatmap)

# Extracción la matriz de datos de metabolitos
metabolite_data <- assay(se, "counts")

# Aplicación de transformación logarítmica para estabilizar la varianza
metabolite_data_log <- log2(metabolite_data + 1)

# Se eliminan las filas con NA
metabolite_data_log_clean <- na.omit(metabolite_data_log)


# Crear un heatmap usando pheatmap
pheatmap(metabolite_data_log, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         main = "Heatmap de Metabolitos",
         scale = "row",  # Escala cada metabolito por filas
         show_rownames = FALSE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50))  # Define los colores



# Crea un factor para los nombres de las muestras en el orden deseado
sample_names <- colnames(metabolite_data_log)
sample_order <- factor(sample_names, levels = c(paste0("a", 1:17), paste0("b", 1:17), paste0("c", 1:17)))

# Reordenar la matriz usando el factor
metabolite_data_log_ordered <- metabolite_data_log[, order(sample_order)]






