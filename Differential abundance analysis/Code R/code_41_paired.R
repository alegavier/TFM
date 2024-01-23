library(dplyr)
library(stats)
library(MASS)
library(ggplot2)

# Cargo los datos ya transformados y con los "missing values" imputados.
table <- read.csv(file = 'ion_abundance_41.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)
reference_channel <- read.csv(file = 'reference_intensity.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Ordenar las muestras para que se analicen pareadas.
orden_personalizado <- c('CPT0063570003', 'CPT0025750003', 'CPT0237840004', 'CPT0124980004', 'CPT0127270003', 'CPT0011520003', 'CPT0183450005', 'CPT0077860003', 'CPT0198420004', 'CPT0203850004', 'CPT0236700005', 'CPT0208690004', 'CPT0170240004', 'CPT0237910004', 'CPT0162760003', 'CPT0091660003', 'CPT0198910005', 'CPT0001730017', 'CPT0162660004', 'CPT0238610009', 'CPT0064090003', 'CPT0183750005', 'CPT0218200004', 'CPT0208840004', 'CPT0239910004', 'CPT0238920004', 'CPT0221240004', 'CPT0166760004', 'CPT0246990003', 'CPT0123650003', 'CPT0241500004', 'CPT0197420004', 'CPT0086190004', 'CPT0236780005', 'CPT0226600004', 'CPT0174700004', 'CPT0248410004', 'CPT0246940003', 'CPT0078000003', 'CPT0108930003', 'CPT0248590004',
                       'CPT0063580003', 'CPT0025760003', 'CPT0237850004', 'CPT0125010003', 'CPT0127290003', 'CPT0011540003', 'C3N-02573-03', 'CPT0077880003', 'CPT0198440004', 'C3N-00709-06', 'CPT0236710004', 'CPT0208710004', 'CPT0170260004', 'CPT0237920004', 'CPT0162780004', 'CPT0091680003', 'CPT0198930003', 'CPT0001750003', 'CPT0162670004', 'CPT0238640004', 'CPT0064110003', 'CPT0183770004', 'CPT0218220004', 'CPT0208850004', 'CPT0239920004', 'CPT0238940004', 'CPT0221250004', 'CPT0166780004', 'CPT0247000003', 'CPT0123670003', 'CPT0241510004', 'CPT0197440004', 'CPT0086210003', 'CPT0236800004', 'CPT0226610004', 'CPT0174800003', 'CPT0248440003', 'CPT0246950003', 'CPT0078030003', 'CPT0108960003', 'CPT0248620003')

# Crear un nuevo dataframe con el orden personalizado
orden_personalizado_df <- data.frame(Sample.ID = orden_personalizado)

# Convertir la columna Sample.ID a factor con el orden personalizado
table$Sample.ID <- factor(table$Sample.ID, levels = orden_personalizado_df$Sample.ID)

# Ordenar el dataframe según la columna Sample.ID
order <- table[order(table$Sample.ID), ]


#-------------------------- LOG 2 FOLD CHANGE --------------------------------------

# Filtrar las muestras por condición (Tumor y Normal)
tumor_samples <- order$Class == "Tumor"
normal_samples <- order$Class == "NAT"

# Seleccionar solo las columnas de proteínas (ignorando las columnas 'Sample.ID' y 'Class')
proteins_data <- order[tumor_samples | normal_samples, -(1:2)]

median_tumor <- apply(proteins_data[tumor_samples, ], 2, median)
median_normal <- apply(proteins_data[normal_samples, ], 2, median)

medianas_df <- data.frame(
  Proteina = colnames(proteins_data),
  MedianaTumor = median_tumor,
  MedianaNormal = median_normal
)

# Añadir el canal de referencia a la tabla de medianas
medianas_df <- merge(medianas_df, reference_channel, by.x = "Proteina", by.y = "Sample.ID")

# Añadir una nueva columna con la diferencia entre la mediana y la abundancia de referencia para cada condición
medianas_df$MedianaTumorReferencia <- medianas_df$MedianaTumor - medianas_df$ReferenceIntensity
medianas_df$MedianaNormalReferencia <- medianas_df$MedianaNormal - medianas_df$ReferenceIntensity

# Calcular la diferencia entre las medianas
medianas_df <- medianas_df %>%
  mutate(Log2FoldChange = MedianaTumorReferencia - MedianaNormalReferencia)


# ----------------------- WILCOXON SIGNED RANK TEST -----------------------------------------

data_for_wilcox <- order[tumor_samples | normal_samples, ]

# Reordenar el dataframe por Sample.ID para asegurarse que las muestras pareadas estén juntas
data_for_wilcox <- data_for_wilcox[order(data_for_wilcox$Sample.ID), ]

# Filtrar las columnas de abundancias (a partir de la 3º columna)
abundance_columns <- colnames(data_for_wilcox)[-(1:2)]

# Normalizar las abundancias por el canal de referencia
for (protein_col in abundance_columns) {
  reference_intensity <- reference_channel$ReferenceIntensity[match(colnames(data_for_wilcox)[which(colnames(data_for_wilcox) == protein_col)], reference_channel$Sample.ID)]
  data_for_wilcox[, protein_col] <- data_for_wilcox[, protein_col] - reference_intensity
}

# Separar las muestras tumorales y normales
tumor_samples <- data_for_wilcox[data_for_wilcox$Class == "Tumor", ]
normal_samples <- data_for_wilcox[data_for_wilcox$Class == "NAT", ]

# Crear una lista para almacenar los resultados
wilcox_results <- list()

# Realizar el test para cada proteína
for (protein_col in abundance_columns) {
  # Obtener las abundancias para la proteína actual
  tumor_abundances <- tumor_samples[, protein_col]
  normal_abundances <- normal_samples[, protein_col]
  
  # Realizar el Wilcoxon signed-rank test
  result <- wilcox.test(tumor_abundances, normal_abundances, paired = TRUE)
  
  # Almacenar el resultado en la lista
  wilcox_results[[protein_col]] <- result
}

# Crear un dataframe con los resultados
wilcox_results_df <- data.frame(do.call(rbind, wilcox_results))

# Ordenar el dataframe por p-value
wilcox_results_df$p.value <- as.numeric(unlist(wilcox_results_df$p.value))
wilcox_results_df <- wilcox_results_df[order(wilcox_results_df$p.value), ]

# Aplicar el ajuste de Benjamini-Hochberg
wilcox_results_df$adjusted_p_value <- p.adjust(wilcox_results_df$p.value, method = "BH")

# Añadir la columna de los p valores ajustados a medianas_df
wilcox_results_df <- tibble::rownames_to_column(wilcox_results_df, var = "Proteina")

medianas_df <- merge(medianas_df, wilcox_results_df[, c("Proteina", "adjusted_p_value")], by = "Proteina", all.x = TRUE)

# Umbral de significancia
umbral_significancia <- 0.01
results_DE <- medianas_df %>%
  filter(adjusted_p_value <= umbral_significancia)

# Se han filtrado 5780 proteínas que superan el umbral de significancia de 12120 iniciales.



# --------------- LISTA DE GENES --------------------------------------------------------

# Obtener los nombres de los genes correspondientes a las proteínas
genes <- read.csv(file = 'genes.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Filtrar el dataframe por las proteínas de interés
genes_interes <- genes[genes$ProteinID %in% results_DE$Proteina, ]

results_DE <- merge(results_DE ,genes_interes[, c("ProteinID", "Gene.Symbol")], by.x = "Proteina", by.y = "ProteinID", all.x = TRUE)
results_DE <- results_DE[, c(1, ncol(results_DE), 2:(ncol(results_DE)-1))]



# --------------- CITOQUINAS -------------------------------------

# Cargar la lista de citoquinas de interés expresadas en ADP
citoquinas <- read.table("CITOQUINAS.txt", header = FALSE, sep = "\t")

# Filtrar cuales de mis proteínas significativas están en la lista de citoquinas
genes_citoquinas <- genes_interes[genes_interes$Gene.Symbol %in% citoquinas$V1, ]
genes_citoquinas <- merge(genes_citoquinas, results_DE[, c("Proteina", "Log2FoldChange", "adjusted_p_value")], by.x = "ProteinID", by.y = "Proteina", all.x = TRUE)

# Hago lo mismo pero incluyendo citoquinas no significativas
genes_medianas <- genes[genes$ProteinID %in% medianas_df$Proteina, ]
cytokines <- genes_medianas[genes_medianas$Gene.Symbol %in% citoquinas$V1, ]
cytokines <- merge(cytokines, medianas_df[, c("Proteina", "Log2FoldChange", "adjusted_p_value")], by.x = "ProteinID", by.y = "Proteina", all.x = TRUE)

# Guardo la información de las citoquinas que están entre el conjunto de proteínas significativas
write.table(genes_citoquinas, "citoquinas_proteomica.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# VOLCANO PLOT CITOQUINAS
label <- ifelse(-log10(cytokines$adjusted_p_value) > 2, cytokines$Gene.Symbol, "")
volcano_plot_cytokines <- ggplot(cytokines, aes(x = Log2FoldChange, y = -log10(adjusted_p_value), label = ifelse(-log10(adjusted_p_value) > 2, Gene.Symbol, ""))) +
  geom_point(aes(color = ifelse(Log2FoldChange > 1 & -log10(adjusted_p_value) > 2, ">2x Up", ifelse(Log2FoldChange > 0 & Log2FoldChange < 1 & -log10(adjusted_p_value) > 2, "Up", ifelse(Log2FoldChange < -1 & -log10(adjusted_p_value) > 2, ">2x Down", ifelse(Log2FoldChange < 0 & Log2FoldChange > -1 & -log10(adjusted_p_value) > 2, "Down", "FDR > 0.01"))))), alpha = 0.7) +
  geom_text_repel(aes(label = label), box.padding = 0.3, segment.color = 'grey50', size = 3) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  scale_color_manual(
    values = c(">2x Down" = "blue", ">2x Up" = "red", "Down" = "lightblue", "Up" = "pink", "FDR > 0.01" = "grey"),
    labels = c(">2x Down", ">2x Up", "Down", "FDR > 0.01", "Up")) +
  labs(title = "DE cytokines Tumor/NAT",
       x = "Median log2 Tumor/NAT FC",
       y = "-log10(Adjusted P-value)",
       color = NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 10)) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )



#------------------------- PLOTS -----------------------------------------------------

# Añadir a la tabla de resultados de medianas una columna de color ya que contiene las proteínas no significativas para que aparezcan en el volcano plot
medianas_df$color <- ifelse(medianas_df$Log2FoldChange > 1 & -log10(medianas_df$adjusted_p_value) > 2, "red",
                            ifelse(medianas_df$Log2FoldChange < -1 & -log10(medianas_df$adjusted_p_value) > 2, "blue",
                                   ifelse(medianas_df$Log2FoldChange >= 0 & medianas_df$Log2FoldChange <= 1 & -log10(medianas_df$adjusted_p_value) > 2, "pink",
                                          ifelse(medianas_df$Log2FoldChange >= -1 & medianas_df$Log2FoldChange < 0 & -log10(medianas_df$adjusted_p_value) > 2, "lightblue", "grey"))))

volcano_plot_2 <- ggplot(medianas_df, aes(x = Log2FoldChange, y = -log10(adjusted_p_value), color = color)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  scale_color_manual(
    values = c("blue", "grey", "lightblue", "pink", "red"),
    labels = c(">2x Down", "FDR > 0.01", "Down", "Up", ">2x Up")) +
  labs(title = "Wilcoxon Signed rank test Tumor/NAT N=41", x = "Median log2 Tumor/NAT FC", y = "-log10(Adjusted P-value)", color = NULL) +
  theme_minimal() +
  scale_x_continuous(limits = c(-3, 3)) + #He puesto esta escala pero hay valores que se encuentran más allá, hay que ver cuales son
  scale_y_continuous(limits = c(0, 10))

# Define el ancho del área del gráfico y su posición
width <- 8  # Ajusta según sea necesario
x_position <- (9 - width) / 2  # Calcula la posición para centrar el gráfico

# Ajusta el margen superior del tema y el tamaño del título
theme_settings <- theme(
  plot.margin = margin(t = 2, r = 0, b = 0, l = 45),
  plot.title = element_text(size = 13, hjust = 0.5)  # Ajusta el tamaño del título y su justificación
)
# Dibuja el gráfico con las etiquetas y ajustes de posición
library(grid)
grid.newpage()
grid.draw(volcano_plot_2 + theme_settings)

# Dibuja los textos fuera del área del gráfico
grid.text("Down\n Sig: 3279 (27%)\n 2-fold: 377", x = unit(0.09, "npc"), y = unit(0.2, "npc"), rot = 360, gp = gpar(col = "black", cex = 0.7))
grid.text("Up\n Sig: 2501 (21%)\n 2-fold: 241", x = unit(0.9, "npc"), y = unit(0.2, "npc"), rot = -360, gp = gpar(col = "black", cex = 0.7))


# ---------------- PROTEIN LOCATION PREDICTION -------------------------------

# Coger las proteínas con 2-fold up para determinar su localización subcelular
for_protein_location <- results_DE$Proteina[results_DE$Log2FoldChange > 1]
write.table(for_protein_location, "ProteinsIDs.txt", sep = "\t", quote = FALSE, row.names = FALSE)
