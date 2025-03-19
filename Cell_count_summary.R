seu_obj = readRDS("~/VisiumHD/Myeloid_Seurat.rds")

meta_data = seu_obj@meta.data

# Summarize cell type information
celltype_summary = meta_data %>%
  group_by(Myeloid_Annotation) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  arrange(desc(cell_count))

# View the summary
head(celltype_summary)

# Calculate total number of cells in the sample
total_cells = nrow(meta_data)

# Compute frequency of each cell type
celltype_summary = celltype_summary %>%
  mutate(frequency = (cell_count / total_cells) * 100)

# Save the results to a CSV file
write.csv(celltype_summary, "~/VisiumHD/celltype_frequencies.csv", row.names = FALSE)

# Print summary
print(celltype_summary)
