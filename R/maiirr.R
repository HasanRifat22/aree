

BiocManager::install("restfulr")
# Install AnnotationHub
BiocManager::install("AnnotationHub")

# Install EnsDb annotation package
BiocManager::install(
  "EnsDb.Hsapiens.v86",
  dependencies = TRUE
)

# Install required dependencies
BiocManager::install(
  c(
    "restfulr",
    "GenomicFeatures",
    "BiocGenerics",
    "BiocParallel",
    "curl",
    "XML",
    "RCurl",
    "httr",
    "jsonlite"
  ),
  ask = FALSE,
  update = TRUE
)