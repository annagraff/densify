# utility to fetch the latest Glottolog taxonomy from the CLDF dataset
library(tidyverse)

# download the Glottolog data from Zenodo
glottolog_zipfile <- tempfile()
on.exit(unlink(glottolog_zipfile))
download.file("https://zenodo.org/records/10804582/files/glottolog/glottolog-cldf-v5.0.zip", glottolog_zipfile)


# the hierarchy data is encoded in the Values table of the CLDF dataset
values <- read_csv(unz(
  glottolog_zipfile,
  filename = "glottolog-glottolog-cldf-4dbf078/cldf/values.csv",
))

# extract the languoid classification
glottolog_languoids <- filter(values, Parameter_ID %in% c("level", "classification")) |>
  select(id = Language_ID, Parameter_ID, Value) |>
  pivot_wider(names_from = Parameter_ID, values_from = Value) |>
  mutate(
    family_id = str_remove(classification, "/.+$"),
    parent_id = str_remove(classification, "^.+/")
  ) |>
  select(id, parent_id, family_id, level)

  save(glottolog_languoids, file='data/glottolog_languoids.rda')
