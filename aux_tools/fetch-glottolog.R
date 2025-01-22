# utility to fetch the latest Glottolog taxonomy from the CLDF dataset
library(tidyverse)

# 1. query API endpoint
o <- jsonlite::fromJSON('https://zenodo.org/api/records/14006636')

# 2. get url
latest_url <- o$files[1,]$links$self

# 3. download
glottolog_zipfile <- tempfile()
download.file(latest_url, glottolog_zipfile, method="curl", extra='-L')

# 4. find values.csv
valuefile <- grep('values.csv$', unzip(glottolog_zipfile, list=TRUE)$Name, ignore.case=TRUE, value=TRUE)

# 5. get file
# the hierarchy data is encoded in the Values table of the CLDF dataset
values <- readr::read_csv(unz(glottolog_zipfile, filename = valuefile))

# extract the languoid classification
glottolog_languoids <- filter(values, Parameter_ID %in% c("level", "classification")) |>
  select(id = Language_ID, Parameter_ID, Value) |>
  pivot_wider(names_from = Parameter_ID, values_from = Value) |>
  mutate(
    family_id = str_remove(classification, "/.+$"),
    parent_id = str_remove(classification, "^.+/")
  ) |>
  select(id, parent_id, family_id, level) |> arrange(id)

  save(glottolog_languoids, file='data/glottolog_languoids.rda')

