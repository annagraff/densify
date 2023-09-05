#' WALS language-feature matrix
#'
#' This dataset is the language-feature matrix generated from WALS Online (retrieved on 5 September 2023).
#'
#' @name wals
#' @docType data
#' @source WALS Online
#' @references
#' - Dryer, Matthew S. & Haspelmath, Martin (eds.) 2013.
#' WALS Online (v2020.3). Zenodo.
#' https://doi.org/10.5281/zenodo.7385533
#' (Available online at https://wals.info, Accessed on 2023-09-05.)

#' @usage data(wals)
#' @import tidyverse

library(tidyverse)

temp <- tempfile()
download.file("https://zenodo.org/record/7385533/files/cldf-datasets/wals-v2020.3.zip",temp)
wals_languages <- read.csv(unz(temp, "cldf-datasets-wals-878ea47/cldf/languages.csv"))
wals_values <- read.csv(unz(temp, "cldf-datasets-wals-878ea47/cldf/values.csv"))
wals_parameters <- read.csv(unz(temp, "cldf-datasets-wals-878ea47/cldf/parameters.csv"))
wals_codes <- read.csv(unz(temp, "cldf-datasets-wals-878ea47/cldf/codes.csv"))
unlink(temp)
rm(temp)

wals_parameters$Feature <- apply(wals_parameters, 1, function(x) paste(x[1]," ",x[2],sep=""))
wals_parameters <- wals_parameters %>% select(c("ID","Feature"))
wals_codes <- wals_codes %>% select(c("ID","Name"))

names(wals_languages)[1]<-"Language_ID"
names(wals_parameters)[1]<-"Parameter_ID"
names(wals_codes)<-c("Code_ID","Code_Name")

wals <- wals_languages %>% merge(wals_values, by = "Language_ID", all = TRUE) %>%
  merge(wals_parameters, by = "Parameter_ID", all = TRUE) %>%
  merge(wals_codes, by = "Code_ID", all = TRUE) %>%
  dplyr::filter(Glottocode!="") %>%
  select(c("Glottocode","Feature","Code_Name")) %>% na.omit()

wals <- as.data.frame(pivot_wider(wals,
                                  names_from = Feature,
                                  values_from = Code_Name,
                                  values_fill = "?",
                                  values_fn = function(x){if (length(unique(x))==1){unique(x)} else ("?")}))
lgs <- wals$Glottocode
rownames(wals) <- lgs

rm(wals_codes, wals_parameters, wals_languages, wals_values, lgs)

"wals" <- select(wals,-Glottocode)
