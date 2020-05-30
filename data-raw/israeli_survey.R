## code to prepare `israeli_survey` dataset goes here
israeli_survey <- as.table(
  matrix(
    c(128, 52, 81, 14, 12,
      118, 28, 32, 6, 7,
      218, 28, 97, 12, 14,
      11, 2, 4, 1, 1,
      104, 22, 61, 8, 5,
      117, 24, 70, 9, 7,
      42, 6, 20, 2, 0,
      48, 16, 104, 14, 9),
    ncol = 5, byrow = TRUE
  )
)

colnames(israeli_survey) <- c("EUAM", "IFEA", "ASAF", "IFAA", "IFI")
rownames(israeli_survey) <- c("OTH", "POL", "MIL", "ECO", "ENR", "SAB", "MTO", "PER")


usethis::use_data(israeli_survey, overwrite = TRUE)
