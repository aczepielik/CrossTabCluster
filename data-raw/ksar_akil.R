ksarakil <- as.table(
  matrix(
    c(2, 12, 6, 12, 4,
      16, 44, 14, 6, 4,
      72, 105, 54, 55, 169,
      111, 87, 114, 148, 115,
      35, 40, 48, 47, 55,
      60, 74, 76, 53, 56,
      62, 51, 206, 127, 66,
      24, 50, 80, 67, 30,
      52, 177, 344, 205, 75,
      21, 81, 138, 31, 22),
    ncol = 5, byrow = FALSE
  )
)

rownames(ksarakil) <- 1:10
colnames(ksarakil) <- c("Partially cortical", "Non cortical", "Flake blades", "Blades", "Bladelets")

usethis::use_data(ksarakil)
