.onAttach  <- function(libname, pkgname){
  packageStartupMessage(
    paste0("Please cite the package as:\n",
           "Dalia Camacho-Garcia-Formenti and Rodrigo Zepeda-Tello ",
           "(2018). bw: Dynamic Body Weight Models for Children and ",
           "Adults. R package version 1.0.0.\n",
           "Feel free to contact us with any questions."), 
    domain = NULL, appendLF = TRUE)
}