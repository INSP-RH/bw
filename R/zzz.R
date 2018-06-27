.onAttach  <- function(libname, pkgname){
  packageStartupMessage(
    paste0("Please cite both our package and the model. Feel free to contact us with any questions."), 
    domain = NULL, appendLF = TRUE)
}