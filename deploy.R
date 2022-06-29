if(!require("devtools"))
  install.packages("devtools")
devtools::install_github("rstudio/rsconnect")

library(rsconnect)
rsconnect::setAccountInfo(name='kasaha1', token='6F76A57070059AF3831C684ACAEFAC4B', secret='i+R1hIDZEpk3HHokrrniVGLEH/fQ3+unnK3lfC/Y')
deployApp()
