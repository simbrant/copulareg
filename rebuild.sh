
# Re-genarate documentation via devtools::document()
Rscript document.R

# Build tar file
R CMD build copulareg

# Check cran-compability
R CMD check copulareg_0.1.0.tar.gz

# re-install package
R CMD INSTALL copulareg_0.1.0.tar.gz