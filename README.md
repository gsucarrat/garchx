# garchx
The R Package *garchx* provides a simple, fast, flexible and robust framework for GARCH-X modelling. CRAN webpage: [https://CRAN.R-project.org/package=garchx]( https://CRAN.R-project.org/package=garchx)

# Installation
The following R command installs the stable version from CRAN:

    install.packages('garchx', dependencies = TRUE)

To install the current development version available here at Github, first download the tarball (i.e. the file named garchx_devel.tar.gz). Next, run:

    system("R CMD INSTALL --build garchx")

Alternatively, you can try installing the development version directly from GitHub with the following code:

    install.packages(
      "https://github.com/gsucarrat/garchx/raw/master/garchx_devel.tar.gz",
      repos = NULL, type = "source"
    )
    
# Resources
* CRAN webpage: [https://CRAN.R-project.org/package=garchx]( https://CRAN.R-project.org/package=garchx)
* Article (in The R Journal): [https://journal.r-project.org/archive/2021/RJ-2021-057/RJ-2021-057.pdf](https://journal.r-project.org/archive/2021/RJ-2021-057/RJ-2021-057.pdf)

# License
This package is free and open source software, licensed under GPL (>= 2)
