FROM r-base:4.5.3 AS base

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    make \
    libuv1-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "\
  pkgs <- c('Rcpp', 'ggplot2', 'gridExtra', 'reshape2', 'survey', 'testthat', 'knitr', 'jsonlite'); \
  install.packages(pkgs, repos='https://cloud.r-project.org', Ncpus=4L); \
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; \
  if (length(missing)) stop('Failed to install: ', paste(missing, collapse = ', ')) \
"

FROM base AS pkg

WORKDIR /pkg
COPY . .
RUN R CMD INSTALL .

CMD ["R"]
