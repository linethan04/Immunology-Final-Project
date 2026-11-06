# One-time setup: downloads the authors' repository at the exact commit this
# analysis used, extracts their nrstats helper package, and installs any
# missing CRAN packages. Run from the repository root:  Rscript wes/setup.R
commit <- '137111c38a32f8626ee6a436e14ac667bd31b1e7'
ext <- 'external'
dir.create(ext, showWarnings = FALSE)
repo <- file.path(ext, 'bms038_analysis')
if (!dir.exists(repo)) {
  stopifnot(system2('git', c('clone', 'https://github.com/riazn/bms038_analysis.git', repo)) == 0)
}
stopifnot(system2('git', c('-C', repo, 'checkout', '--quiet', commit)) == 0)
# nrstats is only used for its survival plotting functions, so extract rather than install
if (!dir.exists(file.path(ext, 'nrstats'))) {
  untar(file.path(repo, 'packages', 'nrstats_0.1.0.tar.gz'), exdir = ext)
}
pkgs <- c('ggplot2', 'survival', 'RColorBrewer', 'ks')
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing, repos = 'https://cloud.r-project.org')
cat('Setup complete. Now run: Rscript wes/run_wes.R\n')
