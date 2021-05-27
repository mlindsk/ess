.PHONY: compile doc check install build_site run_main test clean

all: compile doc check

compile:
	Rscript -e "Rcpp::compileAttributes()"

doc:
	Rscript -e "devtools::document()"

check: 
	Rscript -e "devtools::check()"

check_fast: 
	Rscript -e "devtools::check(build_args = c('--no-build-vignettes'), args = c('--no-build-vignettes'))"

install:
	Rscript -e "devtools::install()"

install2:
	cd ..; R CMD build ess/; \
	R CMD INSTALL ess_1.0.tar.gz

build:
	Rscript -e "devtools::build()"; \
	cd /home/mads/Documents/phd/software/; \
	R CMD check --as-cran ess_1.0.tar.gz

test:
	Rscript -e "devtools::load_all(); tinytest::test_all('~/Documents/phd/software/efs')"

readme:
	Rscript -e "rmarkdown::render('README.Rmd')" #; \

clean:
	rm -f README.html
	cd src/ && rm -f *.o && rm -f *.so
