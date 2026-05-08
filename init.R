#############################
#
#   Program: init.R
#
#   Author: Alvaro Quijano
#
#
#

## install packages if needed
# install.packages(c('devtools','roxygen2'))
options(scipen = 999)

data <- read.csv("../data/samplemiss_10.csv") %>% filter(v_num==1)

# define the sampling design
dsgn <-svydesign(data = data,
                    strata=~strat,
                    ids=~BGid,
                    weights=~bghhsub_s2, nest = TRUE)


# Run survey_mi
fit <- survey_mi(design           = dsgn,
  analysis         = y_gfr ~ x17 + y_bmi + x12 + x14 + x18,
  variables        = c("y_gfr", "y_bmi", "x12", "x14", "x17", "x18"),
  methods          = list(gaussian = "x17", binomial = "x18"),
  analysis_family  = "gaussian",
  B                = 100,
  n_rounds         = 10,
  seed             = 123
)
fit$summary
fit2$summary


