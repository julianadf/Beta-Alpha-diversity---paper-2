# test autocorrelation of model residuals

library(ncf)

# data are in data frame df, with X and Y columns containing longitude and latitude
# data were analysed with model m

ncf.spl.res <- spline.correlog(df$X, df$Y, resid(m, type  = "pearson"), resamp=1000, na.rm = T)
plot(ncf.spl.res)

