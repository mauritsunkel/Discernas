
file <- ""
data <- readRDS(file)

`%notin%` <- Negate(`%in%`)

subset(data, cluster = c(1,2,3))
# or use %notin%
