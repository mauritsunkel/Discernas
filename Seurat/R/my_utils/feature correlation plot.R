file <- ""
data <- readRDS(file)

x = data@assays$SCT@data['XPO7',]
y = data@assays$SCT@data['SETD1A',]
plot(x, y, pch = 19, col = "lightblue", xlab = 'XPO7', ylab = 'SETD1A', main = 'SCT normalised count data')
abline(lm(y ~ x), col = "red", lwd = 3)
text(paste("Correlation:", round(cor(x, y), 2)), x = .3, y = .2)
table(y, x)
table(data@assays$SCT@counts['XPO7',])
table(data@assays$SCT@counts['SETD1A',])
table(data@assays$SCT@counts['XPO7',] > 0, data@assays$SCT@counts['SETD1A',] > 0)
table(data@assays$SCT@counts['ONECUT3',])
