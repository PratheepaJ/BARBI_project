main_factor = "Sample_Type"

if (dim(otu_table(ps))[2] != nsamples(ps)) {
        otu_table(ps) <- t(otu_table(ps))
}

ot <- otu_table(ps) %>% data.frame %>% as.matrix

if (!isTRUE(all(ot == floor(ot)))) {
        # regression for negative binomial data so the
        stop("otu_table entries must be integer")
}

# library size normalization factor.
geo_mean <- function(x) {
        if(all(x == 0)){
                val <- 0
        }else{
                val <- exp(sum(log(x[x > 0]))/length(x))
        }
        return(val)
}

geom_mean_row <- apply(ot, 1, FUN = geo_mean)


sj <- estimateSizeFactorsForMatrix(ot, median, geoMeans = geom_mean_row)

des <- as.formula(paste("otuT", "~", "offset(asinh(allSjT))"))

samdf <- sample_data(ps) %>% data.frame


des_v <- as.formula(paste("~", 1))
mm <- model.matrix(des_v, data = samdf)
v <- asinhVoom(counts = ot, design = mm, sj = sj)
weights.cal <- v$weights


com_beta <- function(taxIndex, sampleDf, otuDf, allSj, weightDf, desingGEE, b) {

        otuT <- as.numeric(otuDf[taxIndex, ])
        allSj <- as.numeric(allSj)
        weightT <- as.numeric(weightDf[taxIndex, ])
        dffT <- cbind(sampleDf, otuT = otuT, allSjT = allSj, weightT = weightT)

        glmft.tx <- tryCatch(MASS::glm.nb(desingGEE, data = dffT, weights = weightT, method = "glm.fit", link = arcsinhLink()),
                error = function(e){
                        dffT$otuT <- dffT$otuT + sample(c(0,1),length(dffT$otuT), replace = TRUE); glm(desingGEE, data = dffT, weights = weightT, method = "glm.fit", family = poisson()) #when count is very small
                })
}


ind <- as.list(c(1:ntaxa(ps)))

df.beta.hat <- lapply(ind, function(x) {
        rt <- com_beta(x, sampleDf = samdf, otuDf = ot, allSj = sj, weightDf = weights.cal, desingGEE = des)
        return(rt)
})

df.beta.hat <- data.frame(do.call("rbind", df.beta.hat))




inv_asinh <- function(x) {
        y <- 0.5*exp(-x)*(exp(2*x)-1)
        return(y)
}

inv_asinh(4.502)

theta = 2.84
