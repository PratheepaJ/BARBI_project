# Find index your desired species:
blk <- 3
desired.taxa.name <- "s__Escherichia_coli"
desired.taxa.index <- which(taxa_names(psPlByBlock[[blk]]) %in% desired.taxa.name)
desired.sample.name <- "P_044"
desired.sample.index <- which(sample_names(psPlByBlock[[blk]]) %in% desired.sample.name)

taxa.post <- taxa.post.all.sam[[desired.sample.index]]
burnIn <- 5000
signal.hist <- taxa.post[[desired.taxa.index]][-(1:burnIn),]
signal.hist <- signal.hist[1:5000]
background.hist <- rgamma(5000,shape=gammPrior[[desired.sample.index]][[1]][desired.taxa.index],rate = gammPrior[[desired.sample.index]][[2]][desired.taxa.index])
hist.df <- data.frame(Type= factor(rep(c("Real Reads", "Contaminant Reads"), each = 5000)), Reads = c(signal.hist, background.hist))

hdi.v <- hdi(taxa.post[[desired.taxa.index]][-(1:burnIn),],credMass = .95)
lower.s[[desired.taxa.index]] <- round(hdi.v[1],digits = 0)
upper.s[[desired.taxa.index]] <- round(hdi.v[2],digits = 0)

hdi.b <- hdi(rgamma(5000,shape=gammPrior[[desired.sample.index]][[1]][desired.taxa.index],rate = gammPrior[[desired.sample.index]][[2]][desired.taxa.index]),credMass = .95)
lower.b[[desired.taxa.index]] <- round(hdi.b[1],digits = 0)
upper.b[[desired.taxa.index]] <- round(hdi.b[2],digits = 0)

xj.raw.reads=as.numeric(gammPrior[[desired.sample.index]][[3]])[desired.taxa.index]

hist.df <- data.frame(Read_Type= factor(rep(c("Real Reads", "Contaminant Reads"), each = 5000)), Reads = c(signal.hist, background.hist))


hist.plot.label = paste(desired.taxa.name, " in ", desired.sample.name,  sep = "")

ggplot(hist.df, aes(x=Reads, fill=Read_Type)) + geom_density(alpha=.3)  + 
        geom_vline(xintercept =lower.s[[desired.taxa.index]], color = "red") + 
        geom_vline(xintercept =upper.s[[desired.taxa.index]], color = "red") + 
        geom_vline(xintercept =lower.b[[desired.taxa.index]], color = "blue") + 
        geom_vline(xintercept =upper.b[[desired.taxa.index]], color = "blue") +
        scale_fill_manual(values = c("blue", "red")) + 
        ggtitle("Result in patient with pos. E. coli blood culture") + theme_classic() +
        geom_vline(xintercept = xj.raw.reads, color = "violet")+
        theme(plot.title = element_text(hjust = 0.5)) 
# + geom_vline(xintercept = xj.raw.reads, color = "purple")

hist.plot.file.save.name <- paste("./Results_Bayesian_non_paired/Blk3/", desired.taxa.name, "_in_", desired.sample.name, ".png",  sep = "")

ggsave(hist.plot.file.save.name, plot = last_plot(), width = 10, height = 5, units = "in")
