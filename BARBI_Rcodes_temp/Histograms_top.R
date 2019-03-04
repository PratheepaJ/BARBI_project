library(tidyr)
blk <- 3
desired.sample.name <- "P_011"
desired.sample.index <- which(sample_names(psPlByBlock[[blk]]) %in% desired.sample.name)
tax_interested <- rownames(sort(otu_table(psPlByBlock[[blk]])[,desired.sample.index],decreasing = TRUE))[c(1:5,504)]
tax_interested_ind <- which(as.character(taxa_names(psPlByBlock[[blk]]))%in%tax_interested)
tax_names <- taxa_names(psPlByBlock[[blk]])[tax_interested_ind]

taxa.post <- taxa.post.all.sam[[desired.sample.index]]
burnIn <- 5001
signal.hist <- taxa.post[tax_interested_ind]
signal.hist <- lapply(signal.hist,function(x){x[-(1:burnIn),]})
signal.df <- data.frame(do.call("cbind",signal.hist))
colnames(signal.df) <- tax_names
signal.df$group <- rep("Real Reads",length=dim(signal.df)[1])

bg <- list()
for(ind in 1:length(tax_interested_ind)){
        bg[[ind]] <- rgamma(5000,shape=gammPrior[[desired.sample.index]][[1]][tax_interested_ind[ind]],rate = gammPrior[[desired.sample.index]][[2]][tax_interested_ind[ind]])
}
bg.df <- data.frame(do.call("cbind",bg))
colnames(bg.df) <- tax_names
bg.df$group <- rep("Contaminant Reads",length=dim(bg.df)[1])

bg.signal <- rbind(signal.df,bg.df)
bg.signal$group <- as.factor(bg.signal$group)
bg_sig_long <- tidyr::gather(bg.signal,key="Taxa",value="Reads",1:6)
bg_sig_long$Taxa <- as.factor(bg_sig_long$Taxa)

ggplot(bg_sig_long,aes(x=Reads,fill=group,color=group))+
        geom_density(alpha=.2)+
        facet_wrap(~Taxa,scales = "free")+
        scale_fill_manual(values=c("blue","red"))+
        scale_color_manual(values=c("blue","red"))+
        ggtitle(desired.sample.name)+
        theme(plot.title = element_text(hjust = 0.5))+
        theme(legend.title=element_blank())


ggsave("./Results_Bayesian_non_paired/Blk3/top62.png", plot = last_plot(), width = 10, height = 5, units = "in")
