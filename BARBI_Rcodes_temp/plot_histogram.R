taxa <- 18
dd <- data.frame(real=taxa_post[[taxa]][-(1:burnIn),], cont=rgamma((itera-burnIn+1), shape= gammPrior[[sam]][[1]][taxa], rate = gammPrior[[sam]][[2]][taxa]))
ddl <- gather(dd)
ggplot(ddl, aes(x=value, col = key,fill=key))+ geom_density()+
        geom_vline(xintercept = as.numeric(gammPrior[[sam]][[3]][taxa]))

