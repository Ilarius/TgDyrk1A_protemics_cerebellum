load("~/OneDrive - CRG - Centre de Regulacio Genomica/Dropbox (CRG)/cerebellum.rda")

library(ggpubr)
###Figure 2


> ggarrange(gl[[1]], gl[[2]], ncol = 2, nrow = 2, labels="AUTO")
> ?ggarrange
widths	
(optional) numerical vector of relative columns widths. For example, in a two-column grid, widths = c(2, 1) would make the first column twice as wide as the second column.

heights	
same as widths but for column heights.

go=clusterProfiler::dotplot(cr_bp_all, showCategory=100, by="rowPercentage",label_format = 100)

#figure2
library("cowplot")
ggdraw() +
  draw_plot(go, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(tE, x = 0, y = 0, width = .5, height = .5) +
  draw_plot(dotplotfirst, x = .5, y = 0, width = .5, height = .5) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0, 0.5), y = c(1, 0.5, 0.5))
##10.5 x 6.5 inches

go2=clusterProfiler::dotplot(cr_bp_all2_s, showCategory=100, by="rowPercentage",label_format = 50)





tgt=unique(c(lista_de_proteome$`TG.EGCG-TG.NT`))
wtt=unique(c(lista_de_proteome$`WT.EGCG-WT.NT`))
int=unique(c(lista_de_proteome$`(TG.EGCG-TG.NT)-(WT.EGCG-WT.NT)`))


a=data_proteomic$pwt_wtegcg
b=data_proteomic$ptg_tgegcg
subsettoEGCG=unique(c(lista_de_proteome$`TG.EGCG-TG.NT`,lista_de_proteome$`WT.EGCG-WT.NT`,lista_de_proteome$`(TG.EGCG-TG.NT)-(WT.EGCG-WT.NT)`))

a=data_proteomic$pwt_wtegcg[subsettoEGCG,]

b=data_proteomic$ptg_tgegcg[subsettoEGCG,]

df=data.frame("TG"=b$log2FC, "WT"=a$log2FC, "significant"="no")
df$significant[which(rownames(a) %in% lista_de_proteome$`(TG.EGCG-TG.NT)-(WT.EGCG-WT.NT)`)]="yes"


gldiff=ggplot(df, aes(x=TG, y=WT,color=significant)) + geom_point()+ scale_color_manual(values=c("black", "red"))+theme_light()
eulerr_options(labels=list(fontsize=11))

vennc=plot(euler(list("TG and WT"=common, "Differently changing"=interactions)), shape = "ellipse", quantities = TRUE,adjust_labels=TRUE)
#figure3
ggdraw() +
  draw_plot(venna, x = 0, y = 0.7, width = 1/3, height = 0.3) +
  draw_plot(vennb, x = 1/3, y = 0.7, width = 1/3, height = 0.3) +
  draw_plot(vennc, x = 2/3, y = 0.7, width = 1/3, height = 0.3) +
  
  draw_plot(gldiff, x = 0, y = 0, width = .35, height = .6) +
  draw_plot(go2, x = .35, y = 0, width = .65, height = .7) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 15,
                  x = c(0, 1/3, 2/3,0,0.35 ), y = c(1, 1, 1, 0.7,0.7))

igraph::degree(graph_joint),main="", ylab="number of proteins", xlab="number of interactions"
h=qplot(igraph::degree(graph_joint),
      geom="histogram",fill=I("gray"), 
      col=I("black"), 
      xlab = "number of interactions",
      ylab="number of proteins",binwidth = 1  )+theme_light()


library("poweRlaw")
m1 = displ$new(igraph::degree(graph_joint))
m1$setPars(estimate_pars(m1))

m2 = dispois$new(igraph::degree(graph_joint))
m2$setPars(estimate_pars(m2))

m3 = disexp$new(igraph::degree(graph_joint))
m3$setPars(estimate_pars(m3))
m4 = dislnorm$new(igraph::degree(graph_joint))
m4$setPars(estimate_pars(m4))


g=plot(m1, draw=FALSE)
g1=lines(m1, draw=FALSE)
g2=lines(m2, draw=FALSE)
g3=lines(m3, draw=FALSE)
g4=lines(m4, draw=FALSE)

g$y1=NA
g$y1[1:13]=g$y[2:14]
g$y1[14]=1-g$y[1]

g1$group="power-law"
g2$group="Poisson"
g3$group="Exponential"
g4$group="log-normal"
gtot=rbind( g1, g2, g3, g4)
colnames(gtot)=c("x", "y", "fit" )


powerlawplot=ggplot(g) + geom_point(aes(x=x, y=y)) + labs(x="interactions", y="CDF")  + theme_bw() + 
  geom_line(data=gtot, aes(x=x, y=y, group=fit, color=fit))+ scale_x_continuous(trans='log2',breaks=c(1,2, 5, 10)) +
  scale_y_continuous(trans='log10', breaks=c(0.5,0.1,0.02,0.005), limits = c(0.0035,1))



#figure3
ggdraw() +
  draw_plot(h, x = 0, y = 0.6, width = .4, height = 0.4) +
  draw_plot(powerlawplot, x = .5, y = 0.6, width = .4, height = 0.4) +
  draw_plot(gl[[1]], x = 0, y = 0, width = .5, height = 0.6) +
  draw_plot(gl[[2]], x = .5, y = 0, width = .5, height = .6) +
  draw_plot_label(label = c("A", "B", "C" ,"D"), size = 15,
                  x = c(0, 0.5, 0,0.5), y = c(1, 1, 0.6, 0.6))



df1=cbind(data_proteomic$pwttg$log2FC*-1, -log10(data_proteomic$pwttg$adj.pvalue))
colnames(df1)=c("log2FC", "significance")
df2=cbind(data_proteomic$pwt_wtegcg$log2FC*-1, -log10(data_proteomic$pwt_wtegcg$adj.pvalue))
colnames(df2)=c("log2FC", "significance")
df3=cbind(data_proteomic$ptg_tgegcg$log2FC*-1, -log10(data_proteomic$ptg_tgegcg$adj.pvalue))
colnames(df3)=c("log2FC", "significance")
df4=cbind(data_proteomic$pint_egcg$log2FC*-1, -log10(data_proteomic$pint_egcg$adj.pvalue))
colnames(df4)=c("log2FC", "significance")



plot(data_proteomic$pint_egcg$log2FC*-1, -log10(data_proteomic$pint_egcg$adj.pvalue),pch=20, xlab="log2(FC)", ylab="significance", main="(TG.EGCG-TG.NT)-(WT.EGCG-WT.NT)")
abline(h=-log10(0.05), col="red", lty="dashed")
abline(v=-0.3, col="blue", lty="dashed")
abline(v=0.3, col="blue", lty="dashed")
v4=grab_grob()