################################################################
## Zhang, et al., Nature, 2025: Companion analysis code ########
## Erythropoietin receptor on cDC1s dictates immune tolerance ##
## 12/05/24; Chris McGinnis, PhD; Stanford, Satpathy Lab #######
################################################################

library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(reshape2)
library(viridis)
library(speckle)
library(RColorBrewer)
library(ComplexHeatmap)


##############
## Figure 3 ##
##############
load("seu_spl.Robj")
load("anno_markers_spl.Robj")
load('anno_freq_spl.Robj')
load('anno_freq_spl_diff.Robj')

# Fig. 3A: Spleen cDC1 UMAP and annotation dotplot ~ xcr vs flox, treat vs untreat 
DimPlot(seu_spl, group.by = 'subtype', cols=c('black','grey80','cadetblue2','dodgerblue','navy','seagreen'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
g <- DotPlot(seu_spl, group.by = 'subtype', features=anno_markers_spl, cols='RdBu', dot.scale = 3, cluster.idents = T) + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size=6, angle = 90, vjust = 0.5, hjust=1))
g$data$id <- factor(g$data$id, levels=rev(c('pre','prolif','imm_early','imm_late','mat_early','mat_late')))
print(g)

# Fig. 3B: Frequency histogram (untreated flox vs treated flox)
ggplot(anno_freq_spl[which(anno_freq_spl$condition %in% c('flox_unt','flox_treat')), ], aes(x=subtype, y=freq, fill=condition)) + 
  geom_col(position=position_dodge(), color='black') +
  theme_classic() + theme(legend.position = 'none') + scale_fill_manual(values=alpha(c('grey50','steelblue3'),0.8))

# Fig. 3C: Sample UMAP (untreated flox vs treated flox), feature plots, violins for Itgae, Lgals3, Apol7c  
DimPlot(seu_spl, cells.highlight = colnames(seu_spl)[which(seu_spl@meta.data$fig_ident == 'flox_unt')], cols.highlight = alpha('grey50',0.8), sizes.highlight = 0.1, pt.size = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_spl, cells.highlight = colnames(seu_spl)[which(seu_spl@meta.data$fig_ident == 'flox_treat')], cols.highlight = alpha('steelblue3',0.8), sizes.highlight = 0.1, pt.size = 0.1) + NoLegend() + NoAxes()
FeaturePlot(seu_spl, 'Itgae', max.cutoff = 'q95', pt.size = 0.1) + NoAxes() + scale_color_viridis()+NoLegend()+theme(plot.title = element_blank())
FeaturePlot(seu_spl, 'Lgals3', max.cutoff = 'q95', pt.size = 0.1) + NoAxes() + scale_color_viridis()+NoLegend()+theme(plot.title = element_blank())
FeaturePlot(seu_spl, 'Apol7c', max.cutoff = 'q95', pt.size = 0.1) + NoAxes() + scale_color_viridis()+NoLegend()+theme(plot.title = element_blank())
VlnPlot(seu_spl, 'Itgae', idents = c('flox_unt','flox_treat'), cols=alpha(c('grey50','steelblue3'),0.8), pt.size = 0)+NoLegend()+theme(plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.line = element_line(linewidth = 0.75), axis.ticks = element_line(linewidth=0.75))
VlnPlot(seu_spl, 'Lgals3', idents = c('flox_unt','flox_treat'), cols=alpha(c('grey50','steelblue3'),0.8), pt.size = 0)+NoLegend()+theme(plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.line = element_line(linewidth = 0.75), axis.ticks = element_line(linewidth=0.75))
VlnPlot(seu_spl, 'Apol7c', idents = c('flox_unt','flox_treat'), cols=alpha(c('grey50','steelblue3'),0.8), pt.size = 0)+NoLegend()+theme(plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.line = element_line(linewidth = 0.75), axis.ticks = element_line(linewidth=0.75))
 
# Fig. 3D: Fold change frequency histogram (xcr vs flox, log2FC TLI vs UNT)
ggplot(anno_freq_spl_diff, aes(x=subtype, y=log2fc, fill=strain)) + 
  geom_col(position=position_dodge(), color='black') +
  theme_classic() + theme(legend.position = 'none') + scale_fill_manual(values=alpha(c('steelblue3','darkred'),0.8))

# Fig. 3F: Sample UMAP (treated flox vs treated xcr), feature plots, violins for Cd83, Rel, and Dnase1l3
DimPlot(seu_spl, cells.highlight = colnames(seu_spl)[which(seu_spl@meta.data$fig_ident == 'flox_treat')], cols.highlight = alpha('steelblue3',0.8), sizes.highlight = 0.1, pt.size = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_spl, cells.highlight = colnames(seu_spl)[which(seu_spl@meta.data$fig_ident == 'xcr_treat')], cols.highlight = alpha('darkred',0.8), sizes.highlight = 0.1, pt.size = 0.1) + NoLegend() + NoAxes()
FeaturePlot(seu_spl, 'Cd83', max.cutoff = 'q95', pt.size = 0.1) + NoAxes() + scale_color_viridis()+NoLegend()+theme(plot.title = element_blank())
FeaturePlot(seu_spl, 'Rel', max.cutoff = 'q95', pt.size = 0.1) + NoAxes() + scale_color_viridis()+NoLegend()+theme(plot.title = element_blank())
FeaturePlot(seu_spl, 'Dnase1l3', max.cutoff = 'q95', pt.size = 0.1) + NoAxes() + scale_color_viridis()+NoLegend()+theme(plot.title = element_blank())
VlnPlot(seu_spl, 'Cd83', idents = c('flox_treat','xcr_treat'), cols=alpha(c('steelblue3','darkred'),0.8), pt.size = 0)+NoLegend()+theme(plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.line = element_line(linewidth = 0.75), axis.ticks = element_line(linewidth=0.75))
VlnPlot(seu_spl, 'Rel', idents = c('flox_treat','xcr_treat'), cols=alpha(c('steelblue3','darkred'),0.8), pt.size = 0)+NoLegend()+theme(plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.line = element_line(linewidth = 0.75), axis.ticks = element_line(linewidth=0.75))
VlnPlot(seu_spl, 'Dnase1l3', idents = c('flox_treat','xcr_treat'), cols=alpha(c('steelblue3','darkred'),0.8), pt.size = 0)+NoLegend()+theme(plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.line = element_line(linewidth = 0.75), axis.ticks = element_line(linewidth=0.75))

# Fig. 3G: Frequency histogram (untreated flox vs untreated xcr)
ggplot(anno_freq_spl[which(anno_freq_spl$condition %in% c('flox_unt','xcr_unt')), ], aes(x=subtype, y=freq, fill=condition)) + 
  geom_col(position=position_dodge(), color='black') +
  theme_classic() + theme(legend.position = 'none') + scale_fill_manual(values=alpha(c('grey50','black'),0.8))

# Fig. 3H: Sample UMAP (untreated flox vs untreated xcr), feature plots, violins for Cd83, Rel, and Dnase1l3
DimPlot(seu_spl, cells.highlight = colnames(seu_spl)[which(seu_spl@meta.data$fig_ident == 'flox_unt')], cols.highlight = alpha('grey50',0.8), sizes.highlight = 0.1, pt.size = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_spl, cells.highlight = colnames(seu_spl)[which(seu_spl@meta.data$fig_ident == 'xcr_unt')], cols.highlight = alpha('black',0.8), sizes.highlight = 0.1, pt.size = 0.1) + NoLegend() + NoAxes()
VlnPlot(seu_spl, 'Dnase1l3', idents = c('flox_unt','xcr_unt'), cols=alpha(c('grey50','black'),0.8), pt.size = 0)+NoLegend()+theme(plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.line = element_line(linewidth = 0.75), axis.ticks = element_line(linewidth=0.75))
VlnPlot(seu_spl, 'Xcr1', idents = c('flox_unt','xcr_unt'), cols=alpha(c('grey50','black'),0.8), pt.size = 0)+NoLegend()+theme(plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.line = element_line(linewidth = 0.75), axis.ticks = element_line(linewidth=0.75))
# seu_spl@meta.data$subtype_ident <- paste(seu_spl@meta.data$subtype, seu_spl@meta.data$fig_ident, sep='_')
# seu_spl <- SetIdent(seu_spl, value=seu_spl@meta.data$subtype_ident)
VlnPlot(seu_spl, 'Cd274', idents = c('mat_early_flox_unt','mat_early_xcr_unt'), cols=alpha(c('grey50','black'),0.8), pt.size = 0)+NoLegend()+theme(plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.line = element_line(linewidth = 0.75), axis.ticks = element_line(linewidth=0.75))
VlnPlot(seu_spl, 'Cd274', idents = c('mat_late_flox_unt','mat_late_xcr_unt'), cols=alpha(c('grey50','black'),0.8), pt.size = 0)+NoLegend()+theme(plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.line = element_line(linewidth = 0.75), axis.ticks = element_line(linewidth=0.75))

##############
## Figure 4 ##
##############
load("seu_spl_epor.Robj")
load('anno_freq_spl_epor.Robj')

# Fig. 4:...

#######################
## Extended Figure 5 ##
#######################
load("seu_spl.Robj")
load('condition_markers_spl.Robj')
load('anno_markers_spl_epor.Robj')
load('markers_epor.Robj')

# Extended Fig. 5A: Sample UMAP (untreated xcr vs treated xcr)
DimPlot(seu_spl, cells.highlight = colnames(seu_spl)[which(seu_spl@meta.data$fig_ident == 'xcr_unt')], cols.highlight = alpha('black',0.8), sizes.highlight = 0.1, pt.size = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_spl, cells.highlight = colnames(seu_spl)[which(seu_spl@meta.data$fig_ident == 'xcr_treat')], cols.highlight = alpha('darkred',0.8), sizes.highlight = 0.1, pt.size = 0.1) + NoLegend() + NoAxes()

# Extended Fig. 5B: Condition marker dotplot 
g <- DotPlot(seu_spl, group.by = 'fig_ident', features=condition_markers, cols='RdBu', dot.scale = 3, cluster.idents = T) + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size=6, angle = 90, vjust = 0.5, hjust=1))
g$data$id <- factor(g$data$id, levels=rev(c('flox_unt','flox_treat','xcr_unt','xcr_treat')))
print(g)

# Extended Fig. 5D: Spleen cDC1 annotation and sample UMAPs ~ EpoR+/-  
DimPlot(seu_spl_epor, group.by = 'subtype', cols=c('cadetblue2','dodgerblue','navy','seagreen','grey80'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
DimPlot(seu_spl_epor, cells.highlight = colnames(seu_spl_epor)[which(seu_spl_epor@meta.data$epor == 'pos')], cols.highlight = alpha('red',0.8), sizes.highlight = 0.1, pt.size = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_spl_epor, cells.highlight = colnames(seu_spl_epor)[which(seu_spl_epor@meta.data$epor == 'neg')], cols.highlight = alpha('navy',0.8), sizes.highlight = 0.1, pt.size = 0.1) + NoLegend() + NoAxes()

# Extended Fig. 5E: Spleen cDC1 annotation dotplot ~ EpoR+/-   
g <- DotPlot(seu_spl_epor, group.by = 'subtype', features=anno_markers_spl_epor, cols='RdBu', dot.scale = 3, cluster.idents = T) + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size=6, angle = 90, vjust = 0.5, hjust=1))
g$data$id <- factor(g$data$id, levels=rev(c('prolif','imm_early','imm_late','mat_early','mat_late')))
print(g)

# Extended Fig. 5F: EpoR +/- marker dotplot
DotPlot(seu_spl_epor, group.by = 'epor', features=markers_epor, cols='RdBu', dot.scale = 3, cluster.idents = T) + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size=6, angle = 90, vjust = 0.5, hjust=1))

# Extended Fig. 5G: Condition marker dotplot for EpoR +/- cDC1s
DotPlot(seu_spl_epor, group.by = 'epor', features=condition_markers, cols='RdBu', dot.scale = 3, cluster.idents = T) + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size=6, angle = 90, vjust = 0.5, hjust=1))

# Extended Fig. 5H: Frequency histogram (EpoR+ vs EpoR- splenic cDC)
ggplot(anno_freq_spl_epor, aes(x=subtype, y=freq, fill=epor)) + geom_col(color='black', position=position_dodge()) + scale_fill_manual(values=alpha(c('red','navy'),0.8)) + theme_classic()


