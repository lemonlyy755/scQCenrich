#' Generic, cross-tissue marker signatures (broad/extended)
#' @param species "mouse" or "human"
#' @param tier "broad" or "extended"
#' @export
generic_signatures <- function(species = c("mouse","human"), tier = c("broad","extended")) {
  species <- match.arg(species)
  tier    <- match.arg(tier)

  sig_broad <- list(
    Tcell        = c("CD3D","CD3E","TRAC","CD2","LCK"),
    Bcell        = c("MS4A1","CD79A","CD79B","CD74"),
    Plasma       = c("MZB1","XBP1","JCHAIN","SDC1"),
    NK           = c("NKG7","GNLY","PRF1","GZMB"),
    Myeloid      = c("LYZ","LST1","S100A8","S100A9","FCGR3A","ITGAM"),
    Dendritic    = c("ITGAX","HLA-DRA","CCR7"),
    Erythroid    = c("HBB","HBA1","ALAS2"),
    Megakaryo    = c("PPBP","PF4","GP9"),
    Endothelial  = c("PECAM1","KDR","VWF"),
    PericyteSMC  = c("ACTA2","TAGLN","MYH11","RGS5"),
    Fibroblast   = c("COL1A1","COL1A2","DCN","LUM","PDGFRA"),
    Epithelial   = c("EPCAM","KRT8","KRT18","KRT19"),
    Neuron       = c("SNAP25","RBFOX3","TUBB3","MAP2"),
    Astrocyte    = c("GFAP","AQP4","SLC1A3"),
    Oligodendro  = c("MBP","MOG","PLP1"),
    Microglia    = c("CX3CR1","P2RY12","TMEM119")
  )

  sig_ext <- list(
    Cardiomyocyte  = c("TNNT2","MYH6","ACTC1"),
    SkelMuscle     = c("MYH1","MYH7","ACTA1"),
    Adipocyte      = c("ADIPOQ","PPARG","FABP4"),
    Hepatocyte     = c("ALB","TTR","APOA1"),
    Cholangiocyte  = c("KRT19","KRT7","SOX9"),
    Melanocyte     = c("PMEL","TYR","DCT"),
    AT2_Lung       = c("SFTPC","SFTPB","SLC34A2"),
    AT1_Lung       = c("PDPN","AGER"),
    Pancreas_Acinar= c("PRSS1","CPA1","CELA3A"),
    BetaCell       = c("INS","IAPP","PDX1"),
    Kidney_PT      = c("SLC34A1","LRP2","ALDOB"),
    Podocyte       = c("NPHS1","PODXL","WT1"),
    Enterocyte     = c("ALPI","AQP8","FABP1"),
    Goblet         = c("MUC2","TFF3","SPDEF")
  )

  out <- sig_broad
  if (tier == "extended") out <- c(out, sig_ext)
  lapply(out, unique)
}
