modTraitCor = cor(mergedMEs, dataExpr, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
corblack=modTraitCor[which(row.names(modTraitCor)=='MEyellow'),]
head(corBlack[order(-corblack)])
#RAD51AP1     HDAC2     FOXM1    NCAPD2      TPI1      NOP2 
#0.9249354 0.9080872 0.8991733 0.8872607 0.8717050 0.8708449


table(colors)
