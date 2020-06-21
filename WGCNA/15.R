

module = "blue"
column = match(module, modNames);
moduleGenes = mergedColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for GBM/LGG",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


table(moduleGenes)

