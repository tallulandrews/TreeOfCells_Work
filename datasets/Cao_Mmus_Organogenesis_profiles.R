source("/nfs/users/nfs_t/ta6/My_R_packages/TreeOfCells/R/Fitting_ZINB.R")
require("Matrix")
require("methods")

require(SingleCellExperiment)
sce <- readRDS("Cao_Mmus_Organogenesis_p1.rds");

set <- as.numeric(commandArgs(trailingOnly=TRUE));

OUT_list <- list();
for (type in levels(factor(sce$cell_type1))[set]) {
        print(type)
        if (sum(sce$cell_type1==type) < 20) {next;}
        out <- fit_ZINB_to_matrix(as.matrix(assays(sce)$counts[,sce$cell_type1==type]));
        OUT_list[[type]] <- out;
}
saveRDS(OUT_list, file=paste("Cao_Mmus_Organogenesis_Profiles_p1_",set,".rds", sep=""));

