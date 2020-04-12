source("/nfs/users/nfs_t/ta6/My_R_packages/TreeOfCells/R/Distances.R")
source("/nfs/users/nfs_t/ta6/My_R_packages/TreeOfCells/R/TM_Benchmark.R")
source("/nfs/users/nfs_t/ta6/My_R_packages/TreeOfCells/R/Normalization.R")
source("/nfs/users/nfs_t/ta6/My_R_packages/TreeOfCells/R/Feature_Selection.R")


gene_lengths <- read.delim("Mus_musculus.GRCm38.79.gene_lengths.txt",header=T, stringsAsFactors=FALSE, sep=" ")
gene_lengths$length <- gene_lengths$average_transcript_length

dist_function_list <- list(
	"euclid"=base_dist, 
	"pearson"=cor_dist, 
	"spearman"=function(x){cor_dist(x, type="spearman")},
	#"kl_poisson"=function(x){dmat_wrapper(x, d_func=KL_Poisson)},
	#"hellinger_poisson"=function(x){dmat_wrapper(x, d_func=Hellinger_poisson)},
	"phi_s"=phi_s,
	"rho_p"=rho_p
)

norm_function_list <- list(
	"none"=function(x, conditions=NULL){x$norm<-x$mus; return(x)},
	"myreg"=my_regression,
	"quantile"=quantile_norm,
	"cpm"=total_norm
#	"logcpm"=function(x, conditions=NULL){if(min(x$mus)>=0){ log_transform(total_norm(x))} else{return(x)}}
	#"glen"=function(x, conditions=NULL){gene_length_norm(x, gene_lengths)},
	#"tpm"=function(x, conditions=NULL){total_norm(gene_length_norm(x, gene_lengths))}
)

fs_function_list <- list(
	"none"=function(x){return(rownames(x$mus))},
	#"d50"=function(x){distribution_olap(x,quantile=0.5)},
	#"d75"=function(x){distribution_olap(x,quantile=0.75)},
	#"d80"=function(x){distribution_olap(x,quantile=0.80)},
	#"d90"=function(x){distribution_olap(x,quantile=0.90)},
	"f25"=function(x){min_zero_frac(x,diff=0.5)},
	"f50"=function(x){min_zero_frac(x,diff=0.75)},
	"f75"=function(x){min_zero_frac(x,diff=0.80)},
	"f90"=function(x){min_zero_frac(x,diff=0.90)},
	"f95"=function(x){min_zero_frac(x,diff=0.95)},
#	"fc"=FS_fold_change,
	"top3000"=function(x){top_diff(x,ntop=3000)},
	"top1000"=function(x){top_diff(x,ntop=1000)},
	"top500"=function(x){top_diff(x,ntop=500)},
	"top100"=function(x){top_diff(x,ntop=100)},
	"top50"=function(x){top_diff(x,ntop=50)},
	"top10"=function(x){top_diff(x,ntop=10)},
	"top5000"=function(x){top_diff(x,ntop=5000)},
	"top10000"=function(x){top_diff(x,ntop=10000)}
)

for (n1 in 1:length(norm_function_list)) {
n2=n1;
print(names(norm_function_list)[n1])
#for (n2 in 1:length(norm_function_list)) {
for (fs in 1:length(fs_function_list)) {
	
outfile<- paste("TM_Distance_Comparison", 
		names(norm_function_list)[n1], 
		names(norm_function_list)[n2], 
		names(fs_function_list)[fs], 
		"out.rds", sep="_")
if (file.exists(outfile)){next;}

print("running benchmark")
out <- list()
for (func in 1:length(dist_function_list)) {
	set.seed(101)
	print(names(dist_function_list)[func]);
	TM <- run_TM_benchmark(dist_function_list[[func]], 
		norm_facs=norm_function_list[[n1]], 
		norm_10x=norm_function_list[[n2]],
		norm_final=norm_function_list[[1]],
		fs_function=fs_function_list[[fs]])
	out[[names(dist_function_list)[func]]] <- TM;
}

print("saving results")
saveRDS(out, paste("TM_Distance_Comparison", 
		names(norm_function_list)[n1], 
		names(norm_function_list)[n2], 
		names(fs_function_list)[fs], 
		"out.rds", sep="_"))

}}#}

exit();

files <- Sys.glob("TM_Distance*_out.rds")
vals <- vector();
metadata <- vector();
for (f in files) {
	d <- unlist(readRDS(f));
	labs <- matrix(unlist(strsplit( names(d), "\\.")), ncol=2, byrow=T)
	labs <- cbind(labs, rep(f, nrow(labs)));
	metadata <- rbind(metadata, labs);
	vals <- c(vals, d);
}
out <- data.frame(distance=metadata[,1], norm=metadata[,3], measure=metadata[,2], val=vals)

saveRDS(out, "all_comparison2.rds")

tmp <- out[out$measure=="prop_recip",]; head(tmp[order(tmp$val, decreasing=TRUE),])


### Plotting

png("log2CPM.png", width=6, height=4, units="in", res=300)
par(mfrow=c(2,3))
par(mar=c(1,1,2,1))
bp_info <- vector(); 
for(n in names(out)) {
	image(out[[n]]$dist_all, main=n, xaxt="n", yaxt="n");
	abline(v=0.5); abline(h=0.5); 
	bp_info<-c(bp_info, out[[n]]$prop_match)
}
par(mar=c(6,4,2,1))
barplot(bp_info, names=names(out), las=2, ylab="%match")
dev.off()

#reciprocal
png("ReciprocalBest.png", width=5, height=5, units="in", res=300)
ylim=c(0,1)
files <- c(
	none="TM_Distance_Comparison_none_none_none_out.rds",
	myreg="TM_Distance_Comparison_myreg_myreg_none_out.rds",
	cpm="TM_Distance_Comparison_cpm_cpm_none_out.rds",
	logcpm="TM_Distance_Comparison_logcpm_logcpm_none_out.rds"
	)
par(mfrow=c(2,2))
for (f in 1:length(files)){
	out <- readRDS(files[f]);
	bp_info <- vector();
	for(n in names(out)) {
	        bp_info<-c(bp_info, out[[n]]$prop_recip)
	}
	par(mar=c(6,4,2,1))
	barplot(bp_info, names=names(out), las=2, ylab="% Reciprocal Best Hit", ylim=ylim, main=names(files)[f])
}
dev.off();

# FS
ylim = c(0,1);
mat = vector();
png("FS_results.png", width=6, height=6, units="in", res=300)
files <- Sys.glob("TM_*myreg_myreg_*.rds")
par(mfrow=c(4,5))
for (f in files) {
	out <- readRDS(f);
	name <- unlist(strsplit(f, "_"));
	name <- name[6]
        bp_info <- vector();
        for(n in names(out)) {
                bp_info<-c(bp_info, out[[n]]$prop_recip)
        }
	par(mar=c(6,4,2,1))
        barplot(bp_info, names=names(out), las=2, ylab="% Reciprocal Best Hit", ylim=ylim, main=name)
	mat = rbind(mat, bp_info)
}
dev.off()

# Tree
require("ape")
require("phangorn")
out <- readRDS("TM_Distance_Comparison_myreg_myreg_none_out.rds")
tree <- NJ(out$pearson$dist_all)
png("tmp_curr_best_tree.png", width=6, height=6, units="in", res=300)
par(mar=c(0,0,0,0))
plot(tree, cex=0.6, type="radial")
dev.off()

