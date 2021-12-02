#!/usr/bin/R

# Eric Minikel
# script to run ExomeDepth
# how to run:
# runExomeDepth.r -b bamlist.txt -o /output/path/ -v > runExomeDepthOutput.txt

start_time = Sys.time()

require(optparse) # http://cran.r-project.org/web/packages/optparse/optparse.pdf
require(ExomeDepth)
require(GenomicRanges)

packageVersion("ExomeDepth")

options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("-b", "--bamlist"), action="store", default='', 
              type='character', help="Path to list of BAMs"),
  make_option(c("-o", "--outdir"), action="store", default='./',
              type='character', help="Output directory [default %default]"),
  make_option(c("-p", "--panel"), action="store", default='', type='character', 
		help="Path to panel.csv"),
  make_option(c("-l", "--cnv_length"), action="store", default=500,
              type='integer', help="Expected CNV length [default 500]"),
  make_option(c("-s", "--sensitivity"), action="store", default=0.1,
              type='double', help="Sensitivity [default 0.1]"),	      
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print verbose output [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Do not print verbose output (this is the default)")
)
opt = parse_args(OptionParser(option_list=option_list))

print(opt,file=stderr())

if (file.exists(opt$panel)) {
	panels.hg19 = read.csv(opt$panel, header=TRUE)
} else {
    cat("You need to specify a valid panel.csv using -p.\n",file=stderr())
    cat(paste("The filename you specified was '",opt$panel,"'.",sep=''),file=stderr())
    stop()
}
panels.hg19.GRanges <- GRanges(seqnames = panels.hg19$chromosome, IRanges(start=panels.hg19$start,end=panels.hg19$end),names=panels.hg19$name)

# read list of BAMs
# to avoid writing a tryCatch I use file.exists, and set default to '', which is
# a file that never exists.
if (file.exists(opt$bamlist)) {
    # read bam list directly into a vector (note use of $V1)
    bams = read.table(opt$bamlist,header=FALSE)$V1
} else {
    cat("You need to specify a valid BAM list using -b.\n",file=stderr())
    cat(paste("The filename you specified was '",opt$bamlist,"'.",sep=''),file=stderr())
    stop()
}

# read output directory
# note right now if not specified, I stop execution. 
# an alternative is to create the dir.
# see http://stackoverflow.com/questions/4216753/check-existence-of-directory-and-create-if-doesnt-exist
if (file.exists(opt$outdir)) {
    setwd(opt$outdir)
} else {
    cat("You need to specify a valid output directory using -o.\n",file=stderr())
    cat(paste("The directory you specified was '",opt$outdir,"'.",sep=''),file=stderr())
    stop()
}

if (opt$verbose) {
    cat(paste("Read BAM list from ",opt$bamlist,"\n",sep=''),file=stdout())
}

counts = getBamCounts(bed.frame = panels.hg19, bam.files = bams, include.chr = TRUE)

if (opt$verbose) {
    cat(paste("Calculated counts\n",sep=''),file=stdout())
}

#####
# If desired, at this point you can save the counts, then have a second script
# which re-loads them. To do that, uncomment this part and split accordingly.

# save(counts,file="counts.rda")

# if (opt$verbose) {
    # cat(paste("Wrote counts to ",getwd(),"\n",sep=''),file=stdout())
# }

# # re-load the counts
# load('counts.rda')

# counts is an S4 object.
# you need to cast it to a data frame (for bin length)
# AND to a matrix (for reference.count)
# and for some reason you can't cast S4 directly to matrix, only via df
countdf = as.data.frame(counts)

countmat <- as.matrix(countdf[,5:dim(countdf)[2]]) # remove cols 1-5 metadata

countdf$chromosome <- gsub(pattern = 'chr', replacement = '', as.character(countdf$chromosome)) #removes the chr letters

# beta version: assume you want CNVs on all samples
for (i in 1:dim(countmat)[2]) {
    sample_name = colnames(countmat)[i]
	print(sample_name)
    reference_list = select.reference.set(test.counts = countmat[,i], 
       	reference.count = countmat[,-i],
        bin.length=(countdf$end-countdf$start)/1000,
        n.bins.reduced = 10000)
    reference_set = apply(
        X = as.matrix(countdf[, reference_list$reference.choice, drop = FALSE]), 
        MAR=1, FUN=sum)
    all_exons = new('ExomeDepth', test=countmat[,i], 
        reference=reference_set,
        formula = 'cbind(test,reference) ~ 1')


    # default expected.CNV.length is 50000 
    #all_exons <- CallCNVs(x = all_exons, transition.probability = opt$sensitivity,
    #    chromosome = countdf$chromosome, start=countdf$start,
    #    end=countdf$end, name=countdf$names, expected.CNV.length = opt$cnv_length)

    all_exons <- CallCNVs(x = all_exons, transition.probability = opt$sensitivity,
        chromosome = countdf$chromosome, start=countdf$start,
        end=countdf$end, name=countdf$exon, expected.CNV.length = opt$cnv_length)

    #all_exons = gsub("chr", "", all_exons@CNV.calls$chromosome)

   #print(nrow(all_exons@CNV.calls))
   # head(all_exons)

   if ((nrow(all_exons@CNV.calls)) > 0) {
	all_exons = AnnotateExtra(x = all_exons, reference.annotation = panels.hg19.GRanges, min.overlap = 0.02, column.name = 'panels.hg19')
   }

    # Specifying gene plots. Make sure it's hg19 coordinate. Sequence is chromosome #. xlim is coordinate range
    png(filename = paste(sample_name,"_MUTYH.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '1', xlim =c(45794900-10000,45806000+10000),count.threshold = 20, main = 'MUTYH gene', cex.lab = 0.8, with.gene = TRUE)   
    dev.off()   
    
    png(filename = paste(sample_name,"_EPCAM.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '2', xlim =c(47596000-10000,47614000+10000),count.threshold = 20, main = 'EPCAM gene', cex.lab = 0.8, with.gene = TRUE)   
    dev.off()   

    png(filename = paste(sample_name,"_MSH2.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '2', xlim =c(47630000-10000,47710100+10000),count.threshold = 20, main = 'MSH2 gene', cex.lab = 0.8, with.gene = TRUE)   
    dev.off()   

    png(filename = paste(sample_name,"_MSH6.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '2', xlim =c(48010357-10000,48034015+10000),count.threshold = 20, main = 'MSH6 gene', cex.lab = 0.8, with.gene = TRUE)   
    dev.off()   
    
    png(filename = paste(sample_name,"_MLH1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '3', xlim =c(37035000-10000,37092200+10000),count.threshold = 20, main = 'MLH1 gene', cex.lab = 0.8, with.gene = TRUE)   
    dev.off()   

    png(filename = paste(sample_name,"_APC.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '5', xlim =c(112000000-10000,112082000+10000), count.threshold = 20, main = 'APC gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_PMS2.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '7', xlim =c(6013000-10000,6048700+10000), count.threshold = 20, main = 'PMS2 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_PTEN.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '10', xlim =c(89624000-10000,89725200+10000), count.threshold = 20, main = 'PTEN gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_ATM.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '11', xlim =c(108098300-10000,108236300+10000), count.threshold = 20, main = 'ATM gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_BRCA2.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '13', xlim =c(32890598-10000,32972299+10000), count.threshold = 20, main = 'BRCA2 gene', cex.lab = 0.8, with.gene = TRUE)   
    dev.off()   

    png(filename = paste(sample_name,"_PALB2.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '16', xlim =c(23614700-10000,23652500+10000), count.threshold = 20, main = 'PALB2 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_TP53.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '17', xlim =c(7572900-10000,7580000+10000), count.threshold = 20, main = 'TP53 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_BRCA1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '17', xlim =c(41197695-10000,41276084+10000), count.threshold = 20, main = 'BRCA1 gene', cex.lab = 0.8, with.gene = TRUE)   
    dev.off()   
    
    png(filename = paste(sample_name,"_STK11.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '19', xlim =c(1206900-10000,1226600+10000), count.threshold = 20, main = 'STK11 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_BARD1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '2', xlim =c(215593400-10000,215675000+10000), count.threshold = 20, main = 'BARD1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_BRIP1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '17', xlim =c(59760000-10000,59940000+10000), count.threshold = 20, main = 'BRIP1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_CDH1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '16', xlim =c(68770000-10000,68870000+10000), count.threshold = 20, main = 'CDH1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_CHEK2.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '22', xlim =c(29083000-10000,29123100+10000), count.threshold = 20, main = 'CHEK2 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_FANCC.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '9', xlim =c(97863000-10000,98011000+10000), count.threshold = 20, main = 'FANCC gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_NBN.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '8', xlim =c(90948000-10000,90996800+10000), count.threshold = 20, main = 'NBN gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_RAD51C.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '17', xlim =c(56770000-10000,56811500+10000), count.threshold = 20, main = 'RAD51C gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_PTEN.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '10', xlim =c(89624200-10000,89725200+10000), count.threshold = 20, main = 'PTEN gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_BMPR1A.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '10', xlim =c(88635600-10000,88683500+10000), count.threshold = 20, main = 'BMPR1A gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_POLD1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '19', xlim =c(50902000-10000,50921200+10000), count.threshold = 20, main = 'POLD1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_POLE.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '12', xlim =c(133201100-10000,133264000+10000), count.threshold = 20, main = 'POLE gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_SMAD4.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '18', xlim =c(48573300-10000,48604900+10000), count.threshold = 20, main = 'SMAD4 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_GREM1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '15', xlim =c(32930000-10000,33026800+10000), count.threshold = 20, main = 'GREM1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_RAD51D.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '17', xlim =c(33428000-10000,33448000+10000), count.threshold = 20, main = 'RAD51D gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   
    
    png(filename = paste(sample_name,"_SDHC.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '1', xlim =c(161284165-5000,161334541+5000), count.threshold = 20, main = 'SDHC gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_SDHB.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '1', xlim =c(17345216-3500,17380527+3500), count.threshold = 20, main = 'SDHB gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_CDC73.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '1', xlim =c(193091087-13000,193223945+13000), count.threshold = 20, main = 'CDC73 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_MITF.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '3', xlim =c(69985750-3000,70017488+13000), count.threshold = 20, main = 'MITF gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_BAP1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '3', xlim =c(52435023-1000,52444024+2000), count.threshold = 20, main = 'BAP1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   
    
    png(filename = paste(sample_name,"_KIT.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '4', xlim =c(55524094-8000,55606881+8000), count.threshold = 20, main = 'KIT gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_PKD2.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '4', xlim =c(88928786-7000,88998929+7000), count.threshold = 20, main = 'PKD2 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_PDGFRA.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '4', xlim =c(55095459-7000,55164412+7000), count.threshold = 20, main = 'PDGFRA gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_CTNNA1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '5', xlim =c(138089113-18000,138270723+18000), count.threshold = 20, main = 'CTNNA1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_MSH3.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '5', xlim =c(79950470-22000,80172634+22000), count.threshold = 20, main = 'MSH3 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_SDHA.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '5', xlim =c(218434-3900,257197+3900), count.threshold = 20, main = 'SDHA gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_PKHD1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '6', xlim =c(51480096-47000,51952411+47000), count.threshold = 20, main = 'PKHD1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_MET.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '7', xlim =c(116312249-13000,116438431+13000), count.threshold = 20, main = 'MET gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_EGFR.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '7', xlim =c(55086709-19000,55279321+19000), count.threshold = 20, main = 'EGFR gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_FH.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '1', xlim =c(241660856-2000,241683085+2000), count.threshold = 20, main = 'FH gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_POT1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '7', xlim =c(124462439-10000,124569879+10000), count.threshold = 20, main = 'POT1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_RPS20.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '8', xlim =c(56979853-700,56987065+700), count.threshold = 20, main = 'RPS20 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_CDKN2A.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '9', xlim =c(21967750-700,21975132+700), count.threshold = 20, main = 'CDKN2A gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_PTCH1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '9', xlim =c(98205261-7000,98271548+7000), count.threshold = 20, main = 'PTCH1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_GALNT12.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '9', xlim =c(101569951-4000,101612363+4000), count.threshold = 20, main = 'GALNT12 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_SUFU.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '10', xlim =c(104263718-13000,104393214+13000), count.threshold = 20, main = 'SUFU gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_RET.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '10', xlim =c(43572516-5000,43625799+5000), count.threshold = 20, main = 'RET gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_SDHD.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '11', xlim =c(111957596-900,111966518+900), count.threshold = 20, main = 'SDHD gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_AIP.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '11', xlim =c(67250496-800,67258574+800), count.threshold = 20, main = 'AIP gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_SDHAF2.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '11', xlim =c(61197596-1700,61214239+1700), count.threshold = 20, main = 'SDHAF2 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_MEN1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '11', xlim =c(64570985-700,64578188+700), count.threshold = 20, main = 'MEN1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_CDK4.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '12', xlim =c(58141509-500,58146093+500), count.threshold = 20, main = 'CDK4 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_CDKN1B.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '12', xlim =c(12870203-500,12875317+500), count.threshold = 20, main = 'CDKN1B gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_RB1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '13', xlim =c(48877882-18000,49056026+18000), count.threshold = 20, main = 'RB1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_MLH3.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '14', xlim =c(75480466-3800,75518235+3800), count.threshold = 20, main = 'MLH3 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_MAX.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '14', xlim =c(65541844-2700,65569235+2700), count.threshold = 20, main = 'MAX gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_DICER1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '14', xlim =c(95552564-7000,95623759+7000), count.threshold = 20, main = 'DICER1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_PKD1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '16', xlim =c(2138708-4700,2185899+4700), count.threshold = 20, main = 'PKD1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_NTHL1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '16', xlim =c(2089820-800,2097835+800), count.threshold = 20, main = 'NTHL1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_PALB2.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '16', xlim =c(23614485-3800,23652631+3800), count.threshold = 20, main = 'PALB2 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_NF1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '17', xlim =c(29421944-20000,29704695+20000), count.threshold = 20, main = 'NF1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_PRKAR1A.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '17', xlim =c(66508519-2000,66529572+2000), count.threshold = 20, main = 'PRKAR1A gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_SMARCE1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '17', xlim =c(38781213-2000,38804070+2000), count.threshold = 20, main = 'SMARCE1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_AXIN2.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '17', xlim =c(63524680-3000,63557766+3000), count.threshold = 20, main = 'AXIN2 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_HOXB13.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '17', xlim =c(46802124-400,46806111+400), count.threshold = 20, main = 'HOXB13 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_RNF43.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '17', xlim =c(56429860-6500,56494943+6500), count.threshold = 20, main = 'RNF43 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_FLCN.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '17', xlim =c(17115525-2500,17140482+2500), count.threshold = 20, main = 'FLCN gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_SMARCA4.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '19', xlim =c(11071705-10000,+1173144+10000), count.threshold = 20, main = 'SMARCA4 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_NF2.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '22', xlim =c(29999544-9500,30094589+9500), count.threshold = 20, main = 'NF2 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_SMARCB1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '22', xlim =c(24129152-5000,24180196+5000), count.threshold = 20, main = 'SMARCB1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_LZTR1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '22', xlim =c(21336585-1700,21353321+1700), count.threshold = 20, main = 'LZTR1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_TSC1.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '9', xlim =c(135759173-1000,135797864+1000), count.threshold = 20, main = 'TSC1 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_TSC2.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '16', xlim =c(2096980-1000,2139621+1000), count.threshold = 20,amain = 'TSC2 gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    png(filename = paste(sample_name,"_VHL.png", sep=""), width = 8, height = 5, units = 'in', res = 300)
    plot(all_exons, sequence = '3', xlim =c(10173347-1000,10205364+1000), count.threshold = 20, main = 'VHL gene', cex.lab = 0.8, with.gene = TRUE)
    dev.off()   

    write.table(all_exons@CNV.calls, file=paste(sample_name,".csv",sep=''), 
        	sep=',', row.names=FALSE, col.names=TRUE, quote=FALSE)
    if (opt$verbose) {
        cat(paste("Wrote CNV calls for ",sample_name,"\n",sep=''),file=stdout())
    }
}

duration = format(Sys.time() - start_time)

if(opt$verbose) {
    cat(paste("Completed execution in ",duration,"\n",sep=''),file=stdout())
}
