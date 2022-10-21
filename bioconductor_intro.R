#_________________BASIC DATA TYPES in Bioconductor________________________#

 #1# Experimental data
 # e.g. gene1 has an expression level of 7.65

   #2# Metadata (more info about the experiment)
   #  e.g. sample1 is a 60 years old woman

     #3# Annotation (data we pick up most often from a big
     # corporate database that gives context to the experiment).
     # e.g., gene1 is highly conserved in yeast.



# A GRanges is a datastructure for storing genomic intervals.

# Every R-user dealing with genomic data needs to master this material.

# Many entities in genomics are intervals or sets of intervals (=integers):
# Promoters, Genes, SNPs, CpG Islands,...Sequencing reads; mapped and processed.


#GenomeInfoDb( Rank = 4 ) 

# packageDescription("Biostrings")$Version   -> Get back the pack's version.


##################1) IRanges - BASIC USAGE    (Rank= 5) ###################################### 

#Foundation of integer range manipulation in Bioconductor.


library(IRanges) 

ir1<- IRanges(start = c(1,3,5), end = c(3, 5,7))
ir1

ir2<- IRanges(start = c(1,3,5), width = 3)
ir2

start(ir1)

width(ir2)<-2
ir2

names(ir1) <- paste("A", 1:3, sep = "" )
ir1

dim(ir1) # = NULL because they all are vectors 
#even though they look a bit like matrices

length(ir1) # it works.

ir1[1]
ir1["A1"]

#concatenate two IRanges
c(ir1, ir2)
length(c(ir1, ir2))

#a specific type of IRanges: a normal IRange
 plotRanges <- function(x, xlim= x, main = deparse(substitute(x)), 
                       col = "black", sep = 0.5, ...) {
        height <- 1
        if(is(xlim, "Ranges"))
xlim <- c(min(start(xlim)), max(end(xlim)))
        bins <- disjointBins(IRanges(start(x), end(x) + 1))
        plot.new()
        plot.window(xlim, c(0, max(bins)*(height + sep)))
ybottom <- bins * (sep + height) - height
rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
title(main)
axis(1)
        }
par(mfrow = c(2,1))
ir <- IRanges(start = c(1,3,7,9), end = c(4,4,8,10))
plotRanges(ir) #plot the IRange
plotRanges(reduce(ir)) #get a normal IRanges = a minimal representation of the original IRanges as a set.
plotRanges(disjoin(ir)) #create a set of disjoint (non-overlapping) intervals

#Manipulating IRanges
#Take the original ranges and produce a single new range for each of the original ranges
resize(ir, width =1, fix = "start")
plotRanges(resize(ir, width =1, fix = "start"))
resize(ir, width =1, fix = "end")
resize(ir, width =1, fix = "center")#More useful: resize them from the center of the intervals
        
ir1<- IRanges( start = c(1, 3, 5), width = 1)
ir2<- IRanges( start = c(4, 5, 6), width = 1)       

union(ir1, ir2)        
plotRanges(union(ir1, ir2))     
reduce(c(ir1, ir2)) 
intersect(ir1, ir2)

ir1<- IRanges( start = c(1, 4, 8), end = c(3, 7, 10))
ir2<- IRanges( start = c(3, 4), width = 3)
ir1
ir2
 
ov <- findOverlaps(query= ir1, subject=ir2)
ov

plotRanges(ir1)
plotRanges(ir2)

countOverlaps(ir1, ir2) #a quick function that summarizes the number of overlaps  

nearest(ir1, ir2) #function ++++ in genomics (e.g. in case of a region of interest, find the nearest gene) # which of these IRanges in ir1 are closer to the ones in ir2





#########################2) GenomicRanges - GRanges   (Rank = 12) #################################

#Representation and manipulation of genomic intervals.
#It's very similar to IRanges with some additional stuff having to do with chromosomes.
#Chromosomes in GRanges are called seqnames.

BiocManager::install("GenomicRanges") 
library(GenomicRanges)


gr <- GRanges(
        seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
        strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        score = 1:10,
        GC = seq(1, 0, length=10))
gr


gr<- GRanges(
        seqnames = Rle(c("chr1", "chr2"), c(2, 1)),
        ranges = IRanges(start= c(10000,11100,20000), end = c(10300, 11500,20030)),
        strand = Rle(strand(c("+", "-", "+") ) ),
        score= c(10, 20, 15)
)
gr


gr= GRanges(seqnames = c("chr1"), 
            strand = c("+", "-", "+"), 
            ranges = IRanges(start = c(1,3,5), width=3))
gr

#strand can be "+" (forward : 5'--> 3'), "-" (reverse : 3' --> 5') or "*"  (unknown or there's an entity present on both strands))


#Splitting GRanges objects
sp <- split(gr, rep(1:2, each=5))
split(gr, names(gr))
split(gr, seqnames(gr))
split(gr, strand(gr))

#Combining GRanges objects
c(sp[[1]], sp[[2]]) # can use append()

#Subsetting GRanges objects
 #GRanges objects act like vectors of ranges
gr[2:3]
gr[2:3, "GC"]

# 2nd row replaced with the 1st row 
singles <- split(gr, names(gr))
grMod <- gr
grMod[2] <- singles[[1]]
head(grMod, n=3)

#repeat, reverse, or select specific portions of GRanges object.

rep(singles[[2]], times = 3)

rev(gr) # flip it downside-up

window(gr, 3,5)

gr[IRanges(start=c(2,7), end=c(3,9))] #(2-3 & 7-9)

g <- gr[1:3]
g

g<- append(g, singles[[10]]) # insert 

start(g)

range(g) # IRanges by chromosome


# Recover regions flanking the set of ranges represented
#by the GRanges object

flank(gr,5) 
#get a GRanges object containing the ranges 
#that include the 5 bases upstream of the ranges.

flank(gr, 5 , start=FALSE)
#get a GRanges object containing the ranges 
#that include the 5 bases downstream of the ranges.

shift(g, 5) #move the ranges by a specific number of base pairs
resize(g, 30) #extend the ranges by a specified width

reduce(g) #align the ranges and merge overlapping ranges to produce a simplified set
# same as range(g).

gaps(g) #the gaps between the ranges 

disjoin(g) #represent the GRanges object as a collection of non-overlapping ranges

coverage(g) #quantify the degree of overlap for all the ranges

g
g2 <- head(gr, n=2)

union(g, g2)
intersect(g, g2)
setdiff(g, g2)

# A list of available methods is discovered with
methods(class="GRanges")


promoters(gr)

seqinfo(gr) # get info about your chromosome
seqlengths(gr) = c("chr1"= 10)
seqinfo(gr)
seqlevels(gr)# see the chromosome's name

#gaps() gives us all the stuff on the X that aren't covered by a range.
gaps(gr)
seqlevels(gr)= c("chr1", "chr2")
seqnames(gr) <- c("chr1", "chr2", "chr1")
gr
sort(gr)
seqlevels(gr)= c("chr2", "chr1")
sort(gr)

genome(gr)= "hg19"
gr
seqinfo(gr)

gr2= gr
genome(gr2)= "hg18"

findOverlaps(gr, gr2)



######################## GenomicRanges - BASIC GRanges Usage ########################

#Dataframe (works like classical R dataframe), it allows many types of objects of arbitrary type.
#It allows storing IRanges inside them.

ir= IRanges(start = 1:3, width = 2)
df <- DataFrame(ir=ir, score=rnorm(3))
df

df[1,1]
df[2,1]
df[3,1]
df[1,2]
df[2,2]
df[3,2]

df$ir
df$score

##Let's return to GRanges.
gr

# add 1 or more metadata columns to the GRanges object.
values(gr)= DataFrame(score=rnorm(3), expression= 1:3)
gr

# N.B: The metadata is like a dataframe.
#      The info to the left of the | is not like a data.frame.
#      e.g.: we cannot do something like gr$seqnames 
#            However we can do sth like gr$score

#Access the metadata columns (as a DataFrame object )
values(gr)
mcols(gr) # do the same.

values(gr)$score
mcols(gr)$expression
gr$score
gr$score[1]
gr$expression

#GRanges extracted without metadata
granges(gr)

#Extract IRanges
ranges(gr)


# add another metadata column
gr$score2 = gr$score/2 
gr

# add again another metadata column
gr$promoter = c("a","b", "c") 

##The main work horse of the GRanges class and the IRanges ecosystem: is the findOverlaps(). 

gr2<- GRanges(seqnames = c("chr1", "chr2", "chr1"), ranges = IRanges(start=c(1,3,5), width=3), strand = "*")

findOverlaps(gr, gr2)
# see one-to-one overlaps between peaks and CpG islands.
#It returns a matrix showing which peak overlaps which CpG island.

countOverlaps(gr, gr2)
# tabulates the number of overlaps for each element in the query.

#ignore.strand = TRUE allows elements on the + strand  to overlap elements on the - strand.
findOverlaps(gr, gr2, ignore.strand = TRUE)

# ignore.strand When set to TRUE, the strand information is ignored in the overlap calculations.

# A BIT MORE OF EXPLANATIONS :

#For set operations: If set to TRUE, then the strand of x and y is set to "*" prior to any computation.

#For parallel set operations: If set to TRUE, the strand information is ignored in the computation and the result has the strand information of x. 



#select only ranges from the GRanges that overlap some other elements.
subsetByOverlaps(gr, gr2)
subsetByOverlaps(gr2, gr)
# returns a subsetted query containing only those 
#elements overlapped by at least one element in subject 

#Convert a classic-R dataframe into GRanges:

df <- data.frame(chr = "chr1", start= 1:3, end= 4:6, score = rnorm(3))
df
makeGRangesFromDataFrame(df) # 'score' column will be missed here
makeGRangesFromDataFrame(df, keep.extra.columns = T) #SOLUTION

###seqinfo
library("GenomicRanges")
gr<-GRanges(seqnames = c("chr1", "chr2"), 
            ranges = IRanges(start = 1:2, end = 4:5))      
gr

# remove a X
dropSeqlevels(gr, "chr2", pruning.mode = "coarse")

# keep a X
keepSeqlevels(gr, "chr2",pruning.mode = "coarse")


#keep only Standard X
gr<-GRanges(seqnames = c("chr1", "chrU2"), 
            ranges = IRanges(start = 1:2, end = 4:5))

keepStandardChromosomes(gr, pruning.mode = "coarse")


#Change the X' names
gr<-GRanges(seqnames = c("chr1", "chr2"), 
            ranges = IRanges(start = 1:2, end = 4:5))      

newStyle <- mapSeqlevels(seqlevels(gr), "NCBI")
newStyle 

gr <- renameSeqlevels(gr, newStyle)
gr


######################3)  AnnotationHub   (Rank=55) ##################################

# AnnotationHub is a wonderful resource for accessing genomic data or 
#querying a large collection of whole genome resources, including ENSEMBL, 
#UCSC, ENCODE, Broad Institute, KEGG, NIH Pathway Interaction Database, etc.


#This package is an interface to a lot of different online resources.
#The idea is to create a hub which is a local database of a lot of different
# online data.
#You take this local database, you query it, and you figure out which data
# do you want and then you go online and you retrieve them.
#This is a scalable way to access tons of different data.

# If working offline, add argument localHub=TRUE.
#   to work with a local, non-updated hub.
#    You gotta Only resources available that have previously been downloaded.


BiocManager::install(c("AnnotationHub"))
library(AnnotationHub)
Sys.setenv(TZ='GMT')
ah= AnnotationHub()
ah
unique(ah$dataprovider)
head(unique(ah$species))
ah[1]
unique(ah$rdataclass) 
unique(ah$genome) 
length(ah)

#To search this database you can use the subset() for example.
# If you work on human data:
ah<- subset(ah, species== "Homo sapiens")
ah

#query(): search a term in all the different components of the database
#e.g. if you wanna query the AnnotationHub to find a specific histonmodification.
query(ah, "H3K4me3")
#Find it within a specific cell line.
query(ah, c("H3K4me3", "Gm12878"))

ah2 <- display(ah)

# get an idea of the different files available
table(ah$sourcetype)


############# Usercase: AnnotationHub and GRanges ##############

library("AnnotationHub")
        ahub = AnnotationHub()



ahub = subset(ahub, species == "Homo sapiens")
qhs = query(ahub, c("H3K4me3", "Gm12878")) 
qhs


# H3K4me3 (H34K trimethylation) is a histone mark found inside
# promoters of genes and associated with gene expression activation.
BiocManager::install("rtracklayer")
library("rtracklayer")

gr1 = qhs[[2]]
gr1            # BroadPeaks
gr2=qhs[[4]]
gr2            # NarrowPeak

#In narrowpeaks (5-20 million reads) the regions bound are pretty much limited. 
#In broadpeaks (20-60 million reads) the regions can be much wider.

summary(width(gr1))
summary(width(gr2))

which(width(gr1)== 10)

gr1[6980] 

table(width(gr2))

peaks<- gr2 # narrowPeak

# ChIP-seq peaks are the sites where DNA-binding proteins interact with DNA.
#  These proteins could be: 
#   1) Transcription factors,
#   2) Chromatin-modifying enzymes,
#   3) Modified histones interacting with genomic DNA,
#   4) Components of the basal transcriptional machinery (e.g., RNA polymerase II).
                  #[beyzal]             


#We're gonna answer if these peaks are enriched in promoters.

#we need to get the promoters.
#Promoter = is a small interval of ~ 2KB
# around the trscrpt° start side (TSS) of genes.

#We're gonna get the TSS of genes.
#The way to do that is to use a transcript database object (TSDB). 

# But instead we're gonna get gene annotation using an AnnotationHub.
#We're gonna go for a specific type of gene annotation called RefSeq.

#RefSeq gives a highly curated (organized, selected) set of validated genes.

qhs = query(ahub, "RefSeq")
qhs

qhs$genome

cache(qhs[1]) <- NULL #remove the cached file to download that sh*t below.

genes= qhs[[1]]
genes

#RefSeq is very conservative.
#how many transcripts do we have per gene.
#how often I see a single transcript, 2 and so on and so forth.
#there's about 1k genes that has multiple transcripts.
table(table(genes$name))

# 2 transcripts from the same gene are gonna have the same name.

length(genes$name) # 50066 transcripts of different genes
length(unique(genes$name)) # only 46340 genes

# looking for genes with several transcripts
tail(table(genes$name), 100) 

# e.g.: how many transcripts do NR_110886 & NR_110999 have?
which(genes$name== "NR_110886") # 3 transcripts = 37341  37350   37355
which(genes$name== "NR_110999") # 2 transcripts =  21384  21414

prom <- promoters(genes)

#how wide are these promoters ?
table(width(prom)) # promoters' size = 2.2 kb

args(promoters) 

peaks
seqnames(peaks)

# we're gonna ask if these histone modification peaks
#are enriched in promoters.

# To do that, we need to know how often they overlap.
ov=findOverlaps(prom, peaks)
ov
# number of promoters that have a peak in them.
length(unique(queryHits(ov)))#9820

# number of peaks that have a promoter in them.
length(unique(subjectHits(ov))) #9850

# out of my peaks, how many overlap a promoter
length(subsetByOverlaps(peaks, prom, ignore.strand = TRUE))
#9850

length(subsetByOverlaps(prom, peaks, ignore.strand = TRUE))
#9820

# Percentage of peaks that overlap a promoter
length(subsetByOverlaps(peaks, prom, ignore.strand = TRUE))/ length(peaks)
#0.132268 ---->  13.2268 %

# H3K4me3 is a histone mark that doesn't just mark actual promoters. 
#It's also sometimes associated with enhancers and other regulatory elements.

# How many promoters have a peak in them.
length(subsetByOverlaps(prom,peaks, ignore.strand = TRUE))/ length(prom)
#0.1708539 ----> 17.08539 %

# Human genome = 3.2 billion base pairs = 3.10^9 base pairs

# How many Megabases do the peaks call ?
sum(width(reduce(peaks, ignore.strand = TRUE)))        #11190440
sum(width(reduce(peaks, ignore.strand = TRUE)))/ 10^6 
# 11.19044 Mb

# How many Megabases do the promoters call ?
sum(width(reduce(prom, ignore.strand = TRUE)))/ 10^6
# 108.3146 Mb

# How big is the overlap ?
sum(width(intersect(peaks, prom, ignore.strand = TRUE)))/ 10^6
#1.356499 Mb 

# N.B: a base can be either in a peak and in a promoter or can be in a peak and not 
#in a promoter and vice versa.



# We're gonna quantify the relationship between the things.
#look for is there some kinda significant enrichment here.
inOut <- matrix(0, ncol= 2, nrow = 2)

colnames(inOut) <- c("in", "out")
rownames(inOut) <- c("in", "out")

inOut[1,1] <- sum(width(intersect(peaks, prom, ignore.strand = TRUE)))
inOut[1,2] <- sum(width(setdiff(peaks, prom, ignore.strand = TRUE)))
inOut[2,1] <- sum(width(setdiff(prom, peaks , ignore.strand = TRUE)))
inOut
?setdiff

colSums(inOut)
rowSums(inOut)
inOut[2,2] = 3*10^9 - sum(inOut)
inOut

#setdiff():  is used to find the elements which are in the first Object
#          but not in the second ( asymmetric difference ).

fisher.test(inOut)$statistic

oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2] )
oddsRatio #3.716635

#oddsRatio is a number between 0 & Inf.
# if  > 1  ----> enrichment

#In this case that means that the overlap between the peaks and the promoters is like
# 3 fold more enriched than we would expect.


################################4) Biostrings   (Rank= 10) ############################

#This is a package that contains functionality for representing and
# manipulating biological strings (DNA strings, RNA strings, or amino acids) 
# and biodata.

library("Biostrings")

dna1<- DNAString("ACGT-T")
dna1

dna2<- DNAStringSet(c("ACG", "ACGT", "ACGTT"))
dna2

nchar(dna1 )  # count the number of nucltd. in a string.
nchar(dna2 )

IUPAC_CODE_MAP

dna1[2:4]
dna2[1:2]
dna2[[3]]

(dna2[[1]])[3] # Access inside the sequence 

unlist(dna2) # paste all in 1 ( convert to a DNAString object )

letter(dna1, 3) # extract a letter by its position

matchPattern(DNAString("T"),dna1)  

subseq(dna1, end=3)

subseq(dna1, start=3)

subseq(dna1, 5,5 ) # extract the substring located in 5th position


names(dna2) <- paste0("sequence", 1:3)

width(dna2)

sort(dna2, decreasing=T) #reverse the sequences' order 
rev(dna2) #do the same thing

rev(dna1) #Reverse clean the sequence: Biological reversion !!!
reverse(dna1) # same.

reverse(dna2) # Biological reversion for DNAStringSet !

complement(dna1)
reverseComplement(dna1)

complement(dna2) 
reverseComplement(dna2)

translate(dna2)
translate(dna1)

AMINO_ACID_CODE


# The frequency of occurrence of each nucleotide 
alphabetFrequency(dna1)
alphabetFrequency(dna2)

# Get the GC content for ex.
letterFrequency(dna2, letters = "GC")

dinucleotideFrequency(dna2)
trinucleotideFrequency(dna2)
oligonucleotideFrequency(dna2, 4)

consensusMatrix(dna2) # tells how many strings have a specific nucleotide
# of that persistance.
# This is useful when you build up precision weight matrcies or
# motifs for transcription factor binding sites.


#__Biostrings__Matching functionality for finding subsequences & other sequences:

library(BSgenome.Scerevisiae.UCSC.sacCer2)
dnaseq<- DNAString("ACGTACGT")
dnaseq

#If you have millions of short reads it's recommended 
# to use a dedicated short read aligner such as bowtie and MAQ.

#However you gotta be able to match a small set of sequences
# to a given genome. 

# Matching a single string to a single string.
     matchPattern(dnaseq, Scerevisiae$chrV)


# Counting how many matches we have. 
    countPattern(dnaseq, Scerevisiae$chrV)

#Most of the time,
# We're not interested in matching up against a single chromosome.
# But rather in matching up against a set of chromosomes. 
# For this we use:
              vmatchPattern(dnaseq, Scerevisiae)
# N.B: Contrary to matchPattern(), vmatchPattern() do search both
#the + and - strands.

vcountPattern(dnaseq, Scerevisiae) # searching both strands across all chromosomes

sum(vcountPattern(dnaseq, Scerevisiae)["count"]) # 170 counts= 170 sequences

strand(vmatchPattern(dnaseq,Scerevisiae))

#The DNA sequence we are searching for
# turns out to be its own reverse complement.
dnaseq== reverseComplement(dnaseq)

complement(dnaseq)
reverseComplement(dnaseq)  

# matchPDict() takes a set of sequences such as
# short reads of the same length and it builds a dictionary
# on them, and then it matches them against the full genome.


# matchPDict(), matchpattern() and vmatchPattern() all
# allow a small set of mismatches and indels,
# and is a very fast and efficient way of searching the 
# genome for a small set of sequences.

# matchPWM() (precision weight matrix or sequence logo) 
# allows us to search the genome for example
# for binding site for a given transcription factor. 


# pairwiseAlignment() implements a classic pair-wise Alignment algorithm,
# either a global pair-wise Alignment or a local pair-wise Alignment
# like a Smith Waterman or Needleman Bush, a type algorithm.
# And it allows to map millions of reads against a short sequence
# such as a gene.

# It turns out that these local and global alignment route genes using dynamic
# programming are impossible to use when you map up against the entire genome.

# But they are still very useful even if you have a million of reads,
# As long as you align them up to a very small section of the genome (e.g.,gene). 


# bowtie2 and bwa : The most commonly used programs for mapping in bioinformatics.

# Mapping: is the process of comparing each one of the reads with
#the reference genome.

# A READ: is a sequence fragment produced by a sequencing machine,
#they can be long, short, paired, single, etc

# OR MORE FORMALLY:

# A READ is a raw sequence that comes off a sequencing machine. A read may
#consist of multiple segments. 
#For sequencing data,reads are indexed by the order in which they are sequenced.


# trimLRPatterns()      
# trimming off specific patterns on the left and the right of a DNA string set.
# The use case here is trimming off sequence adapters.
# trimLRPatterns has a very rich set of functionality allowing indels
# and mismatches in the sequence adapters. 

# Sequence Alignment Softwares (ELAND, MAQ, and Bowtie) perform the bulk
#of short read mappings to a target genome.  



##--------------------------#5) BSgenome  (Rank=51/2140) #----------------------------

# the full genome sequences for living organisms..

BiocManager::install(c("BSgenome"))
library("BSgenome")

available.genomes() # list all the genomes you can download directly from 
# the Bioconductor website


# load a genome:
BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer2")
library(BSgenome.Scerevisiae.UCSC.sacCer2)

Scerevisiae
BSgenome.Scerevisiae.UCSC.sacCer2 # Same thing.

isCircular(Scerevisiae) 

seqnames(Scerevisiae)
seqlengths(Scerevisiae)

# 2micron: is a plasmid of Saccharomyces cerevisiae that is a relatively small 
#multi-copy selfish DNA element that resides in the yeast nucleus at a 
#copy number of 40-60 per haploid cell.

Scerevisiae$chrI 

# Compute the GC content
letterFrequency(Scerevisiae$chrI, "GC")

# GC content as percentage  :  90411 / 230208 = 0.39 * 100
letterFrequency(Scerevisiae$chrI, "GC", as.prob= TRUE)*100 # ~ 40%

# GC content as percentage
letterFrequency(Scerevisiae$chrI, "AT", as.prob= TRUE)*100 # ~ 60%


dinucleotideFrequency(Scerevisiae$chrI)
# % of each dinucltd.
sum(dinucleotideFrequency(Scerevisiae$chrI, as.prob= TRUE )*100 )

alphabetFrequency(Scerevisiae$chrI)
# %A %T %C %G :
alphabetFrequency(Scerevisiae$chrI, as.prob= TRUE)*100


#Computing the GC content for the entire genome.
# we can use lapply()

# But there's a new type of apply called bsapply().
# To run bsapply you need to set up BSParams.

param <- new("BSParams", X = Scerevisiae, FUN = letterFrequency) 
  
    bsapply(param, "GC")

unlist(bsapply(param, "GC"))

#Gettting the GC content across the entire genome
sum(unlist(bsapply(param, "GC")) / sum(seqlengths(Scerevisiae)))
# ~ 0.38 = 38 %

#Checking the GC % for each single X
unlist(bsapply(param, "GC", as.prob = TRUE))

#______________________BSgenome - Views_____________________________________

library(BSgenome.Scerevisiae.UCSC.sacCer2)

dnaseq<- DNAString("ACGTACGT")
vi<- matchPattern(dnaseq, Scerevisiae$chrI)
vi

#underneath the hood, a views object is the same as an IRanges object.
ranges(vi)

#Get out the entire seq
Scerevisiae$chrI[57932:57939]

alphabetFrequency(vi)

# We can Shift the view by 10 paces for ex.
shift(vi, 10)
shift(vi, -10)

#Getting the coordinate of each IRanges object
# associated with its DNAstring.
# That allows us to easily represent many 
# subsequences, such as promoters or exons.

gr <- vmatchPattern(dnaseq, Scerevisiae)
gr
vi2 <- Views(Scerevisiae, gr)
vi2


#Computing the GC content of the
#promoters in the yeast genome,

#in order to get the promoters,
# let's load up AnnotationHub:

library(AnnotationHub)
ahub <- AnnotationHub()
ahub

# Let's do a query on our AnnotationHub for
#this specific version of the genome.

query(ahub, c("Homo sapiens"))


qh <- query(ahub, c("sacCer2", "genes"))
qh

genes<- ahub[["AH7048"]]
genes

prom <- promoters(genes)

#Lots of error & warnings 'cause G ranges 
# complain  when it gets genome
# indices that are less than zero.
# There are some of these sequences that are
# right at the boundary of the genome. 

prom 

#PROMOTER= 2.2 kB
# 2k paces upstream and 200 paces downstream
# of the transcription stop site.

#Trim off the promoters by cutting off
# anything outside the sequence length
# of the genome.

prom <- trim(promoters(genes))
#Some promoters are gonna be less than 2.2 kb long. 
prom
 
promViews <- Views(Scerevisiae, prom)
promViews

gcProm<- letterFrequency(promViews, "GC",as.prob = T)
                     
plot(density(gcProm)) # a density plot of the GC content in the promoters
abline(v = 0.38)

#The promoters don't seem to have different GC content
# than the rest of genome.

#That may be surprising, but
#remember that the yeast genome is quite small.
#It's filled up with coding sequences and
#regulatory regions unlike the human genome. 


##-------- * Run-length encoding (RLE) - GRanges -------------------
#Rle is a way of representing very long vectors (compression)

library(GenomicRanges)

rl <- Rle(c(1, 1, 1,1,1,1,2,2,2,2,2,4,4,2))
rl

runLength(rl)
runValue(rl)
#This is a way of compressing a vector. 'Cause
# We have taken 14 numbers and compressed them
# into 8 numbers.

# Convert Rle to a normal vector.
as.numeric(rl)

#Coverage is defined as the number of sample 
# nucleotide bases sequence
# aligned to a specific
# locus in a reference genome.

#E.G. computing average signal
# across a set of pre-specified
# genomic regions.

# Coverage of signal across inside genome
 rl
 
# a set of genomic regions 
 
ir <- IRanges(start = c(2,8), width = 4)

# we wanna for ex. compute the average 
# signal across a set of pre-specified genomic
# regions of a coverage vector. 

aggregate(rl, ir, FUN = mean)

#[1] 1 2

vec<- as.numeric(rl)
mean(vec[2:5])  #[1] 1
mean(vec[8:11]) #[1] 2


#Construct a coverage vector out of IRanges.

ir <- IRanges(start = 1:5, width = 3)
ir

rl


#
# If you wanna visualize that sh*t, run this f**king code !
plotRanges <- function(x, xlim= x, main = deparse(substitute(x)), 
                       col = "black", sep = 0.5, ...) {
        height <- 1
        if(is(xlim, "Ranges"))
                xlim <- c(min(start(xlim)), max(end(xlim)))
        bins <- disjointBins(IRanges(start(x), end(x) + 1))
        plot.new()
        plot.window(xlim, c(0, max(bins)*(height + sep)))
        ybottom <- bins * (sep + height) - height
        rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
        title(main)
        axis(1)
}
par(mfrow = c(2,1))
ir <- IRanges(start = 1:5, width = 3)
plotRanges(ir)
#

#So I have these overlapping intervals,
# what's the coverage of them?
# I mean, how often for integer number 1 ?
# How many of these ranges is a part of integer number 2 ?
# How many of these ranges is a part of integer number 3 ?
# And so forth.

 
# coverage() computes the coverage across a set of ranges.
coverage(ir) # gives a Rle 

#counts, for each integer, how many ranges overlap the integer."

# For each position in the space underlying a set of ranges, counts
#the number of ranges that cover it."

as.numeric(coverage(ir))


# to figure out areas where the vector 
# is big use slice() which creates views objects.
# N.B: 1st arg. by default is lower.

rl

slice(rl, upper =1)

slice(rl, lower =2)

slice(rl, lower= 4)


vi <- Views(rl, IRanges(c(2,8), width = 2))
vi

# We can compute functions on these views, 
mean(vi)
# Recall, that's the same as what we did with the
# aggregate() above.


# Let's do the same thing with GenomicRanges.

library(GenomicRanges)

gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:10, width=3))
rl <- coverage(gr)       
rl        
# we have a RleList 'cause we gotta  1 Rle for each X.

vi <- Views(rl, GRanges("chr1", ranges= IRanges(3,7)))
vi



##--------* Lists - GRanges ------------------------------------------

# A GRangesList is a list where each element is a GRanges object.

# A use case for GRanges: 
 #a single GRanges list could be describing the exons of a transcript.

#e.g., we gotta gene which has multiple transcripts and each transcript
# has multiple exons.
# It's very natural to think about a GRanges list to encode that gene.

gr1 <- GRanges(seqnames = "chr1", ranges = IRanges(start= 1:4, width = 3))
gr2 <- GRanges(seqnames = "chr2", ranges = IRanges(start= 1:4, width = 3))
gr1
gr2

gL <- GRangesList(gr1 = gr1, gr2 = gr2)
gL

options(warn=-1) #Get rid of the f--king warnings.

gL[[1]]
gL$gr1

start(gL)
seqnames(gL)

elementNROWS(gL) # how long are each element.
sapply(gL, length) # same.

gL

shift(gL, 10) # Shift everything by +10. 

findOverlaps(gL, gr2)
findOverlaps(gL, gL)


#
#From the vignette
#

gr1 <- GRanges(
        seqnames = "chr2",
        ranges = IRanges(103, 106),
        strand = "+",
        score = 5L, GC = 0.45)

gr2 <- GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges(c(107, 113), width = 3),
        strand = c("+", "-"),
        score = 3:4, GC = c(0.3, 0.5))

gr3 <- GRangesList("txA" = gr1, "txB" = gr2)

gr3

elementNROWS(gr3) # a faster alternative to: sapply(gr3, NROW)
                                        
isEmpty(gr3)


# Metadata at the list level differs from the one 
#that is at the level of the individual GRanges objects.

mcols(gr3) <- c("Transcript A","Transcript B")

# List-level metadata:
mcols(gr3)
values(gr3)

mcols(gr1)
mcols(gr2)

mcols(unlist(gr3)) # Element-level metadata 

coverage(gr3)

sapply(gr3, length)

grl2 <- shift(gr3, 10)

names(grl2) <- c("shiftTxA", "shiftTxB")
names(grl2)







##---------#6) GenomicFeatures (Rank=31/2140) #----------------------------

# It retrieves and manages transcript-related features 
#from the UCSC and BioMart databases.

# This package contains functionality and provides support 
# for Transcript Data Base or TxDb objects.

# This database stores information about the transcripts, 
#such as the location of a transcript and its exons, 
#its CDS start and end position, and the transcript 
#sequence.  (CDS= coding sequence start)

#The CDS contains start & stop codon and does not include any UTR and introns. 
#Therefore, CDS does not correspond to the actual mRNA sequence.

# The package consists of a set of tools and methods for
#making and manipulating transcript centric annotations.
#It allows downloading the genomic locations of the transcripts,
#exons and CDS of a given organism, from either
#the UCSC Genome Browser or a BioMart database.


# The package is useful for ChIP-chip, ChIP-seq, and RNA-seq analyses.

BiocManager::install("GenomicFeatures")
library(GenomicFeatures)

#Housekeeping genes are typically constitutive genes that are required
#for the maintenance of basic cellular function.
# They are examples of regions in a genome that tend to be highly
#conserved and evolve slower than other genes such as the tissue-specific genes .

 #Get hold of a transcript Data Base
BiocManager:: install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

txdb   # create a txdb object  

seqlevels(txdb)

seqlevels(txdb) <- "chr1"      # only set Chromosome 1 to be active

genes(txdb)      # Only genes on chromosome 1 

seqlevels(txdb) <- seqlevels0(txdb)     # reset back to the original seqlevels 
                      
     
# Let's examining a very specific loci
# Picking up a GRanges that selects a small region of the X 1
# Then looking at all the exons, the transcripts and the genes 
# that are on this particular loci.
gr <- GRanges(seqnames = "chr1", strand = "+", 
              ranges = IRanges( start = 11874, end = 14409))

genes(txdb)

# let's look at which genes actually
#overlap this little genomic range I have. 
subsetByOverlaps(genes(txdb), gr)

# Actually it turns out that there's another gene that 
#overlaps the same loci on the reverse strand.
subsetByOverlaps(genes(txdb), gr, ignore.strand = T)

# let's look at which transcripts overlap our genomic
# range :

subsetByOverlaps(transcripts(txdb), gr)

# These three different transcripts are different 'cause 
# they've different exons.

# let's look at what exons do we have here:

subsetByOverlaps(exons(txdb), gr)

# These 6 exons are combined in different ways to form these
# 3 different transcripts.


# let's figure out how these 6 exons are combined together to 
#form the 3 transcripts

subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)


# CDS :
# Not each transcript can have a coding sequence.
# Having a CDS implies by the name that it gets
#translated into a protein.

# A given transcript may have multiple ORF (open reading
# frame)
# In some databases and some organisms, you look for all ORFs
# inside the transcript, and you just pick the longest ORF
# as the potential ORF or CDS for that particular transcript.
# That turns out to not always give you the correct result.


subsetByOverlaps(cds(txdb), gr)
#These are not necessarily ORF 'cause if it was
# an ORF we should have some splicing information here
#These are the coding sequence in the different
#splice transcript intersected with the exons.
# This is a weird thing.
#I find it weird that it's called cds in the package.
#But we are really interested in is cdsBy.

subsetByOverlaps(cdsBy(txdb, by = "tx"), gr) 
#Now we get the actual CDS for the three different transcripts.

#To see which transcript contains a real CDS:
subset(transcriptLengths(txdb, with.cds_len = T), gene_id == "100287102") 

# transcriptLength ()
#gives us transcript lengths not in terms of their pre-MRA, 
#but in terms of the spliced RNA. 

# Here we see the three different transcripts.
#We see that each of them have three exons.
#We see they've different lengths 'cause they're spliced differently.
#And we see that the transcripts 1 & 3 have a CDS length of zero, 
#which basically says that there's no coding sequence inside those
#two transcripts.And the transcript 2 has a cds length of 402. 

# useful to have tx and gene together. 
transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))

# Retrieve a subset of the transcripts available such as those on the strand of chromosome.
transcripts(txdb, filter=list(tx_chrom = "chr1", tx_strand = "+"))


# The intronsByTranscript, fiveUTRsByTranscript and threeUTRsByTranscript are convenience 
#functions that provide behavior equivalent to the grouping functions, but in prespecified form.
#These functions return a GRangesList object grouped by transcript for introns, 5' UTR's, and 3' UTR's,


# Translate transcripts into proteins :

#   (that requires an appropriate BSgenome package)
 # for example:

library(BSgenome.Hsapiens.UCSC.hg19)

tx_seqs1 <- extractTranscriptSeqs(Hsapiens, TxDb.Hsapiens.UCSC.hg19.knownGene,
                                  use.names=TRUE)

suppressWarnings(translate(tx_seqs1))
# Note: But of course this is not a meaningful translation, because the call to extractTranscriptSeqs will have extracted all the transcribed regions of the genome regardless of whether or not they are translated.

# Translate only the coding regions :

cds_seqs <- extractTranscriptSeqs(Hsapiens,
                                  cdsBy(txdb, by="tx", use.names=TRUE))
translate(cds_seqs)


##---------#7) rtracklayer -Data import (Rank=24/2140) #-------------------------

# The rtracklayer package is an interface (or layer ) between R and genome
#browsers. 

#Its main purpose is the visualization of genomic annotation tracks.

# And it's a bidirectional interface: we can both take
# data in R and put it onto the genome browser, and
# we can take data from the genome browser and put it into R. 

#This section is about the latter: the import of certain type of data into R.

library(rtracklayer)
library(AnnotationHub)
ahub <- AnnotationHub()
ahub

#The different types of data that we can connect to.
table(ahub$rdataclass)

ahub.bw <- subset(ahub, rdataclass == "BigWigFile" & species == "Homo sapiens")
bw <- ahub.bw[[1]] 
bw

# bigWig files contain alignment data to the genome. It is similar to BED file which contains location of the genome where the read aligned to. 


#Read in  a part of the BigWigFile

gr.chr22 <- import(bw, which= GRanges("chr22", ranges = IRanges(1, 10^8)))

rle.chr22 <- import(bw, which = GRanges("chr22", ranges = IRanges(1, 10^8)), as= "Rle")
  
rle.chr22$chr22 


#liftOver() is a tool from UCSC allowing to convert between
# different genome versions.

# e.g., Convert genome position from one genome assembly to another genome assembly.

# it requires a ChainFile which contains info about converting
# one specific genome to another.

ahub.chain <- subset(ahub, rdataclass == "ChainFile")
ahub.chain

# Keep only human data and cut out the data about other species
ahub.chain <- subset(ahub.chain, species== "Homo sapiens")
ahub.chain

#let's get the lift-over chain from hg19 to hg18

chain <- query(ahub.chain, c("hg18", "hg19"))[[1]]
chain

gr.chr22 <- import(bw, which=GRanges("chr22", ranges 
                                     = IRanges(1, 10^8)))
gr.chr22 

#We wanna convert that from hg19 to hg18
gr.hg18 <- liftOver(gr.chr22, chain)
gr.hg18

class(gr.hg18)
length(gr.hg18)
length(gr.chr22)


##---------# * ExpressionSet data 'Class' #-------------------------

# (!!!) ExpressionSet is generally used for array-based experiments, 
#     where the rows are features (=genes).

# This is an old and very important data container in Bioconductor.

# ExpressionSet is an ExperimentData Package,
# which contains a data set.
# This kind of package re-package existing data
#from publications into a nice easy-to-use data set.

BiocManager::install("ALL") 

# ALL = Data "Acute Lymphoblastic Leukemia"

library(ALL)

data(ALL) 
ALL  # feature = gene.

experimentData(ALL) #Get info about the experiment.

?ALL #Get additional info.

# The most important thing in the data is 
#the expression measure = the quantificatiion
#of the expression of different genes.

# That can be accessed using exprs().
head(exprs(ALL))

class(ALL)
dim(ALL)    

exprs(ALL)[1:4, 1:4]

sampleNames(ALL) 
featureNames(ALL) 

# Get the phenotype data (data about the samples)
pData(ALL)  #Biobase (Rank=6/2140)

colnames(pData(ALL)) # all covariates we have.

dim(pData(ALL))  
 
summary( pData(ALL)) # useful !

str(pData(ALL) )

# Acess a specific covariate
pData(ALL)$sex
ALL$sex    # same.

# Find out how many F & M you have.
table(pData(ALL)$sex)  

range(ALL$age, na.rm = T)

# An expression set satisfies two dimensional subsetting.
#The first dimension gives you features, and
#the second dimension gives you samples

exprs(ALL)[1, 1:5]
exprs(ALL[1, 1:5]) # same but with feature's name.

#Get info about the genes.
featureData(ALL)
# the author didn't provide any data.

ids= featureNames(ALL)[1:5]
ids

# Let's find out which genes are these 5 identifiers

# We need an annotation package called "hgu95av2"
#the information is inside it.

BiocManager::install("hgu95av2.db")
library(hgu95av2.db)

#it's possible to map these ids using various
#maps inside this package.

as.list(hgu95av2ENTREZID)

?hgu95av2ENTREZID
# hgu95av2ENTREZID is an R object that provides mappings
#between manufacturer identifiers and Entrez Gene 
#identifiers.

as.list(hgu95av2ENTREZID)["1001_at"] 

#    $`1001_at`
#    [1] "7075"


# From https://www.ncbi.nlm.nih.gov/gene :

# Official Symbol: TIE1
# Official Full Name tyrosine kinase with
#immunoglobulin like and EGF like domains 1



##-----#8) SummarizedExperiment -data container   Rank= 13/2140 #-------------------------

# (!!!) SummarizedExperiment is generally used for sequencing-based experiments,
#    where the rows are GenomicRanges. 

# This is a more modern or a new version of the ExpressionSet. 

# We will use a data package called airway to illustrate it.

library(GenomicRanges)

BiocManager::install("GenomicFiles") #Rank= 188/2140
# This package provides parallel computations distributed 'by file' or 'by range.

library(GenomicFiles)         

BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)        


BiocManager::install("airway") #Experiment Package: Rank = 4/410
library(airway)
data(airway) 
airway
??airway

# lots of things are much the same as what we find in ExpressionSet.
# But here the syntax and the output are slightly different.

colData(airway) 

airway$cell # get back the different cell lines.

metadata(airway) # Data about the Experiment 

colnames(airway) 
rownames(airway)

# Get back the data expression
assayNames(airway)
assay(airway, "counts")[1:4, 1:4] 

length(rowRanges(airway))

dim(assay(airway, "counts"))

# How many exons per gene ?
rowRanges(airway) 

elementNROWS(rowRanges(airway))

# The 1st gene listed has 17 exons located on chromosome X.

# How many exons in all ?
sum(elementNROWS(rowRanges(airway)))

# We've got almost: 750,000 exons in 64,000 genes.

start(airway) # Get back the start coordinate of each exon.
 
# Let's get the genes inside "airway" that overlap a specific genomic
#interval. 

gr<- GRanges("1", range = IRanges(start = 1, end = 10^7))

subsetByOverlaps(airway, gr) 
#There's 329 genes inside this 10 Mbp genomic interval.



#**** Creating a SummarizedExperiment object ________________

# simulate an RNA-seq read counts table
nrows <- 200
ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)

# create gene locations
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, T),
                     feature_id=paste0("gene", 1:200))

# create table for the columns
colData <- DataFrame(timepoint=1:6,
                     row.names=LETTERS[1:6])


# create SummarizedExperiment object:

se = SummarizedExperiment( assays= list(counts=counts),
                        rowRanges= rowRanges, colData=colData )

se

rowRanges(se)

assay(se, "counts")








##--------------------- #9) GEOquery  (Rank= 48/2140) #----------------------------

# GEOquery is the bridge between GEO and BioConductor. 
#It allows interfacing with NCBI GEO.

# GEO is a public functional genomics data repository.

#The NCBI Gene Expression Omnibus (GEO) is the largest public
#repository of high-throughput microarray experimental data.

#  GEO (The Gene Expression Omnibus (GEO) database) is
# a widely-used repository of high-throughput gene expression 
# and other functional genomics data sets.

#And unlike what the name suggests,
#it contains other data types than gene-expression data.
#There's also a lot of epigenetic data in there. 


#  N.B. (!!!)
#Illumina array studies (Affymetrix arrays) have no celfiles.


# GEO2R : is an interactive web tool that allows users to compare two
#or more groups of Samples in a GEO Series in order to identify genes 
#that are differentially expressed across experimental conditions.

BiocManager::install("GEOquery")
library("GEOquery")


eList <- getGEO("GSE11675") 

# Chronic myelogenous leukemia hematopoietic stem cells
# Expression profiling by array

# It downloads the data into a list 'cause there might be one 
#data set associated with RNA sequencing, and
#another data set associated with ChIP sequencing...for ex. 

length(eList)
names(eList)
eData <- eList[[1]]
eData

exprs(eData )
names(pData(eData)) 
# Fore more about this dataset:
# SHOWUP_Bioconductor.R / 06 Showing Up: Differential Expression Analysis with limma 


  #***Acessing the RAW DATA from GEO

# NCBI GEO accepts raw data such as .CEL files, .CDF files, images, etc.
# getGEOSuppFiles() download all the raw data associated with a given accession number. 
# the function will create a directory in the current working directory to store the raw data.

eList2 <- getGEOSuppFiles("GSE11675")
eList2 




##--------------------- #10) biomaRt  (Rank= 20/2140) #----------------------------

# This package is an Interface to BioMart databases:
#(Ensembl, Uniprot, HGNC, Reactome, InterPro, HapMap...)

# The BioMart project provides free software and data services
#to the international scientific community in order to foster
#scientific collaboration and facilitate the scientific discovery
#process.
#The project adheres to the open source philosophy that 
#promotes collaboration and code reuse.


BiocManager::install("biomaRt")
library("biomaRt")

# In biomaRt lingo: a Mart = a DataBase

listMarts()

# Establish a connection to a database 
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart

 # -> This database provides access to gene annotation information.


# Get the datasets inside the selected database 
listDatasets(mart) 

ensembl<- useDataset("hsapiens_gene_ensembl", mart)
ensembl


#Let's get back the gene name associated with these probe IDs
values <- c("202763_at", "209310_s_at", "207500_at")


# getBM() is the workhorse function for querying in biomaRt.

# attributes = whatever you wanna retrieve from the database.

# filters is a way of selecting what you wanna retrieve.

# "affy_hg_u133_plus_2" is the microarray's name (Affymetrix micro array)


getBM(attributes = c("ensembl_gene_id", "affy_hg_u133_plus_2"), 
   filters = "affy_hg_u133_plus_2", values = values, mart = ensembl)

# We gotta a dataframe wih 2 columns ( gene ID & Affymetrix ID )
# we need both of these to link up the probe ID to the gene.

attributes <- listAttributes(ensembl)
nrow(attributes)
head(attributes)


filters <- listFilters(ensembl)
head(filters)
nrow(filters)

which(filters$name == "affy_hg_u133_plus_2")

# Hint:
library("dplyr")
filters %>% arrange(name) %>% head

#  Setting up a query is like going through all of this and
# understanding it, and it can get rather daunting. 

attributePages(ensembl)
unique(attributes$page) #Same.

attributes <- listAttributes(ensembl, page = "feature_page")
attributes


#---------------# *) R S4 Classes  #-------------------------------

# Bioconductor Project is a heavy user of the S4 system,
# Normal R packages available from CRAN usually don't use
#this S4 system very much.

# in R, there're 3 different ways of doing Object-oriented
#Programming:

# [ Object-oriented Programming (OOP) is a computer programming
#model that organizes software design
#around data, or objects, rather than functions and logic. ]

# There's the classic system, also known as S3.
#There's a more advanced system, known as S4, and
#there's a recently introduced system known as S5, or 
#reference classes. 

# In the S4 system, S4 classes and S4 methods are two separate things.
#You could have a package that uses S4 classes without using S4 methods,
#and vice versa. 

# S4 classes is a way of representing complicated data structures.
#We have used this a lot already.
#We have seen things like expression sets and summarized experiments.
#And we have seen the wealth of information and
#links between different data structures like they have.
#And this is possible through the S4 system.

# The S4 classes have really proven their worth in the bioconductor 
#project, and we are better off using them. 


# We're interested in getting an example dataset.
library(ALL)
library(GenomicRanges)

# We can make any object into any kind of class:

df <- data.frame(y = rnorm(10), x = rnorm(10))
df

lm.object <- lm(y ~ x, data = df)
lm.object

class(lm.object) # "lm"

names(lm.object)

xx <- list(a = letters[1:3], b= rnorm(4))
xx

# Look here, we're able to turn this list into a class of object "lm"
#without any kind of error message or warning. 
class(xx) <- "lm"
xx


# This is really useful for complicated data structures,
#because as a programmer you know what you deliver to people,
#and what you get from people.

#When you say you have an object,
#you know it's an expression set, you know exactly what that means.

#But let's look a little at an expression set. 

data(ALL)
ALL

class(ALL) # ExpressionSet Class defined in the Biobase package.

isS4(ALL) # check if ALL is a S4 object.

# How do I get help on that class:
class?ExpressionSet
?"ExpressionSet-class" # Same.

# Get back the definition of the class.
getClass("ExpressionSet")

slot(ALL, "annotation")
annotation(ALL)


#---------------# *) R S4 methods #-------------------------------

# a method is a function that allows you to
#run different sets of code based on different values of
#the argument. 
 
library(GenomicRanges)

showMethods("findOverlaps")

##--------------#11) ShortRead  (Rank= 71/2140) #----------------------------

# It provides functionality for working with FASTQ files from high-
#throughput sequence analysis.

# This package contains two different functionalities for
#reading in "raw sequencing reads" (in the form of FASTQ files).
# And it contains also a functionality for reading in "aligned reads".

# Although a lot of its functionality is now outdated and
#you should be using the Rsamtools or GenomicAlignments package instead. 

# The ShortRead package was one of the first Bioconductor packages to deal with
#low-level analysis of high-throughput sequencing data. Some of its functionality 
#has now been superseded by other packages,
#but there is still relevant functionality left.

#* The FASTQ file format is the standard way of representing raw (unaligned) 
#next generation sequencing reads, particular for the Illumina platform.

# FASTQ files are often used for basic quality assessment.

#The most important parameters to check for raw sequencing data quality are
# the base quality, the nucleotide distribution, GC content distribution and the duplication rate.

# c.f. "FastqCleaner" Bioconductor package :
# A Shiny Application for Quality Control, Filtering and Trimming of FASTQ Files.

BiocManager::install("FastqCleaner") #Already installed.

#In the GEO page, the link to the SRA (Sequence Read Archive) project
#contains the raw FASTQ files.



BiocManager::install("ShortRead")
library(ShortRead)

data(package="ShortRead")  # No data sets found

# system.file() provides the path to the files
# It outputs the full path to the file that is installed with the data package.

# Creating the path to the FASTQ file associated with the package
fastqDir <- system.file("extdata", "E-MTAB-1147", package = "ShortRead")

# The directory path (dirpath) of the FASTQ file
fastqPath <- list.files(fastqDir, pattern = ".fastq.gz$", full= T)[1]

# The FASTQ files are read by readFastq() which produces an object of class ShortReadQ
reads<- readFastq(fastqPath)
reads

args(readFastq)

# We can also read the file by FastqFile ().
# we need this sometimes in order to read in smaller chunks (sections) 
#of the file.
fqFile <- FastqFile(fastqPath)
fqFile
reads <- readFastq(fqFile)

#Acess the info inside of FastqFile
sread(reads)[1:2]
sread(reads)[1:5]
#ShortReadQ class is very similar to a DNAStringSet 


# Get back the quality values (in some weird letters).
quality(reads)[1:2]
# in order to do any real computation on qualities
#you have to convert them into integers.

id(reads)[1:2] 

# To do that, you gotta coerce the quality scores into a matrix.
as(quality(reads), "matrix")[1:2,1:10]
# Here we show only the 10 pages of the 1st two reads.
# 'H' gets translated to 39 and 'I' to 40.

# Each letter is matched to an integer between 0 and 40. This matching
#is known as the "encoding" of the quality scores and there has been different
#ways to do this encoding.

# These numbers are supposed to related to the probability that the reported base 
#is different from the template fragment (ie. a sequence error). One should be aware
#that this probabilistic interpretation is not always true; methods such as "quality-remapping"
#helps to ensure this.



##--------------#12) Rsamtools  (Rank= 22/2140) #----------------------------

# SAM: (Sequence Alignment/Map) 

#The main functionality of the package is support for reading BAM files.

# This package provides an interface to BAM files.

# It is an interface to the widely-used samtools/htslib library.

# It contains functionality for reading and examining 
#aligned reads in the BAM format.

# BAM files are produced by samtools and other software,
#and represent a flexible format for storing 'short'reads'
#aligned to reference genomes. 

# BAM files typically contain sequence and
#base qualities, and alignment coordinates and quality measures.

# Samtools library is laid up and replaced by HGS or high sequencing library.
# This is a set of collections for dealing with files in SAM and BAM format.
# SAM is a text file format for representing aligned reads.

# BAM is a binary version of SAM.

# For all practical purposes,everybody should exclusively work with BAM files.
# 'cause they are smaller, much faster and  more convenient. 

BiocManager::install("Rsamtools")
library(Rsamtools)
data(package="Rsamtools") # no data sets found

# How to read a BAM file.

# First set up a BamFile object:
#A pointer to the file is created by the BamFile() constructor:
bamPath <- system.file("extdata", "ex1.bam", package="Rsamtools")
bamFile <- BamFile(bamPath)
bamFile
Bam
# Some high-level information can be accessed here.
# Which sequences or chromosomes are the reads aligned to.
seqinfo(bamFile)

# The possibility of selecting reads using:
# ScanBamParams()

seqlevels(bamFile)

# Read in a Bam File using scanBam()
aln <- scanBam(bamFile) 
aln[[1]][1] # show the 1st element
length(aln) 
class(aln) 
 # We get back a list of length 1; this is because scanBam() can return
#output from multiple genomic regions, and here we have only one (everything).

names(aln) # >>> NUL
 # The reason why it doesn't have a name,
#is 'cause that this is the entire files.

# So let's subset it to the first element in the list .
aln <- aln[[1]]

# That gives us a list and we show the information from the first element.
names(aln)
aln$qname

# These names are the names used in the BAM specification.
#Here is a quick list of some important ones:

# qname: The name of the read.
# rname: The name of the chromosome / sequence / contig it was aligned to.
# strand: The strand of the alignment.
# pos: The coordinate of the left-most part of the alignment.
# qwidth: The length of the read.
# mapq: The mapping quality of the alignment.
# seq: The actual sequence of the alignment.
# qual: The quality string of the alignment.
# cigar: The CIGAR string (below).
# flag: The flag (below).


# Let's have a look at the first element:
lapply(aln, function(xx) xx[1]) 
 aln[[1]][1]   # same
                       

# BAM files can be extremely big and it is often necessary to read
#in parts of the file.

#Let's talk about reading in BAM files in small chunks.

# You can do this in different ways:
#Read a set number of records (alignments).
#Only read alignments satisfying certain criteria.
#Only read alignments in certain genomic regions.


# Let's start with the first of this.

#By specifying yieldSize when you use BamFile(),
yieldSize(bamFile) <- 1

#every invocation of scanBam() will only read yieldSize number
#of alignments.

#You can then invoke scanBam() again to get the next set of
#alignments; 
scanBam(bamFile)[[1]]$seq

#This requires you to open() the file first.
open(bamFile)

#(otherwise you will keep read the same alignments).

## Cleanup;
close(bamFile)
yieldSize(bamFile) <- NA

# The other two ways of reading in parts of a BAM file is to use ScanBamParams(), 
#specifically the what and which arguments.

gr <- GRanges(seqnames = "seq2",
              ranges = IRanges(start = c(100, 1000), end = c(1500,2000)))


 # ScanBamParam is like a class that encodes our query to the BAM file.
  # It's similar to bsapply().

params <- ScanBamParam(which = gr, what = scanBamWhat())

  # which = "which genomic region we're gonna read"
 
  # what = which pieces of the BAM files we wanna read.
        #scanBamWhat() tells the function to read everything. 

  
#This useful, because sometimes you are only interesting in the position,
#you don't wanna retrieve the actual sequence, or the actual base qualities of
#the read, because they're rather big they take up a lot of space. 

aln <- scanBam(bamFile, param = params)

names(aln)

head(aln[[1]]$pos)

# Notice how the pos is less than what is specified in the which argument; this is because 
#the alignments overlap the which argument.


#Other functionality from Rsamtools:
#BamViews: reading multiple BAM files

# Instead of reading a single file, it is possible to construct something called a BamViews,
#a link to multiple files. 

# In many ways, it has the same Views functionality as other views. 

#A quick example should suffice, first we read everything;

bamView <- BamViews(bamPath)
aln <- scanBam(bamView)
names(aln)

BamV
# We can also set bamRanges() on the BamViews to specify that only certain ranges are read;
#this is similar to setting a which argument to ScanBamParams().

bamRanges(bamView) <- gr
aln <- scanBam(bamView)
names(aln) 

names(aln[[1]])

#quick summary of what is in the file
quickBamFlagSummary(bamFile)

#* Sometimes, all you wanna do is count. use countBam() instead of scanBam().

countBam(bamFile) # nb of nucleotides, nb of records...



##------------# 13) oligo (Rank= 133/2140) #----------------------------

# Analyze oligonucleotide arrays (expression/SNP/tiling/exon) at probe-level.

# Preprocessing tools for oligonucleotide arrays

#  This package has supplanted an earlier package
#  called Affy package.

# It's used for handling Affymetrix and nimbleGen
# microarrays (especially gene expression, exon expression
#and SNP arrays).


devtools::install_url("https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.2.tar.gz")

BiocManager::install("oligo")
library(oligo)

# We will use the dataset deposited as a GEO accession number :
#  GSE38792. 

# In this dataset, the experimenters profiled fat
#biopsies from two different conditions: 10 patients with 
#obstructive sleep apnea (OSA) and 8 healthy controls.

# The profiling was done using the Affymetrix Human Gene ST 1.0 array.

# First we need to get the raw data; this will be a set of binary 
#files in CEL format.

# There will be one file per sample.
#The CEL files are accessible as supplementary information 
#from GEO; we get the files using GEOquery.

# CEL file = is a binary format in which microarrays and raw data 
#are stored for Affymetrix .

# These are not the standard files we get from GEO.

# These are always submitted as supplementary files (suppFiles) from GEO.

# So the way we get these CEL files is we know what the system number
#of the experiment is, and then we get the supplemental files. 

library(GEOquery)

getGEOSuppFiles("GSE38792")

list.files("GSE38792") 

# TAR archive is like ZIP archive and it contains CEL files.

untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")
list.files("GSE38792/CEL")


# we construct a vector of filenames.
# And we feed it to read.celfiles().

celfiles <- list.files("GSE38792/CEL", full = T)
celfiles

# read.celfiles() : Convenience function for reading in 
#many files at once into a data container:
rawData <- read.celfiles(celfiles)
rawData 
# This data is in the form of an GeneFeatureSet;
#which is an ExpressionSet-like container.

# There's more than 1 million different probes.

getClass("GeneFeatureSet")
# There's some stuff that looks very much like an
#expression Set (eSet).

# Look at some data
exprs(rawData)[1:4,1:3]

# We can see that we get integers that are pretty large in this case here,
#between roughly 200 and 10,000, and that indicates that
#this expression data is raw intensity measurement from the scanner.

# the unit of measure is integer measurements on a 16 bit scanner,

#   A micro-array scanner is typically a 16 bit scanner,
# which means when you scan a probe you get a number between 0 and
# 2^16 (two to the sixteenth ) which is 65536.

# We can confirm that by looking at the highest value inside
# the expression of data.
max(exprs(rawData))

  #  N.B: Research has shown that this is not a very good scale to work
#on with microarray data.

#Usually we log transform these datas here.
#When we log transform them, If we use a log with base 2,
#we get a number basically between zero and 16. 


# Note the large number of features in this dataset, more than 1 million. 
#Because of the manufacturing technology, Affymetrix can only make very
#short oligos (around 25bp) but can make them cheaply and at high quality.
#The short oligos means that the binding specificity of the oligo is not 
#very good. To compensate for this, Affymetrix uses a design where a gene
#is being measured by many different probes simultaneously; this is called
#a probeset. As part of the preprocessing step for Affymetrix arrays,
#the measurements for all probes in a probeset needs to be combined into
#one expression measure.


# Let us clean up the phenotype information for rawData.
filename <- sampleNames(rawData)
filename

pData(rawData)$filename <- filename
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames
sampleNames(rawData) <- sampleNames


#/a/learning/curve/
args(sub)

v <- c("abcdeabfgh")
v
sub("ab", " ", v) # replaces only the first occurrence.
gsub("ab", " ", v) # replaces all occurences found.
sub(".*^", "XXX", v)# put sth. on the beggining
sub(".*", "XXX", v) # replaces the whole string

vv<- c("ijklm_pqrst")
vv
sub(".*_","", vv) # replaces what's on the left of "_"
sub("_.*","", vv) # replaces what's on the right of "_"
sub("rst$", "", vv) # replaces what's on the end.

vec<-"I was born on July 23,1996"
#eliminating the numeric values
gsub( "[0-9]*" , "" , vec )
gsub('[0-6]*','',vec)
#eliminating letters
gsub( "[a-c]*" , "" , vec) # abc
gsub( "[an]*" , "" , vec)  # a & n

#_


# Create a  which is gonna tell us
#the different experimental groups:

pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)),
                               "OSA", "Control")
pData(rawData)


# Normalization

#Let us look at the probe intensities across the samples, using
#boxplot().

boxplot(rawData, target="core")

# Boxplots are great for comparing many samples because 
#it is easy to display many box plots side by side.
#We see there is a large difference 
# in both location and spread
#between samples.

# There are three samples with very low intensities;
#almost all probes have intensities less than 7
#on the log2 scale.

# The Y axis in this box plot here is on the log scale.
#So, a difference of 1 or 2 is a quite massive difference.

# This is an extremely low
#intensity. 

# One hypothesis is that nothing was really hybridized for
#these samples here (the array hybridization failed
# for these arrays).

# Really understanding that, or assessing that, or deciding
#that really it's gonna require a lot more of exploratory data 
#analysis for this particular experiment.


# We have a dataset like this, most often we wanna start off by 
#normalizing it.

# A classic and powerful method for preprocessing Affymetrix gene 
#expression microarrays is the RMA method.

# [Data preprocessing = transforming raw data into an 
#understandable format].

#Robust Multichip Average preprocessing methodology. This strategy allows background subtraction, quantile normalization and summarization 
normData <- rma(rawData)
normData

# Data Normalization aims at reducing data redundancy and inconsistency in order to achieve data integrity.
# Data Normalization is a best practice for processing and utilizing stored data.
# Data Inconsistency : When the same data exists in different formats in multiple tables.
# Data redundancy occurs when the same piece of data exists in multiple places.

# Note how normData has on the order of 33k features which 
#is closer to the number of genes in the human genome.

# We can check the performance of RMA by looking at
#boxplots again.
boxplot(normData)

# Here, it is important to remember that the first set of boxplots
#is at the probe level (~1M probes) whereas the second set of boxplots
#is at the probeset level (~33k probesets).

#The data is now ready for differential expression analysis.
exprs(normData)[1:3]
length(exprs(normData) )

# Performing DEA :
 pData(normData)

 comp <- factor(pData(normData)$group)
 comp
 
 design <- model.matrix(~   comp)
 design
 
 fit <- lmFit( normData , design)
 
 mean(fit$Amean) #7
 
 hist(fit$Amean)
 plotSA(fit)
 
 keep <- (fit$Amean >= 7)
 
 fit2 <- eBayes(fit[keep,], trend=TRUE)
 
 topTable(fit2)
 
 fit <- eBayes(fit)          
 topTable(fit) 
 
 allg <- topTable(fit, number= 33297 )
 
 length(which ( allg $adj.P.Val < 0.05 )) # 0  It's strange here !!
 
 
 
##------------# 14) limma (Rank= 15/2140) #----------------------------

# "limma" stands for "linear models for microarray data".

# It allows performing Transcriptomic Data Analysis :
#     Differential gene expression for:
#                              Microarray and RNA-seq Data.

# Linear models is a general class of models that is being used for
# continuous data.

BiocManager::install("limma")
library(limma)

 
 
BiocManager::install("leukemiasEset")
library(leukemiasEset)
data(leukemiasEset) 

?leukemiasEset # info about this dataset.

# raw data  is available at GEO, under accession number "GSE13159".

head(exprs( leukemiasEset))


summary(pData( leukemiasEset ))

sampleNames(leukemiasEset ) 

experimentData(leukemiasEset)

unique(leukemiasEset$LeukemiaTypeFullName)
table(leukemiasEset$LeukemiaTypeFullName)
table(leukemiasEset$LeukemiaType)

# Our question is :
# which genes are differentially expressed between the "ALL" type and the "NoL" type.

# Let's Subsetting and cleaning up the data:
ourData <- leukemiasEset[, leukemiasEset$LeukemiaType %in% c("ALL", "NoL")]

leukemiasEset$LeukemiaType

ourData$LeukemiaType

ourData$LeukemiaType = factor(ourData$LeukemiaType)
ourData$LeukemiaType

# -> By convention, the control level should be the first level.

# A linear model:

# Now we do a standard limma model fit.

# The design of the experiment (The comparison of interest) is 
#defined through a design matrix:

#N.B.The approaches for two groups extend easily to any number of groups.

design <- model.matrix(~ ourData$LeukemiaType)
head(design)
# Each row corresponds to an array in the experiment.
# Each column corresponds to a coefficient that is used
# to describe the RNA sources in the experiment.


       #To Remove the intercept, we use : 
        model.matrix(~ ourData$LeukemiaType - 1) # or
        model.matrix(~ ourData$LeukemiaType + 0)
        model.matrix(~ 0 + ourData$LeukemiaType)# same

#___________________________________________________________________
 # Using the second approach: extract the difference as a contrast:
designII <- model.matrix(~ ourData$LeukemiaType + 0)

colnames(designII) <- c("ALL", "NoL")
designII

fit <- lmFit(ourData, designII)

contrast.matrix <- makeContrasts(NoL-ALL, levels = designII)
                                  
fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

topTable(fit2, adjust="BH")
#____________________________________________________________________


# Why does the intercept column in model.matrix 
# replace the first level? ------------------->
# In statistics, when we have a factor variable with k levels,
#we need to convert it to k-1 indicator variables.
# We choose one level as the baseline, and then have
#an indicator variable for each of the remaining levels.
#  In this case we have one factor variable (LeukemiaType)
# with two levels (ALL and NoL).
# The baseline level is "ALL".
# --> The design matrix shows an intercept & one indicator variable.

# N.B. If there's only one factor variable, we can specify not
# to produce the intercept and keep the two dummy variables produced
# (i.e., an indicator variable for each single level).
# This encodes actually the same information.


# lmFit() fits a linear model to each gene separately.

fit <- lmFit(ourData, design)

# Empirical Bayes (eBayes) techniques to borrow information across features (genes) :
fit <- eBayes(fit) 

topTable(fit, adjust="BH" ,sort.by= "logFC") #listing the 10 top differentially expressed genes.

# Why is there only one list of different expressed genes?
# Because this design only has two groups so there is only one comparison
#that's worth making: which genes differs between the two groups ?

#* We have here a log2 fold-change (logFC) from Nol to ALL (Nol - ALL).
# # LogFC is the log difference between our groups.

#So here we notice a down-regulation of these 10 genes in cancer (ALL). 

# Multidimensional scaling plot of distances between gene expression profiles:
plotMDS(ourData, dim.plot = c( 3,4) , col=  c( rep(" blue", 12 ), rep(" red", 12)), 
           labels = c( rep(" ALL", 12 ), rep(" NoL", 12)))
          # That makes sense man !!

# A volcano plot displays log2 fold changes on the x-axis versus
#a measure of statistical significance on the y-axis.
volcanoplot(fit, coef = 2, highlight = 10,  names = fit$genes$NAME) 
#In a volcano plot, the most upregulated genes are towards the right, the most downregulated genes are towards the left, and the most statistically significant genes are towards the top.

args(volcanoplot )


toptable <- topTable(fit, n = Inf)
names(toptable)

#*We can also use the "EnhancedVolcano" Bioconductor package to get a very nice visualization.

#Hints to INTERPRET MDS PLOT:
# Objects that are closer together on the plot are more alike than those further apart.
# The distances among each pair of points correlates as best as possible to the dissimilarity between those two samples.
# The values on the two axes tell you nothing about the variables for a given sample - the plot is just a two dimensional space to arrange the points.
# That is not a statistical test, just a visualisation that helps you seek patterns.

  # How many genes are DE here ? 
Total <- topTable(fit, number=20172 )
length(which ( Total$adj.P.Val < 0.05 )) #  Only 5437 ( 27 %) 


# we can compute logFC by hand to confirm:

topTable(fit, n=1) # get back the 1st gene

genename <- rownames(topTable(fit, n=1))
genename

#Compute the mean of expression inside each group

typeMean <- tapply(exprs(ourData)[genename,], ourData$LeukemiaType, mean)
typeMean

# Here's the logFC :
typeMean["NoL"] - typeMean["ALL"]


# Testing that the two groups have the same expression
#level is done by testing whether the second parameter 
# (equal to the difference in expression between
#the two groups) is equal to zero.

design2 <- model.matrix(~ ourData$LeukemiaType - 1)
head(design2)

colnames(design2) <- c("ALL", "NoL")

fit2 <- lmFit(ourData, design2)
contrast.matrix <- makeContrasts("ALL-NoL", levels=design2)
contrast.matrix
# Contrast is gonna test whether ALL - NoL = 0
fit2C <- contrasts.fit(fit2, contrast.matrix)
fit2C <- eBayes(fit2C)
topTable(fit2C) 
# ALL - NoL =/= 0
#We've got the reverse of what we got for NoL-ALL. 


##------------# 15) minfi  (Rank= 118/2140) #----------------------------

# minfi: Analyzing Illumina Infinium DNA methylation arrays

# Handling data from DNA methylation microarrays.
#(Illumina Infinium methylation arrays).

# The human genome contains a 30k CpG islands,
#and  30 million CpG dinucleotides.

# IDAT files = Raw scanning files from the Illumina platform. 
#They're available as supplementary material.

BiocManager::install("minfi")
library(minfi)
library(GEOquery)

# We're gonna play with a file (450K data set) from a study
#that have tried to understand how whether or
#not DNA methylation changes are associated with acute mania.

# Illumina 450k array has been recently (in 2016) replaced by
#The Infinium MethylationEPIC Array which covers over 850k CpGs
#(the double number of Illumina 450K) distributed genome-wide.

getGEOSuppFiles("GSE68777") 
# This is a f--king big download, it takes ages !

list.files("GSE68777") 

# There are 2 supplementary files.

# We are interested in this file called "GSE6877_RAW.tar",
#that is an archive containing the idat files. 

#  Unpack the archive:
untar("GSE68777/GSE68777_RAW.tar", exdir ="GSE68777/idat")

# inside the directory that has been created, we have a list of
# idat files.

head(list.files("GSE68777/idat", pattern = "idat"))

rgSet <- read.450k.exp("GSE68777/idat") #(!!!!!!!!!!) deprecated f()
rgSet 
# *Defunct functions in minfi package:

#read.450k:       >>>> Use read.metharray
#read.450k.sheet: >>>> Use read.metharray.sheet
#read.450k.exp:   >>>> Use read.metharray.exp

### >>>>>>>>>> USE THAT instead:
rgSet <- read.metharray.exp("GSE68777/idat")
rgSet

pData(rgSet) # No pData

head(sampleNames(rgSet))

# Download the original processed GEO data to get pData.

geoMat <- getGEO("GSE68777")#That's another f**king big  download!

# Get the phenotype data from the original GEO matrix.
pD.all <- pData(geoMat[[1]])

table(pD.all )

# We're only interested in 4 columns. 
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1",
                 "characteristics_ch1.2")]

head(pD) # Here's the 4 columns 

# Cleaning the pData up :

names(pD)[c(3,4)] <- c("group", "sex")
pD$group <- sub("^diagnosis:","", pD$group)
pD$sex <- sub("^Sex:", "", pD$sex)
head(pD)

sampleNames(rgSet) <- sub(".*_5", "5", sampleNames(rgSet) )
rownames(pD) <- pD$title
pD <- pD[sampleNames(rgSet),]
sampleNames(rgSet)


head(sampleNames(rgSet))
head(pD)

pData(rgSet) <- pD
rgSet

# Normalizing the microarray raw data
grSet <- preprocessQuantile(rgSet)   # It takes ages !!
# This f() does a couple of things :
#Normalizing and Mapping to genome.


# Which normalization should we apply ????????
#===> RULE OF THUMB______________________________________

#If there exist global biological methylation differences 
#between your samples, e.g.,
#Dataset with cancer and normal samples, 
#or a dataset with different tissues/cell types,use 
# --> the preprocessFunnorm().

#If you do not expect global differences between your 
#samples, e.g.,
#a blood dataset, or one-tissue dataset, use
#--> the preprocessQuantile().



# This is a process of assigning each probe to a given
#location in genome where that particular CpG is located.

grSet

# Get the location of different CpGs
granges(grSet)

# *CpG shore: just by the CpG island (inside an area close
#                                      to the island)
# *CpG shelf: a little bit further away the CpG shore.
# *CpG open sea: far away from CpG island.

head(getIslandStatus(grSet))   

#  Beta value = % methylation :
getBeta(grSet)[1:3, 1:3]


# Now we can perform some Differential Methylation Analyses 
 #using limma package.







##------------# 16) edgeR & DESeq2   (Rank= 25 ; 29 /2140) #----------------------------
  
    # Count-based RNA-seq analysis


library(airway)
data(airway)

??airway   # RNA-Sequencing based-experiment

# GEO --> GSE52778 

assay(airway, "counts")[1:3, 1:3]

# this is unormalized data

# The variable of interest for this particular 
#comparison is the "dex" covariate (dexamethasone)
#which tells us whether the sample belongs
#to the treatment group or
#the untreatment group. 

colData(airway) #Phenotype Data

airway$dex

# The reference level is treatment.
#So if we start fitting models to the data here,
#we're gonna find differences from treated to untreated.
#That's a little unnatural,

# so we use the relevel function and
#we reset the reference level

airway$dex <- relevel(airway$dex, "untrt")
airway$dex

granges(airway)

# Construct count matrix

BiocManager::install("edgeR")  #(similar to limma)
library(edgeR)
# (Rank: 25/2140)
# Differential expression analysis 
#of RNA-seq expression profiles
#with biological replication.


# Unlike limma, we can't run edgeR directly
#on summarized experiments.

# we gotta take our data from the summarized 
#experiment and convert it into a limma class.


#______ # NOTE __________________________________________________________________________

        #When it comes down to fitting common types of models:

#1) limma is useful for "continuous data" such as microarray data.

#2) edgeR / DESeq2 are useful for count data ("discrete data")
#such as high-throughput sequencing (RNA-seq).
 
 # --> N.B:   
# count data can also be transformed into continuous data and then limma can be used.

#~~~~~~~~~~~~ A bit more of statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1/*** CONTINUOUS DATA : (~ Normal distribution = Gaussian distribution : The most well-known continuous distribution)
# It can take on any numeric value.
#There are an infinite number of possible values between any two values.
# We often measure a continuous variable on a scale (height, weight, temperature )

# 2/*** COUNT DATA :  (~ Poisson distribution )
# It counts the frequency of occurrences as integers,
# and can have only non-negative integers. 
# (e.g., 0, 1, 2, etc.)

# 3/ BINARY DATA : (~ Binomial distribution )
# It  only counts two states
# (e.g., Yes/No, Pass/Fail, Accept/Reject)

#*N.B. -->  2/ & 3/ are discrete data with two discrete Distributions.

 # Discrete data ->
 # If there are only a set array of possible outcomes 
 # (e.g. only zero or one, or only integers), then the data are discrete.

 # A key difference:
 # Discrete data is countable while continuous is measurable.

 #*** Negative Binomial Distribution ( Pascal Distribution ) 
 # = Discrete Distribution
 # with a regular binomial distribution, you're looking at the number of successes.
 #*** With a negative binomial distribution, it's the number of failures that counts.

#_________________________________________________________________________________________________


# Give the DGEList a count matrix and a group variable
#which is the group of interest.

dge <- DGEList(counts = assay(airway, "counts"),
                   group = airway$dex) 


#  Take phenodata into the DGEList: 

dge$samples <- merge(dge$samples, 
                     as.data.frame(colData(airway)),
                     by= 0)

#  Take info about the genes into the DGEList: 

dge$genes <- data.frame(name = names(rowRanges(airway)),
  stringsAsFactors = F)   

head(dge$genes)                  


# Normalization of the library sizes
# and getting the effective of library size.
dge <- calcNormFactors(dge)                         
dge
                     
# Estimate the dispersion of the variability
#of the data.

# estimate a common dispersion
dge <- estimateGLMCommonDisp(dgel) 
# estimate a checkwise dispersion
dge <- estimateGLMTagwiseDisp(dge)

# Model fitting
design <- model.matrix(~dge$samples$group)
head(design)

fit <- glmFit(dge, design) # fit glm for each gene.

# figure out for our coefficients, our contrast of interest,
#which genes are most significant there.

#To do that, we can do this in various ways, good example here,
#I do it using a test, using the GLM test.

#The coefficient equal 2 means that I'm testing the second coefficient
#in the sign matrix which is equal to the coefficient
# representing differences between
#the treatment and the untreatment group.

lrt <- glmLRT(fit, coef = 2)

# Get a list of my top tags, the ones that are most
#differentially expressed between cases treatment and non-treatment

topTags(lrt)
# This is how we do it with edgeR.


#Let's do it with DESeq2.

BiocManager::install("DESeq2")
library(DESeq2)

#Rank: 29/2140
# Differential gene expression analysis 
#based on the negative binomial distribution


# like edgeR, in order to use DESeq2,
#we gotta put our data into a DESeq2 container.

# Unlike edgeR, it's quite easy to do because we can convert it
#directly from a summarized experiment. 

dds <- DESeqDataSet(airway, design = ~ dex)
# Notice that we need to store the design matrix inside the data object.

dds <- DESeq(dds)
#That performs lots of edgeR commands wrapped into one.

res <- results(dds)
res 








#----------------------------------____________________________---------
# plyranges ; Rank= 1931/2140
library(plyranges)

# A dplyr-like interface for interacting with the common 
#Bioconductor classes Ranges and GenomicRanges. 

# all 'chr21' ranges
 
i %>% 
        filter(seqnames == 'chr12')

# filter by one region (stringent, i.e., fully contained in region)
data_GR %>% 
        filter(seqnames == 'chr21', start >= 2e7L, end <= 3e7L)

# filter by one region (permissive, i.e., any overlap with region)
data_GR %>% 
        filter_by_overlaps(as('chr21:20000000-30000000', 'GRanges'))

# filter by multiple regions
data_GR %>% 
        filter_by_overlaps(as(c('chr1:1-10000000', 'chr21:20000000-30000000'), 'GRanges'))


#______________________________
BiocManager::install("genomation") #Rank 246/ 2140
# Summary, annotation and visualization of genomic data.
library(genomation ) 










