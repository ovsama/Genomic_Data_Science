#----------------------### 01 SHOWING UP : Biostrings exercises###----------------------------------------------------


###_https://www.r-exercises.com/2017/01/22/bioinformatics-lessons-in-r-1-biostrings-exercises/

               

myDNA <- c("ATGTTGCATTCATTAATTAAGAACGACCCAATACA")
myDNA

mySequencing <- c("CTGATTT-GATGGTC-NAT")
mySequencing

length(mySequencing)

myDNASeq <- DNAString("CTGATTT-GATGGTC-NAT")
myDNASeq

length(myDNASeq)

class(myDNA)

class(myDNASeq)


pastedDNA<- paste(myDNA, myDNASeq)
class(pastedDNA)
 
DNAString(pastedDNA) # it doesn't work

#But here :
pastedDNA<- paste(myDNA, myDNASeq, sep = "")
DNAString(pastedDNA) # it works

myDNASeq
complement(myDNASeq)
reverse(myDNASeq) 
reverseComplement(myDNASeq) 

alphabetFrequency(myDNASeq) 

translate(myDNASeq)


#__https://www.r-exercises.com/2017/05/28/manipulate-biological-data-using-biostrings-package-exercises-part-2/


# Sequence Alignment Techniques:

#Sequence Alignment is comparing the similarity between the sequences to check how much the DNA,RNA or Protein are related to each other.

 #There are three types of Sequence Alignment:

#1. Global Alignment
#2. Local Alignment
#3. Overlap Alignment

lsf.str("package:Biostrings")

myd<- DNAString("CTGATTTCA")
complement(myd)

myr<- RNAString("CUGAUUUCA")
myr

reverseComplement(myd)


translate(RNAString("CUGAUUUCA"))

abc <- translate(RNAString("CUGAUUUCA"))

#Print the three letter of each amino acid:
Amino<-AMINO_ACID_CODE[strsplit(as.character(abc), NULL)[[1]]]
Amino

#N.B: strsplit() splits the sequences in a set of sequences according to a
#pattern.

DNAString(myDNA)
translate(DNAString(myDNA))
AMINO_ACID_CODE[strsplit(as.character(translate(DNAString(myDNA))), NULL)[[1]]]
# YOU'RE KILLING IT MAN ! KEEP GRINDING ! IT'LL BE WORTH IT !

# Frequency of each amino acid:
table(AMINO_ACID_CODE[strsplit(as.character(translate(DNAString(myDNA))), NULL)[[1]]] ) 
# Better visualization of that :
plot( table(AMINO_ACID_CODE[strsplit(as.character(translate(DNAString(myDNA))), NULL)[[1]]] ) )


# Global Alignment technique:
d1 <- DNAString("ATCGGAAGTT")
d2 <- DNAString("AGCTGAACCG")
PairwiseAlignmentsSingleSubject(d1, d2)

# Global Alignment technique after specifying your own score for match and mismatch among the sequence:
d1 <- DNAString("ATCGGAAGTT")
d2 <- DNAString("AGCTGAACCG")
matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
pairwiseAlignment(d1, d2,substitutionMatrix = matrix)

#  Local Alignment technique:
d1 <- DNAString("ATCGGAAGTT")
d2 <- DNAString("AGCTGAACCG")
matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
pairwiseAlignment(d1, d2,substitutionMatrix = matrix, type="local")


# Overlap Alignment technique
d1 <- DNAString("ATCGGAAGTT")
d2 <- DNAString("AGCTGAACCG")
matrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
pairwiseAlignment(d1, d2,substitutionMatrix = matrix,type="overlap")


#-------------------------------------------------------------------------------------

# Exercise from --> https://compgenomr.github.io/book/exercises-4.html

library(rtracklayer)
library(AnnotationHub)
ahub <- AnnotationHub()
ahub

# Get hold of CpG islands from chr12 mouse mm9

fil <- subset(ahub, species == "Mus musculus")
fil
fjl <- query(fil, c("mm9", "UCSC", "CpG islands")) 
fjl
CpG <- fjl[["AH6117"]] #Granges 
CpG

seqlevels(CpG)

CpG_chr12 <- keepSeqlevels(CpG, "chr12",pruning.mode = "coarse")
CpG_chr12 

mean(width(CpG_chr12))
range(width(CpG_chr12))

# Get hold of RefSeq genes from chr12 mouse mm9

fkl <- query(fil, c("mm9", "UCSC", "refGene")) 
fkl
Refseq <- fkl[["AH6090"]]
Refseq 

seqlevels(Refseq)

Rfsq_chr12 <- keepSeqlevels(Refseq, "chr12",pruning.mode = "coarse")
Rfsq_chr12 

args(promoters)
 
prrfsq <- promoters(Rfsq_chr12, 1000, 1000 )
prrfsq

subsetByOverlaps(prrfsq, CpG_chr12)
# 2271 promoters that overlap CpG

length((subsetByOverlaps(prrfsq, CpG_chr12))$name) 

length (prrfsq) # 5198 promoters
length (CpG_chr12) # 627 CpG 

# Percentage of promoters which overlap with CpG islands
2271/5198 #= 0.44 = 44 % 

subsetByOverlaps(CpG_chr12, prrfsq)
# Percentage of CpG islands overlap with promoters
a <- 425/627
b<- round(a, 2)
c <- b*100
c # 68 % 

CpGoverlap <- subsetByOverlaps(CpG_chr12, prrfsq)
CpGoverlap

width(CpGoverlap)
 
plot( width(CpGoverlap)) 


#---- 02 SHOWING UP : Enrichment CpG, Promotors, TSS...---------------------------------

BiocManager::install("genomation")
library(genomation)


library("AnnotationHub")
ahub = AnnotationHub()

#______How transcription factor binding sites can be annotated using CpG islands

qhub = query(ahub, c("wgEncodeHaibTfbsGm12878Sp1Pcr1xPkRep1.broadPeak.gz")) 
qhub

peaks <- qhub[["AH24939"]]

peaks  # =  TF binding sites /  = transcription start sites (TSS)

# How far are promoters from TSS: Promoters are upstream of TSS.

# The majority of mammalian gene promoters are encompassed within regions
#of the genome called CpG islands that have an elevated level of
#non-methylated CpG dinucleotides.

# In vertebrates, DNA is pervasively methylated at CG dinucleotides, a modification 
#that is repressive to transcription. However, approximately 70% of vertebrate gene 
#promoters are associated with CpG islands (CGIs) that are refractory 
#to DNA methylation.


peaks21 <- keepSeqlevels(peaks, "chr21",pruning.mode = "coarse")
peaks21
#  620 peaks (TSS)


qhub2 = query(ahub, c("cpgi", "hg19"))  # Get CpG
qhub2

cpg <- qhub2[["AH5086"]]  
cpg

cpg21 <- keepSeqlevels(cpg, "chr21",pruning.mode = "coarse")
cpg21
#   365 CpG

subsetByOverlaps(peaks21, cpg21)

subsetByOverlaps(cpg21, peaks21) 

# 55 :
#       15% of CpG 
#                    ~9 % of peaks

length(subsetByOverlaps(cpg21, peaks21, ignore.strand = TRUE))/ length(cpg21)
# 15% 

length(subsetByOverlaps(peaks21, cpg21, ignore.strand = TRUE))/ length(peaks21)
# ~ 9%


# Check for each peak if it overlaps with a CpG:
countOverlaps(peaks21, cpg21) 
#N.B. tabulates the number of overlaps for each element in the query.

# Doing it the other way around
countOverlaps(cpg21, peaks21) 


# see one-to-one overlaps between peaks and CpG islands.
#It returns a matrix showing which peak overlaps which CpG island.
Overl <- findOverlaps(peaks21, cpg21)
Overl 

mean(width(peaks21)) # 685 
range(width(peaks21)) # 216 4319

mean(width(cpg21)) # 734 
range(width(cpg21)) # 202 5698

which( (countOverlaps(peaks21, cpg21)) >= 1) # queryHits in overl

which( (countOverlaps(cpg21, peaks21)) >= 1) # subjectHits in overl

sum(countOverlaps(cpg21, peaks21)) #56 overlaps

sum(countOverlaps(peaks21, cpg21)) #56


which( (countOverlaps(peaks21,cpg21)) == 2) #462th position

# look at [9] & [10] : the same peak (n°462) overlaps with 2 different CpG.
Overl[15:25] 


peaks21[462]
width(peaks21[462]) # 1062 bp

cpg21[93:94]
width(cpg21[93:94]) # 381 bp ; 969 bp

# " peak14901 " overlap with " CpG:_33 " and  " CpG:_105 ".


# There's 56 overlaps but why subsetByOverlaps() shows only 55 ?

 # Is because :

#subsetByOverlaps() returns a subsetted query containing only those 
#elements overlapped by at least one element in subject 

cpg21[12]

which(cpg21$name== "CpG:_33")  #58  (93) 320
which(cpg21$name== "CpG:_105")  # 94

length(unique(cpg21$name)) # Only 141 CpG  !!!!  /365
cpg21[c(58, 93, 320)]   

length(unique(peaks21$name))  # 620 /620 

cpg21[12]  # why do we have 2 overlaps here??

Overl[4:20] # CpG in the 12th position overlaps with 2 different peaks

which(countOverlaps(cpg21, peaks21) == 2 ) # 12

width(cpg21[12]) # 487 bp

width(peaks[250:251]) # 397 bp ; 707 bp



#*** find nearest CpG to each peak (TSS):

nearest(peaks21, cpg21)

distanceToNearest(peaks21, cpg21)

distanceToNearest(peaks21, cpg21)[462] # distance O !




###__________________ TSS and CpG for autosomes :

TSS <- peaks

CpG <- cpg

granges(TSS)$chr

TSS_aut <- keepSeqlevels( TSS, (paste0("chr", 1:22)), pruning.mode = "coarse"   )
TSS_aut

table(seqnames(TSS_aut))


CpG_aut <- keepSeqlevels( CpG, (paste0("chr", 1:22)), pruning.mode = "coarse"   )
CpG_aut

fovr <- findOverlaps(TSS_aut, CpG_aut)
fovr  # 5207 overlaps in total


sovr <- subsetByOverlaps(TSS_aut,CpG_aut , ignore.strand = TRUE)
sovr # 5156 of TSS &  4938 of CpG

which (countOverlaps(TSS_aut, CpG_aut) == 2) # 51 ( max overlaps = 2 )
# 5156 + 51 = 5207 

length(subsetByOverlaps(TSS_aut, CpG_aut, ignore.strand = TRUE))/ length(TSS_aut)
# 20 % of TSS have a CpG

length(subsetByOverlaps(CpG_aut, TSS_aut,  ignore.strand = TRUE))/ length(CpG_aut)
# 18 % of CpG are located inside TSS


# Is there any significant relationship between the two :

inOut <- matrix(0, ncol= 2, nrow = 2)

colnames(inOut) <- c("in", "out")
rownames(inOut) <- c("in", "out")

inOut

inOut[1,1] <- sum(width(intersect(TSS_aut, CpG_aut, ignore.strand = TRUE)))
inOut[1,2] <- sum(width(setdiff(TSS_aut,CpG_aut , ignore.strand = TRUE)))
inOut[2,1] <- sum(width(setdiff(CpG_aut, TSS_aut , ignore.strand = TRUE)))
inOut

inOut[2,2] = 3*10^9 - sum(inOut)
inOut

fisher.test(inOut)$statistic #ERROR

oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2] )
oddsRatio # = 26    : (>1) there's a strong relationship = significant enrichment.

# FROM LITERATURE:
# CpG islands typically occur at or near the transcription start site of genes.
# CpG islands often include transcription initiation sites and display 'active' histone marks,
#notably histone H3 lysine 4 methylation ( H3K4me3 ).


###_____________ How promoters of refseq genes are enriched in these CpG islands ?

#1)

cpg

qhub3 = query(ahub, "RefSeq")
qhub3

genes= qhub3[[1]]  
genes                 # SRY gene:  NM_003140 (human) & NM_011564 (mouse) from UCSC
 
which(genes$name==  "NM_003140" ) # WOOW ---> 24792

genes[24792 ]

width(genes[24792 ]) # 887 bp

seqlevels(cpg)

seqlevels(genes)

seqlevels(cpg) == seqlevels(genes) # We have the same seqlevels.

cpgstX <- keepStandardChromosomes(cpg, pruning.mode = "coarse")
cpgstX

genesstX <- keepStandardChromosomes(genes, pruning.mode = "coarse")
genesstX


seqlengths(cpgstX)/10^6

seqlengths(genesstX)/10^6

# Remove mitDNA
cpg_1Y <- dropSeqlevels(cpgstX, "chrM", pruning.mode = "coarse"  )
cpg_1Y 
seqlengths(cpg_1Y)

genes_1Y <- dropSeqlevels(genesstX, "chrM", pruning.mode = "coarse"  )
genes_1Y 
seqlengths(genes_1Y)


# N.B: 

# chromosome M (chrM) = Mitochondrial DNA = 16500 base pairs

# ChrUn contains clone contigs that cannot
#be confidently placed on a specific chromosome.


width(cpg_1Y)
range(width(cpg_1Y)) # 201 bp  ; 45712 bp
mean(width(cpg_1Y))  # 763 bp


width(genes_1Y)
range(width(genes_1Y)) # 20 bp ; 2320934 bp
mean(width(genes_1Y))  # 58336 bp

which(width(genes_1Y) == 20)
(genes_1Y)[c(17970, 18654, 18655, 29283, 32890)] #genes with only 20 bp


# How many transcripts we've got here ?

length(genes_1Y$name)  # 47460 transcripts
length(unique(genes_1Y$name))  # Only 46308 genes  ( - 1152 ) 

# Find out which genes have several transcripts:

# how many transcripts in all genes ? 
table(table(genes_1Y$name ))

# Here are the 500 genes who have 2 different transcripts: 
length(which(table(genes_1Y$name ) == 2 ) )

# e.g., this gene has 2 transcripts located in 13191  and 13199 
which((genes_1Y$name) == "NM_000344" )


# We're gonna ask if these Promoters are enriched in CpG

length( promoters(genes_1Y) ) # 47,460 total promoters
length( cpg_1Y )              # 27,718 total CpG 

prom_gen1Y <- promoters(genes_1Y)
prom_gen1Y

fov <- findOverlaps(prom_gen1Y, cpg_1Y)
fov  # 30,973 overlaps
 
cov <- countOverlaps(prom_gen1Y, cpg_1Y)
cov
sum(cov) # 30,973 overlaps

(queryHits(fov))[ 561] 

cov[561] 

prom_gen1Y[762]


sov <- subsetByOverlaps(prom_gen1Y, cpg_1Y )
sov #29534 


# Nbr of promoters that have a CpG in them :
length(unique(queryHits(fov))) # 29534 promoters
29534/ 47460*100                     
                   # 62 % of promoters have a CpG in them.

#** That makes sense 'cause:

        #In humans, about 70% of promoters located near the TSS
                      # of a gene (proximal promoters) contain a CpG island.

# number of CpG that have a promoter in them.
length(unique(subjectHits(fov))) # 14077 CpG
14077/ 27718*100 # 50 % of CpG have a promoter in them.

# How many Megabases do the promoters call ?
sum(width(reduce(prom_gen1Y, ignore.strand = TRUE))) 
# 59.44159 Mbp

# How many Megabases do the CpGs call ?
sum(width(reduce(cpg_1Y, ignore.strand = TRUE)))
# 21.15239 Mbp

# How big is the overlap ?
sum(width(intersect(prom_gen1Y, cpg_1Y, ignore.strand = TRUE)))/ 10^6
# 8.392124 Mbp

distanceToNearest(prom_gen1Y, cpg_1Y)


#explain that sh*t : 

sum(width(prom_gen1Y)) /10^6 # 104.412 Mbp !!!!

sum(width(reduce(prom_gen1Y)))/10^6 #  61.23916 Mbp !!!!
# is 'cause some genes have many transcripts.

###   THE RIGHT WAY TO DO IT :
sum(width(reduce(prom_gen1Y, ignore.strand = TRUE)))/10^6 #  59.44159 Mbp !!!!
#N.B.:
# ignore.strand = TRUE allows elements on the + strand 
#to overlap elements on the - strand.

#ignore.strand When set to TRUE, the strand information is ignored in the overlap calculations.

# For set operations: If set to TRUE, then the strand of x and y is set to "*" prior to any computation.

# For parallel set operations: If set to TRUE, the strand information is ignored in the computation and the result has the strand information of x. 


#Compute the Odd ratio: 
# Is a measure of association. It quantifies the relationship. 

inOut <- matrix(0, ncol= 2, nrow = 2)

colnames(inOut) <- c("in", "out")
rownames(inOut) <- c("in", "out")

inOut

inOut[1,1] <- sum(width(intersect(prom_gen1Y, cpg_1Y, ignore.strand = TRUE)))
inOut[1,2] <- sum(width(setdiff(prom_gen1Y,cpg_1Y , ignore.strand = TRUE)))
inOut[2,1] <- sum(width(setdiff(cpg_1Y, prom_gen1Y , ignore.strand = TRUE)))
inOut

inOut[2,2] = 3*10^9 - sum(inOut)
inOut

fisher.test(inOut)$statistic #ERROR

oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2] )
oddsRatio # = 37  ; (>1) there's a significant enrichment.



##_____________________  what about CpG islands & H3K27me3  ??? 


#
# H3K4me1 and H3K27ac are associated with enhancer regions;
# H3K4me3 is associated with promoter regions;
# H3K36me3 is associated with transcribed regions in gene bodies;
# H3K27me3 is associated with Polycomb repression;
# H3K9me3 is associated with heterochromatin.
#


cpg_1Y  

ahub = subset(ahub, species == "Homo sapiens")
qhub4 = query(ahub, c("H3K27me3"))  # Get H3K27me3
qhub4

H3K27me3 <- qhub4[[  "AH23223" ]]
H3K27me3

table(seqlevels(H3K27me3))

seqlevels(H3K27me3) 

seqlengths(H3K27me3)

hkm27st <- keepStandardChromosomes(H3K27me3, pruning.mode = "coarse"   )
hkm27st 

length(table(seqlevels(hkm27st )))

hkm27.1Y <-dropSeqlevels(hkm27st, "chrM", pruning.mode = "coarse" )
hkm27.1Y

width(hkm27.1Y)
range(width( hkm27.1Y ))  # min 100  max 4611210


length(table(seqlevels(hkm27.1Y)))

length( cpg_1Y )              # 27 718 total CpGs 
length( hkm27.1Y )            # 40 693 total histone marks

ovf<- findOverlaps(cpg_1Y, hkm27.1Y   )
ovf   # 7996 overlaps

ovc <- countOverlaps( cpg_1Y, hkm27.1Y    )
ovc

ovs <- subsetByOverlaps(cpg_1Y, hkm27.1Y , ignore.strand = TRUE )
ovs


length(unique(queryHits(ovf)))   # 7588  cpgs have that histone mark
length(unique(subjectHits(ovf))) # 3810  histone marks are within CpG

length(subsetByOverlaps(cpg_1Y, hkm27.1Y, ignore.strand = TRUE))/ length(cpg_1Y)
# 27 % of these CpG have a histone mark within them.

length(subsetByOverlaps(hkm27.1Y, cpg_1Y , ignore.strand = TRUE))/ length(hkm27.1Y)
# 9 % of these histone marks are located within CpG.

# How many Megabases do the histone marks call ?
sum(width(reduce(hkm27.1Y, ignore.strand = TRUE)))/ 10^6
# 861.3961 Mbp

# How many Megabases do the CpGs call ?
sum(width(reduce(cpg_1Y, ignore.strand = TRUE)))/10^6
# 21.15239 Mbp

# How big is the overlap ?
sum(width(intersect(hkm27.1Y, cpg_1Y, ignore.strand = TRUE)))/ 10^6
#  5.86855 Mbp


#Compute the Odd ratio: 
# Is a measure of association. It quantifies the relationship. 

inOut <- matrix(0, ncol= 2, nrow = 2)

colnames(inOut) <- c("in", "out")
rownames(inOut) <- c("in", "out")

inOut

inOut[1,1] <- sum(width(intersect(hkm27.1Y, cpg_1Y, ignore.strand = TRUE)))
inOut[1,2] <- sum(width(setdiff(hkm27.1Y,cpg_1Y , ignore.strand = TRUE)))
inOut[2,1] <- sum(width(setdiff(cpg_1Y , hkm27.1Y, ignore.strand = TRUE)))
inOut

inOut[2,2] = 3*10^9 - sum(inOut)
inOut

fisher.test(inOut)$statistic #ERROR

oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2] )
oddsRatio # = 0.95  ; (< 1) weak association = no significant enrichment.

# H3K27me3-only promoters are CpG poor (from literature). 


##_____________________  How about prom_gen1Y & hkm27.1Y ???

length(subsetByOverlaps(prom_gen1Y, hkm27.1Y, ignore.strand = TRUE))/ length(prom_gen1Y)
# 28 % 

length(subsetByOverlaps(hkm27.1Y, prom_gen1Y, ignore.strand = TRUE))/ length(hkm27.1Y)
# 12 % 

inOut <- matrix(0, ncol= 2, nrow = 2)

colnames(inOut) <- c("in", "out")
rownames(inOut) <- c("in", "out")

inOut

inOut[1,1] <- sum(width(intersect(hkm27.1Y, prom_gen1Y, ignore.strand = TRUE)))
inOut[1,2] <- sum(width(setdiff(hkm27.1Y,prom_gen1Y , ignore.strand = TRUE)))
inOut[2,1] <- sum(width(setdiff(prom_gen1Y , hkm27.1Y, ignore.strand = TRUE)))
inOut

inOut[2,2] = 3*10^9 - sum(inOut)
inOut

oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2] )
oddsRatio # = 0.93  ; (< 1) weak association = no significant enrichment.


##_____________________  What about prom_gen1Y & H3K4me3 ???

prom_gen1Y 

ahub = subset(ahub, species == "Homo sapiens")
qhub5 = query(ahub, c("H3K4me3", "Gm12878")) 
qhub5

H3K4me3 <- qhub5[[1]] 
H3K4me3

table(seqlevels(H3K4me3 ))

H3K4me3.aut <- keepSeqlevels(H3K4me3, paste0("chr", 1:22), pruning.mode = "coarse" )
H3K4me3.aut

length(subsetByOverlaps(prom_gen1Y, H3K4me3.aut, ignore.strand = TRUE))/ length(prom_gen1Y)
# 61 % 

length(subsetByOverlaps(H3K4me3.aut, prom_gen1Y, ignore.strand = TRUE))/ length(H3K4me3.aut)
# 26 % 

inOut <- matrix(0, ncol= 2, nrow = 2)

colnames(inOut) <- c("in", "out")
rownames(inOut) <- c("in", "out")

inOut

inOut[1,1] <- sum(width(intersect(H3K4me3.aut, prom_gen1Y, ignore.strand = TRUE)))
inOut[1,2] <- sum(width(setdiff(H3K4me3.aut,prom_gen1Y , ignore.strand = TRUE)))
inOut[2,1] <- sum(width(setdiff(prom_gen1Y , H3K4me3.aut, ignore.strand = TRUE)))
inOut

inOut[2,2] = 3*10^9 - sum(inOut)
inOut

oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2] )
oddsRatio # = 14.5  ; (> 1) Strong association = significant enrichment.



##_____________________  what about CpG islands  & H3K4me3 ???


length(subsetByOverlaps(cpg_1Y, H3K4me3.aut, ignore.strand = TRUE))/ length(cpg_1Y)
# 51 % 

length(subsetByOverlaps(H3K4me3.aut, cpg_1Y, ignore.strand = TRUE))/ length(H3K4me3.aut)
# 26 % 

inOut <- matrix(0, ncol= 2, nrow = 2)

colnames(inOut) <- c("in", "out")
rownames(inOut) <- c("in", "out")

inOut

inOut[1,1] <- sum(width(intersect(H3K4me3.aut, cpg_1Y, ignore.strand = TRUE)))
inOut[1,2] <- sum(width(setdiff(H3K4me3.aut,cpg_1Y , ignore.strand = TRUE)))
inOut[2,1] <- sum(width(setdiff(cpg_1Y , H3K4me3.aut, ignore.strand = TRUE)))
inOut

inOut[2,2] = 3*10^9 - sum(inOut)
inOut

oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2] )
oddsRatio # = 22  ; ( >1) Strong association = significant enrichment.



#----03 SHOWING UP : Sequence Matching and Aligning using Biostrings------------

# BioConductor contains a number of string matching/pairwise
#alignment tools in the Biostrings package that can be invaluable
#in answering complex scientific questions.


dnaseq <-   DNAString("ACGTACGT")

dnaseq2 <-  DNAString("AGGTATCGT")


pairwiseAlignment(dnaseq,   dnaseq2) # by default: type= "global"

pairwiseAlignment(dnaseq,   dnaseq2, type= "local")

pairwiseAlignment(dnaseq,   dnaseq2, scoreOnly= T ) #Get only the score.

pairwiseAlignment(dnaseq,Scerevisiae$chrV, type= "local")

pairwiseAlignment(dnaseq,Scerevisiae$chrI)


# generate a set of 16 sequences in order to map them to a sequence in Scerevesiae genome.

which(LETTERS== "A") #1

which(LETTERS== "C") # 3

which(LETTERS== "G") #7

which(LETTERS== "T") # 20

nuc <- LETTERS[c(1,3,7,20 )]   

sample(nuc, 8, replace= T  )

paste( sample(nuc, 8, replace= T  ) , collapse="") 


    # A much quicker way to do it:

sample (DNA_BASES, 8,  replace= T ) # create a generator of sequences

paste( sample (DNA_BASES, 8,  replace= T ) , collapse="") 


# BAM ! Here are my sequences :

setseq <- DNAStringSet(  )

for (i in 1:16) {
        
        setseq[i] <- paste( sample (DNA_BASES, 8,  replace= T ) , collapse="") 
}

setseq

# Alignment of these sequences to ScScerevisiae genome:

pairwiseAlignment( c( "ACAB","ABAC" ), "ACBA ")    # If u wanna play with it !

pairwiseAlignment( setseq, Scerevisiae$chrI ) 

   #About pairwiseAlignment():

# It aligns one or more strings specified in the pattern argument with a single 
#  string specified in the subject argument.
# It returns only the first occurrence of the best scoring alignment.
# If more than one pairwise alignment has the maximum alignment score exists,
#  the first alignment along the subject is returned.

pairwiseAlignment( setseq, Scerevisiae$chrM, type= "local" )

vmatchPattern(setseq, Scerevisiae$chrI)


           find.package("GenomicAlignments")
         library("GenomicAlignments") # Rank 26 / 2140 

           BiocManager::install("msa") # Multiple Sequence Alignment
           library("msa") #  Rank 161 / 2140
     
    
# Load amino acid sequences from one of the example files that are supplied 
# with the msa package: 
           
mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
mySequences <- readAAStringSet(mySequenceFile) 
mySequences           
           
           ?readAAStringSet() # it contains a filepath

myFirstAlignment <- msa(mySequences) # runs ClustalW with default parameters   

print(myFirstAlignment, show="complete")

args(msa)

# Alignment of our 16 seq. :

msa(setseq) # (ClustalW) # woow, that's awesome !

msa(setseq, "ClustalOmega") # a newer version of CLUSTALW

msa(setseq, "Muscle") # MUltiple Sequence Comparison by Log- Expectation

# we've got 3 different consensus sequence.

# Consensus Sequence: In a Nucleotide or an amino acid sequence,
#each base pair (an amino acid or a nucleotide) may occur more 
#frequently at a particular region in different sequences of nature.      
     
seqcer <- Scerevisiae$chrI[1:23 ]     

# hint:      
DNAStringSet(seqcer)

# Finally I managed to align my sequences to cerevisiae genome:
msa(c(setseq, DNAStringSet(seqcer) ) ) 

msa(c(setseq, DNAStringSet(Scerevisiae$chrI) ) ) # It takes ages !

# Get back consensus sequences
msaConsensusSequence( msa(c(setseq, DNAStringSet(seqcer) ) ) ) 
consensusString( msa(c(setseq, DNAStringSet(seqcer) ) ) )

# The main purpose of consensus sequences is to get an impression of 
#conservation at individual positionsof a multiple alignment.
          
msaConsensusSequence(msa(c(setseq, DNAStringSet(Scerevisiae$chrI) ) ))
           

msa(setseq) 

msaConsensusSequence(msa(setseq) ) 

consensusString(msa(setseq) ) # returns some arcane letters. 


###--------04 SHOWING UP :   Biostrings + Alignment     ------------------------------


     # Saccharomyces cerevisiae

letterFrequency(Scerevisiae$chrI, "GC", as.prob= TRUE)*100

sum(seqlengths(Scerevisiae)) #  12.163 Mbp = genome size 

# Note: 
#The term "megabase" (Mb) is commonly used inter-changeably with megabase pair (Mbp), 
#although strictly this would refer to a single-stranded nucleic acid.

seqnames(Scerevisiae) # 16 chromosomes 

seqlengths( Scerevisiae ) /10^6

alphabetFrequency(Scerevisiae$chrI)/10^6

#Gettting the GC content across the entire genome :

param <- new("BSParams", X = Scerevisiae, FUN = letterFrequency) 

   bsapply(param, "GC")

unlist(bsapply(param, "GC")

sum(unlist(bsapply(param, "GC")) / sum(seqlengths(Scerevisiae)))
# ~ 0.38 = 38 %

#Checking the GC % for each single X :
unlist(bsapply(param, "GC"))

            available.genomes()


   #Chlamydomonas reinhardtii genome
BiocManager::install("BSgenome.Creinhardtii.JGI.v5.6")
library(BSgenome.Creinhardtii.JGI.v5.6 )

Creinhardtii

seqnames(Creinhardtii)  # 17 chromosomes 

seqlengths( Creinhardtii ) /10^6

sum(seqlengths( Creinhardtii )) /10^6  # 111.1007 Mbp = genome size

alphabetFrequency(Creinhardtii$chromosome_13)/10^6

letterFrequency(Creinhardtii$chromosome_13, "GC", as.prob= TRUE)*100
# 62 %

sort(dinucleotideFrequency(Creinhardtii$chromosome_13)/10^3)

# % of each dinucltd.
sort (dinucleotideFrequency(Creinhardtii$chromosome_13, as.prob= TRUE )*100)

# %A %T %C %G :
alphabetFrequency(Creinhardtii$chromosome_13, as.prob= TRUE)*100


# GC content across the entire genome
param <- new("BSParams", X = Creinhardtii, FUN = letterFrequency) 

bsapply(param, "GC")

unlist(bsapply(param, "GC")

# % GC across the whole genome:
(sum(unlist(bsapply(param, "GC")) / sum(seqlengths(Creinhardtii)))) *100
# 62 %

#Checking the GC % for each single X
unlist(bsapply(param, "GC", as.prob = TRUE))


# Let's play for a while  with these genomes :

seqalgae <-    Creinhardtii$chromosome_1[1:16]
seqyeast <- Scerevisiae$chrI[49:64]

   #Note : An optimal alignment is an alignment giving the highest score.

pairwiseAlignment(seqalgae,   seqyeast)

pairwiseAlignment(seqalgae,   seqyeast, type= "local") # That makes sense.



# Automatic generators of sequences from these 2 genomes:
    sample ( Creinhardtii$chromosome_1  , 8,  replace= T )
    sample ( Scerevisiae$chrI  , 8,  replace= F )
    
pairwiseAlignment(sample ( Scerevisiae$chrI  , 8,  replace= T ), Scerevisiae$chrI , type= "local")

   
    #1st generator
algset <- list( ) 
for (i in 1:6) {

        algset[i] <-  sample ( Creinhardtii$chromosome_1  , 8,  replace= T )
}
algset

DNAStringSet(algset) # I'm killing it man !!!

   #2nd generator
yeaset <- list( ) 

for (i in 1:6) {

        yeaset[i] <-   sample ( Scerevisiae$chrI  , 8,  replace= T )
}
yeaset

DNAStringSet(yeaset) # I f--king kill it man every single solitary day !!! 

#'cause u know what ? I'm a bad motherf--ker !
# and one way or another I'm gonna make it f--king happen sooner or later !!!

 algseq <-   DNAStringSet(algset)
 yeaseq <-   DNAStringSet(yeaset)
 
names(algseq) <- paste0("algae", 1:6)
names(yeaseq) <- paste0("yeast", 1:6)


# Woow it works wonders ! 

msa( c(algseq,yeaseq) )                     # TG?A????
 
msa( c(algseq,yeaseq), "ClustalOmega")      # TGCA??C

msa( c(algseq,yeaseq),  "Muscle")           # ?T?A?CA?
      

pairwiseAlignment(algseq,yeaseq) # returns the best score

pairwiseAlignment(algseq,yeaseq[3], type= "local")

msaConsensusSequence(  msa( c(algseq,yeaseq) )  )


library(AnnotationHub)
ahub <- AnnotationHub()
ahub


           # Compute the GC content of the
#promoters in the yeast genome:

qhub <- query(ahub,   "sacCer2"  )
qhub 

Scgenes <-  qhub[["AH7048"]]
Scgenes 

proscgen <- promoters( Scgenes )
proscgen                          #  Out-of-bound Ranges !!!!


proscgen <- trim(proscgen )
proscgen

width( proscgen  ) 

# Get back DNA sequence of each promoter: 
proView <- Views(Scerevisiae,  proscgen      )
proView 

proView       

# DNA seq. of promoters
DNAStringSet(proView)

# GC % 
mean ( letterFrequency (proView, "GC" , as.prob= T) *100 )
     # 38 % GC 


GCprom <- letterFrequency (proView, "GC" , as.prob= T) *100 

plot(density(GCprom)) # a density plot of the GC content in the promoters
abline(v = 38)
 


###------------------05 SHOWING UP:  TxDB / transcripts, exons, cds...   ----------------------------------------------------

BiocManager::install("GenomicFeatures")
library(GenomicFeatures)

BiocManager:: install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

tdb <- TxDb.Hsapiens.UCSC.hg19.knownGene # creating a txdb object
tdb

genes(tdb)

#______________A slight detour to look in on my lovely gene ________

       # Finding SRY gene:
which(genes(tdb)$gene_id== "6736")  # 6736 =  gene id for SRY from NCBI 
 # GRanges object in position : 17998 

SRYtxdb <- genes(tdb)[17998] # Here's my lovely SRY gene !!

width(genes(tdb)[17998]) # 897 bp 

#__________________________Go back ! fool go back !_______________________


grSRY <- GRanges(seqnames = "chrY", strand = "-", 
              ranges = IRanges( start = 2654896, end = 2655792))
grSRY


# Which genes overlap our genomic region:
subsetByOverlaps( genes(tdb),  grSRY ) 

# 2 another genes on the reverse strand (+) overlap the same loci.  
subsetByOverlaps(genes(tdb), grSRY , ignore.strand = T)



 #_____If u wanna visualize a bit that sh*t___________

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

ir <- IRanges(start = c(2620337 , 2620337, 2654896), end = c( 2655675, 2734997, 2655792))

names(ir) <- c(" gene1 ", "gene2 ", "SRY" )

plotRanges(ir) #plot the IRange

#_________________________________________________


# Which transcripts overlap our genomic region:
subsetByOverlaps( transcripts(tdb),  grSRY ) # 1 transcript ( 1 exon)
     # That doesn't make much sense.


###_________ let us take another genomic region :

genr <- GRanges(seqnames = "chr1", strand = "+", 
                 ranges = IRanges( start =  11874, end = 794826 ))

width(genr)/10^3  # 782.953 Kbp


# we find out 3 genes in this genomic region: 
 subsetByOverlaps( genes(tdb),  genr ) 
 
# 19 transcripts WTF !!! :
subsetByOverlaps( transcripts(tdb), genr )

 # NOTE :
# More transcripts from a gene means higher expression 
#and fewer means lower expression.

# get a GRangesList of transcripts grouped by gene:
transcriptsBy(tdb)
transcriptsBy(tdb, by = "gene") # same.


# group exons by transcript.
exonsBy(tdb, by = "tx")

# elementNROWS : how long are each element of the list.
(elementNROWS( transcriptsBy(tdb)))[ "79501"] # 1
(elementNROWS( transcriptsBy(tdb)))[ "643837"] # 6
(elementNROWS( transcriptsBy(tdb)))[ "100287102"] # 3         

length(elementNROWS( transcriptsBy(tdb))) # 23459

transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))

ttt <- transcripts(genes(tdb), columns=c("tx_name", "gene_id"))
   
which(ttt$tx_name == c("uc001aaa.3", "uc010nxq.1", "uc010nxr.1", "uc001aal.1", "uc001aaq.2",
              "uc001aar.2", "uc009vjk.2", "uc001aau.3", "uc021oeh.1", "uc021oei.1",
              "uc010nxu.2", "uc001aax.1", "uc001abb.3", "uc031pjj.1", "uc001abp.2",
              "uc009vjn.2", "uc009vjo.2", "uc021oem.2", "uc031pjk.1") 
                                                                                       
                                                                                                                                                                          
      ttt$tx_name[1:19]
      
     shit <-  ttt[1:19]
     
(shit$gene_id)[10:19]

transcripts(genr, columns=c("tx_name", "gene_id"))

dft <- as.data.frame(ttt)

dim(dft)

which( dft[ ,8]== "79501") # 1
which( dft[ ,8]== "643837") # 6
which( dft[ ,8]== "100287102") # 3


# look at exons: My gosh 31 exons !!  
 subsetByOverlaps( exons(tdb),  genr) 

# figure out how these 31 exons are combined together to 
#form 19 different transcripts of 3 genes :
subsetByOverlaps(exonsBy(tdb, by = "tx"), genr) 

# coding sequence start
subsetByOverlaps(cds(tdb), genr) # 11 cds

# actual cds
subsetByOverlaps(cdsBy(tdb, by = "tx"), genr) 

#To see which transcript contains a real CDS:
subset(transcriptLengths(tdb, with.cds_len = T), gene_id == "100287102")# 3transcripts
subset(transcriptLengths(tdb, with.cds_len = T), gene_id == "643837")# 6transcripts
subset(transcriptLengths(tdb, with.cds_len = T), gene_id == "79501")# 1transcript

# Only gene_id ="100287102" $  gene_id = "79501" have cds.

  ### the f--king question that jacks my f-king mind up:
      # WHERE ARE THE 9 f--king other TRANSCRIPTS ???
          # SOME GENES DON'T HAVE an "id".
   
genor <- GRanges(seqnames = "chr22", strand = "+", 
                 ranges = IRanges( start =  19023795, end = 19543795 ))


genor

width(genor)/10^6 # =  0.5 Mbp

subsetByOverlaps( genes(tdb),  genor )  # 3 genes

subsetByOverlaps( transcripts(tdb), genor ) # 8 transcripts

subsetByOverlaps( exons(tdb),  genor) # WTF ! 29 exons !

subset(transcriptLengths(tdb, with.cds_len = T), gene_id == "23617")
# 1 transcript ; 1 cds = 1077

subset(transcriptLengths(tdb, with.cds_len = T), gene_id == "64976")
# 1 transcript ; 1 cds = 621

subset(transcriptLengths(tdb, with.cds_len = T), gene_id == "8318")
# 5 transcripts ; 4 cds 


# Translate the transcripts who have one or more cds :

cds_seqs <- extractTranscriptSeqs(Hsapiens,
                                  cdsBy(tdb, by="tx", use.names=TRUE))
cds_seqs

myaaseq <- translate(cds_seqs)
myaaseq

# Get amino acid seq of the 1st protein :
seq.aa1 <- AMINO_ACID_CODE[strsplit( as.character(myaaseq), NULL)[[1]]]

as.vector(seq.aa1) # Better visualization  ( 133 aa )

table(as.vector(seq.aa1)) # Get the frequency of each aa

inplot <- table(as.vector(seq.aa1))

prop.table(inplot)*100  # % of each aa 

# useful plot to notice the frequency of each aa:
plot(inplot, ylab="Frequency", main="Amino acids Frequency" ) 


# Get amino acid seq of the 3rd protein :

seq.aa3 <- AMINO_ACID_CODE[strsplit( as.character( translate (cds_seqs )), NULL)[[3]]]

as.vector(seq.aa3) 

toplot <- table(as.vector(seq.aa3))

prop.table(toplot)*100

plot(toplot, ylab="Frequency", main="Amino acids Frequency" ) 

# Get amino acid seq of the last protein :
seq.aa63691 <- AMINO_ACID_CODE[strsplit( as.character(myaaseq), NULL)[[63691]]]
seq.aa63691 
 
as.vector(seq.aa63691) 

plot3 <- table(as.vector(seq.aa63691))
 
plot(plot3, ylab="Frequency", main="Amino acids Frequency" ) 

prop.table(plot3 )*100

boxplot(plot3 , target="core")


#-----------06 SHOWING UP: Differential Gene Expression Analysis with limma -------------------------------

# Gene Set Enrichment Analysis (GSEA):
# is a computational method that determines whether an a priori defined set
# of genes shows statistically significant,concordant differences between two biological states (e.g., phenotypes).

# This is a powerful tool to associate a disease phenotype to a group of genes/proteins


    # Let us perform "Differential Expression Analysis" (DEA) 

# = The measurement of differences in genes expression between experimental groups.

#  The objective of an DEA is frequently to discover which genes
# are differentially expressed between 2 groups (e.g., cases & controls).

# Knowledge of DE genes facilitates the discovery of causative genes
#and gene pathways for a disease of interest.


library("GEOquery")

eList <- getGEO("GSE11675")

# Chronic myelogenous leukemia hematopoietic stem cells
# Expression profiling by array

eList

eList[[1]] #  access the data at the same level as "ALL".

exprs(eList[[1]]) # access the expression data 

eList[[1]]

dim(eList[[1]])

class(eList[[1]]) # ExpressionSet like "ALL".  

sampleNames(eList[[1]])[1:3]  
featureNames(eList[[1]])[1:3]

likeALL <- eList[[1]]

names(pData(likeALL ))
colnames(pData(likeALL )) # same

summary(pData(likeALL ))

colnames(pData(likeALL ))

featureData( likeALL )

pData(likeALL )

pData(likeALL )$title  # CML : chronic myelogenous leukemia
                       # = is an uncommon type of cancer of the bone marrow   


#  Differential Expression Analysis between CML $ no CML (Normal) 
 
 mydata <-  pData(likeALL)$title  
  

 Mydata <- ifelse(grepl("^CML", pData(likeALL)$title),
                  "CML", "Normal")
 
#   GSM296630          Lin-CD34-
#   GSM296635      CML Lin-CD34- 
#   GSM296636      CML Lin-CD34+ 
#   GSM296637          Lin-CD34+
#   GSM296638          Lin+CD34+
#   GSM296639      CML Lin+CD34+ 
 
 Mydat <-  factor(Mydata)
 
 
 # design matrix
 
 design <- model.matrix(~ Mydat )
 design
 
 likeALL

 
 # Fitting a linear model
 
 # 1) Estimate the fold changes and standard errors 
 #by fitting a linear model for each gene:
 
 fit <- lmFit(likeALL, design)
 
 # 2) Apply empirical Bayes smoothing to the standard errors:
 fit <- eBayes(fit)
 
 # 3) Show statistics for the top differentially expressed genes.
 
 topTable(fit) 
 
 topTable(fit, n=1) # Show statistics only for the 1st gene.
 
 
 plotMDS( likeALL)
 plotMDS( likeALL, dim.plot =c(3,4) , col= c ( rep("blue",1), rep("red",2 ), rep("blue",2) ,rep ("red",1 )),
          labels= c(rep("Normal",1), rep("CML" ,2), rep("Normal",2) , rep("CML",1)))
 ?plotMDS
              #INTERPRETATION :
 
  #* We have here a log2 fold-change (logFC) from "Normal" to "CML" (Normal - CML).
         # We've logFC of Normal/CML.
 
 # 32064_at (logFC is negative) has higher expression in Cancer (CML):
  # Up-regulation in cancer.
 
 # 464_s_at (logFC is positive) has lower expression in Cancer (CML).
  # Down-regulation in cancer.

 # 40095_at (logFC is negative) has higher expression in Cancer (CML):
 # Up-regulation in cancer.
 
 # How many genes are DE here ? 
 
 allg <- topTable(fit, number=12625 )
 
 length(which ( allg $adj.P.Val < 0.05 )) #  Only 3 ( 0.02 %) 
 
 
 
 fit$genes[11271,] 
 
 
# AveExpr: Average expression across all samples, in log2 CPM.
 
# t: the t-statistic used to assess differential expression.
   # = logFC divided by its standard error.
 
# P.Value: Raw p-value (based on t) 
 #the p-value for differential expression; this value is not adjusted for 
 # multiple testing
 
# adj.P.Val: the p-value adjusted for multiple testing. Different adjustment
 #methods are available, the default is Benjamini-Horchberg. 
 
# B: log-odds that gene is DE (arguably less useful than the other columns). 
 
 

  

  
 # Let us do it the other way around ( i.e.,  logFC of CML/Normal)
 
  Mydat2 <- relevel( Mydat, "Normal" )
  
  
  design <- model.matrix(~  Mydat2 )
  design
  
  
  fit <- lmFit( likeALL, design)
  
  fit <- eBayes(fit)          
  topTable(fit)      # Same intrepretation.
  
  
 

    
  #   playing a bit  with  leukemiasEset: 
  
  leukemiasEset$LeukemiaType
  
mytypes <- leukemiasEset[, leukemiasEset$LeukemiaType %in% c("AML", "CML")]  
mytypes   
  
mytypes$LeukemiaType

Mytyp <-  factor ( mytypes$LeukemiaType )
Mytyp

  design <- model.matrix(~   Mytyp )
  design
  
  fit <- lmFit( mytypes, design)
  
  fit <- eBayes(fit)          
  topTable(fit)      #  Up-reguation in CML and Down-regulation in AML.
  
     # 1e-6 = 0.000001
     # 1e3 =  1000
 
  
  plotMDS( mytypes, top= 3952,  dim.plot =c(3,4) , col= c ( rep("blue",12 ), rep ("red",12 ) ),
           labels= c(rep("AML" ,12), rep("CML",12) ) ) 
  
  
  # How many genes are DE here ? 
  
  allg <- topTable(fit, number=20172 )
  
  length(which ( allg $adj.P.Val < 0.05 )) #  Only 3952 ( 19 %) 
  

  
  
     # Let us play with another dataset :  
  
  BiocManager::install(c("breastCancerVDX"))
 #The vdx dataset from the breastCancerVDX package is an ExpressionSet of women
#diagnosed with breast cancer. It contains information about 22,000 genes and
#several variables including estrogen receptor (ER) status.  

  library(breastCancerVDX)  
  data(vdx)
  vdx

pData(vdx )
summary(  pData(vdx )) 

exprs(vdx)

# When breast cancer cells test positive for estrogen receptors, it's called estrogen
#receptor-positive (ER-positive) breast cancer. It means that estrogen is fueling
#the growth of the cancer.

erdata <- factor(pData(vdx)$er)

vdx

design <- model.matrix(~   erdata )
design

sort( erdata)

fit <- lmFit( vdx, design)

mean(fit$Amean)
range(fit$Amean)

hist(fit$Amean)
plotSA(fit)
mean(fit$Amean) #  = 7

keep <- (fit$Amean > 7)

fit2 <- eBayes(fit[keep,], trend=TRUE)

topTable(fit2, coef=2)

args(topTable)


fit <- eBayes(fit)          
topTable(fit)
    # we've a logFC of 1/0 (ER+ vs ER-) --> (1-0). 
  # These 10 genes are Up-regulated in breast cancer cells ER+.
args(topTable)

plotMDS( vdx, dim.plot =c(3,4) , col= c ( ifelse(grepl("^0", erdata), "red", "blue") ),
         labels= c( ifelse(grepl("^0", erdata), "NO", "YES")       ) )   
# That's the sweetest victory of ALL !!
# That builds confidence that's unbreakable man !
# All the efforts done today, even the failure felt today
# are all measures to create the moment later !!

plotMDS( vdx, dim.plot =c(3,4))

# Check if everything went right in the plot:
which( colnames(vdx)== "VDX_1661" ) # 326
erdata[326] # 0 ( NO / ER-)


plotMDS( vdx, dim.plot =c(3,4), top=11280  , col= c ( ifelse(grepl("^0", erdata), "red", "blue") ),
         labels= c( ifelse(grepl("^0", erdata), "NO", "YES")       ) )   


volcanoplot(fit,coef=2,  highlight= 50, names = fit$genes$NAME)
# towards the Right = the most upregulated genes 
# towards the  Left  = the most downregulated genes
# towards the Top   =  the most statistically significant genes


# 2 genes downregulated among the first 50 genes.
# 22270 
# 4278 

# Mapping probe IDs to gene symbols

fit$genes[ 4278, ] # "desmocollin 2"   
fit$genes[ 22270, ] #  No info available.

# OH MAN I'M KILLING IT !!

args(volcanoplot )
?volcanoplot


# Hint to compute how many genes are DE (DEGs) : 

allgen <- topTable(fit, number = 22283)

length(which ( allgen$adj.P.Val < 0.05 )) # 11280 genes (50% ) are DE.

# Notice there are genes which have a negative logFC:

topTable(fit, number = 50)$logFC 

  # e.g. gene "204751_x_at" in row 42 has a negative logFC.
   #  i.e. down-regulated in breast cancer cells ER+.


# If we're interested in a particular gene for example:

which( allgen==  "212195_at"  ) # 9
allgen[9,] 




#-----------07 SHOWING UP: GEOquery & GEO-NCBI -------------------------------

# Reading some data downloaded in my WD from GEO ncbi website : 

getwd()

setwd("C:/Users/OS/Desktop/RLang")

list.files()


# Reading in some data by getGEO()

fil <- getGEO("GSE205670")

fil[[1]] #  0 features  !!

# the submitters did not submit the data into GEO.  Instead, 
#the data were actually submitted as .txt files, one per sample.


# to deal with that, run :

supf <- getGEOSuppFiles( "GSE205670")

#This will download GSE205670_Raw.tar which you can untar, resulting
#in one file per sample.  You can then read these files into R one-at-a-time
#and attach the data to the existing ExpressionSet.


list.files("GSE205670") 

untar("GSE205670/GSE205670_RAW.tar", exdir = "GSE205670/CEL")
list.files("GSE205670/CEL")


celfils <- list.files("GSE205670/CEL", full = T)
celfils


library(oligo)

setwd("C:/Users/OS/Documents")
getwd()

file.exists("GSE68849/CEL/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz")
 



      


#_______
library(GEOquery)

#  Impact of "influenza A" on human plasmacytoid dendritic cells (pDC) gene expression.

 # Expression profiling by array ( DNA microarray analysis ).

 # Over 2600 genes were differentially expressed in pDCs exposed to influenza A
#compared with the controls (no virus). Multiple functional group clusters of genes 
#were impacted, including those involving antiviral responses and also metabolism, 
#such as the glycolysis, oxidative phosphorylation, glucose metabolic process, 
#and positive regulation of macromolecule biosynthesis clusters.

#******** Download all entities and....
library("GEOquery")


gds <- getGEO("GDS6063") 

gse <- getGEO("GSE68849") # it downloads also " gpl "
#*Processed data included within Sample table
#*Raw data are available on Series record



gse2 <- getGEO("GSE68849", GSEMatrix=F)

gpl <- getGEO("GPL10558") # USING LOCALLY CACHED VERSION...

gsm96 <- getGEO("GSM1684096")  

supp <- getGEOSuppFiles( "GSE68849") # 2 Supplementary files
     # Processed data included within Sample table
     # Raw data are available on Series record

soft <- getGEOfile("GSE68849") # USING LOCALLY CACHED VERSION
# a soft format file associated
#   with the GEO accession number.
# Simple Omnibus Format in Text (SOFT) is designed for batch 
#submission (and download) of data.


# The message "Using locally cached version..."
#means that getGEO() is not fetching the file from GEO, 
#but trying to use a local file on your machine.
#************


gds    # GEO DataSets   

    # 1) Curated sets of GEO Sample data.

  # 2) Information reflecting experimental
            #design is provided through GDS subsets.    

Columns(gds)

Columns(gds)[, 2 ]

Table(gds)

Meta(gds)


gpl   # Platform record 

 # 1) The list of elements on the array 
 #(e.g., cDNAs, oligonucleotide probesets, ORFs, antibodies).

 # 2) or The list of elements that may be detected and quantified in
# in that experiment (e.g., SAGE tags, peptides).

Columns ( gpl)
Meta( gpl)



gsm96  # Sample record

    # 1) The conditions under which an individual Sample was handled,
# the manipulations it underwent, and the abundance measurement of 
#each element derived from it.

Columns( gsm96)

Table(gsm96)

head(Table(gsm96))

dim( Table(gsm96))

Meta( gsm96)



gse    # Series record

  # 1) A set of related Samples considered to be part of a group.
  # 2) A focal point and description of the experiment as a whole.
  # 3) Tables describing extracted data.

esetgse <-  gse[[1]]

Meta(gse2)

gse2 # GSEMatrix=F : No expression data

GSMList(gse2)

names( GSMList(gse2)) # names of all the GSM objects contained in the GSE.

GSMList(gse2)[[1]]

head(exprs(esetgse) )

pData(esetgse) # phenotypic data

phenoData(esetgse)
# access the phenotypic data ( covariates )
# and meta-data ( descriptions of covariates ).
phenoData(esetgse)[[1]]
phenoData(esetgse)[[2]]
phenoData(esetgse)[[3]]
phenoData(esetgse)[[4]]
phenoData(esetgse)[[5]]
phenoData(esetgse)[[6]]
phenoData(esetgse)[[12]]

pData(phenoData(esetgse)) 

#  Converting GDS to an ExpressionSet: 
esetgds <- GDS2eSet( gds, do.log2=T )


           #** Let us perform DEA ;
          
          
library(limma)

pData(esetgse)$agent

cle <-   ifelse(grepl("^No", pData(esetgse)$agent ),"ctrl", "virus"  ) 
cle

# (!)
grep("^No", pData(esetgse)$agent ) 
grepl("^No", pData(esetgse)$agent )


groups <- factor(cle)
groups

design <- model.matrix(~   groups)
design

fit <- lmFit( esetgse , design)

fit <- eBayes(fit)          
topTable(fit)

# looking at the log2FC:
plotMDS( esetgse, dim.plot =c(3,4), top=50  , col= c ( ifelse(grepl("^ctrl", groups), "blue", "red") ),
         labels= c( ifelse(grepl("^ctrl", groups), "0", "V") ) )   


volcanoplot(fit, coef = 2,  highlight = 50 ,  names = fit$genes$NAME)
# We notice lots of genes upregulated in pDCs exposed to influenza A.
# coef: index indicating which coefficient of the linear model is to be plotted.

total <- topTable(fit, number= 47321 )

length(which ( total$adj.P.Val < 0.05 )) 
# 274 GDE     vs  2600 GDE brought up in the paper.
     #   what's wrong ???

length(which ( total$adj.P.Val < 0.01 )) 
# 63 GDE
     

### taking another shot:
covar <- pData(esetgse)$title

groups2 <-   factor(ifelse(grepl("No virus control", covar)," ctrl", "virus" ) )
groups2

design <- model.matrix(~   groups2)
design

fit <- lmFit( esetgse , design)

fit <- eBayes(fit)          
topTable(fit)

total <- topTable(fit, number= 47321 )

length(which ( total$adj.P.Val < 0.05 )) 
# Same results !!! 



###RAW DATA of: gse <- getGEO("GSE68849")

#  Impact of "influenza A" on human pDC gene expression.


library(GEOquery)
library(limma)

getGEOSuppFiles("GSE68849")

getwd()
setwd("C:/Users/OS/Documents/Datasets_GEO-ncbi")

#GSE68849 contains Illumina Beadchip arrays and not Affymetrix arrays
rawd <- read.ilmn("GSE68849_non-normalized.txt.gz")



### EList-class {limma}

#A list-based S4 classes for storing expression values (E-values), for example for a set of one-channel microarrays or a set of RNA-seq samples. EListRaw holds expression values on the raw scale. EList holds expression values on the log scale, usually after background correction and normalization. 
#it inherits directly from class list so any operation appropriate for lists will work on objects of this class. In addition, EList objects can be subsetted and combined. EList objects will return dimensions and hence functions such as dim, nrow and ncol are defined.
#an EList object can be created directly by new("EList",x), where x is a list. 

class(rawd  )

nrow(rawd )

summary(rawd)
dim(rawd)

summary(rawd$other)
str(rawd$other)

rawd$source

rawd$other

rawd$E
length(rawd$E)
dim( rawd$E )

summary(rawd$E )

head(rawd)

dim(exprs(esetgse)) # 47321 genes
length(rawd$E) # 473230 probes 'cause we use a probeset for every single gene.

summary(rawd$E) 
range(rawd$E) #            60.65833       30903.71000
range( exprs(esetgse)) #   67.16471      26775.34000

rawd[[2]][1,]
exprs(esetgse)[1,]

sampleNames(esetgse)

pData(esetgse)[1]

colnames(rawd)  
sampNames <-  colnames(rawd) 

boxplot(rawd$E, target="core")

# creating an ExpressionSet for raw data

 #NOTE(!!!!!)
# You should not be converting to an ExpressionSet,
# before the normalization is complete.

esetraw <- new("ExpressionSet",exprs=rawd$E)
class(esetraw)

exprs(esetraw )
dim (esetraw)

# Adding phenotypic data to esetraw:
samplz <- c("GSM1684101", "GSM1684102", "GSM1684103", "GSM1684104", "GSM1684095", "GSM1684096", "GSM1684097",
            "GSM1684098", "GSM1684099", "GSM1684100")
pData(esetraw) <- data.frame( sample= samplz, group= c( rep(c("ctrl","virus" ) )   )  ) 
pData(esetraw)
esetraw


pData(esetraw)$sample
pData(esetraw)$group

length(featureNames(esetraw)) # Only 47 323 from 473 230 what happened ?
    #                           only 10 % of what we had.
sampleNames(esetraw)

head(exprs(esetraw ))

# *** let's perform DEA on the raw data ***

# figuring out to which group belongs each sample:

rr<- read.delim("GSE68849_non-normalized.txt", header=F) 

class(rr) # dataframe

head(rr,10)

pData(esetgse)[1]

  #GSM1684095      5455178010_A    ctrl
  #GSM1684096      5455178010_B         virus
  #GSM1684097      5455178010_E    ctrl
  #GSM1684098      5455178010_F         virus
  #GSM1684099      5455178010_I    ctrl
  #GSM1684100      5455178010_J         virus
  #GSM1684101      5522887032_E    ctrl
  #GSM1684102      5522887032_F         virus
  #GSM1684103      5522887032_I    ctrl
  #GSM1684104      5522887032_J         virus


# In the raw data : colnames(rawd)

# "5522887032_E"    ctrl        GSM1684101 
# "5522887032_F"        virus   GSM1684102 
# "5522887032_I"    ctrl        GSM1684103
# "5522887032_J"        virus   GSM1684104
# "5455178010_A"    ctrl        GSM1684095  
# "5455178010_B"        virus   GSM1684096 
# "5455178010_E"    ctrl        GSM1684097
# "5455178010_F"        virus   GSM1684098 
# "5455178010_I"    ctrl        GSM1684099
# "5455178010_J"        virus   GSM1684100


pData(esetraw)$group

groop <- factor(pData(esetraw)$group)
groop 

design <- model.matrix(~   groop )
design

fit <- lmFit( esetraw , design)

fit <- eBayes(fit)          
topTable(fit)

total <- topTable(fit, number= 47323 )

length(which ( total$adj.P.Val < 0.05 )) 
# 113 GDE
# lots of genes are upregulated here.
volcanoplot(fit, coef = 2,  highlight = 200 ,  names = fit$genes$NAME)



# Performing RMA Normalization

#Preprocessing expression data
library(oligo) 
normD <- rma(rawd$E)
normD

args(rma)

BiocManager::available("affy")

class(esetraw)

help.search("rma") # find out which package contains this function()

GFset <- new("GeneFeatureSet",exprs=rawd$E)
GFset

args(new)

exprs(GFset )
dim(GFset)

# Another way of creating an ExpressionSet object.
 library(Biobase)
 datatable <- matrix(rnorm(473230 * 10), ncol = 10)
 esetdata <- ExpressionSet(datatable)
 dim(esetdata)

 EsetR <- ExpressionSet(rawd$E )
 class(EsetR) 
 
rma(esetraw)
rma(GFset)
 

 gse <- getGEO("GSE68849")
 gse

 esetgse <-  gse[[1]]
 esetgse 
 
pData(esetgse)[1]

length(featureNames(esetgse) )
length(featureNames(esetraw))
 
head(exprs(esetraw))        

# PROBLEM rma() DON'T READ NEITHER esetraw NOR GFset !!
 #NOTE:
 #(!!!) 
# rma() must not be used here 'cause these are not Affymetrix arrays.


#*___THE RIGHT WAY OF DEALING WITH THE RAW DATA of getGEO("GSE68849")

library(limma)
 # reading in the raw data
#Reading Illumina BeadChip Data

x <- read.ilmn("GSE68849_non-normalized.txt.gz", probeid="ID_REF")
x

getwd()
setwd("C:/Users/OS/Documents/Datasets_GEO-ncbi")

# the normalization's method for Illumina BeadChips arrays:
y <- neqc(x)
y             # it also automatically removes the control probes,
             #leaving only the regular probes in y.

dim(y) #   47323    10

dim(y$E) # 47323    10

length(y$E) # 473230

class(y$E) "matrix" "array"

head(y$E)
head(rawd$E)

#______________________________________________________________________________
# Negative controls are probes designed for DNA sequences that should not be present 
#in the sample.

# Positive controls are replicate probes for sequences that should occur;
#often these are in fact abundant.
#______________________________________________________________________________

#Now that the Normalization is done, let's proceed to a DEA:

#The ftp link allows getting the metadata for the file to know which 
#sample belongs to which conditions by clicking on 
# the "Series Matrix File(s)" link.

getwd()

setwd()

gse <- getGEO("GSE68849")

#-> read sample annotation from a GEO series matrix file into data.frames
targets  <- readSampleInfoFromGEO("C:/Users/OS/AppData/Local/Temp/RtmpAtlwL9/GSE68849_series_matrix.txt.gz")
targets 

# How I managed to do it????
# 1) ran getGEO("GSE68849")
# 2) copied the 1st path showed and supplied it to readSampleInfoFromGEO().


#If u throw sh*t against the wall eventually it sticks !
# Keep on keeping on, keep going motherf*cker...
#My glory doesn't happen in front of the crowd, it happens in
# that dark room during these dark moments, in solitude, alone,
# where I try and I try and I try again to be everything that possibly
# I can be, when life knocks me down, drops me to my knees, then I find
# the intestinal fortitude to pull myself together, pick myself up, and
# try that same endeavor and find myself victorious in that same endeavor
# that the greatest glory of all.
# One day's gonna be my day, and I'm gonna fly !  # (11:00 p.m ~ 6/20/2022)

summary(targets) 

class(targets$CharacteristicsCh1)


# Proceeding with DEA

# N.B. here we have Paired Samples.

# Paired samples occur whenever we compare two treatments and each independent 
#subject in the experiment receives both treatments.


Donor <- factor(targets$CharacteristicsCh1[,"donor"])
Agent <- factor(targets$CharacteristicsCh1[,"agent"])

data.frame(Donor,Agent)

design <- model.matrix(~Donor+Agent)

fit <- lmFit(y,design)

#__________________________________________________________________________________
# ***Removing Unexpressed genes using Amean :

# Amean is supposed to hold the mean log-intensity 
#for each probe, an indicator of the overall
# expression level of the corresponding transcript.

#However, if the data values input to lmFit() are log-ratios rather
#than log-intensities, then there is no way to compute the mean
#log-intensities. If the first argument to lmFit() is a matrix,then there is
#no way for the function to determine whether the input 
#values are log-ratios or log-intensities.
#If they are log-ratios, then calculating Amean would be meaningless.

#If the first argument to lmFit() is an MAList or exprSet or
#ExpressionSet then Amean is computed. It is actually not clear that
#this is right behaviour in the case of an ExpressionSet, but that's
#how it is for now.

#  The average log-expression values (AveExpr) of the probes is computed bu lmFit
# and stored as Amean.

# To identify a cutoff below which Amean values
#can be filtered use:
hist(fit$Amean)
plot(fit$sigma, fit$Amean)
plotSA(fit)

#________________________________________________________________________________

keep <- (fit$Amean > 5)

fit2 <- eBayes(fit[keep,], trend=TRUE)

topTable(fit2, coef=6)
topTable(fit2, coef="AgentNo virus control") # same

# coef: column number or column name specifying which coefficient or contrast of the
#linear model is of interest.

volcanoplot(fit2, coef = 6,  highlight = 20 ,  names = fit$genes$NAME)

class(fit)

total <- topTable(fit2, number= 47323 )

length(which ( total$adj.P.Val < 0.05 )) 
# 8105 ???

expressed <- rowSums(y$other$Detection < 0.05) >= 3

y <- y[expressed,]
 dim(y)


###\\\\_______________________________________

getwd()
setwd("C:/Users/OS/Documents/Datasets_GEO-ncbi")

getGEOSuppFiles("GSE206103")

list.files("GSE206103") 

untar("GSE206103/GSE206103_RAW.tar", exdir = "GSE206103/CEL")
list.files("GSE206103/CEL")

celfilz <- list.files("GSE206103/CEL", full = T)
celfilz

library(limma)
library(oligo)

rawDada <- read.celfiles(celfilz)
rawDada   # ExpressionFeatureSet
                    # 300304 features
pData(rawDada )

sampleNames(rawDada)
colnames(rawDada) # same

head(exprs(rawDada))

# Cleaning up the pdata
sampleNames(rawDada)
sampleNames(rawDada)  <- sub("_.*","", sampleNames(rawDada))
sampleNames(rawDada)

# Creating a group variable
pData(rawDada)$group <- c(sort(rep(c("CTL ","HOXC13-AS"),3)))
pData(rawDada)

dim(exprs(rawDada))
range(exprs(rawDada))

boxplot(rawDada, target="core")

normDada <- rma(rawDada)
normDada   # 27189 features

boxplot(normDada)

# DEA

pData(normDada)

gdea <- factor(pData(normDada)$group)
gdea

design <- model.matrix(~  gdea)
design

fit <- lmFit( normDada , design)

fit <- eBayes(fit)          
topTable(fit) 

total <- topTable(fit, number= 27189 )

length(which ( total$adj.P.Val < 0.05 ))
# 101 GDE
 
volcanoplot(fit, coef = 2,  highlight = 100 ,  names = fit$genes$NAME)
# --> many genes significantly up or down regulated.


#-----------08 SHOWING UP: BiomaRt & Biomart Project-------------------------------

#  BioMart est développé dans le but de permettre d'interroger, 
#avec un seul outil, un grand nombre de bases de données
#internationales.


library("biomaRt")

# In biomaRt lingo --> a Mart is a DataBase.

listMarts()

# Establish a connection to a database 
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart
      

# Get the datasets inside the selected database 
listDatasets(mart) 

ensembl<- useDataset("hsapiens_gene_ensembl", mart)
ensembl


#**** All that in a single command
#*select a both the database and dataset in one step:

ensembl2 <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl")
ensembl2
#****


listEnsembl() # All ENSEMBL's databases in Biomart project.

 # Looking for a dataset in a given database
searchDatasets(mart = mart , pattern = "hsapiens")
searchDatasets(mart = mart , pattern = "zebr")


#Using archived versions of Ensembl
listEnsemblArchives()

listEnsembl(version = 105)

ensembl105 <- useEnsembl(biomart = 'genes', 
                        dataset = 'hsapiens_gene_ensembl',
                        version = 105)

#Using Ensembl Genomes
listEnsemblGenomes()

ensembl_plants <- useEnsemblGenomes(biomart = "plants_mart")

searchDatasets(ensembl_plants, pattern = "Arabidopsis")


#*** Builiding a biomaRt query via getBM()

 #Filters and values are used to define restrictions on the query.

 head(listFilters(ensembl))
   # -> displays all available filters in the selected dataset.

 #Attributes define the data we are interested in retrieving
 
 head(listAttributes(ensembl))
   # -> displays all available attributes in the selected dataset.

 
 
 # Mapping probes id's to gene symbols
 
#   We have a list of Affymetrix identifiers from the u133plus2
# platform and we want to retrieve the corresponding EntrezGene
# identifiers using the Ensembl mappings:
 
 affyids <- c("202763_at","209310_s_at","207500_at") 
 
 
which(listAttributes(ensembl)== "affy_hg_u133_plus_2")
 listAttributes(ensembl)[107,]
 
which(listAttributes(ensembl)== "entrezgene_id")
 listAttributes(ensembl)[77,] 
 

 getBM(attributes = c('affy_hg_u133_plus_2', 'entrezgene_id'),
       filters = 'affy_hg_u133_plus_2',
       values = affyids, 
       mart = ensembl) 
 
# we can use
 searchAttributes() and searchFilters() 
       # the same way we used searchDatasets()

 
#****__ Using Biomart for non-Ensembl databases ______
 
  #Demonstrating the setting required to query Wormbase ParaSite and Phytozome
 
#1) Connecting to Wormbase BioMart
 
 listMarts(host = "parasite.wormbase.org")
 
 wormbase <- useMart(biomart = "parasite_mart", 
                     host = "https://parasite.wormbase.org", 
                     port = 443)
 wormbase
 
 listDatasets( wormbase) 
 
 ds_wormbase <- useDataset("wbps_gene", wormbase)
 
 head(listFilters(ds_wormbase))
 head(listAttributes(ds_wormbase))
 
 
 #2) Connecting to Phyotozome version 13 
 
 phytozome_v13 <- useMart(biomart = "phytozome_mart", 
                          dataset = "phytozome", 
                          host = "https://phytozome-next.jgi.doe.gov")
 phytozome_v13
 
 

 
 
 
 
 
#_________________

#_Working with the processed data of GEOSuppFiles("GSE38792") c.f.oligo section___

directa <-   getGEO("GSE38792")[[1]]
# microarray

pData(directa )

experimentData( directa)
# OSA is associated with alterations in visceral fat gene expression.

# !!! !!!!!!
# A total of 368 DEGs (176 upregulated and 192 downregulated)
# were identified in OSA samples.

# !!!!  The genes with p value < 0.01 were deemed to be
  #differentially expressed genes (DEGs).

# Differential expression analysis was conducted using limma package,

summary( pData(directa ))

length(exprs( directa))

# We're gonna perform DEA on this data.

pData(directa )$title

pData(directa )$characteristics_ch1

pData(directa )$subject # that seems easy to work with.

group <- factor(pData(directa )$subject)
group


design <- model.matrix(~   group)
design

fit <- lmFit( directa , design)

fit <- eBayes(fit)          
topTable(fit) 

allg <- topTable(fit, number= 33297 )

length(which ( allg $adj.P.Val < 0.05 )) # 0  It's strange here !!

# Let's do it with pData(directa )$title


grup <- ifelse(grepl("^Visceral fat_Control_", pData(directa)$title),
               "Control", "OSA")

grup

grup1 <- factor( grup)
grup1 


design <- model.matrix(~   grup1)
design

fit <- lmFit( directa , design)

fit <- eBayes(fit)          
topTable(fit) 

allg <- topTable(fit, number= 33297 )

length(which ( allg $adj.P.Val < 0.05 )) # 0 again what the f-ck is wrong !


# if u keep on throwing sh*t against the wall eventually it will stick!

# 3rd shot:


design <- model.matrix(~   grup1 )
design

fit <- lmFit( directa , design)
fit <- eBayes(fit)          
topTable(fit) 

design <- model.matrix(~   grup1 )
design


fit <- lmFit(directa, design)
fit <- eBayes(fit)  

plotMDS( directa, dim.plot =c(5,4))

plotMDS( directa , col= c ( rep("blue",8 ), rep ("red",10 ) ),
        labels= c(rep("c" ,8), rep("t",10) ) )   

# Generating a volcano plot to visualize differential expression
# Highlighting the top 2 genes
volcanoplot(fit, coef = 2,  highlight = 100 ,  names = fit$genes$NAME)
    # It took ages to highlight all genes : 33297 
   # with nonsensical plot.

  # Anyway that seems to make sense...
  # Oh man ! that seems to confirm what's brought up in the paper.

# We notice much more genes located on the top of the left side.

  # REMEMBER :
# !!! !!!!!!
# A total of 368 DEGs (176 upregulated and 192 downregulated)
# were identified in OSA samples.

# !!!!  The genes with p value < 0.01 were deemed to be
#differentially expressed genes (DEGs).


#_________ 09 SHOWING UP:  LINEAR MODELS    ___________________________________________

# 1) Single-Channel Designs ( One-color microarray)


# Three RNA sources to be compared

#Suppose that the first three arrays
#are hybridized with RNA1, the next two with RNA2 and
#the next three with RNA3. 

#all pair-wise comparisons between the RNA sources are of interest. We assume that the data has
#been normalized and stored in an ExpressionSet object

# 3 (rna1)  + 2(rna2) + 3(rna3) = 8 arrays = 8 samples


esetdata <- ExpressionSet(datatable)


val <- matrix(runif(32323* 8, 1, 10),32323 )
val 

dim(val)

esetd <- ExpressionSet(val)

pData(esetd) <- data.frame( array= paste0("array",1:8), 
                            RNA_source=(rep( (paste0("RNA", 1:3)), c(3,2,3)) ))

pData(esetd)

head(exprs(esetd))

colnames(esetd) <- paste0("ar",1:8)
sampleNames(esetd)


# all pair-wise comparisons between the RNA sources are of interest

design <- model.matrix(~ 0+factor(c(1,1,1,2,2,3,3,3)))

colnames(design) <- c("group1", "group2", "group3")
design

#______________________________________________________
#The usual way of creating a design matrix :
pData(esetd)$RNA_source
rna <- factor(pData(esetd)$RNA_source)
design0 <- model.matrix(~ rna +0)
design0 
#___________________________________________________

fit <- lmFit(esetd, design)

# To make all pair-wise comparisons between the three groups
#the appropriate contrast matrix can be created by

contrast.matrix <- makeContrasts(group2-group1, 
                                 group3-group2, group3-group1, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

# A list of top genes differential expressed in 
# group2 vs group1 can be obtained from

topTable(fit2, coef=1, adjust="BH")
# N.B. : switch the coef (to 2, 3)  to view the other comparisons.

# The BH-adjusted p-values control the false discovery rate. 
#This ensures that the expected proportion of false positives
# in your set of significant DE genes is 
#below a certain threshold (usually 5%).


#The outcome of each hypothesis test can be assigned using

#Identify which genes are significantly DE for each contrast
results <- decideTests(fit2)

#A Venn diagram showing numbers of genes significant 
#in each comparison (contrast)
vennDiagram(results)


# Taking a concrete case

library(leukemiasEset)
data(leukemiasEset)  

leukemiasEset

head(exprs( leukemiasEset))

pData( leukemiasEset )

ctrst <- pData(leukemiasEset)$LeukemiaType

design <- model.matrix(~ ctrst +0)
design 

colnames(design) <- c("ALL", "AML", "CLL", "CML", "NoL")
design

fit <- lmFit(leukemiasEset, design)

contrast.matrix <- makeContrasts(NoL-ALL, 
                                 NoL-AML, NoL-CLL, NoL-CML ,levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

topTable(fit2, coef=1, adjust="BH", sort.by= "logFC")

results <- decideTests(fit2)

vennDiagram(results) 



#-----------09 SHOWING UP: minfi -----------------------------------------------

     # Data from DNA methylation microarrays

library(minfi)
library(GEOquery)

getGEOSuppFiles("GSE68777") 

getwd()

setwd("C:/Users/OS/Documents/Datasets_GEO-ncbi")

list.files("GSE68777") 

untar("GSE68777/GSE68777_RAW.tar", exdir ="GSE68777/idat")

BiocManager::install("crossmeta")
# This package automates common tasks such as downloading,
#normalizing, and annotating raw GEO data.
library(crossmeta)


#  >>>>>>>>>>>>>>>

 #Unzip a tar.gz file with untar:
filename <- "C:/Users/OS/Documents/Datasets_GEO-ncbi/GSE68777_RAW.tar"
untar(filename, list= T ) # works
untar(filename, exdir ="GSE68777/idat") # works


idatfiles <- list.files("GSE68777/idat", pattern = "idat.gz")

# We must gunzip each idatfile such that we can read the files into the R environment.
for(i in 1:length(idatfiles)){
        gunzip(filename = idatfiles[i], destname = gsub("[.]gz$", "", idatfiles[i]))
}


# We've removed ".gz" in the destination file : 
idatfiles[1]
gsub("[.]gz$", "", idatfiles[1])
                                                

#  Look at the current wd to find the path to the gunzipped IDAT files then feed it to
# the read.metharray.exp():

# READING THE IDAT FILES IN:
rgSet <- read.metharray.exp("C:/Users/OS/Documents/Datasets_GEO-ncbi/GSE68777/idat")

# Oh man ! Eventually it works. My perseverance's paid off.

rgSet

# Class '"RGChannelSet"'
#These classes represents raw (unprocessed) data from a two
#color micro array; 
#specifically an Illumina methylation array.

# CONSTRUCTOR :
#    RGChannelSet(Green = new("matrix"), Red = new("matrix"),annotation = "", ...)
     
sampleNames(rgSet)
featureNames(rgSet)

     
pData(rgSet) # Not available let's get it through GSE68777_series_matrix

geoMat <- getGEO("GSE68777")
geoMat

egse <- geoMat[[1]] 

head(exprs(egse))

pData(egse)


library(Biobase)
 
# Save the data for the next session: 

saveRDS(egse, 'MyExpressionSet.Rds')
saveRDS(rgSet, 'MyRGChannelSet.Rds')


# Let's go back, Jack !

getwd()
setwd("C:/Users/OS/Documents/Datasets_GEO-ncbi")

# Read back the processed data in
egse <- readRDS('MyExpressionSet.Rds')
egse

# Read back the raw data in
rgSet <- readRDS('MyRGChannelSet.Rds')
rgSet

#  We're only interested in 4 columns. 
pD <- pData(egse)[, c("title", "geo_accession", "characteristics_ch1.1",
                      "characteristics_ch1.2")]

colnames(pData(egse))
colnames(pD)
head(pD)

# Cleaning the "pD" up

names(pD)[c(3,4)] <- c("group", "sex")
pD$group <- sub("^diagnosis:","", pD$group)
pD$sex <- sub("^Sex:", "", pD$sex)
head(pD)

sampleNames(rgSet) <- sub(".*_5", "5", sampleNames(rgSet))

rownames(pD) <- pD$title

pD <- pD[sampleNames(rgSet),]

sampleNames(rgSet)

pData(rgSet) <- pD #ERROR
rgSet

library(methods)
pdata <- as(pD, "DataFrame")
pdata

pData(rgSet) <- pdata
pData(rgSet)   # Man, I've killed it !
# That's the sweetest victory of all !
# That's the moment where my glory happens.

pData(rgSet)$group

# Microarray raw data after Normalization 
# & Mapping to genome 
GRSet <- preprocessQuantile(rgSet) 
GRSet <- readRDS('GenomicRatioSet.Rds')
GRSet

sampleNames <- sampleNames(GRSet)
probeNames <- featureNames(GRSet)
pheno <- pData(GRSet)

# Returns the probe locations:
granges(GRSet)

# Type of each CpG
head(getIslandStatus(GRSet)) 


#USEFUL: 
# Methylation of cytosine bases in DNA CpG islands is an important
# epigenetic regulation mechanism in the organ development, 
#aging and different disease statuses.


### Measure the methylation level:


MethM <- getM(GRSet) # --> A matrix of normalized M values.

# + M value -->  more molecules are methylated than unmethylated,
# - M value -->  the opposite.


# Note :
#1) Both Beta-value and M-value statistics are used 
#as metrics to measure methylation levels.
#2) The Beta-value has a more intuitive biological interpretation,
#but the M-value is more statistically valid for the differential 
#analysis of methylation levels.


MethBeta <- getBeta(GRSet) # --> A matrix of normalized B values
MethBeta <- readRDS('BetaValueGRSet.Rds')
MethBeta

range(MethBeta)
#  0 -> CpG completely unmethylated 
#  1 -> CpG completely methylated



###*** Proceeding to a Differential Methylation Analysis :

#**1) Using the processed data 'egse' 
grp <-  factor(pData(rgSet)$group) 
grp  

design <- model.matrix(~  grp )
design
library(limma)
fit <- lmFit( egse, design)
fit <- eBayes(fit)          
topTable(fit, number= 485512) 
    volcanoplot(fit, highlight = 5,  names = fit$genes$ID)
 topp <- topTable(fit, number=485512 )
    length(which ( tot$adj.P.Val < 0.05 ))# 262 CpG sites Differentially Methylated
    
    
  #*Looking for CYP11A1 gene found
    #hypomethylated in Mania according to the paper :
   
   colnames(fit$genes) 
    
  fit$genes$UCSC_RefGene_Name
  
  which(fit$genes$UCSC_RefGene_Name== 'CYP11A1')
   
(fit$genes$UCSC_RefGene_Name)[c(117074, 324310, 328906, 365165, 429679, 435847)]
         #  We found 6 genes 'CYP11A1'
  dim(fit$genes)
  
  fit$genes[c(117074, 324310, 328906, 365165, 429679, 435847), ] 
  # UCSC_RefGene_Accession = NM_000781     #chr15
  
  #cg05939495   S_Shore  
  #cg17790333   S_Shore   
  #cg18068537   S_Shore 
  #cg20330630   S_Shore 
  #cg24482024   S_Shore 
  #cg24849517   S_Shelf 
  
  # looking at their logFC:
  topp[ c("cg05939495","cg17790333", "cg18068537", "cg20330630", 
               "cg24482024", "cg24849517"), ]
  
  tot2[ c("cg05939495","cg17790333", "cg18068537", "cg20330630", 
           "cg24482024", "cg24849517"), ]
  
  # no significant p-values
  
#*** 2) Using normalized Beta-Values from the raw data  
design2 <- model.matrix(~  grp )
design2
fit2 <- lmFit(MethBeta, design2)
fit2 <- eBayes(fit2)          
topM <- topTable(fit2, number= 50) 

range(topM$logFC)  
range(topM$adj.P.Val)

#From the paper:
#We identified a methylation locus in the CYP11A1 gene,
#which is regulated by corticotropin, that is hypo-methylated
#in individuals hospitalized for mania compared with unaffected controls.


write.table(topM, 'TopmethylatedCpGSites.txt', quote=FALSE, sep='\t', row.names=FALSE)
write.csv(topM, file = 'TopmCpGSites.csv')

tot2 <- topTable(fit2, number=485512 )
length(which ( tot2$adj.P.Val < 0.05 )) # 262  # 137 ???

volcanoplot(fit2, highlight = 5,  names = fit2$genes$ID)

# Checking for potential outliers by 
#analyzing a multidimensional scaling (MDS) plot. 

# Informative method to see if any sample groups
#cluster separately based on their DNA methylation.

plotMDS(MethBeta, sampGroups = pData(rgSet)$group) 
         
         
# Analysis of methylation data

# 'Scatter plots' are used to correlate the methylation data; 
# 'Bar plots' to visualize relative levels of methylation at each site tested;
# 'Heat maps' to cluster the data to compare the methylation profile at the sites tested.



  #** A simple analysis of Methylation level 


 #1) Based on Beta-values

summary(MethBeta)

methmania <- which( pData(rgSet)$group  == " Mania" )

MethBeta[,methmania]

ncol(MethBeta[,methmania]) # 20 colums Mania

metmania <- (MethBeta[,methmania])

range(metmania) # 0.004226364 -> 0.994707754
mean(metmania)    # 0.4879465
summary(metmania)

methctr <- which( pData(rgSet)$group  == " Ctr" )

MethBeta[,methctr]

metctr <- (MethBeta[, methctr])

range(metctr ) # 0.005213475  --> 0.993948831
mean(metctr)      # 0.4882371

heatmap(MethBeta) # Error: cannot allocate vector of size 878.1 Gb

barplot(MethBeta)
 
range(MethBeta)
range(exprs(egse))
range(GRSet)


 # 2) Based on M-values

range(MethM) # -7.8  -->  7.5

methmania <- which( pData(rgSet)$group  == " Mania" )

MethM[,methmania]

metmania <- (MethM[,methmania])

range(metmania) # -7.8 -->   7.5
mean(metmania)  # -0.3


methctr <- which( pData(rgSet)$group  == " Ctr" )

MethM[,methctr]

metctr <- (MethM[,methctr])

range(metctr) # -7.5 -->  7.3
mean(metctr)  # -0.3

range(metctr) # -7.8  -->   7.5
mean(metctr)  # -0.3


# For annotation
# Probe sequences for microarrays of type IlluminaHumanMethylation450
annotation(rgSet)
BiocManager::install("IlluminaHumanMethylation450kprobe")
library(IlluminaHumanMethylation450kprobe)
data(IlluminaHumanMethylation450kprobe)


# Note:
#A p-value of 5% means that 5% of all tests will result in false positives.
#A q-value of 5% means that 5% of significant results will result in false positives. 
# When you report your results, you should set your FDR to 0.1 or 0.05.


# Identifying differentially methylated regions (DMRs) and 
#differentially methylated positions (DMPs)

#1)____ Finding DMPs by dmpFinder() 

group <- pData(rgSet)$group 
dmp <- dmpFinder(MethBeta, pheno = group  , type = "categorical")
head(dmp)
# **dmpFinder() is a wrapper around lmFit from limma.
# When sample size is small (< 10) use the argument 'shrinkVar=TRUE'.

#2)_____ Finding DMRs by  bumphunter() 

# This f() in minfi is a version of the bump hunting algorithm
# adapted to the 450k array.

# Instead of looking for association between a single genomic location
#and a phenotype of interest, bumphunter looks for genomic regions that
#are differentially methylated between two conditions.


   #1 - Define your phenotype of interest
pheno <- pData(rgSet)$group 
designMatrix <- model.matrix(~ pheno)

   #2 - Run the algorithm with B=0 permutation on the Beta-values, 
#with a medium difference cutoff, say 0.2 (which corresponds to 20\%
#difference on the Beta-values):
dmrs <- bumphunter(MethBeta[c(1:10),], design = designMatrix, 
                   cutoff = 0.2, B=0, type="Beta")

   #3 - If the number of candidate bumps is large, say >30000, increase the cutoff.
   #4 - Run the algorithm with a large number of permutations, say B=1000:
dmrs <- bumphunter(MethBeta[c(1:10),], design = designMatrix, 
                   cutoff = 0.2, B=1000, type="Beta")

#Error in .local(object, ...) : 
       # If object is a matrix, pos must be specified


#Another Package for Differnetial Methylation Analysis:
# COHCAP (pronounced "co-cap") provides a pipeline to analyze single-nucleotide resolution methylation data (Illumina 450k/EPIC methylation array, targeted BS-Seq, etc.). It provides differential methylation for CpG Sites, differential methylation for CpG Islands, integration with gene expression data, with visualizaton options. 



#-----------10 SHOWING UP:  egdeR / DESeq2  -----------------------------------------------
     
            # Count-based RNA-seq analysis  

# Both packages implement a statistical methodology
#based on the negative binomial distributions.



  ## RNA-Sequencing based-experiment: GSE52778 = single-cell dataset
library(SummarizedExperiment)

library(airway)
data(airway) 
airway  # Single Cell RNA-Seq (human airway smooth muscle cell lines)

#Single-cell RNA sequencing (scRNA-seq) is a widely used technique for 
#profiling gene expression in individual cells. This allows molecular
#biology to be studied at a resolution that cannot be matched by bulk
#sequencing of cell populations. 



# class: RangedSummarizedExperiment 
# The "Ranged" part refers to the fact that
#the rows of the assay data (here, the counts)
#can be associated with genomic ranges (the exons of genes). 
#This association facilitates downstream exploration of results,

dim(airway)

head(assay(airway, "counts"))
head(assay(airway))

lsf.str("package:SummarizedExperiment") # list all the f()s of the package.

metadata(airway)

colData(airway) 
#This package provides a RangedSummarizedExperiment object of read counts in 
#genes for an RNA-Seq experiment on four human airway smooth muscle 
#cell lines treated with dexamethasone. 
#(that's why we've 'untrt' in "albut" colum).

summary(colSums(assay(airway)))

granges(airway)
rowRanges(airway) #same

as.vector(width(granges(airway)[1])) 

length(rowRanges(airway)) # 64102 genes

# Get genomeInfo, among others, nb of transcripts, exons, cds. 
metadata(rowRanges(airway))

# Retrieve element-level metadata 
unlGRl <- unlist(rowRanges(airway))
unlGRl
mcols(unlGRl)

assay(airway, "counts")# UNormalized data

kounts <- assay(airway, "counts")

# Here we've 2 treatment groups, each with 4 replicates:
airway$dex
colData(airway)$dex # same

airway$dex <- relevel(airway$dex, "untrt")
airway$dex

granges(airway)

# The counts represent the total number of reads aligning to each gene
#(or other genomic locus)

# A replicate is an exact copy of a sample that is being analyzed.


#----------------------* edgeR -----------------------------------------------------------------------------

library(edgeR)
# Differential analysis of sequence read count data.

#edgeR can be applied to differential expression at 
#the gene, exon, transcript or tag level.

edgeRUsersGuide() #open the User's Guide in a pdf viewer.


# There are two testing methods under the generalized linear models framework
# >>> Testing for Differentially Expressed Genes in more complex experiments 
#     (= with multiple factors ): 

# 1) Quasi-likelihood F-tests : glmQLFit() / glmQLFTest
#       --> Bulk RNA-seq

# 2) Likelihood ratio tests   : glmFit() / glmLRT() 
#       --> Single cell RNA-seq, datasets with no replicates.


#*** Generalized linear models are an extension of classical linear models
     #to nonnormally distributed response data.

#~~~~ Typical edgeR analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  x <- read.delim("TableOfCounts.txt",row.names="Symbol")
  group <- factor(c(1,1,2,2))
  y <- DGEList(counts=x,group=group)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~group)
  y <- estimateDisp(y,design)

  # to perform Quasi-Likelihood F-Tests:
   fit <- glmQLFit(y,design)
   qlf <- glmQLFTest(fit,coef=2)
   topTags(qlf)

  # to perform Likelihood Ratio Tests:
   fit <- glmFit(y,design)
   lrt <- glmLRT(fit,coef=2)
   topTags(lrt)
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   
# Filtering out genes with very low counts ***
# As a rule of thumb, a gene is required to have
#a count of 5-10 in a library to be considered 
#expressed in that library.

colSums(kounts)/10^6 # library sizes

colMins(kounts)
colMaxs(kounts)/10^3
colMeans(kounts) 

length(which(kounts <5  )) #  369525
sum(kounts <5 ) #same

length(kounts) # 512816 total number of counts
 
kounts[which(kounts <5  )] # 72 % of low counts ( < 5 reads )

airway # 64102 * 8 

which(kounts[, 1] < 5)

# We find 46194 genes lowly expressed in the 1st library
length(which(kounts[, 1] < 5))
sum(kounts[, 1] < 5) #same

# Here're the values
kounts[which(kounts[, 1] < 5)]


#__ _ __ _ _# Let us perform our analysis #__ _#


######### I/ Classic approach (Pairwise comparisons between the groups)


#1) make a DGEList object
grp <- airway$dex

dgel <- DGEList(counts=kounts, group=grp)

dgel$samples

#2) keep only rows (genes) that have worthwhile counts :
keep <- filterByExpr(dgel)
dgel <- dgel[keep, , keep.lib.sizes=F]


    #How many genes have been filtered out in each library:
    colSums(kounts) - dgel$samples$lib.size

    
    #Take phenotypic data into the DGEList:
    dgel$samples <- merge(dgel$samples, 
                       as.data.frame(colData(airway)),
                                                       by= 0) 
                
    #Take info about the genes into the DGEList: 
    dgel$genes <- data.frame(name = rownames(dgel$counts),
                            stringsAsFactors = F)   
    

#3) Normalization of the library sizes 
dgel <- calcNormFactors(dgel)                         
dgel  
 
#4) Estimating dispersions
#estimate common dispersion and tagwise dispersions in one run:
dgel <- estimateDisp(dgel)


#~~~~Note:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#For experiments with single factor
# (Pairwise comparisons between 2 ou + groups) we use :
#--> estimateDisp(y)
# * qCML method :  It fails to take into account the
#effects from multiple factors in a more complicated experiment.

#For experiments with multiple factors we use:
#--> estimateDisp(y, design)
# * CR method :  It takes care of multiple factors by
# fitting generalized linear models (GLM) with a design matrix
#________________________________________________________________________

 # Testing for DE genes using the exact test
# The exact test is only applicable to experiments with a single factor,
# it allows both common dispersion and tagwise dispersion approaches:
 et <- exactTest(dgel)
 DE1 <- topTags(et)

 write.table( DE1, 'edgeR_RNA-seq_DEGmethode1.txt', quote=F, sep='\t', row.names=T)
 getwd()
 setwd("C:/Users/OS/Desktop/RLang")
 rDE1 <- read.delim('edgeR_RNA-seq_DEGmethode1.txt',  header=T)

 
 
##########    II/ GLM approach:
 
 
# similar to the classic approach, but permits more general comparisons to be made.
                                
# It allows an infinite variety of contrasts to be tested between the groups. 
 
 
 # A design matrix is required to describe the treatment conditions :
 
 
design <- model.matrix(~grp)

colnames(design)[2] <- c("Dexamethasone")


#estimate common dispersion, trended dispersions and tagwise dispersions in one run:
dgel <- estimateDisp(dgel ,design)

#likelihood ratio tests (scRNA-seq)
fit <- glmFit(dgel ,design)
lrt <- glmLRT(fit,coef=2)
DE2 <- topTags(lrt)
                     # -->> the same result !
rownames(rDE1)  
rownames(DE2)

write.table( DE2, 'edgeR_RNA-seq_DEGmethode2.txt', quote=F, sep='\t', row.names=T)
rDE2 <- read.delim('edgeR_RNA-seq_DEGmethode2.txt',  header=T)

dim(lrt$table)

 tot <- topTags(lrt, n= 15926) 
 tot
 length(which(tot$table$PValue< 0.05)) # 3804 DEG 

 
# CPM is an appropriate unit for inter-sample comparison. 
# When comparing feature expression within samples, TPM should be used
#instead of RPKM/FPKM.  
 
 
# LogCPM can be understood as measuring expression level.
 
# CPM (counts per million) 
# RPKM (Reads Per Kilobase Million)  
# FPKM (Fragments Per Kilobase Million)
# TPM (Transcripts Per Kilobase Million)  
# --> 3 metrics attempt to normalize for sequencing depth and gene length.
# Those methods aim at adjusting for differences in overall read counts among 
# libraries. 
 
#(!)N.B.: Neither edgeR nor DESeq2 uses those methods.

# "Sequencing Depth" (also known as "read depth"): 
# is defined as the ratio of the total number of bases obtained by sequencing
# to the size of the genome or the average number of times each base is
# measured in the genome. 
 
 
 
# Practical statistical analysis of RNA-Seq data - edgeR 
# >>>>>> https://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html#differential-expression-analysis


 
 
#----------------------* DESeq2 -----------------------------------------------------------------------------
 
library(DESeq2)
# It provides methods to test for differential expression by use of
#negative binomial generalized linear models; the estimates of dispersion 
#and logarithmic fold changes.

# The DESeq2 model internally corrects for library size.
# A DESeqDataSet object should take only un-normalized counts as input.

# DESeqDataSet class extends the RangedSummarizedExperiment class
#of the SummarizedExperiment package .
 
 
 # Constructing a DESeqDataSet from a SummarizedExperiment object:
 
 ddsSE <- DESeqDataSet(airway, design = ~ cell + dex)
 ddsSE
 
 design(ddsSE)  
 
 # Here we use two factor variables.
 # Note: 
 #In order to benefit from the default settings of the package,
 #you should put the variable of interest at the end of the design formula 
 #and make sure the control level is the first level. 
 
 counts(ddsSE) # access the count data.
 assay(ddsSE)  # (equivalent command)
 
 # If we wish to pre-filter low count genes (but it's not necessary):
 keep <- rowSums(counts(ddsSE)) >= 10
 ddsSE <- ddsSE[keep,]
 
 
## --> Differential expression analysis ~~~~:
  
#The standard steps are wrapped into a single function: DESeq().
  ddsSE2 <- DESeq(ddsSE)
  
# results() generates results tables.
  
 #With no additional arguments to results(), the log2 fold change and Wald test 
 #p value will be for the last variable in the design formula, and
 #if this is a factor, the comparison will be the last level of this 
 #variable over the reference level.
  
 #the order of the variables of the design do not matter so long as we
 #specify the comparison to build a results table for, using the name or
 #contrast arguments of results.
  
  tab <-  results(ddsSE2)
  # we have a logFc from untrt to trt, which seems kinda unnatural !
  
  # let's switch it to trt - untrt :
  tab2 <- results(ddsSE2, contrast=c("dex","trt","untrt"))
  
  
  #* We could have defined 'untrt' as the reference level :
  airway$dex <- relevel(airway$dex, "untrt")
  ddsSE <- DESeqDataSet(airway, design = ~ cell + dex) 
  ddsSE$dex # verifying our setting
  ddsSE2 <- DESeq(ddsSE)
  tab3 <-  results(ddsSE2)  # we got the same result of tab2.
 
  
#  baseMean: mean of normalized counts for all samples.
#  log2FoldChange: log2 fold change.
#  lfcSE: standard error.
#  stat: Wald statistic.
#  pvalue: Wald test p-value.
#  padj: BH adjusted p-values.
 
  # Subsetting the results based on an adjusted p value threshold: 
  thresh <- subset(tab3, padj < 0.1)
  dim(thresh) # We got 4819 genes DE / 64102 genes.
  
  length(which(thresh$log2FoldChange > 0)) # 2604 Up
  length(which(thresh$log2FoldChange < 0)) # 2215 Down
  
  # All these infos and more can be retrieved through :
  summary(tab3)
 
  #(N.B):
  #If there are very many outliers (e.g. many hundreds or thousands)
  #reported by summary(res), one might consider further exploration 
  #to see if a single sample or a few samples should be removed due
  #to low quality (c.f. vignette).
  
  # Ordering our results table by the smallest p value:
  permut <- order(tab3$pvalue) 
  reord <- tab3[permut,] 
  
  
  # The significance level alpha = 0.1 (by default) : the significance cutoff. 
  
  # If the adjusted p value cutoff will be a value other than 0.1,
  #alpha should be set to that value.
  
  # In the majority of analyses, an alpha of 0.05 is used as the cutoff 
  #for significance.
  res05 <- results(ddsSE2, alpha = 0.05)
  summary(res05) 
  

  # Exploring results 
  
 #MA-plot : 
  
#1) This plot allows us to evaluate the magnitude of fold changes and how they 
#are distributed relative to mean expression. Generally, we would expect to
#see significant genes across the full range of expression levels.  
  
#2) This is also a great way to illustrate the effect of logFC shrinkage:  
 
  plotMA(tab3, ylim=c(-2,2), colSig = "red") # (the unshrunken results).
 
#The genes that are significantly DE are colored to be easily identified.
#Points will be colored red if the adjusted p value is less than 0.1.   
#Points which fall out of the window are plotted as open triangles 
# pointing either up or down.  
# grey = not significant. 

  #After calling plotMA, one can use the function identify() to interactively 
  # detect the row number of individual genes by clicking on the plot.
  #One can then recover the gene identifiers by saving the resulting indices:
  idx <- identify(tab3$baseMean, tab3$log2FoldChange)
  rownames(tab3)[idx]   
  
# >>>>>>> LogFC shrinkage for visualization and ranking :<<<<<<<<<<
 
 # The shrunken logFC are useful for ranking genes by effect size
 #and for visualization.  
  
 # Shrinking improves accuracy and helps to avoid false positives.
  
 # LogFC shrinkage removes the noise associated with log2 fold changes 
 #from low count genes without requiring arbitrary filtering thresholds.
 
 # (N.B.): 
 # Statistical noise is unexplained variability within a data sample. 
 #Noisy data is data that is rendered meaningless by the existence
 # of too much variation
  
# To generate more accurate logFC estimates, 
#DESeq2 allows for the shrinkage of the logFC estimates toward
  #zero when the information for a gene is low.
  
  # We use shrinkage in two cases :
    #- Low counts ;
    #- High dispersion values.
  
# As with the shrinkage of dispersion estimates, logFC shrinkage uses
#information from all genes to generate more accurate estimates.
#Specifically, the distribution of LFC estimates for all genes is 
#used to shrink the LFC estimates of genes with little
#information or high dispersion toward more likely (lower) logFC estimates. 
  

   resultsNames(ddsSE2)
   
  BiocManager::install("apeglm")
  library(apeglm)
  # (' Approximate posterior estimation for GLM coefficients')
  # --> Effect size estimation. 
  
# "apeglm" provides Bayesian shrinkage estimators for effect sizes 
#for a variety of GLM models.  
  
resLFC <- lfcShrink(ddsSE2, coef="dex_trt_vs_untrt", type= "apeglm")
resLFC
 
 
 plotMA(resLFC, ylim=c(-2,2), colSig = "red") # (the shrunken results)


 
#__________________________  A bit of STATISTICS   __________________________

# While "statistical significance" shows that an effect exists in a study,
# "practical significance" shows that the effect is large enough to be meaningful
#in the real world.  
 
#Statistical significance is denoted "by p-values", whereas practical significance 
#is represented by "effect sizes". 

# To declare statistical significance, we need a criterion : The alpha (also known 
# as the Type I error probability or the significance level).
 
#Statistical significance alone can be misleading because
#it's influenced by the sample size. 
#Increasing the sample size always makes it more likely to find
# a statistically significant effect, no matter how small 
#the effect truly is in the real world.

#In contrast, effect sizes are independent of the sample size. 
#Only the data is used to calculate effect sizes.

#That's why it's necessary to report effect sizes in research papers to indicate 
#the practical significance of a finding. 
#The APA guidelines require reporting of effect sizes
# and confidence intervals wherever possible.
#__________________________________________________________________________

 
###~~~~~~~~~~ The Real Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 ### Following recommendations for single-cell analysis (vignette):
   

#Default values for DESeq2 were designed for bulk data and will not be
#appropriate for single-cell datasets (dominance of low and zero counts).
 
#For the best use of DESeq2 for single-cell datasets, we need some 
#additional settings.

 #1) For significance testing, 
 # we should use ( test="LRT" ) over the Wald test.
 
 #2) Set the following DESeq arguments to these values:
 #useT=TRUE, minmu=1e-6, and minReplicatesForReplace=Inf.
 
 
 
 #3) The default size factors are not optimal for single cell count matrices, 
 #instead consider setting sizeFactors from scran::computeSumFactors.
  BiocManager::install("scran")
  library(scran)
   
   #size factors for all cells:
   SFscran <- calculateSumFactors(airway, assay.type = "counts" )
  
  
  
   
   #DESeq() performs a default analysis through the steps:
     #- estimation of size factors: estimateSizeFactors.
     #- estimation of dispersion: estimateDispersions.
     #- Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest.
   
   
   
 #4) Set fitType = "glmGamPoi" to address the speed concerns.
   # One important concern for single-cell data analysis is the size of
   # the datasets and associated processing time.
   
  BiocManager::install("glmGamPoi")
  library(glmGamPoi)
  # Note that glmGamPoi's interface in DESeq2 requires use of test="LRT" 
  #and specification of a reduced design.
  
  
  airway$dex <- relevel(airway$dex, "untrt")
  ddsr <- DESeqDataSet(airway, design = ~ cell + dex) 
  
   recd <- DESeq(ddsr, fitType = "glmGamPoi", test="LRT", useT = T
        ,minmu=1e-6, minReplicatesForReplace=Inf, reduced=  ~ cell)

 #The reduced model is the full model with the variable(s) of interest 
 #removed.
 #If there is only the variable of interest in the full model, then
 # the reduced model will be as follows : reduced = ~ 1
 
   # Setting sizeFactors from scran:
  colData(recd)$sizeFactor <- SFscran 
  colData(recd)

   realres <-  results(recd)
   summary(realres)
   
   
   # MA-plot
   plotMA(realres, ylim=c(-2,2), colSig = "red")
   
   
  # These are the coefficients from the model,
  # we can specify them using 'coef' by name or number below;   
  resultsNames(recd)
  
  # LogFC Shrinkage :
  reaLFC <- lfcShrink(recd, coef= 2, type= "apeglm")
  reaLFC #ERROR !
  
  
  
  
  
 
  
  
  
  
  
  
  
  
  
