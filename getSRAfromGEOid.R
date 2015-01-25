library(GEOquery)
library(SRAdb)
datadir <- "/PATH/SRA"
GEOfile <- "GEOid.txt"

setwd(datadir)
SRA <- c()

#Download SRA database
sqlfile <- 'SRAmetadb.sqlite'
sqlfile <- getSRAdbFile()

#Open connection to local SRA db
sra_con <- dbConnect(SQLite(),sqlfile)

#Import GEO id from text file
GEOids <- read.table(GEOfile, header = FALSE, sep = "\t", quote = "")

#Queary GEO to get SRX information
for (i in as.character(GEOids[,1])){
  gsm <- getGEO(as.character(i))
  links <- gsm@header$relation
  SRAlink <- links[grep("SRA", links)]
  SRX <- unlist(lapply(strsplit(SRAlink, "term="), function (x) x[2]))
  SRA <- c(SRA, SRX)

}
SRA <- sort(unique(SRA))

#Get a conversion table
conversion <- sraConvert( SRA, sra_con = sra_con )

#Get the run ids
SRR<- conversion$run

#Download fastq files
getSRAfile( SRR, sra_con, fileType = "fastq")
getSRAfile( SRR, sra_con, fileType = "sra")

