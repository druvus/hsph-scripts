library(GEOquery)
library(SRAdb)

datadir <- “PATH”
GEOfile <- "GEOid.txt"

#Set working directory
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


#Get a conversion table
conversion <- sraConvert( SRA, sra_con = sra_con )
metadata <- data.frame(GEO=GEOids[,1], conversion)


#Export metadata to file
write.table(metadata, file = "metadata.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = F )

#Get the run ids
SRR<- unique(conversion$run)

#Download fastq files
getSRAfile( SRR, sra_con, fileType = "fastq")
#getSRAfile( SRR, sra_con, fileType = "sra")

