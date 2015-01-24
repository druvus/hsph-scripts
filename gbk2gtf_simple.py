#!/usr/bin/env

from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import sys, re

# Functions
def grep(l,s):
    return [i for i in l if s in i]

# Program
source = sys.argv[3]
frame = "."
feattype = "exon"
score = "."
feats=set()
exonnumber = 1

f = open(sys.argv[2],'w')

for rec in SeqIO.parse( sys.argv[1], "genbank"):
    if rec.features:
        for feature in rec.features:
            if feature.type == "CDS":
                feats.add(feature.type)
                seqname = rec.id
                if feature.strand == 1:
                    strand = "+"
                elif feature.strand == -1:
                    strand = "-"
                else:
                    strand = "."
                start = feature.location.start.position+1
                end = feature.location.end.position
                if 'gene' in feature.qualifiers:
                    genename = feature.qualifiers['gene'][0]
                    transcriptname = genename + "-1"
                if 'product' in feature.qualifiers:
                    annotation = feature.qualifiers['product'][0]
                if 'db_xref' in feature.qualifiers:
                    db_xref= feature.qualifiers['db_xref']
                    geneid =  grep(db_xref,'GeneID')[0].replace('GeneID:', '')
                    transcriptid = "T" + geneid
                else:  
                    continue
                output = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%s\"; gene_name \"%s\"; transcript_name \"%s\"\n" % (seqname, source, feattype, start, end, score, strand, frame, geneid, transcriptid, exonnumber, genename, transcriptname)
                f.write(output) 

f.close()

print ("Converted : " + rec.id)
 
