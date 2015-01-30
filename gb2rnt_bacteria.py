#!/usr/bin/env python
"""
Convert genbank to rnt.

USAGE:
gb2rnt_bacteria.py yourfile.gb > file.rnt

NOTE:
It's designed to work with gb files coming from GenBank. gene is used as gene_id and transcript_id (locus_tag if gene not present).
Only entries having types in allowedTypes = ['tRNA','tmRNA','rRNA','ncRNA'] are stored in rnt file. 

Based on a script by Leszek Pryszcz (lpryszcz@crg.eu)

AUTHOR
Andreas Sjodin

Version 0.2

"""

import os, sys
from datetime import datetime
from Bio      import SeqIO
import sys, re

def grep(l,s):
    return [i for i in l if s in i]


def gb2gtf( source='gb2gtf',allowedTypes=set(['tRNA','rRNA', 'ncRNA', 'tmRNA']) ):
  """
  """

  header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %  ("Location", "Strand", "Length", "PID", "Gene", "Synonym", "Code", "COG", "Product")
  print header
  
  for gb in SeqIO.parse( sys.argv[1],'gb' ):
    acc     = gb.id #gb.name #gb.description # # 
    skipped = 0
    skippedTypes = set()
    for f in gb.features:
    
      #process only gene and CDS entries
      if f.type not in allowedTypes:
        skipped += 1
        skippedTypes.add( f.type )
        continue
      
      #generate comments field
      if 'locus_tag' in f.qualifiers:
        #use locul tag as gene_id/transcript_id
        gene_id = f.qualifiers['locus_tag'][0]
        transcript_id = 'T' + gene_id
      elif 'gene' in f.qualifiers:
        gene_id = f.qualifiers['gene'][0]
        transcript_id = 'T' + gene_id

      elif 'label' in f.qualifiers:
        gene_id = f.qualifiers['label'][0]
        transcript_id = 'T' + gene_id
    
      comments = 'gene_id "%s"; transcript_id "%s"' % ( gene_id,transcript_id )
        
      if 'gene' in f.qualifiers:
        gene_name = f.qualifiers['gene'][0]
        comments += '; gene_name "%s"' % f.qualifiers['gene'][0]
        comments += '; transcript_name "%s-1"' % f.qualifiers['gene'][0]
      if 'protein_id' in f.qualifiers:
        comments += '; protein_id "%s"' % f.qualifiers['protein_id'][0]

      if 'product' in f.qualifiers:
        annotation = f.qualifiers['product'][0]
      else:
        annotation = ''

      
      #add external IDs  
      if 'db_xref' in f.qualifiers:
        db_xref= f.qualifiers['db_xref']
        if len(grep(db_xref,'GI')) != 0 :
          PID =  grep(db_xref,'GI')[0].replace('GI:', '')
        else:
          PID = ''
      
      #code strand as +/- (in genbank 1 or -1)
      if int(f.strand)>0: strand = '+'
      else:               strand = '-'
      
      #define gb
      """
      Location - start and stop of feature
      strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
      Length - length of feature
      PID - protein identifier
      Gene - name of gene
      Synonym - gene identifier
      Code -
      COG -
      Product - Annotation
      """
      rnalength = ((f.location.end.position-(f.location.start.position)) ) 
      COG = '-'
      code = '-'

      
      gtf = '%s..%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % ( f.location.start.position+1, f.location.end.position, strand, rnalength, PID, gene_name, gene_id, code, COG, annotation ) #f.frame,
      print gtf
      
      
    sys.stderr.write( "%s\tSkipped %s entries having types: %s.\n" % ( gb.id,skipped,', '.join(skippedTypes) ) )

if __name__=='__main__': 
  t0=datetime.now()
  gb2gtf()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
