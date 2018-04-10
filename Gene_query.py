"""
Gene Query of Keeney Sequence

Works require annotation.py in the same directory!

This script accepts four inputs from commandlines and outputs a fasta file of 
extracted sequences:
    -s: sk1 genome sequence filename
    -g: the gff filename
    -i/-n: gene id/gene name
    -r: the range on either side of the gene
    
This scripts runs by the commandline:
    e.g. python Gene_query.py -s SK1_MvO_V1___GENOME/sk1_MvO_V1.fasta -g MvO_v1_annotation.gff -n hop1 -r 100
"""

##################################################################################
# Modules

from Bio import SeqIO
from annotation import GFFfile
import optparse

##################################################################################
# Functions

def get_reverse_complement (seq):
    """
    This function takes a sequence (str) as input and returns
    the reverse of its complementary sequence (str).
    """
    # imports string module
    import string
    # changes input into upper case
    seq = seq.upper()
    # creates a table of complementary base pairs
    compbase = string.maketrans('ATCG', 'TAGC')
    # returns the complementary bases according to the compbase
    # table, and put into reversed order
    return seq.translate(compbase)[::-1]



def gene_query(sk1_filename, gff_filename, range_input, gene_id=None, gene_name=None):
    # reads the sk1 fasta file    
    sk1 = list(SeqIO.parse(sk1_filename, 'fasta'))
    
    # reads the GFF file
    g = GFFfile(gff_filename)
    
    # finds all the genes in GFF file
    genes = [i for i in g.annotations if i.feature == 'sacCer2_sgdGene']
    
    # if input is the gene id
    if gene_id:
        gene_id = gene_id.upper()
        for i in genes:
            if i.attributes['ID'] == gene_id:
                chromosome = i.seqname
                start = i.start
                end = i.end
                strand = i.strand
                name = i.attributes['ID']+'-' +i.attributes['Name']
        # extract the sequence with a range option
        for i in sk1:
            if i.id == chromosome:
                seq = str(i.seq[start-1-range_input:end+range_input])
        # if strand is '-', outputs the reverse complementary sequence
        if strand == '-':
            seq = get_reverse_complement(seq)
        
        # saves the sequence into a auto named fasta file
        f=open(name+'.fasta','w')
        f.write('>'+name+':'+chromosome+'-'+str(start)+'-'+str(end)+':'+strand+'\n'+seq+'\n')
        f.close()
    
    # if input is a gene name
    if gene_name:
        gene_name = gene_name.upper()
        for i in genes:
            if i.attributes['Name'] == gene_name:
                chromosome = i.seqname
                start = i.start
                end = i.end
                strand = i.strand
                name = i.attributes['ID']+'-' +i.attributes['Name']
                gene_id = i.attributes['ID']
        # extract the sequence with a range option
        for i in sk1:
            if i.id == chromosome:
                seq = str(i.seq[start-1-range_input:end+range_input])
        # if strand is '-', outputs the reverse complementary sequence
        if strand == '-':
            seq = get_reverse_complement(seq) 
            
        # saves the sequence into a auto named fasta file
        f=open(name+'.fasta','w')
        f.write('>'+name+':'+chromosome+'-'+str(start)+'-'+str(end)+':'+strand+'\n'+seq+'\n')
        f.close()
            
 def calculate_length(gtf_filename):
    f=open(gtf_filename,'r')
    lines=f.readlines().strip().split('\t')
    f.close()   
###############################################################################
# Main

# parse object for managing input options.      
parser = optparse.OptionParser()

# essential data, defines a commanline option "-i"
parser.add_option('-s', dest = 'sk1_filename', default = '', help = 'This input\
 is the fasta file of the sk1 genome sequences')
parser.add_option('-g', dest = 'gff_filename', default = '', help = 'This input\
 is the gff file of the sk1 genes')
parser.add_option('-i', dest = 'gene_id', default = '', help = 'This input\
 is the id of the gene')
parser.add_option('-n', dest = 'gene_name', default = '', help = 'This input\
 is the name of the gene') 
parser.add_option('-r', dest = 'range_input', default = '0', help = 'This input\
 is the range on either side of the gene sequence') 

# loads the inputs
(options, args) = parser.parse_args()

# reads the inputs from command lines
sk1_filename = options.sk1_filename
gff_filename = options.gff_filename
gene_id = options.gene_id
gene_name = options.gene_name
range_input = int(options.range_input)

# runs the function
gene_query(sk1_filename, gff_filename, range_input, gene_id, gene_name)
    
    
    
    
    
    
