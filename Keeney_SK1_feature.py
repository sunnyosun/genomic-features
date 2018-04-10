"""
Promoters in Keeney SK1
"""
##################################################################################
# Modules

from Bio import SeqIO
from annotation import GFFfile
import optparse

##################################################################################
# Functions

def promoter_region(gff_filename, range_input):
    # reads the GFF file
    g = GFFfile(gff_filename)
    
    # finds all the genes in GFF file
    genes = [i for i in g.annotations if i.feature == 'sacCer2_sgdGene']

    output = []
    for i in genes:
        if i.strand == '+':
            start = i.start - range_input
            end = i.start
            if start < 0:
                start = 1
        elif i.strand == '-':
            start = i.end
            end = i.end + range_input
        output.append(i.seqname+'\t'+i.source+'\t'+'GenePromoter'+'\t'+str(start)+'\t'+str(end)+'\t'+i.frame+'\t'+i.strand+'\t'+i.score+'\t'+'ID='+i.attributes['ID']+';Name='+i.attributes['Name'])
    f=open('SK1_Promoters.gff','w')
    f.write('\n'.join(output))
    f.close()




def three_end_region(gff_filename, range_input):
    # reads the GFF file
    g = GFFfile(gff_filename)
    
    # finds all the genes in GFF file
    genes = [i for i in g.annotations if i.feature == 'sacCer2_sgdGene']

    output = []
    for i in genes:
        if i.strand == '+':
            start = i.end - range_input
            end = i.end
            if start < 0:
                start=1
        elif i.strand == '-':
            start = i.start
            end = i.start + range_input
        output.append(i.seqname+'\t'+i.source+'\t'+'Gene3end'+'\t'+str(start)+'\t'+str(end)+'\t'+i.frame+'\t'+i.strand+'\t'+i.score+'\t'+'ID='+i.attributes['ID']+';Name='+i.attributes['Name'])
    f=open('SK1_3_end.gff','w')
    f.write('\n'.join(output))
    f.close()



def three_end_region(gff_filename, range_input):
    # reads the GFF file
    g = GFFfile(gff_filename)    
    # finds all the genes in GFF file
    genes = [i for i in g.annotations if i.feature == 'sacCer2_sgdGene']
    output = []
    for i in genes:
        if i.strand == '+':
            start = i.end - range_input
            end = i.end + range_input
            if start < 0:
                start=1
        elif i.strand == '-':
            start = i.start - range_input
            end = i.start + range_input
        output.append(i.seqname+'\t'+i.source+'\t'+'Gene3end'+'\t'+str(start)+'\t'+str(end)+'\t'+i.frame+'\t'+i.strand+'\t'+i.score+'\t'+'ID='+i.attributes['ID']+';Name='+i.attributes['Name'])
    f=open('SK1_3_end.gff','w')
    f.write('\n'.join(output))
    f.close()




def sk1_clean(gff_filename):
    # reads the GFF file
    g = GFFfile(gff_filename)
    
    # remove gaps or patches
    genes = [i for i in g.annotations if i.feature == 'sacCer2_sgdGene' or i.feature == 'sacCer2_sgdOther']

    # unique ones
    uniq = [i for i in genes if len(i.attributes) == 4 and i.feature == 'sacCer2_sgdGene']
    output1 = []
    for i in genes:
        if len(i.attributes) == 4:
            output1.append(i.seqname+'\t'+i.source+'\t'+i.feature+'\t'+str(i.start)+'\t'+str(i.end)+'\t'+i.frame+'\t'+i.strand+'\t'+i.score+'\t'+'ID='+i.attributes['ID']+';Name='+i.attributes['Name'])
    f=open('SK1_uniq.gff','w')
    f.write('\n'.join(output1))
    f.close()

            
    # gene id starts with Y
    dup_Y = [i for i in genes if len(i.attributes) == 5 and i.attributes['ID'][0] == 'Y']
    temp=[]
    for i in dup_Y:
        temp.append(i.seqname+'\t'+i.source+'\t'+i.feature+'\t'+str(i.start)+'\t'+str(i.end)+'\t'+i.frame+'\t'+i.strand+'\t'+i.score+'\t'+'ID='+i.attributes['ID']+';Name='+i.attributes['Name'])
    f=open('SK1_Y.gff','w')
    f.write('\n'.join(temp))
    f.close()


    # not starts with Y
    others = [i for i in genes if len(i.attributes) == 5 and i.attributes['ID'][0] != 'Y']
    temp=[]
    for i in others:
        temp.append(i.seqname+'\t'+i.source+'\t'+i.feature+'\t'+str(i.start)+'\t'+str(i.end)+'\t'+i.frame+'\t'+i.strand+'\t'+i.score+'\t'+'ID='+i.attributes['ID']+';Name='+i.attributes['Name'])
    f=open('SK1_others.gff','w')
    f.write('\n'.join(temp))
    f.close()
                

    
