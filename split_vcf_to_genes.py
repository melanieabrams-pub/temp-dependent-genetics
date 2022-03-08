import vcf
import os

#usage: #python split_vcf_to_genes.py
#go from split chrom to gzipped gene VCFs

###PARAMETERS###
gff='/usr2/people/mabrams/Amended_Genomes/S288C/saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'

##split_chrom_dir='/usr2/people/mabrams/data/all_1011_popgen/split_chrom/'
##infile_prefix='_'
##out_dir='/usr2/people/mabrams/data/all_1011_popgen/split_chrom/split_gene/'
##outfile_extension='_1011.vcf.gz'


##split_chrom_dir='/usr2/people/mabrams/data/1WineEuropean_1011genomes/split_chrom/'
##infile_prefix='merged_1WineEuropean_'
##out_dir='/usr2/people/mabrams/data/1WineEuropean_1011genomes/split_chrom/split_gene/'
##outfile_extension='_1Wine.vcf'

##split_chrom_dir='/usr2/people/mabrams/data/M3_MosaicRegion3/split_chrom/'
##infile_prefix='merged_M3_MosaicRegion3_'
##out_dir='/usr2/people/mabrams/data/M3_MosaicRegion3/split_chrom/split_gene/'
##outfile_extension='_M3.vcf'

##split_chrom_dir='/usr2/people/mabrams/data/3BrazilianBioethanol_1011genomes/split_chrom/'
##infile_prefix='merged_3BrazilianBioethanol_'
##out_dir='/usr2/people/mabrams/data/3BrazilianBioethanol_1011genomes/split_chrom/split_gene/'
##outfile_extension='_3Brazil.vcf'

##split_chrom_dir='/usr2/people/mabrams/data/25Sake_1011genomes/split_chrom/'
##infile_prefix='merged_25Sake_'
##out_dir='/usr2/people/mabrams/data/25Sake_1011genomes/split_chrom/split_gene/'
##outfile_extension='_25Sake.vcf'

split_chrom_dir='/usr2/people/mabrams/data/8MixedOrigin_1011genomes/split_chrom/'
infile_prefix='merged_8MixedOrigin_'
out_dir='/usr2/people/mabrams/data/8MixedOrigin_1011genomes/split_chrom/split_gene/'
outfile_extension='_8Mixed.vcf'


roman_numerals_in_gff=True

###FUNCTIONS###


def ParseFromGFF(gfffile):
    '''
    Parses gff
    Output: dict of {gene:[start,stop,chrom]}
    '''
    roman_to_numerals={
        'chrI':'chromosome1','chrII':'chromosome2','chrIII':'chromosome3','chrIV':'chromosome4',
        'chrV':'chromosome5','chrVI':'chromosome6','chrVII':'chromosome7','chrVIII':'chromosome8',
        'chrIX':'chromosome9','chrX':'chromosome10','chrXI':'chromosome11','chrXII':'chromosome12',
        'chrXIII':'chromosome13','chrXIV':'chromosome14','chrXV':'chromosome15','chrXVI':'chromosome16',
        'chrMito':'chromosomeMito','2-micron':'chromosome2u'}
    
    ann_dict={}
    
    f = open(gfffile)
    lines=[]
    for line in f:
        if line[0]!='#': #skip header rows
            row_data=line.split("\t")
            #print(row_data)
            if roman_numerals_in_gff==True: #convert roman numerals if needed
                chrom=roman_to_numerals[row_data[0]]
            else:
                chrom=row_data[0]
            start=int(row_data[3])
            stop=int(row_data[4])
            info=row_data[8].split(";")
            yName=info[0].split('=')[1]
            if yName[0]=='Y' and len(yName)>5:
                ann_dict[yName]=start,stop,chrom
    f.close()

    #ann_dict={'YLR397C':ann_dict['YLR397C']} #for debug
    
    return ann_dict

def getGeneVCF(gene,start,stop,chrom):
    '''searches by pos in chrom'''
    vcf_in=split_chrom_dir+infile_prefix+chrom+'.vcf.gz'
    vcf_out=out_dir+gene+outfile_extension

    print(gene,start,stop)

    cmd='tabix -h '+vcf_in+' '+chrom+':'+str(start)+'-'+str(stop)+'>'+vcf_out #thanks to https://www.biostars.org/p/46331/
    #print(cmd)
    os.system(cmd)
    

##    with open(vcf_in, 'r') as f:  #The commented out secion below ran very slowly so replaced it
##        vcf_reader = vcf.Reader(f, 'r')
##        record_counter=0
##        with open(vcf_out, 'w') as wf:
##            vcf_writer = vcf.Writer(wf, vcf_reader)
##
##            
##            for record in vcf_reader:
##                record_counter+=1
##                if record_counter%1000==0:
##                    print(record_counter)
##                if start<=record.POS<=stop:
##                    vcf_writer.write_record(record)
##                if stop<record.POS:
##                    os.system('gzip vcf_out')
##                    return
##    os.system('gzip vcf_out')
    return


###RUN###
            
ann_dict=ParseFromGFF(gff)
print('imported annotations')

##for gene in ann_dict: for debug
##    if gene=='YLR397C':
##        getGeneVCF(gene, ann_dict[gene][0],ann_dict[gene][1],ann_dict[gene][2])

for gene in ann_dict:
    getGeneVCF(gene, ann_dict[gene][0],ann_dict[gene][1],ann_dict[gene][2])




