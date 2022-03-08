import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



#modified from Claire Dubin's script, found at https://github.com/clairedubin/thermotolerance/blob/master/MK_analysis.ipynb

##PARAMETERS###
gff='saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'
infile='TajimaD_by-gene_25Sake.tsv'

essential_file='essential.csv' #essential genes - from Winzeler et al., 1999


##figname='25Sake_TD_resample_barseq_tophits.png'
##goi=['YGR198W', 'YMR207C', 'YGL082W', 'YNL049C', 'YDL035C', 'YDR508C', 'YBR136W',
## 'YML099C', 'YPL254W', 'YIL152W', 'YKL017C', 'YGR140W', 'YJR127C', 'YDR375C',
## 'YOR091W', 'YLR397C', 'YNL132W', 'YMR078C', 'YLR422W', 'YMR125W', 'YOR371C',
## 'YMR094W', 'YMR167W', 'YDR103W', 'YDR318W', 'YAL026C', 'YDR180W', 'YOR092W',
## 'YDR235W', 'YER151C', 'YMR275C', 'YKL114C', 'YOL081W', 'YPR049C', 'YGL095C',
## 'YDR456W', 'YKL197C', 'YIL068C', 'YOR326W', 'YNR045W', 'YJR107W', 'YPL268W',
## 'YJL062W', 'YCR042C'] #barseq_tophits

##figname='25Sake_TD_resample_popgen_hits.png'
##goi=['YGR098C','YMR168C','YKR054C', 'YPL174C','YPR164W','YCR042C','YMR016C',
##     'YJR135C','YJL025W','YDR443C','YLR397C','YHR023W','YDR180W','YKL134C'] #popgen_hits

##figname='25Sake_TD_resample_Weiss_goi.png'
##goi=['YLR397C','YGR098C','YKR054C','YHR023W',
##     'YDR180W','YCR042C','YNL172W','YMR168C'] #thermotolerance loci from Weiss et al 2018

figname='25Sake_TD_resample_shortcomm_goi.png'
goi=['YLR397C','YGR098C','YHR023W','YMR168C'] #top 4 thermotolerance loci from Weiss et al 2018


##FUNCTIONS###

def ParseFromGFF(gff):
    '''
    Parses gff
    Output: dict of {yName:{geneName}}
    '''
    ann_dict={}
    
    f = open(gff)
    lines=[]
    for line in f:
        if line[0]!='#': #skip header rows
            row_data=line.split("\t")
            if row_data[2]=='gene':
                info=row_data[8].split(";")
                yName=info[0].split('=')[1]
                geneName=info[2].split('=')[1]
                ann_dict[yName]=geneName
    f.close()
    return ann_dict

def resample_med(df1, df2, col, essential_genes, direction, n=10000, graph=False, noisy=True,figName='resampling_distribution.png'):

    # inputs: df1 = DataFrame of candidate genes, needs column of gene names named 'gene'
    #         df2 = DataFrame of all genes to resample from, needs column of gene names named 'gene'
    #         col = name of df column with statistic to resample
    #         essential = list of essential genes
    #         direction = ['less_than', 'greater']: test whether statistic is less than or greater than random samples
    #         n = number of random samples to generate
    #         graph = show KDE plot of resampling distribution
    #         noisy = print statistics

    
    # output: p value (proportion of random samples with statistic ['less_than', 'greater_than'] or equal to candidates


    actual_med = df1[col].median()
    
    essential_count = len([i for i in df1['gene'] if i in essential_genes])
    nonessential_count = len(df1[col]) - essential_count

    essential_df = df2[df2['gene'].isin(essential_genes)]
    nonessential_df =  df2[~df2['gene'].isin(essential_genes)]

    sample_meds = []
    results = 0

    for _ in range(n):
        s = random_sample_med(essential_df, essential_count, nonessential_df, nonessential_count, direction, actual_med, col)
        results += s[0]
        sample_meds += [s[1]]
    
    if graph:  
        #sns.displot(sample_meds, kind='kde') #new seaborn
        sns.distplot(sample_meds, kde=True) #old seaborn
        plt.axvline(x=actual_med, color='r', label = 'True median')
        plt.xlabel('Sample medians')
        plt.ylabel('Frequency')
        plt.title(col + ' resampling distribution, p= ' + str(results/n))
        plt.savefig(figName)

    if noisy:
        print('candidate gene median {}: {}'.format(col, actual_med))
        print('essential count: {}; nonessential_count: {}'.format(essential_count, nonessential_count))
        print('resampling pool size: {}'.format(df2.shape[0]))
        print('p = {}'.format(results/n))
    
    return results/n

def random_sample_med(essential_df, essential_count, nonessential_df, nonessential_count, direction, actual_med, col):

        sample = essential_df.sample(n=essential_count)
        sample = sample.append(nonessential_df.sample(n=nonessential_count))
        sample_med = sample[col].median()

        rv = 0
        if direction == 'less_than':
            if sample_med <= actual_med:
                rv += 1

        elif direction == 'greater_than':
            if sample_med >= actual_med:
                rv += 1

        return rv, sample_med

###RUN###

#parse essential file
essential = pd.read_csv(essential_file, header=None)
essential[1] = essential[1].str.strip('\t')
essential_genes = [i.split(' ')[0] for i in essential[1]]


#build gene dict
ann_dict=ParseFromGFF(gff)

gene_dict = {}
for yName in goi:
    geneName=ann_dict[yName]
    gene_dict[yName]=geneName

#read TD input
df = pd.read_csv(infile,sep='\t')
df.dropna(inplace=True)
print(df)
candidates = df[df['gene'].isin(gene_dict.keys())]
candidates['name'] = [gene_dict[g] for g in candidates['gene']]
candidates.sort_values(by=['TajimaD'],inplace=True)

##Resample TD
np.random.seed(7777)

#print('missing: ', [gene_dict[i] for i in gene_dict.keys() if i not in candidates['gene'].tolist()])
print('TD p={}'.format(resample_med(candidates, df, 'TajimaD', essential_genes, direction='less_than', graph=False,n=10000,figName=figname)))

#plot TD by gene

fig,ax=plt.subplots(1,1,figsize=(10, 5))
plt.bar([gene_dict[gene] for gene in candidates.gene], candidates['TajimaD'])
plt.ylabel("Tajima's D")
plt.tick_params(axis='x', labelbottom=True, labeltop=False, bottom=False)
#plt.xticks(rotation=90)
plt.axhline(y=candidates['TajimaD'].median(),linewidth=1, color='orange', label='Candidate median')
plt.axhline(y=df['TajimaD'].median(),linewidth=1, color='purple', label='Genome wide median')

ax.xaxis.tick_top()
plt.setp(ax.get_xticklabels(), ha="right", rotation=90)

plt.legend(loc='upper right', bbox_to_anchor=(1, 0),
      ncol=3, fancybox=False, shadow=False)
plt.savefig(figname)
                                                    

                                                      
