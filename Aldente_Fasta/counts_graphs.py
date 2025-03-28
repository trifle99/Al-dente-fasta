"""
MODULE INFO
"""


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import gzip

###################################
###################################
###################################

def counts_graph(file, histo_or_kde):
    #TODO: Add error exceptions...
    #^ todo: file detection for correct file format
    """
    Takes file containing average quality scores of each read and visualize into graphs
    :param file: <aldente_fasta.avg_quals>
    :param histo_or_kde: <0,1,2> Shows different graphs according to number
    :return: Graph plot
    """

    data=np.loadtxt(file) #Load data from <avg_quals>

    if histo_or_kde==2:
        kde_TF=True
        histo_title='Histogram + KDE comparisson plot'
    else:
        kde_TF=False
        histo_title='Histogram plot'

    if histo_or_kde in (0, 2): #Checks if histogram plot should have KDE line as well

        plt.figure(figsize=(20,10))

        sns.histplot(data,
                     #bins=40, #How many bars
                     kde=kde_TF, #<histo_or_kde>
                     edgecolor='white', #Better visualization for histogram+KDE lines
                     linewidth=0.5) #Better visualization of KDE line

        #plt.xlim(0,40) #Set limits on x-axis
        plt.title(histo_title)
        plt.xlabel('Quality Score')

        plt.show()

    elif histo_or_kde==1:

        plt.figure(figsize=(20, 10))

        sns.kdeplot(data)

        #plt.xlim(0, 40)
        plt.title('Density estimation of quality scores among reads')
        plt.show()
    else: print('Enter a value of 0 or 1 for <histo_or_kde>. 0 for histoplot, 1 for kde plot. (2 for custom histoplot with added kde overlay)')
    return

###################################
###################################
###################################

def gc_per_read(file,out):


    line_no=0
    gc=['G','C']
    gc_count=0

    with gzip.open(file,'r') as fastq: #gzip
        with open(out,'w') as outfile:

            for line in fastq:
                line=line.rstrip()
                line=line.decode() #decode because file was written in bytes, but if im writing in 'wb' i dont need to decode

                line_no+=1

                if (line_no % 2 == 0) and (line_no % 4 == 2):  # If statement to check every 2nd sequence line of each reads (fastqfile has a 4 line repeating format containing read information)

                    for bases in line:
                        if bases in gc:
                            gc_count+=1
                    outfile.write((str(gc_count))+'\n') #string to write in
                    gc_count=0


    return


# file=('C:\\Users\\gaura\\Documents\\Python Learning\\Projects\\fastq parse\\SRR2093871_1.fastq.gz')
# file2=('C:\\Users\\gaura\\Documents\\Python Learning\\Projects\\fastq parse\\out\\srr2093871_1_subset.txt')
#
# out=('C:\\Users\\gaura\\Documents\\Python Learning\\Projects\\fastq parse\\out\\gc_per_read.txt')
# out1=('C:\\Users\\gaura\\Documents\\Python Learning\\Projects\\fastq parse\\out\\avg_quals')
# import time
#
# start=time.time()
# counts_graph(out1,0)
# end=time.time()
#
# print(end-start)

#ADD IF STATEMENTS TO <COUNTS_GRAPH> AND ANOTHER PARAM ARG TO DISTINGUISH USING COUNTS FOR AVG_QUALITY_SCORES OR GC_MEAN_CONTENT>
#THEN TAILOR IF STATEMENTS FOR PLOT ARGUEMENTS (BINS, AXIS-LIMS, AXIS LABELS, PLOT TITLE ETC)