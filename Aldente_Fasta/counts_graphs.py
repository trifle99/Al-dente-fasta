"""
Aldente Fasta: Tile sequences Module

This module contains functions for visualizing quality scores and gc content into graph plots.

"""


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import gzip

###################################
###################################
###################################

def avg_quals(file, out):
    """
    Creates a 2d np array from file and outputs a file containing average quality score of each read
    :param file: Input file containing quality scores (numerical value after <qual_convert>) of each base in each read
    :param out: Output file containing average quality score of each read
    :return: param out
    """
    nparray=np.loadtxt(file, delimiter=',') #//np.genfromtxt (from numpy import genfromtxt)
    #^ Loads all data into nparray, this can take a while depending file size
    reads_avg=(np.sum(nparray, axis=1))/(len(nparray[0])) #Calculate average of each read TODO: BUG ALERT: IMPORTANT: This assumes each reads are the same length>divide by length of that specific read>add a for loop to iterate over nparray, and sum that nparray[x]

    with open(out,'w') as out:
        #Write each read average into a file <out>
        for x in reads_avg:
            out.write(str(x)+'\n') #'\n' needed to add new line #ALSO> Change <str> support to write as <float> or <int> straight in, (To support int/float write in, I might have to write into csv files)
            #IF writing csv files, add a out.write('') as a header to make it empty, so the values can index correctly from 1,2,3...etc

    del nparray #Deleting just to clear memory space

    return

###################################
###################################
###################################

def gc_per_read(file,out):
    """
    Takes a fastq file and writes the GC percent of each read into <out> file
    :param file: Fastq gz compressed file
    :param out: File containing GC percent of each read
    :return: <param out>
    """
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
                    lenseq=len(line)
                    for bases in line:
                        if bases in gc:
                            gc_count+=1 #increment by one for each base (g or c) found in read
                    outfile.write((str(100*(gc_count/lenseq)))+'\n') #string to write in
                    gc_count=0


    return

###################################
###################################
###################################

def counts_graph(file, histo_or_kde, gc_or_qual):
    #TODO: Add error exceptions...
    #^ todo: file detection for correct file format
    """
    Takes file containing average quality scores of each read and visualize into graphs
    :param file: <aldente_fasta.avg_quals> or <gc_per_read>
    :param histo_or_kde: <0,1,2> Shows different graphs according to number
    :return: Graph plot
    """

    if gc_or_qual=='qual':
        xlimplot=40
        xlabelplot='Quality Score'
    elif gc_or_qual=='gc':
        xlimplot=100
        xlabelplot='GC content'

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

        plt.xlim(0,xlimplot) #Set limits on x-axis
        plt.title(histo_title)
        plt.xlabel(xlabelplot)

        plt.show()

    elif histo_or_kde==1:

        plt.figure(figsize=(20, 10))

        sns.kdeplot(data)

        plt.xlim(0, xlimplot)
        plt.title('Density estimation of '+xlabelplot+' of every read')
        plt.show()
    else: print('Enter a value of 0 or 1 for <histo_or_kde>. 0 for histoplot, 1 for kde plot. (2 for custom histoplot with added kde overlay)')
    return

###################################
###################################

