"""
Aldente Fasta Module:

This module contains functions related around analysing sequencing data (fasta/fastq files).

"""
test=123456789

# **********************************************************************************************************************
# **********************************************************************************************************************
# **********************************************************************************************************************
# *Create a working pipeline for implenetations>run <summary> first and then <...> and then <...> etc...
# *Share on github
# *Turn module into a package where I add another project of protein/gene annotation analyzer OR I just create a new module/web app for that
# **********************************************************************************************************************
# **********************************************************************************************************************
# **********************************************************************************************************************
# IMPORTANT !!! IMPORTANT !!! NEED TO INSTALL NUMPY AND MATPLOTLIB FOR FIRST TIME PROJECT INSTALLATIONS !!! IMPORTANT !!! IMPORTANT
# **********************************************************************************************************************
# **********************************************************************************************************************
# **********************************************************************************************************************

import csv #Needed for opening and writing files
import gzip #Needed for opening compressed gz files
import numpy as np #Needed for creating 2d arrays and processing it
import matplotlib.pyplot as plt #Needed for plotting graphs for data analysis

#TODO: Might include line_no as a global variable instead of nested within each function (if i do this, i have to remember to reset line_no to 0 so theres no errors when calling other functions including line_no
#TODO: Add error exceptions (try, except errors)

def summary(file):
    """
    Gives a basic overview of your fastqfile
    :param file: FASTQ file (gz compression) TODO: Add file support/detection for other compression types
    :return: summary (dict-containing basic stats)
    """

    line_no=0 #Variable to act as index by counting line number in fastqfile

    #Variables to store nucleotide base information
    a=0
    c=0
    t=0
    g=0
    n=0

    with gzip.open(file, 'r') as fastqfile: #Open fastqfile #Fastqfiles should be written in binary so I could technically specifyy 'rb' to read bytes but 'r' should read all types anyway > i specify byte in counts anyway @38
        for line in fastqfile: #Read all lines in fastqfile by iterating over all lines> TODO: might be a way to speed this up or make it more memory efficient> using yield and generators?
            line_no+=1 #Count line number

            if (line_no%2==0) and (line_no%4==2): #If statement to check every 2nd sequence line of each reads (fastqfile has a 4 line repeating format containing read information)
                #TODO: Make more efficient by calling count only once >perhaps create a list which contains my counts (A,C,T,G) and iterate over this list to increment by 1 if line detects list contents
                #Simple counter for nucleotide bases
                a+=line.count(b'A') #Need to specifiy b'X' because it has to read in bytes
                c+=line.count(b'C')
                t+=line.count(b'T')
                g+=line.count(b'G')
                n+=line.count(b'N')

    #Basic operations to extract information:
    reads=(line_no/4)
    gc=g+c
    total_seq_len=a+c+t+g+n
    gc_percent=round(((gc/total_seq_len)*100),2)

    summary_dict={'reads':reads,
                  'A count': a,
                  'C count': c,
                  'T count': t,
                  'G count': g,
                  'N count': n,
                  'GC%': gc_percent,
                  'Total Bases': total_seq_len
                  }

    return summary_dict

#res=summary('C:\\Users\\gaura\\Documents\\Python Learning\\Projects\\fastq parse\\SRR2093871_1.fastq.gz')

#print(res) #{'reads': 7022387.0, 'A count': 222500881, 'C count': 128502570, 'T count': 224049496, 'G count': 127141349, 'N count': 44404, 'GC%': 36.4, 'Total Bases': 702238700}

#25-27seconds

#SRR2093871_1.fastq.gz

def qual_out(file, out):
    """
    Reads quality scores of each read from fastq file and outputs a file containing only these quality scores
    :param file: Fast q file (gz file format) TODO: Add file detect and other compression file support
    :param out: Text file containing quality scores (4th line of each read, in ASCII format>not decoded)
    :return: param out
    """

    line_no=0 #Variable to index by counting each lines

    with gzip.open(file, 'r') as fastq:
        with open(out, 'wb') as outfile: #Writing as byte because original gz file is written and loaded in bytes

            for line in fastq:
                line_no+=1 #Index counter
                if (line_no%4==0): #If statement to check 4th line (quality score) of every read
                    outfile.write(line)

#19seconds for 7mill reads

def qual_convert(file, out, enc):
    """
    Converts quality scores into numerical values (Phred Quality Score, depending on the specified encoding) into an output file
    :param file: File containing quality scores (output from <qual_out> func)
    :param out: File
    :param enc: Encoding used to decode quality score (should only accept 33 or 64 as those are the only encoding used) TODO: Enforce this encoding to only accept values of 33 or 64
    :return: param out
    """

    with open(file,'r') as filein:
        with open(out, 'w') as fileout:

            csv_writer=csv.writer(fileout, delimiter=',')

            for line in filein:
                line=line.rstrip() #Removing trailing artefacts at end of each line
                numbers=[(int(ord(x)-enc)) for x in line] #Converts the quality score characters into numerical values depending on the encoding used (33 or 64)
                csv_writer.writerow(numbers)

#118seconds to do 646.7mb into 2gigs file (for 7mil reads of SRR2093871_1.fastq.gz)

#TODO: might have to incorporate these two <avg_quals> and <quals_boxplot> into one for efficiency's sake, because they both need to load nparray which takes a long time depending on file sizes
#TODO: do some extra analysis with these average quals of each read>maybe do another boxplot or another graph idk
def avg_quals(file, out):
    """
    Creates a 2d np array from file and outputs a file containing average quality score of each read
    :param file: Input file containing quality scores (numerical value after <qual_convert>) of each base in each read
    :param out: Output file containing average quality score of each read
    :return: param out
    """
    nparray=np.loadtxt(file, delimiter=',') #//np.genfromtxt (from numpy import genfromtxt)
    #^ Loads all data into nparray, this can take a while depending file size
    reads_avg=(np.sum(nparray, axis=1))/(len(nparray[0])) #Calculate average of each read TODO: IMPORTANT: This assumes each reads are the same length>divide by length of that specific read>add a for loop to iterate over nparray, and sum that nparray[x]

    with open(out,'w') as out:
        #Write each read average into a file <out>
        for x in reads_avg:
            out.write(str(x)+'\n') #'\n' needed to add new line #ALSO> TODO: Change <str> support to write as <float> or <int> straight in, (To support int/float write in, I might have to write into csv files)
            #TODO: IF writing csv files, add a out.write('') as a header to make it empty, so the values can index correctly from 1,2,3...etc

    del nparray #Deleting just to clear memory space

#Benchmarking to measure time taken to complete this func for large datasets todo
#Might take a long while as I have to load huge 2d nparray AND then write average of each qual into file <out>

#TODO: Histogram plots for outliers? Theres an analysis on boxplot outliers which recommends doing histogram plot for outliers>look more into this data analysis for outliers
def quals_boxplot(file):
    #TODO: Error exceptions to check for read seq lengths, as for this boxplot to run>each read has to be same seq length
    """
    Takes nparray from file containing quality scores, and returns a quality score boxplot of each base position in each read
    :param file: File containing quality scores of each base in each read, read in as a 2d np array
    :return: NONE
    """
    nparray=np.loadtxt(file, delimiter=',') #//np.genfromtxt (from numpy import genfromtxt)
    #^ Loading all qual score data into 2d nparray, can take a while depending on file size

    plt.figure(figsize=(20,5)) #TODO: customize figsizes according to how large datasets are, can check for how many bases in each read to make figsizes wider(stretched) or narrower(slim)

    plt.boxplot(nparray, showfliers=False) #Boxplot qual scores per base for each read>computationaly intensive and long process depending on dataset size: !!! showfliers=False to remove outliers from making boxplots
    #unreadable, as there can be a LOT of outliers depending on dataset sizes/quality which can clutter the graph with dots, making it harder to read

    plt.title('Quality scores box plot of each base position per read')

    #Fixing xtick labels (base sequence lengths) so it starts correctly from 1 (because python counts from 0)
    #TODO: This assumes again that each read seq is same length...
    xlabels=0
    xlist=[]
    for x in range(len(nparray[0])):
        xlabels+=1
        xlist.append(xlabels)

    plt.xticks(rotation=90) #Make xtick values readable, as if they are lined up along the horizon, the values overlap each other making it unreadable, so I rotate it sideways to prevent overlapping texts
    plt.ylabel('Quality Score')
    plt.show()

    del nparray #clear mem space

#Benchmarking to measure time taken to complete this func for large datasets todo


def test(x):
    t=x+1
    print(t)
    return(t)
