"""
MODULE INFO
"""


#WHAT ABOUT SEQUENCES OF DIFFERENT LENGTH?
#WILL MY BASE SEQUENCE QUALITY SCORE STILL WORK IF READS ARE DIFFERENT LENGTH?...

########################
########################
########################

import gzip
import matplotlib.pyplot as plt

def base_count(file):
    ###takes a fastq file and returns dict of base content in all reads

    base_dict={'A':0,'C':0,'T':0,'G':0,'N':0}
    bases=['A','C','T','G','N']
    line_no=0
    with gzip.open(file, 'r') as fastq:

        for line in fastq:
            line_no+=1

            if (line_no % 2 == 0) and (line_no % 4 == 2):  # If statement to check every 2nd sequence line of each reads (fastqfile has a 4 line repeating format containing read information)
                for base in line:
                    if base in bases:
                        base_dict[base]+=1

    return base_dict

########################
########################
########################

def base_pos_count(file):
    """
    Reads all bases in each position of every read sequences and returns statistics on them
    :param file: Fastq file
    :return: Dictionary containing stats of bases
    """

    #ASSUMES EACH READ IS SAME LENGTH...
    #TODO: Different file support?
    #TODO: Add support for when read sequences are different lengths
    #TODO: Improve efficiency and speed of this code...try minimize the loops

    import aldente_fasta as af #calling af.head function for intializing first part:

    #Initializing variables for this func
    a=[]
    c=[]
    t=[]
    g=[]
    n=[]
    line_no=0

    #Reading sequence line from fastq file
    seq_prev=af.head(file,2,True) #TODO: CHANGE IF IM DOING ADDING DIFFERENT FILE SUPPORT
    seq=(seq_prev[1].decode()) #Decode as it is in bytes (utf-8) format: b'string format example'

    #Setting up base lists and filling up its values with 0 for each position #TODO: ASSUMING EACH READ IS SAME LENGTH...
    #Probably a better way to do all this TODO: IMPROVE EFFICIENCY AND SPEED
    for x in seq:
        a.append(0)
        c.append(0)
        t.append(0)
        g.append(0)
        n.append(0)


    with gzip.open(file,'r') as fastq:
        #Iterating over each line
        for line in fastq:
            line_no+=1 #Used to count line number and use if loop to check for sequence line
            #Reformat byte into string to run loops over and check for bases in the string sequence
            line=line.rstrip()
            line=line.decode()
            pos=0 #Used to count position for while loop to increment base lists by 1 for each position the base is counted in our sequence
            if (line_no % 2 == 0) and (line_no % 4 == 2):  # If statement to check every 2nd sequence line of each reads (fastqfile has a 4 line repeating format containing read information)
                while pos<len(seq): #Run loop over all positions of seq
                    #Increments base list by 1 if it finds(==) the base in that position
                    #TODO: Improve efficiency and speed of this by implementing other methods...maybe using a dictionary would be faster?
                    if line[pos]=='A':
                        a[pos]+=1
                    elif line[pos]=='C':
                        c[pos]+=1
                    elif line[pos]=='T':
                        t[pos]+=1
                    elif line[pos]=='G':
                        g[pos]+=1
                    elif line[pos]=='N':
                        n[pos]+=1
                    pos+=1

    #Calculating total bases for each position and the GC totals
    #Has to be a better efficient way to do this: TODO: Increase eficiency and speed....
    pos=0
    total=[]
    gc=[]
    for x in a:
        total.append(a[pos]+c[pos]+t[pos]+g[pos]+n[pos])
        gc.append(g[pos]+c[pos])
        pos+=1

    #Add base list results into a dictionary for cleaner format to return data
    res={'A':a,'C':c,'T':t,'G':g,'N':n,'GC':gc,'Total':total}

    #Converting to percentages
    pos=0
    a_perc=[]
    c_perc=[]
    t_perc=[]
    g_perc=[]
    n_perc=[]
    gc_perc=[]
    for x in a:
        a_perc.append(100*(a[pos]/total[pos]))
        c_perc.append(100*(c[pos]/total[pos]))
        t_perc.append(100*(t[pos]/total[pos]))
        g_perc.append(100*(g[pos]/total[pos]))
        n_perc.append(100*(n[pos]/total[pos]))
        gc_perc.append(100*(gc[pos]/total[pos]))
        pos+=1
    res_perc={'A%':a_perc,'C%':c_perc,'T%':t_perc,'G%':g_perc,'N%':n_perc,'GC%':gc_perc,'Total':total}

    graph_res={'A%':a_perc,'C%':c_perc,'T%':t_perc,'G%':g_perc,'N%':n_perc}
    graph_gc={'GC%':gc_perc}
    return res, graph_res, graph_gc

########################
########################
########################

def base_graphs(file):

    first,res,res_gc=base_pos_count(file)

    #Plotting each base's data as a line
    plt.figure(figsize=(20,10))
    plt.plot(res['A%'],label='A%')
    plt.plot(res['C%'],label='C%')
    plt.plot(res['T%'],label='T%')
    plt.plot(res['G%'],label='G%')
    plt.plot(res_gc['GC%'],label='GC%')

    # Creating true xlabels since python counts from 0, but we want to show base seq position which starts from 1 to 100.
    xlabel = []
    xlabelold=[]
    count = 0

    for x in res['A%']:
        xlabelold.append(count)
        count += 1
        xlabel.append(count)


    plt.xlabel('Base position')
    plt.ylabel('%')
    plt.title('Sequence % in each base position across all reads')
    plt.xticks(ticks=xlabelold, labels=xlabel, rotation=90)
    plt.ylim(0,100)
    plt.legend()
    plt.grid(True)
    plt.show()

    return

########################
########################
########################

