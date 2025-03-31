"""
Aldente Fasta: Tile sequences Module

This module contains functions for illumina sequenced files to process tile seq id's and quality scores.

"""
################################
################################
################################

### ONLY FOR ILLUMINA SEQUENCES AS THEY CONTAIN SEQ ID'S ###
### ONLY FOR ILLUMINA SEQUENCES AS THEY CONTAIN SEQ ID'S ###
### ONLY FOR ILLUMINA SEQUENCES AS THEY CONTAIN SEQ ID'S ###


import re #for regex
import gzip #read gz compressed files
import numpy as np #process array matrices
import pandas as pd #to load data from raw file into pd dataframe
import seaborn as sns #for plotting graphs
import matplotlib.pyplot as plt #for plotting graphs


################################
################################
################################

def tile_num(file, out):
    """
    Takes a fastq file containing illumina seq id's and creates a text file output of tile number of each read sequence
    :param file: fastq file TODO: File support/detection for other file type/compressions AND error detections for inappropriate file types
    :param out: text file containing tile number (str format)
    :return: param out
    """

    line_no=3 #adjusting to 3 to relate if statement to seq id lines of fastq files
    tile_num_regex=rb':(\d{4}):' #rb for binary or 'r' for normal TODO: if adding other file support, need to add appropriate reading of str/byte >improve efficiency and speed of writing and reading data

    with gzip.open(file,'r') as fastq:  #gzip.open for fastq gz file or just open for normal TODO: if adding file support, need to add appropriate changes
        with open(out,'wb') as outfile: #wb for writing bytes or just 'w' for normal TODO: if adding file support, need to add appropriate changes

            for line in fastq:
                line_no+=1
                if (line_no%4==0): #line_no starting from 3 so, this detects seq id line of every read
                    tile_n=re.search(tile_num_regex, line).group(1) #Searches line for our regex pattern to write into <outfile>
                    outfile.write(tile_n+b'\n') #b'\n' for adding new line when writing in byte TODO: if adding file support, need to add appropriate changes

    return

################################
################################
################################

def raw_qual_to_num(raw_qual, enc):
    """
    Decodes the raw quality scores using the ASCII 33/64 code
    :param raw_qual: Raw quality score of fastq file
    :param enc: Decoding number
    :return: Returns a list of numbers converted from the raw quality scores
    """

    raw_qual=raw_qual.rstrip()  #rstrip needed to remove trailing artefacts (eg: '\n') when reading from files
    convert=[(int(ord(x) - enc)) for x in raw_qual]  #Ord returns number value of the raw_qual characters and subtract the enc value (33/64) to recieve Phred Qualty Score

    return convert

################################
################################
################################

def tile_link(tile_num, qual_out, enc):
    """
    Takes raw files of tile numbers and quality scores of fastq files and the encoding value to return the averages of the quality score of each base position in each tile number
    :param tile_num: File containing tile numbers of your fastq file (<tile_num> func)
    :param qual_out: File containing raw quality scores of your fastq file (<aldente_fasta.qual_out>
    :param enc: A number value of 33 or 64 to decode raw quality scores (64 for Illumina 1.0-1.8, 33 for Sanger and Illumina 1.8+)
    :return: Dictionary containing average scores of each base positions quality score grouped by its tile number
    """

    with open(tile_num, 'r') as raw_tiles:
        with open(qual_out,'r') as raw_quals:
            #Setting dictionaries for tile_num:quals_score and tile_num:counts/freq
            tile_dict={}
            tile_counts={}

            for tile,qual in zip(raw_tiles, raw_quals): #Iterating over both files, they SHOULD contain same length lines as they come from same fastq file

                tile=tile.rstrip() #Removing formatting/trailing artefacts from raw_tiles
                qual=raw_qual_to_num(qual,enc) #Func to decode and convert raw_quals to phred quality score

                if tile in tile_dict: #If tile number is already a key in tile_dictionary then:
                    #Add values into its tile_number key:
                    tile_dict[tile]=tile_dict[tile]+np.array(qual) #Takes the values in tile_dict[tile] and sums/adds the decoded quality score
                    tile_counts[tile]+=1 #Increments tile_counts[tile] values by 1 to count each time a specific tile_number is seen>frequency of tile_number

                else:
                    #If tile number is NOT a key or NOT found in tile_dict:
                    tile_dict[tile]=np.array(qual) #Creates a new key of tile_number in tile_dict and adds the values as the decoded quality score
                    tile_counts[tile]=1 #Creates a new key of tile number in tile_counts and sets its frequency/counts value as 1

    #Calculating the average scores since we have sums of quality scores in tile_dict and the counts of tile numbers in dict_counts
    avg_tiles={}

    for key in tile_dict:
        avg_tiles[key]=tile_dict[key]/tile_counts[key]

    return avg_tiles

################################
################################
################################

def heat_map(avg_dict):
    """
    Heatmap to show distribution of avg quality scores for each base seq. position in every tile number
    :param avg_dict: Dictionary containing tile number for key and values containing a list of avg quality scores (<tile_link>)
    :return: NONE (plt.show() a heatmap)
    """

    #Creating true xlabels since python counts from 0, but we want to show base seq position which starts from 1 to 100.
    xlabel=[]
    count=0

    for x in list(avg_dict.values())[0]:
        count+=1
        xlabel.append(count)


    df=pd.DataFrame(avg_dict, index=xlabel) #Creating a pandas dataframe of our dictionary
    df=df.T #Transposing dataframe to get x values as base seq positions and y values as tile numbers
    plt.figure(figsize=(20,10)) #plot figure size
    #Using seaborn to generate a heatmap from our data,
    sns.heatmap(df,
                annot=False, #False because we don't want to clutter our heatmap with annotations of quality score values
                fmt='.2f', #Formatting for floats because our average quality scores contains decimal values
                cmap='coolwarm', #Colour scheme for heatmap
                cbar_kws={'label':'Quality Score'}) #Adding a label for the colour bar

    plt.ylabel('Tile Number')
    plt.xlabel('Base Seq. Position')
    plt.yticks(rotation=0) #Setting rotation to 0 so we can read tile numbers (decluttering it)
    plt.title('Heat map of avg quality scores for each base sequence position in all tile numbers')
    plt.show()

    return

################################
################################
################################
