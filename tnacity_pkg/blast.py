"""
Description
"""
# Imports --------------------------------------------------------------------------------------------------
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
#Turn off the SettingWithCopy warning, b/c danger is my middle name
pd.options.mode.chained_assignment = None
import subprocess
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET
from tnacity_pkg import fasta_lib as FU

# Functions ------------------------------------------------------------------------------------------------

def run_blast(fasta):
    """
    subroutine to run blast
    """
    file = Path(fasta)
    outfile = file.with_suffix('.blastn')
    subprocess.run(f'blastn -query {file} -subject {file} -word_size 4 -outfmt "6 std frames" -max_hsps 100 -evalue 100 > {outfile}',
                  shell = True,
                  check = True)
    
def count_num_unique_ign_pairs(df):
    """
    to do
    """
    #create a column that joins (Seq1 - Seq2 or Seq2 - Seq1) to Seq1Seq2
    df["pairs"] = df.apply(lambda row: ''.join(sorted([row["qaccver"],row["q_ign"], row["subj_ign"]])), axis = 1)
    
    #remove rows with duplicate pairs
    df = (df.sort_values(by = ['qaccver'])
               .drop_duplicates(subset = 'pairs')
               .drop(columns = ['pairs'])
         )
    
    #count the number of unique pairs per query
    grouped = df.groupby(['qaccver'])
    df_size = grouped.agg(unique_ign_pairs=pd.NamedAgg(column='qaccver', 
                                                   aggfunc='size'
                                                  )
                          )

    df_size2 = df_size.reset_index()
    try:
        num_unique_pairs = df_size2.unique_ign_pairs.tolist()[0]
    except IndexError:
        num_unique_pairs = 0
    return num_unique_pairs

def format_blast_output(blastoutfmt6):
    df = pd.read_csv(blastoutfmt6, 
                     header = None, 
                     sep = '\t',
                     names = ["qaccver", 
                             "saccver", 
                             "pident", 
                             "length", 
                             "mismatch",
                             "gapopen", 
                             "qstart",
                             "qend",
                             "sstart",
                             "send",
                             "evalue",
                             "bitscore",
                            "frames",
                            ]
                    )

    #Change query/subject to reflect actual coordinates
    ##Split the query/subject names
    df[["saccver2", "subj_ign", "subj_ign_start", "subj_ign_end"]] = df.saccver.str.split("|", expand = True)
    df[["qaccver2", "q_ign", "q_ign_start", "q_ign_end"]] = df.qaccver.str.split("|", expand = True)
    
    ##Change dtype to int64
    df[["q_ign_start", "q_ign_end", "subj_ign_start", "subj_ign_end"]] = df.loc[:,["q_ign_start", 
                                                                                   "q_ign_end", 
                                                                                   "subj_ign_start", 
                                                                                   "subj_ign_end"
                                                                                  ]
                                                                               ].astype('int64')
    
    #Keep the original ids 
    df["qaccver3"] = df["qaccver"]
    df["saccver3"] = df["saccver"]

    ##Update the query/subject coordinates
    df["saccver"] = df["saccver2"]
    df["qaccver"] = df["qaccver2"]
    df["sstart"] = df["sstart"] + df["subj_ign_start"] - 1
    df["send"] = df["send"] + df["subj_ign_start"] - 1
    df["qstart"] = df["qstart"] + df["q_ign_start"] - 1
    df["qend"] = df["qend"] + df["q_ign_start"] - 1

    #remove extra columns
    df.drop(columns = ["saccver2", 
                       "subj_ign_start", 
                       "subj_ign_end",
                       "qaccver2", 
                       "q_ign_start", 
                       "q_ign_end"
                      ], 
            inplace = True
           )
    
    return df

def filter_blast_for_IRs(df, min_len, max_len, max_gaps, max_mismatches, max_ir_dist):
    """
    Selects IRs from a blast output. IRs must be:
        On different query/subjects (i.e., no self IRs)
        <= max_gaps/mismatches/len
        <= max_ir_dist apart
    Returns a pandas dataframe
    """ 
    
    #Remove self-hits
    df1 = df.query('q_ign != subj_ign')
    
    #distance between IRs
    df1['IR_dist'] = (df1['qstart'] - df1['sstart']).abs()    
    
    #Select only inverted repeats that fit the provided criteria
    df2 = df1.query('length < @max_len and length > @min_len and frames == "1/-1" and gapopen <= @max_gaps and mismatch <= @max_mismatches and IR_dist <= @max_ir_dist')

    return df2
                                
def filter_blast_HSPs_near_coord(df, 
                                coordinate,
                                max_dist):
    """
    Filters the blast output for:
        qstart distance to coordinate <= max_dist
        The query/subject pair encloses the coordinate (e.g., qstart < coordinate < send)

    Returns a pandas dataframe
    """ 
    
    #Make a column calculating the distance of the query to coordinate
    df['dist'] = (df['qstart'] - coordinate).abs()
    
    #select queries less than @max_dist from TnsB
    df1 = df.query('dist <= @max_dist')
    
    #select candidate ends that contain TnsB
    ##all the hits are +/-, so sstart > send.
    ##For some reason, I get an error if I combine the two queries into one, so just split and recombine
    df2 = df1.query('qstart <= @coordinate and sstart >= @coordinate')
    df3 = df1.query('send <= @coordinate and qend >= @coordinate')
    df4 = pd.concat([df2, df3])
    
    #sort by distance from coordinate, provide new indexes
    df5 = df4.sort_values(by = 'dist').reset_index(drop = True)
    
    return df5
 
def filter_blast_for_multi_HSPs(df, min_num_HSPs=2):
    """
    Counts the number of HSPs between a query/subject pair
    """
    #select IGN pairs that have multiple HSPs
    df["num_HSPs"] = df.groupby(['q_ign', 'subj_ign'])['qaccver'].transform('count')
    df1 = df.query('num_HSPs >= @min_num_HSPs')

    return df1

def run_blast_filtering(blast_output, tnsB_coordinate, min_len, max_len, max_gaps, max_mismatches, max_ir_dist, max_anchor_dist):
    """
    Accepts path/to/blastt output, the names of intergenic loci and a tnsB coordinate
    runs the filtering pipeline
    At each step, counts the total number of unique IGN pairs
    """
    
    #Format the blast output, check/count ends
    df = format_blast_output(blast_output)
    assert not df.empty, "The blast output is empty. Check the fasta file input to see what is going wrong"
    count1 = count_num_unique_ign_pairs(df)
    
    #Filter for IRs
    df_IR = filter_blast_for_IRs(df, min_len, max_len, max_gaps, max_mismatches, max_ir_dist)
    count2 = count_num_unique_ign_pairs(df_IR)
    
    #filter by TnsB location
    df_tnsb = filter_blast_HSPs_near_coord(df_IR, tnsB_coordinate, max_anchor_dist)
    count3 = count_num_unique_ign_pairs(df_tnsb)
    
    #filter for multiple IRs
    df_multi_IR = filter_blast_for_multi_HSPs(df_tnsb)
    count4 = count_num_unique_ign_pairs(df_multi_IR)

    count_d = {'blast' : count1,
             'IR' : count2,
             'TnsB_dist' : count3,
             'multi_IR' : count4,
             }
    
    accn = df.qaccver.unique().tolist()[0]
    df_count = pd.DataFrame(count_d, index = [accn])
    
    return df_multi_IR, df_count

def deduplicate_reciprocal_pairs(df):
    """
    Accepts a blastoutput dataframe from the function "run_blast_filtering"
    Deduplicates reciprocal pairs
    Returns a df
    """
    
    df2 = (df.loc[:,["qaccver", 
                    'subj_ign', 
                    'pairs', 
                    "q_ign", 
                    'qaccver3', 
                    'saccver3']
                ]
            .sort_values(by = ['qaccver'])
            .drop_duplicates(subset = 'pairs')
            .drop(columns = ['pairs'])
         )
    
    return df2
    
   