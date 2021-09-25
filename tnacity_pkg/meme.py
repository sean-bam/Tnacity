"""
Description
"""
# Imports --------------------------------------------------------------------------------------------------
import pandas as pd
#Turn off the SettingWithCopy warning, b/c danger is my middle name
pd.options.mode.chained_assignment = None
import subprocess
from pathlib import Path
import xml.etree.ElementTree as ET

# Functions ------------------------------------------------------------------------------------------------
def make_meme_background_file(fasta, output, order=1):
    """
    Subroutine to run meme's fasta-get-markov script
    """
    if Path(output).exists():
        Path(output).unlink()
        
    p1 = subprocess.run(f'fasta-get-markov -m {order} {fasta} > {output}',
                  shell = True,
                  check = True,
                  stdout=subprocess.DEVNULL,
                  stderr=subprocess.DEVNULL)
    
def run_meme(fastafile, 
             output_dir,
             backgroundfile,
             num_motifs = 10, 
             minwidth = 15, 
             maxwidth = 20, 
             minsites = 4, 
             maxsites = 6,
             threads = 2,
             markov_order = 0,
             evt = 0.1):
    """
    Subroutine to run meme
    """
    #background_input = ""
    #if backgroundfile:
    #    background_input = f"-bfile {backgroundfile}"
 
    subprocess.run(f'meme {fastafile} -bfile {backgroundfile} -dna -oc {output_dir} -mod anr -nmotifs {num_motifs} -minw {minwidth} -maxw {maxwidth} -minsites {minsites} -maxsites {maxsites} -revcomp -p {threads} -markov_order {markov_order} -evt {evt}',
                  shell = True,
                  check = True,
                  stdout=subprocess.DEVNULL,
                  stderr=subprocess.DEVNULL)
    
def get_meme_results_as_df(meme_xml):
    """
    Converts the XML results of meme to a dataframe
    """   
    tree = ET.parse(meme_xml)
    root = tree.getroot()
    
    df_list = []
    seq_id_dict = {}
    
    #initialize an empty dataframe
    df = pd.DataFrame(columns = ['motif_id', 
                                'motif_name',
                                'width',
                                'sites',
                                'ic',
                                're',
                                'llr',
                                'motif_pvalue',
                                'motif_evalue',
                                'bayes_threshold',
                                'seq_id',
                                'motif_start',
                                'strand',
                                'pvalue'],
                     )
    
    #make a dictionary of the meme seqIDs (e.g., sequence_0) to their real accession (e.g., NC123_01)
    for sequence_tree in tree.find("training_set").findall("sequence"):
        seq_id = sequence_tree.get("id")
        seq_name = sequence_tree.get("name")
        seq_id_dict[seq_id] = seq_name
        #print(seq_id, seq_name)
    
    #Add motif hits to the dataframe
    i = 0
    for motif in root.find("motifs").findall("motif"):
        for hit in motif.find("contributing_sites").findall("contributing_site"):
            #motif info
            df.loc[i, "motif_id"] = motif.attrib.get('id')
            df.loc[i, "motif_name"] = motif.attrib.get('name')
            df.loc[i, "width"] = motif.attrib.get('width')
            df.loc[i, "sites"] = motif.attrib.get('sites')
            df.loc[i, "ic"] = motif.attrib.get('ic')
            df.loc[i,"re"] = motif.attrib.get('re')
            df.loc[i, "llr"] = motif.attrib.get('llr')
            df.loc[i, "motif_pvalue"] = motif.attrib.get('p_value')
            df.loc[i, "motif_evalue"] = motif.attrib.get('e_value')
            df.loc[i,"bayes_threshold"] = motif.attrib.get('bayes_threshold')
            
            #motif hits 
            hit_seq_raw = hit.attrib.get('sequence_id') #get the fake sequence ID
            df.loc[i, "seq_id"] = seq_id_dict.get(hit_seq_raw) #get the real sequence ID
            df.loc[i, "motif_start"] = hit.attrib.get('position')
            df.loc[i, "strand"] = hit.attrib.get('strand')
            df.loc[i, "pvalue"] = hit.attrib.get('pvalue')
            
            
            #df_list.append(df2)
            i += 1
    
    #If there are motif hits, update the coordinates
    if df.shape[0] > 0:
        df2 = df.astype({'width' : 'int64',
                      'sites' : 'int64',
                      'ic' : 'float',
                      're' : 'float',
                      'llr' : 'float',
                      'motif_pvalue' : 'float',
                      'motif_evalue' : 'float',
                      'bayes_threshold' : 'float',
                      'motif_start' : 'int64',
                      'pvalue' : 'float'
                     }
                    )
        df3 = update_meme_coordinates(df2)
    else:
        df3 = df.copy()
    
    return df3

def add_counts_of_motifs(df):
    """
    Adds the following counts to a motif dataframe:
    1. The total number of sequences where the motif is found 
    2. The total number of occurrences of the motif per sequence
    3. The total number of strands where the motif is found
    """
    df["num_seqs_motif_present"] = (df.groupby(by = ['motif_id'])['seq_id']
                                      .transform('nunique')
                                     )
    
    df["motif_count_per_seq_strand"] = (df.groupby(by = ['motif_id', 'seq_id', 'strand'])['seq_id']
                                    .transform('count')
                                )
    
    #using .transform('count') doesn't work with .nunique, so I have to do it this way
    df["num_strands_motif_present"] = (df.groupby('motif_id')['strand']
                                         .transform('nunique')
                                      )
    
    return df

def update_meme_coordinates(df):
    """
    Meme reports the start of motifs in zero-based coordinates and doesn't report the stop location
    This function converts the motif_start to the real start and adds the stop coordinate
    
    Expect a dataframe from get_meme_results_as_df() and the seqid to look like this:
    accession|something|start|stop  
    """
    
    #Change query/subject to reflect actual coordinates
    ##Split the query/subject names
    df[["accession", "ign", "ign_start", "ign_end"]] = df.seq_id.str.split("|", expand = True)

    ##Change dtype to int64
    df[["ign_start", "ign_end"]] = df.loc[:,["ign_start", "ign_end"]].astype('int64')

    ##Update the query/subject coordinates (MEME is zero based)
    df["motif_start"] = df["motif_start"] + df["ign_start"] - 1
    df["motif_end"] = df["motif_start"] + df["width"] - 1
    
    ##drop the extra column
    #df.drop(columns = ["accession"], inplace = True)

    return df

def filter_meme_results(df, min_num_seqs=2, min_num_strands=2, min_num_hits=2, evalue = 0.1, pvalue = 0.00001):
    """
    Filters meme output for motifs 
        <= evalue 
        <= max_motif_count
        Present on more than one sequence
        At least one copy of the motif is on the opposite strand of the pair
    """
    
    #remove bad motifs/hits
    df2 = df.query('motif_evalue <= @evalue and pvalue <= @pvalue')
    
    #count the good hits
    df3 = add_counts_of_motifs(df2)
    
    #remove motif hits that only occur once per seq/strand
    df4 = df3.query('motif_count_per_seq_strand >= @min_num_hits')
    
    #update the counts
    df5 = add_counts_of_motifs(df4)
    
    #Select motifs that fit the criteria
    #If the motif is present on BOTH input sequences, 
    #simply counting the number of unique elements in the strand column will tell me if
    #if there is at least one copy of the motif present on a different strand between the two seqs
    df6 = df5.query('num_seqs_motif_present >= @min_num_seqs and num_strands_motif_present >= @min_num_strands')
    
    #If only one sequence in the pair passes the thresholds, clear the dataframe
    if df6.seq_id.nunique() < 2:
        df6 = df6.iloc[0:0]
        
    return df6

def select_ends_from_meme(df):
    """
    Selects the top motif by evalue
    """
    #Sort by motif_evalue then groupby motif WITHOUT RE-SORTING
    group = (df.sort_values(by = ['accession', 
                                   "motif_evalue"
                                  ]
                            )
                 .groupby('motif_name', 
                          sort = False
                         )
            )

    #Get only the first group
    df2 = group.get_group((list(group.groups)[0]))

    #remove unnecessary columns
    df2.drop(columns = ['seq_id',
                       'motif_id',
                       'num_seqs_motif_present',
                       'num_strands_motif_present',
                      ],
             inplace = True
            )
    
    return df2


def get_motifs_from_meme(motif_name, meme_xml, output):
    """
    accepts a motif_name (e.g., "CCHSTCCCCT")
    Exports to a new meme-formatted output file
    """
    #i = 0
    #num_motifs_to_get = len(motif_names)
    #num_motifs_to_get = 1
    
    p = Path(output)
    #if p.exists():
    #    p.unlink()
        
    subprocess.run(f'meme-get-motif -id {motif_name} {meme_xml} >> {p}',
                  shell = True,
                  check = True)
        
    #with p.open() as f: 
    #    for line in f:
    #        if line.startswith("MOTIF"):
    #           i += 1
    
    #assert num_motifs_to_get == i, f"You passed {num_motifs_to_get} motifs, but only {i} were retrieved"
    
def run_mcast(fastafile, backgroundfile, output, memefile, gapsize=75, maxwidth=250):
    """
    Runs mcast
    Note: I disable check = True, because Mcast throws and error if the pvalue is really high.
    note: I use a high evalue threshold for reporting, so that I don't have empty datatables to parse.
    It doesn't seem to increase the search time
    """
    
    subprocess.run(f'mcast -oc {output} --bfile {backgroundfile} --max-gap {gapsize} --output-ethresh 100 --synth --max-total-width {maxwidth} {memefile} {fastafile} ',
                  #check = True,
                  shell = True,
                  stdout=subprocess.DEVNULL,
                  stderr=subprocess.DEVNULL)

def get_mcast_output_as_df(mcast_tsv):
    """
    to do
    """
    df = pd.read_csv(mcast_tsv, sep = "\t", comment = "#")
    
    #mcast spits out an empty table if the results aren't significant
    #If this is the case, update_mcast_coordinates throws an error
    #just return the empty df if so
    #if df.shape[0] < 1:
    #    return df
    
    #Fix the coordinates
    df2 = update_mcast_coordinates(df)

    #Fix the P-value and Q-value columns
    df2.rename(columns = {"p-value" : "p_value", 
                         "q-value" : "q_value",
                         "E-value" : "e_value"}, inplace = True)
    df2[["p_value", "q_value", "e_value"]] = df2.loc[:,["p_value", "q_value", "e_value"]].astype('float')

    #rearrange the columns
    df2 = df2.loc[:,["sequence_name",
                   "seq_id",
                   "ign",
                   #"pattern_name",
                   "cluster_start",
                   "cluster_stop",
                   "score",
                   "e_value",
                   "p_value",
                   "q_value",
                   "matched_sequence"
                  ]
               ]

    return df2

def update_mcast_coordinates(df):
    """
    Fimo reports the start of motifs in zero-based coordinates
    This function converts the motif_start/stops to the real start/stops
    
    Expect a dataframe from get_meme_results_as_df() and the seqid to look like this:
    accession|something|start|stop
    
    **Note**: FIMO will do this automatically if my headers are in the format 
    >name:start-stop
    But, at this point, I'll just use this function
    """
    
    #Change query/subject to reflect actual coordinates
    ##Split the query/subject names
    df[["seq_id", "ign", "ign_start", "ign_end"]] = df.sequence_name.str.split("|", expand = True)

    ##Change dtype to int64
    df[["ign_start", "ign_end"]] = df.loc[:, ["ign_start", "ign_end"]].astype('int64')

    ##Update the query/subject coordinates (MEME is zero based)
    df["cluster_start"] = df["start"] + df["ign_start"] - 1
    df["cluster_stop"] = df["stop"] + df["ign_start"] - 1
    
    df[["cluster_start", "cluster_stop"]] = df.loc[:, ["cluster_start", "cluster_stop"]].astype('int64')
    
    return df

def filter_mcast(df, max_q=0.01, max_igns=4):
    """
    Expects an mcast dataframe with only two sequences under consideration
    """
    #remove bad hits
    df1 = df.query('q_value <= @max_q')
    
    #if the motif is too common
    if df1.ign.nunique() > 4:
        df1 = df1.iloc[0:0]
    
    #If only one sequence in the pair passes the thresholds, clear the dataframe
    if df1.sequence_name.nunique() < 2:
        df1 = df1.iloc[0:0]
        
    return df1

def select_ends_from_mcast_results(df, coordinate):
    """
    TBD
    """
    #Calculate the product of the pvalues between each pair with a shared motif
    df2 = df.groupby('motif_file').apply(calculate_pair_qvalue)
    df3 = df2.reset_index(drop = True)
    
    #select the top pair
    df4 = df3.sort_values(by = 'pair_qvalue').head(2)
    
    #set proximal/distal columns
    df4["dist"] = (df4["cluster_start"] - int(coordinate))
    df4["proximal_end"] = df4.set_index('ign').loc[:,"dist"].abs().idxmin()
    df4["distal_end"] = df4.set_index('ign').loc[:,"dist"].abs().idxmax()
    #df5 = df4.query('pair_qvalue <= @max_q')
    
    
    return df4

def calculate_pair_qvalue(group):
    #Pick the best cluster per sequence
    group2 = group.sort_values(by = 'q_value').drop_duplicates(subset = 'sequence_name')
    #compute the product of the pvalues
    group2['pair_qvalue'] = group2['q_value'].median()
    return group2

def get_good_motifs_from_mast(mast_xml):
    """
    Parses the output of MAST to select uncorrelated motifs 
    Returns as list
    """
    tree = ET.parse(mast_xml)
    root = tree.getroot()

    good_motifs = []
    for motif in root.findall("./motifs/"):
        
        #Get the name of the sequence
        motif_id = motif.get('id')
        bad_quality = motif.get('bad')
        
        if bad_quality != 'y':
            good_motifs.append(motif_id)
            
    return good_motifs


def count_motifs_in_meme_xml(meme_xml):
    
    i = 0
    motif_ids = []
    with Path(meme_xml).open() as f: 
        for line in f:
            if line.startswith("MOTIF"):
                i += 1
                motif_id = line.strip().split()[1]
                motif_ids.append(motif_id)
                
    motif_ids_unique = set(motif_ids)
    if len(motif_ids_unique) < len(motif_ids):
        return False
    return True
    
    
def run_mast(fastafile, memefile, output):
    """
    Runs MAST, scoring each strand separately based on the given motif
    """
    
    subprocess.run(f'mast -oc {output} -remcorr {memefile} {fastafile}',
                  check = True,
                  shell = True)
    
def dereplicate_motifs(meme_xml, output):
    """
    Accepts a single file of meme motifs,
    runs a fake MAST search
    Parses the output of MAST to return uncorrelated motifs
    """
    assert Path(output).exists() == False, f"{output} already exists and we don't want to overwrite it"
       
    #Make a fake fasta file    
    with open('fake_fasta_file_for_mast.fna', 'w') as f:
        print(">fake_seq", file = f)
        print("ATGCATGCATGCATGCATGCATGC", file = f)
    
    #Search with mast
    run_mast('fake_fasta_file_for_mast.fna', meme_xml, 'tmp_mast_dir/')
    
    #Parse the mast output to select the good motif IDs
    good_motifs_list = get_good_motifs_from_mast('tmp_mast_dir/mast.xml')
    
    #Get the motifs
    num_motifs_to_get = len(good_motifs_list)
    for motif in good_motifs_list:
        get_motifs_from_meme(motif, meme_xml, output)
    
    #make sure then number of motifs in the output agrees
    num_dereplicated_motifs = count_motifs_in_meme_xml(output)
    if len(good_motifs_list) != num_dereplicated_motifs:
        print(f"""
        hmm, there were {num_motifs_to_get} non-redundant motifs, but there are
        {num_dereplicated_motifs} in {output}
        This is probably because {meme_xml} has identical motifs in it with the same
        name. Remove those and hopefully this message will go away :)
        """)
    
    #cleanup
    Path('fake_fasta_file_for_mast.fna').unlink()
    for file in Path('tmp_mast_dir/').glob('**/*'):
        file.unlink()
    Path('tmp_mast_dir/').rmdir()
    
def run_fimo(fastafile, memefile, output, backgroundfile, q_thresh):
    """
    Runs FIMO
    """
    #check to make sure there aren't duplicate motifs
    assert count_motifs_in_meme_xml(memefile), f'{memefile} has duplicate motifs. Make it non-redundant'
    
    subprocess.run(f'fimo -oc {output} --bfile {backgroundfile} --qv-thresh --thresh {q_thresh} {memefile} {fastafile}',
                  check = True,
                  shell = True,
                  stdout=subprocess.DEVNULL,
                  stderr=subprocess.DEVNULL)
    
def get_fimo_output_as_df(fimo_tsv):
    """
    to do
    """
    df = pd.read_csv(fimo_tsv, sep = "\t", comment = "#")
    
    #Fix the coordinates
    df2 = update_fimo_coordinates(df)
    
    #Fix the P-value and Q-value columns
    df.rename(columns = {"p-value" : "p_value", 
                         "q-value" : "q_value"}, inplace = True)
    df[["score", "p_value", "q_value"]] = df.loc[:,["score", "p_value", "q_value"]].astype('float')
    
    #rearrange the columns
    df = df.loc[:,["sequence_name",
                   "seq_id",
                   "ign", 
                   #"length",
                   "motif_start",
                   "motif_end",
                   "strand",
                   "score",
                   "p_value",
                   "q_value",
                   "matched_sequence",
                   "motif_id",
                   #"motif_count",
                   #"ign_start",
                   #"ign_end",
                  ]
               ]
    
    return df

def update_fimo_coordinates(df):
    """
    Fimo reports the start of motifs in zero-based coordinates
    This function converts the motif_start/stops to the real start/stops
    
    Expect a dataframe from get_meme_results_as_df() and the seqid to look like this:
    accession|something|start|stop
    
    **Note**: FIMO will do this automatically if my headers are in the format 
    >name:start-stop
    But, at this point, I'll just use this function
    """
    
    #Change query/subject to reflect actual coordinates
    ##Split the query/subject names
    df[["seq_id", "ign", "ign_start", "ign_end"]] = df.sequence_name.str.split("|", expand = True)

    ##Change dtype to int64
    df[["ign_start", "ign_end",]] = df.loc[:, ["ign_start", "ign_end"]].astype('int64')

    ##Update the query/subject coordinates (MEME is zero based)
    df["motif_start"] = df["start"] + df["ign_start"] - 1
    df["motif_end"] = df["stop"] + df["ign_start"] - 1
    
    return df


def calculate_ign_qvalue(group):
    #compute the product of the qvalues
    group['ign_motif_qvalue'] = group['q_value'].product()
    return group


def filter_fimo_output_across_ends(df, max_q=0.1, min_ubiquity=0.8):
    """
    """
    num_seqids = len(df.seq_id.unique().tolist())
    
    #Select only hits that pass the Q-value threshold on the contig of interest
    df2 = df.query('q_value <= @max_q') #and seq_id == @accession')
    
    #Count how many accessions have the motif
    df2["num_seqids_motif_present"] = (df2.groupby(by = ['motif_id'])['seq_id']
                                      .transform('nunique')
                                     )
    
    #Count the occurrences of each motif on each ign
    df2["motif_count_per_ign"] = (df2.groupby(by = ['motif_id', 'sequence_name'])['sequence_name']
                                      .transform('count')
                                     )
    
    #select sequences that meet the criteria
    min_num_seqids = num_seqids * min_ubiquity
    df3 = df2.query('motif_count_per_ign > 1 and num_seqids_motif_present > @min_num_seqids')
    
    #Get only motifs that are present on more than two sequences and in multiple copies
    #df4 = df3.query('num_seqids_motif_present > @min_num_seqids')

    return df3

def calc_motif_spacing(group):
    #group2 = group.sort_values(by = 'motif_start')
    group["dist_to_next_motif1"] = (group["motif_start"].shift(periods = -1, axis = 'index') - group['motif_end']).abs()
    group["dist_to_next_motif2"] = (group["motif_start"] - group['motif_end'].shift(periods = 1, axis = 'index')).abs()
    #this doesn't work for some reason
    #group["inter_motif_dist"] = group[["dist_to_next_motif1", "dist_to_next_motif2"]].min(axis='index', skipna=True)
    return group

def calc_motif_footprint(group):
    group["motif_width"] = group["motif_end"] - group["motif_start"]
    group["footprint"] = group.motif_width.sum()
    return group

def filter_fimo_output(df, max_q, max_spacing, max_footprint):
    """
    For each sequence/motif combination, filters fimo output for motifs that:
    1. Have a combined total width <= @max_footprint 
    2. The spacing in between instances of the same motif <= @max_spacing
    3. The combined qvalues < @max_q
    Note: if a motif instance is a singleton, it is dropped, the spacing will be NaN.
    It is dropped as a result.
    This is the desired behavior, b/c we want at least two motifs per seq.
    """
    
    #Calc total width of motif hits
    df2 = df.groupby(['sequence_name', 'motif_id']).apply(calc_motif_footprint)
    
    #calculate motif spacing
    df3 = (df2.sort_values(by = 'motif_start')
             .groupby(['sequence_name', 'motif_id'])
             .apply(calc_motif_spacing)
            )
    df3["inter_motif_dist"] = df3[["dist_to_next_motif1", "dist_to_next_motif2"]].min(axis=1, skipna=True)
    df4 = df3.drop(columns = ["dist_to_next_motif1", "dist_to_next_motif2"])
    
    #calculate a single qvalue for each motif on each sequence
    df5 = df4.groupby(by = ['motif_id', 'sequence_name']).apply(calculate_ign_qvalue)
    
    #drop sequences that don't pass the criteria
    df6 = df5.query('ign_motif_qvalue <= @max_q and inter_motif_dist <= @max_spacing and footprint <= @max_footprint')
    
    return df6

def select_ends_from_fimo(df, coordinate, max_dist, total_dist):
    """
    Filters a fimo dataframe to select the ends by:
     2. The motif with the best ign_motif_qvalue is selected
     3. The result is checked to make sure there are seqs up/downstream,
        and one of the seqs is within @max_dist from coordinate, and all are <= @total_dist. 
     4. If so, returns a dataframe; otherwise, repeat steps 2-3. If no motifs match the criteria, 
        returns an empty dataframe
    """
    #Make a column calculating the distance of the query to coordinate
    df['dist'] = (df['motif_start'] - coordinate)
    
    #Get the a list of the motifs, in order from best to worst
    motifs = (df.sort_values(by='ign_motif_qvalue')
                     .drop_duplicates('motif_id')
                     .motif_id
                     .unique()
                     .tolist()
                 )
    
    for motif in motifs:
        print(f'starting on {motif}')
        df2 = df.query('motif_id == @motif')
        
        num_igns = df2.ign.nunique()
        strands = df2.strand.nunique()
        down = df2.query('dist < 0 and dist.abs() <= @total_dist')
        up = df2.query('dist > 0 and dist.abs() <= @total_dist')
        adjacent = df2.query('dist.abs() <= @max_dist')

        
        if not down.empty and not up.empty and not adjacent.empty and num_igns <= 4 and strands > 1:
            #motifs pass all thresholds. 
            
            #get the seqids of the best hits
            if adjacent.sort_values('ign_motif_qvalue').iloc[0]['dist'] < 0:
                #best adjacent hit is downstream
                down_seqid = adjacent.sort_values('ign_motif_qvalue').iloc[0]['sequence_name']
                up_seqid = up.sort_values('ign_motif_qvalue').iloc[0]['sequence_name']
            else:
                #best adjacent hit is upstream
                up_seqid = adjacent.sort_values('ign_motif_qvalue').iloc[0]['sequence_name']
                down_seqid = down.sort_values('ign_motif_qvalue').iloc[0]['sequence_name']
            
            #return a dataframe of the ends
            df_ends = df2.query('sequence_name == @down_seqid or sequence_name == @up_seqid')
            
            #double check the motifs are on opposite strands
            #sometimes >2 IGNs are recovered, and can contribute to the 'strand' column. 
            #now that we picked two IGNs, assert that they are indeed opposite
            strands = df_ends.strand.nunique()
            if strands > 1:
            
                #Set proximal/distal columns
                df_ends["proximal_end"] = df_ends.set_index('ign').loc[:,"dist"].abs().idxmin()
                df_ends["distal_end"] = df_ends.set_index('ign').loc[:,"dist"].abs().idxmax()
                df_ends["tn_start"] = df_ends.motif_start.min()
                df_ends["tn_end"] = df_ends.motif_end.max()
                df_ends["tn_length"] = df_ends["tn_end"] - df_ends["tn_start"]
                return df_ends
        
    #No motifs pass the critera, return an empty dataframe
    df_ends = df.iloc[0:0]
    
    return df_ends

def run_streme(fasta, output, min_width=15, max_width=20):
    """
    """
    subprocess.run(f'streme --p {fasta} --oc {output} --minw {min_width} --maxw {max_width}',
                   shell = True,
                   check = True,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    
def run_tomtom(query, subject, output, distance_metric='pearson', min_overlap = 10,):
    """
    """
    subprocess.run(f'tomtom -oc {output} -dist {distance_metric} -min-overlap {min_overlap} {query} {subject}',
                    shell = True,
                    check = True,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    
def get_tomtom_results_as_df(tomtomtsv):
    """
    """
    df = (pd.read_csv(tomtomtsv, sep = '\t', comment = "#")
        .rename(columns = {"Query_ID" : "motif_id",
                           "Target_ID" : "target_motif",
                          "p-value" : "p_value",
                          "E-value" : "e_value",
                          "q-value" : "q_value"})
    )
    return df

def filter_tomtom(df, min_fraction=0.5):
    """
    Given a self vs. self tomtom output, 
    removes motifs that don't align with >= min_fraction of the input sequences
    If the input is just a single motif,
    returns the dataframe
    """
    num_motifs = len(df.motif_id.unique().tolist())
    if num_motifs < 2:
        return df
    #remove self-hits
    df2 = df.query('motif_id != target_motif')
    
    #count number of hits per query motif
    df2["num_hits"] = df2.groupby('motif_id')['target_motif'].transform('nunique')
    
    #remove motifs that hit fewer than min_fraction of target_seqs
    min_num_hits = num_motifs * min_fraction
    df3 = df2.query('num_hits >= @min_num_hits')
    #dropped_motifs = len(df2.query('motif_id not in @df3.motif_id').motif_id.tolist())
    #print(f'removed {dropped_motifs} motifs')
    return df3