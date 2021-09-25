"""
Description
"""
# Imports --------------------------------------------------------------------------------------------------
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
#Turn off the SettingWithCopy warning, b/c danger is my middle name
pd.options.mode.chained_assignment = None
import subprocess
from pathlib import Path

# Functions ------------------------------------------------------------------------------------------------

def get_seqs_from_fasta(fastafile, iterable, output):
    """
    Accepts a python iterable object, gets seqs from a fasta file by name
    """

    outfile = Path(output)
    if outfile.exists():
        outfile.unlink()
    #outfile.mkdir(parents = True)
    
    seqs_to_get = set(iterable)
    assert len(seqs_to_get) > 0, "You did not pass any accessions"
    
    count = 0
    
    with open(fastafile) as f:
        with open(outfile, 'a') as outfile:
            for title,seq in SimpleFastaParser(f):
                if title in seqs_to_get:
                    count += 1
                    outfile.write(">%s\n%s\n" % (title, seq))
                    
    if count < len(seqs_to_get):
        print(f'You passed {len(seqs_to_get)} unique accession, but only {count} were written to {output}')
        
def get_interregions_gbk(genbank_path, output, output_bed, intergene_length=50, min_orf_size=300):
    """
    Records the start/stop coordinates of all annotated genes (CDS, tRNA, rRNA, ncRNA),
    extracts regions between them, including contig termini
    """
    #If the outputs exist, remove
    if Path(output).exists():
        Path(output).unlink()
    if Path(output_bed).exists():
        Path(output_bed).unlink()
    
    for seq_record in SeqIO.parse(genbank_path, "genbank"):
        
        #Check to make sure there is sequence data
        n_count = seq_record.seq.count('N')
        seq_len = len(seq_record.seq)
        assert n_count < seq_len and seq_len > 0, f"There doesn't seem to be any sequence data in {genbank_path}"
        
        #Get the name of the contig
        contig = seq_record.id.split("|")[0]
        
        cds_list = []
        intergenic_records = []

        # Loop over the genome file, get the gene features on each of the strands
        for feature in seq_record.features:
            if feature.type == "tRNA" or feature.type == "ncRNA" or feature.type == "rRNA" or feature.type == "RNA":
                mystart = feature.location.start.position
                myend = feature.location.end.position
                cds_list.append((mystart, myend))
            if feature.type == "CDS" or feature.type == "misc_feature":
                mystart = feature.location.start.position
                myend = feature.location.end.position
                if myend - mystart >= min_orf_size:
                    cds_list.append((mystart, myend))
                    
        if len(cds_list) == 0:
            print(f"There are no features annotated as CDS/tRNA/ncRNA/rRNA in {seq_record.id}, skipping")

        #Get the 5' end sequence
        if len(cds_list) > 0:
            first_feature_start = cds_list[0][0]
            if first_feature_start >= intergene_length:
                intergene_seq = seq_record.seq[0:first_feature_start]
                ign_num = "ign_" + "0"
                ign_start = 1
                ign_stop = first_feature_start
                intergenic_records.append(
                        SeqRecord(
                            intergene_seq,
                            id = f"{contig}|{ign_num}|{ign_start}|{ign_stop}",
                            description = ""
                        )
                    )

        
        #Get all intergenic loci
        if len(cds_list) > 1:
            for i, pospair in enumerate(cds_list[1:]):
                # Compare current start position to previous end position
                last_end = cds_list[i][1]
                this_start = pospair[0]
                if this_start - last_end >= intergene_length:
                    intergene_seq = seq_record.seq[last_end:this_start]
                    ign_num = "ign_" + str(i)
                    ign_start = last_end + 1
                    ign_stop = this_start
                    
                    intergenic_records.append(
                        SeqRecord(
                            intergene_seq,
                            id = f"{contig}|{ign_num}|{ign_start}|{ign_stop}",
                            description = ""
                        )
                    )

        #Get the 3' end sequence
        if len(cds_list) > 0:
            last_feature_end = cds_list.pop()[1]
            if len(seq_record.seq) - last_feature_end >= intergene_length:
                intergene_seq = seq_record.seq[last_feature_end:]
                ign_num = "ign_" + str(len(cds_list))
                ign_start = last_feature_end + 1
                ign_stop = len(seq_record.seq)
                intergenic_records.append(
                        SeqRecord(
                            intergene_seq,
                            id = f"{contig}|{ign_num}|{ign_start}|{ign_stop}",
                            description = ""
                        )
                    )
        #Write the outputs            
        with open(output, 'a') as f:
            SeqIO.write(intergenic_records, f, "fasta")
            
        with open(output_bed, 'a') as f:
            for seq_record in intergenic_records:
                contig,ign,start,stop = seq_record.id.split("|")
                print(contig,ign,start,stop, sep = ",", file = f)
                
def calculate_ign_coords(group):
    """
    Accepts a pandas groupby object from "get_ign_df" and adds the intergenic coordinates
    """
    group["ign_start"] = group["stop"] + 1
    group["ign_end"] = group["start"].shift(periods = -1, axis = 'index') - 1
    group["ign_length"] = group["ign_end"] - group["ign_start"]
    
    return group

def get_ign_df(feature_table, intergene_length=50, min_orf_size=300):
    """
    Expects a three-column CSV-delimited input with column names
    "accession", "start", "stop", "feature"

    Returns a dataframe with the intergenic coordinates/length
    that fit the criteria
    """
    df = pd.read_csv(feature_table)
    assert 'accession' and 'start' and 'stop' and 'feature' in df.columns.tolist(), f"""
    Your feature_table file {feature_table} is missing one of the following columns:
    'accession' 'start' 'stop' 'feature'
    """
    
    #The start/stop may be given as start < stop OR start > stop.
    #Force start < stop
    df["new_start"] = df.loc[:,["start", "stop"]].min(axis = 1)
    df["new_stop"] = df.loc[:,["start", "stop"]].max(axis = 1)
    df.drop(columns = ['start', 'stop'], inplace = True)
    df = df.rename(columns = {"new_start" : "start", "new_stop" : "stop"})
    
    #calculate the length
    df["length"] = df["stop"] - df["start"]
    
    #Remove short ORFs, add an extra row at the beginning
    df2 = df.query('feature != "CDS" or (feature == "CDS" and length >= @min_orf_size)')

    #Make columns of the intergenic coordinates
    #grouping by accession is necessary b/c 
    #there may be more than two accessions in the dataframe
    df3 = df2.groupby('accession').apply(calculate_ign_coords)
    
    #remove small intergenic loci
    df4 = df3.query('ign_length >= @intergene_length').reset_index(drop = True)
    
    return df4

def get_interregions_ft(feature_table, in_fasta, output, output_bed):
    """
    Expects a feature table with three-columns (CSV-delimited) labelled:
    "accession", "start", "stop", "feature"
    (Start can be > stop)
    
    Outputs the intergenic sequences fitting the above criteria
    """
    #If the outputs exist, remove
    if Path(output).exists():
        Path(output).unlink()
    if Path(output_bed).exists():
        Path(output_bed).unlink()
    
    #Make a dataframe of intergenic regions
    ign_df = get_ign_df(feature_table)
    accessions = ign_df.accession.tolist()
        
    for seq_record in SeqIO.parse(in_fasta, "fasta"):
        
        #Check to make sure there is sequence data
        n_count = seq_record.seq.count('N')
        seq_len = len(seq_record.seq)
        assert n_count < seq_len and seq_len > 0, f"There doesn't seem to be any sequence data in {in_fasta}"
        
        #Get the name of the contig
        contig = seq_record.id.split("|")[0]
        
        #Check if the contig has intergenic regions we want. 
        #(reset the indexes for labelling the ign_num below)
        ign_df2 = ign_df.query('accession == @contig').reset_index(drop=True)
        
        #Get the subsequences
        if ign_df2.shape[0] > 0:
            intergenic_records = []
            
            #Get the 5' end
            first_feature_start = ign_df2.iloc[0,3]
            intergene_seq = seq_record.seq[0:first_feature_start]
            ign_num = "ign_00"
            ign_start = 1
            ign_stop = first_feature_start
            intergenic_records.append(
                                    SeqRecord(
                                        intergene_seq,
                                        id = f"{contig}|{ign_num}|{ign_start}|{ign_stop}",
                                        description = ""
                                            )
                                    )
                
            #Get the intergenic seqs
            for index, row in ign_df2.iterrows():
                accession = row[0]
                start = row[3]
                stop = row[4]
                ign_start = int(row[6]) - 1    #Adjusting for 0-based
                ign_stop = int(row[7]) + 1     #Adjusting for 0-based
                
                intergene_seq = seq_record.seq[ign_start:ign_stop]
                ign_num = "ign_" + str(index)
                intergenic_records.append(
                        SeqRecord(
                            intergene_seq,
                            id = f"{contig}|{ign_num}|{ign_start}|{ign_stop}",
                            description = ""
                        )
                    )
            
            #Get the 3' end
            last_feature_end = ign_df2.iloc[-1,4]
            intergene_seq = seq_record.seq[last_feature_end:]
            ign_num = "ign_" + str(index + 1)
            ign_start = last_feature_end + 1
            ign_stop = len(seq_record.seq)
            intergenic_records.append(
                    SeqRecord(
                        intergene_seq,
                        id = f"{contig}|{ign_num}|{ign_start}|{ign_stop}",
                        description = ""
                    )
                )
            #Write the outputs            
            with open(output, 'a') as f:
                SeqIO.write(intergenic_records, f, "fasta")

            with open(output_bed, 'a') as f:
                for seq_record in intergenic_records:
                    contig,ign,start,stop = seq_record.id.split("|")
                    print(contig,ign,start,stop, sep = ",", file = f)
                    
def get_interregions_gff(fasta, gff, output_fasta, output_csv, intergene_length = 50):
    """
    """
    in_fasta = Path(fasta)
    out_fasta = Path(output_fasta)
    out_csv = Path(output_csv)
    
    if out_fasta.exists():
        out_fasta.unlink()
    
    #make a "genome" file for bedtools
    tmp_genome = out_csv.with_suffix('.tab')
    if tmp_genome.exists():
        tmp_genome.unlink()
        
    with open(tmp_genome, 'w') as o:
        for seq_record in SeqIO.parse(in_fasta, "fasta"):
        
            #Check to make sure there is sequence data
            n_count = seq_record.seq.count('N')
            seq_len = len(seq_record.seq)
            assert n_count < seq_len and seq_len > 0, f"""
            There doesn't seem to be any sequence data in {in_fasta}
            """

            #Get the name of the contig
            contig = seq_record.id.split("|")[0]
            
            #make a two-column file
            print(contig, seq_len, sep = '\t', file = o)
    
    #get a bedfile of the intergenic regions
    bedfile = out_csv.with_suffix('.bed')
    if bedfile.exists():
        bedfile.unlink()
        
    subprocess.run(f'bedtools complement -i {gff} -g {tmp_genome} > {bedfile}',
                       shell = True,
                       #capture_output = True,
                       check = True)
    
    assert bedfile.exists(), f"""
    something went wrong trying to extract intergenic regions from the fasta file {in_fasta}.
    Maybe check the file {tmp_genome}? This should be a two-column tab-delimited file with
    accession --> seq_length
    """
    
    #convert the output to a CSV file for reasons
    ign_df = pd.read_csv(bedfile, sep = '\t', names = ["accession", "start", "stop"])
    #calculate the length
    ign_df["length"] = ign_df["stop"] - ign_df["start"]
    
    
    #get the intergenic sequences
    #note: I do this so that I can name the headers. 
    #Otherwise I could just use bedtools getfasta
    for seq_record in SeqIO.parse(in_fasta, "fasta"):

        #Check to make sure there is sequence data
        n_count = seq_record.seq.count('N')
        seq_len = len(seq_record.seq)
        assert n_count < seq_len and seq_len > 0, f"There doesn't seem to be any sequence data in {in_fasta}"

        #Get the name of the contig
        contig = seq_record.id.split("|")[0]

        #Check if the contig has intergenic regions we want. 
        #(reset the indexes for labelling the ign_num below)
        ign_df2 = ign_df.query('accession == @contig and length >= @intergene_length').reset_index(drop=True)

        #Get the subsequences
        if not ign_df2.empty:
            intergenic_records = []

            for index, row in ign_df2.iterrows():
                accession = row[0]
                ign_start = row[1]
                ign_stop = row[2]

                intergene_seq = seq_record.seq[ign_start:ign_stop]
                ign_num = "ign_" + str(index)
                intergenic_records.append(
                        SeqRecord(
                            intergene_seq,
                            id = f"{contig}|{ign_num}|{ign_start}|{ign_stop}",
                            description = ""
                        )
                    )
            #Write the outputs            
            with open(out_fasta, 'a') as f:
                SeqIO.write(intergenic_records, f, "fasta")
                
            with open(out_csv, 'a') as f:
                for seq_record in intergenic_records:
                    contig,ign,start,stop = seq_record.id.split("|")
                    print(contig,ign,start,stop, sep = ",", file = f)
                
    #convert the output to a CSV file for reasons
    #ign_df.drop(columns = 'length').to_csv(out_csv, index = False)
    
    #remove the intermediate files
    tmp_genome.unlink()
    bedfile.unlink()
                
def run_dust(fastafile, dusted_fastafile, cutoff):
    """
    Runs the Dust script included in meme/bin/)
    """
    p1 = subprocess.run(f'dust {fastafile} {cutoff} > {dusted_fastafile}', 
                    shell = True,
                    check = True,
                    capture_output = True,
                    text = True)
                
def get_igns_near_tnsB_as_df(tnsB_contig_accession, tnsB_coordinate, ign_coordinates, windowsize = 250000):
    """
    Accepts a csv file of the intergenomic start/stop locations, and
            a csv file of the tnsB location
    outputs a dataframe of the intergenomic loci within a given windowsize
    Handles multiple TnsB's and coordinates given as start > stop
    """
    
    #make a dataframe of the intergenomic coordinates
    df_ign = pd.read_csv(ign_coordinates, 
                         header = None, 
                         names = ['accession', 
                                  'ign',
                                  'ign_start',
                                  'ign_stop'
                                 ]
                        )
   
    window_start = int(tnsB_coordinate) - windowsize/2
    window_stop = int(tnsB_coordinate) + windowsize/2
        
    #extract rows where the intergenic loci are within the window
    df_ign2 = df_ign.query('accession == @tnsB_contig_accession and ign_start >= @window_start and ign_stop <= @window_stop')

    return df_ign2

def get_igns_in_window(df_ign_window, ign_fasta, output):
    """
    Accepts a dataframe and python paths to a fasta file and output directory
    Constructs an accession from a dataframe
    extracts the subsequences from the fasta file, labelled by window #
    """

    #construct the fasta header
    df = df_ign_window.astype('str')
    df["fasta_header"] = df["accession"] + "|" + df["ign"] + "|" + df["ign_start"] + "|" + df["ign_stop"]
    
    #make a list of seqs to get
    fasta_seqs_to_get = df.fasta_header.tolist()
    
    #Get the seqs
    get_seqs_from_fasta(ign_fasta, fasta_seqs_to_get, output)

        
