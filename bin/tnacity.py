#!/usr/bin/env python
"""
Identifies the boundaries of Tn7-like transposons
"""
# Imports --------------------------------------------------------------------------------------------------
import argparse
from pathlib import Path
import shutil
import pandas as pd

from tnacity_pkg import fasta_lib as FU
from tnacity_pkg import blast as B
from tnacity_pkg import meme as M

# Args -----------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(
        description='')
parser.add_argument('-i',
                    '--input',
                    help="/path/to/fasta or /path/to/genbank",
                    type = str,
                    required=True)

parser.add_argument('-f',
                    '--features',
                    help="/path/to/featureTable",
                    type = str,
                    required=False)

#parser.add_argument('-g',
#                    '--gff',
#                    help="/path/to/GFF file",
#                    type = str,
#                    required=False)

parser.add_argument('-infmt',
                    '--in_format',
                    help="specify 'fasta' or 'genbank'. Default: fasta",
                    type = str,
                    required=False)


parser.add_argument('-o',
                    '--output',
                    help="/path/to/outputDirectory/ (must be unique for each input)",
                    type = str,
                    required=True)

parser.add_argument('-c',
                    '--coordinate',
                    help="Nucleotide coordinate to serve as the anchor for the search",
                    type = int,
                    required=True)

parser.add_argument('-a',
                    '--accession',
                    help="accession of the contig to be searched",
                    type = str,
                    required=True)

parser.add_argument('-k',
                    '--keep',
                    help="keep all intermediate files",
                    action="store_true")

parser.add_argument('-t',
                    '--task',
                    type = str,
                    help="Identify motifs 'de_novo' or 'search' with motifs provided using the flag -m, default 'de_novo'",
                    required=False)

parser.add_argument('-m',
                    '--motifs',
                    help="/path/to/motifs (meme XML format)",
                    type = str,
                    required=False)

parser.add_argument('-bf_minlength',
                    '--blast_prefilter_min_length',
                    help="minimum alignment length of an inverted repeat detected with blast. Default = 15",
                    type = int,
                    required=False)

parser.add_argument('-bf_maxlength',
                    '--blast_prefilter_max_length',
                    help="max alignment length of an inverted repeat detected with blast. Default = 30",
                    type = int,
                    required=False)

parser.add_argument('-bf_maxgaps',
                    '--blast_prefilter_max_gaps',
                    help="maximum number of gaps in an inverted repeat detected with blast. Default = 2",
                    type = int,
                    required=False)

parser.add_argument('-bf_maxmm',
                    '--blast_prefilter_max_mismatches',
                    help="maximum number of mismatches in an inverted repeat detected with blast. Default = 6",
                    type = int,
                    required=False)

parser.add_argument('-bf_leash',
                    '--blast_prefilter_leash_len',
                    help="One of the inverted repeats detected with blast must be less than X bp from the provided coordinate. Default = 20kb",
                    type = int,
                    required=False)

parser.add_argument('-dust',
                    '--dust_cutoff',
                    help="Value for DUST low complexity filtering. Lower = more masking of low complexity regions. Default = 30",
                    type = int,
                    required=False)

parser.add_argument('-m_leash',
                    '--meme_leash_len',
                    help="One of the motifs detected with Meme must be less than X bp from the provided coordinate. Default = 20kb",
                    type = int,
                    required=False)

parser.add_argument('-m_FDR',
                    '--meme_max_FDR',
                    help="Maximum FDR rates (Q-values) to accept a motif hit. Default 0.1",
                    type = float,
                    required=False)

parser.add_argument('-m_cFDR',
                    '--meme_max_combined_FDR',
                    help="Maximum combined motif FDR rates (Q-values) for one input sequence. Default 0.01",
                    type = float,
                    required=False)

parser.add_argument('-m_space',
                    '--meme_inter_motif_space_max',
                    help="Maximum spacing in between motif hits. Default = 75bp",
                    type = int,
                    required=False)

parser.add_argument('-m_foot',
                    '--meme_max_footprint',
                    help="Maximum footprint size of all motif hits on a sequence. Default = 120bp",
                    type = int,
                    required=False)
args = parser.parse_args()

# Constants ------------------------------------------------------------------------------------------------

## file type, task
input_type = "fasta"
mode = "de_novo"

#fasta prefiltering
dust_cutoff = 30

## blast_prefilter
blast_min_len = 15 
blast_max_len = 30 
blast_max_gaps = 2 
blast_max_mismatches = 6
blast_max_ir_dist = 250000
blast_max_leash_len = 20000

#MEME options
meme_leash_len = 20000 
meme_total_dist = 125000
meme_max_q = 0.01
meme_max_combined_q = 0.01
meme_max_spacing = 75
meme_max_footprint = 120


# Functions ------------------------------------------------------------------------------------------------

def run_motif_searching(fasta, motifs, output_dir, background):
    
    pass

def get_ends(fimo_tsv, motif_file, out_table, out_motif):
    """
    """
    try:
        df_fimo1 = M.get_fimo_output_as_df(fimo_tsv)
        
        df_fimo2 = M.filter_fimo_output(df_fimo1, 
                                        meme_max_combined_q, 
                                        meme_max_spacing,
                                        meme_max_footprint)
        
        df_fimo3 = M.select_ends_from_fimo(df_fimo2, 
                                           args.coordinate, 
                                           meme_leash_len, 
                                           meme_total_dist)

        if not df_fimo3.empty:
            #get a table of the motifs
            df_fimo3["motif_file"] = out_motif
            df_fimo3.to_csv(out_table, index = False)
            
            #get the motif
            if Path(out_motif).exists():
                Path(out_motif).unlink()
            motif = df_fimo3.motif_id.unique().tolist()[0]
            M.get_motifs_from_meme(motif, motif_file, out_motif)
            
            #get the gff file
            #TBD
    except pd.errors.EmptyDataError:
        pass
    
                        


#-----------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    
    #check the inputs
    file = Path(args.input)
    
    #figure out what to do
    if args.task:
        mode = args.task
        assert mode == "de_novo" or mode == "search", f"""
        the flag --task must either be 'de_novo' or 'search',
        but you put '{args.task}'
        """
    if mode == "search":
        assert args.motifs, f"""
        You are trying to search fasta files with motifs (--task = 'search'), but you didn't set the flag
        --motifs
        """

        assert Path(args.motifs).is_file(), f"""
        Cannot find the file {args.motifs}
        """
    
    if args.in_format:
        input_type = args.in_format
        assert input_type == "fasta" or input_type == "genbank", f"""
    The input format must be 'fasta' or 'genbank', but you put '{args.in_format}'
    """

    if input_type == "fasta":
        assert args.features or args.gff, f"""
        You are trying to provide a fasta file as input, so you must also provide a feature table/GFF file
        """

    if args.dust_cutoff:
        dust_cutoff = args.dust_cutoff    
        
    if args.blast_prefilter_min_length:
        blast_min_len = args.blast_prefilter_min_length
    if args.blast_prefilter_max_length:
        blast_max_len = args.blast_prefilter_max_length
    if args.blast_prefilter_max_gaps:
        blast_max_gaps = args.blast_prefilter_max_gaps
    if args.blast_prefilter_max_mismatches:
        blast_max_mismatches = args.blast_prefilter_max_mismatches
    if args.blast_prefilter_leash_len:
        blast_max_leash_len = args.blast_prefilter_leash_len
    
    if args.meme_leash_len:
        meme_leash_len = args.meme_leash_len
    if args.meme_max_FDR:
        meme_max_q = args.meme_max_FDR
    if args.meme_max_combined_FDR:
        meme_max_combined_q = args.meme_max_combined_FDR
    if args.meme_inter_motif_space_max:
        meme_max_spacing = args.meme_inter_motif_space_max
    if args.meme_max_footprint:
        meme_max_footprint = args.meme_max_footprint
    
        
    basename = file.stem
    output_dir = Path(args.output)
    if not output_dir.exists():
        output_dir.mkdir(parents = True)

    #extract the intergenic regions
    output_ign_dir = output_dir / Path('ign/')
    if not output_ign_dir.exists():
        output_ign_dir.mkdir()
    ign_fasta = output_ign_dir / Path(basename).with_suffix('.fna')
    ign_csv = output_ign_dir /  Path(basename).with_suffix('.csv')
    if input_type == "genbank":
        FU.get_interregions_gbk(file, ign_fasta, ign_csv)
    elif args.gff:
        FU.get_interregions_gff(file, args.gff, ign_fasta, ign_csv)
    else:
        FU.get_interregions_ft(args.features, file, ign_fasta, ign_csv)
    
    #Make a file of intergenic sequences around TnsB
    output_blast_dir = output_dir / Path('ign_blast/')
    if not output_blast_dir.exists():
        output_blast_dir.mkdir()
    ign_window = output_blast_dir / Path(f'{args.accession}_window.fna')
    tn_neighborhood = FU.get_igns_near_tnsB_as_df(args.accession, args.coordinate, ign_csv)
    FU.get_igns_in_window(tn_neighborhood, ign_fasta, ign_window)
    
    #Dust the window of intergenic sequences
    ign_window_dust = output_blast_dir / Path(f'{args.accession}_window.dust.fna')
    FU.run_dust(ign_window, ign_window_dust, dust_cutoff)
    ign_window_dust.rename(ign_window)
    
    #setup the output directory
    output_meme_dir = output_dir / Path('meme/')
    if not output_meme_dir.exists():
        output_meme_dir.mkdir(parents = True)
        
    #make a background file
    backgroundfile = output_meme_dir / Path(f'{basename}.background')
    if not backgroundfile.exists():
        print('making backgroundfile')
        M.make_meme_background_file(ign_fasta, backgroundfile)
    
    if mode == 'search':
        #search for ends with the given motifs
        fimo_results = output_meme_dir / Path('fimo.tsv')
       
        M.run_fimo(ign_window,
                 args.motifs, 
                 output_meme_dir,
                 backgroundfile,
                 meme_max_q)
        
        if fimo_results.exists():
            out_csv = output_dir / Path('predicted_ends.csv')
            out_motif = output_dir / Path('tns_motifs.xml')
            get_ends(fimo_results, args.motifs, out_csv, out_motif)
            
        
    else:
        #predict ends de novo
        #Blast each TnsB window against itself
        B.run_blast(ign_window)

        #filter the blast output
        ign_blast_output = output_blast_dir / Path(f'{args.accession}_window.blastn')
        blast_df, count_df = B.run_blast_filtering(ign_blast_output, 
                                                   args.coordinate, 
                                                   blast_min_len, 
                                                   blast_max_len, 
                                                   blast_max_gaps, 
                                                   blast_max_mismatches,
                                                   blast_max_ir_dist,
                                                   blast_max_leash_len)
        ### Save the outputs
        if blast_df.empty:
            print("No IRs could be found with blast, try relaxing the blast prefiltering options")
        else:
            filtered_blast = output_blast_dir / Path(f'{args.accession}_IRs.csv')
            filtered_counts = output_blast_dir / Path(f'{args.accession}_counts.csv')
            blast_df.to_csv(filtered_blast, index = False)
            count_df.to_csv(filtered_counts, index = False)


        #Run meme
        if not blast_df.empty:
            blast_df2 = B.deduplicate_reciprocal_pairs(blast_df)

            #get input seqs
            for index, row in blast_df2.iterrows():
                accession = row[0]
                ign1 = row[1]
                ign2 = row[2]
                ign_seqids = {row[3], row[4]}

                #output paths/files expected
                p = output_meme_dir / Path(f'{ign1}_vs_{ign2}')
                out_fasta = p / "meme_in.fna"

                if out_fasta.exists():
                    out_fasta.unlink()

                if not p.exists():
                    p.mkdir(parents=True)

                FU.get_seqs_from_fasta(ign_window, ign_seqids, out_fasta)

            df_list = []
            print(f'running motif detection using MEME')
            for meme_infile in output_meme_dir.rglob('*.fna'):
                #print(f'running meme on {meme_infile}')
                M.run_meme(meme_infile, meme_infile.parent, backgroundfile)
                df = M.get_meme_results_as_df(meme_infile.parent / Path('meme.xml'))
                df2 = M.filter_meme_results(df)

                #Add a column pointing towards the path of the meme output
                df2["motif_file"] = meme_infile.parent / Path('meme.xml')
                df_list.append(df2)

            #Extract the motif(s)
            df3 = pd.concat(df_list)
            for index, row in df3.drop_duplicates(subset = 'motif_name').loc[:,['motif_name', 'motif_file']].iterrows():
                motif_name = row[0]
                motif_path = row[1]

                motif_output = output_meme_dir / Path(f'selected_motifs.xml')
                M.get_motifs_from_meme(motif_name, motif_path, motif_output)

            #Search with the motifs
            motif_output = output_meme_dir / Path(f'selected_motifs.xml')
            if motif_output.exists():
                fimo_results = output_meme_dir / Path('fimo.tsv')
            
                M.run_fimo(ign_window,
                         motif_output, 
                         output_meme_dir,
                         backgroundfile,
                         meme_max_q)

                if fimo_results.exists():
                    #get at table of the motifs
                    out_csv = output_dir / Path('predicted_ends.csv')
                    out_motif = output_dir / Path('tns_motifs.xml')
                    get_ends(fimo_results, motif_output, out_csv, out_motif)

            
    #clean up
    if not args.keep:
        shutil.rmtree(output_ign_dir)
        shutil.rmtree(output_blast_dir)
        shutil.rmtree(output_meme_dir)

        #if out_csv.exists():
            #Get the directory of the chosen motif
            #motif_file = Path(df_mc5.motif_file.tolist()[0]) 
            #motif_dir = motif_file.parent

            #remove all other directories
            #for x in output_meme_dir.iterdir():
            #    if x.is_dir():
            #        if Path(x) != motif_dir: 
             #           shutil.rmtree(x)
            
            
        #remove empty folders
        try:
            output_dir.rmdir()
        except OSError:
            pass