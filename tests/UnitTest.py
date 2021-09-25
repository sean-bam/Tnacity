import unittest
from pathlib import Path
import subprocess

from tnacity_pkg import blast as B
from tnacity_pkg import fasta_lib as FU
from tnacity_pkg import meme as M



class FixturesTest(unittest.TestCase):

    def setUp(self):
        #these variables are initiated for every test function written below
        print('setting up tests..')
        self.input = Path('tests/NZ_CP067307.1.gb')
        self.ign_fasta = Path('tests/Tn7/ign/NZ_CP067307.fna')
        self.ign_csv = Path('tests/Tn7/ign/NZ_CP067307.csv')
        self.ign_window = Path('tests/Tn7/ign_blast/NZ_CP067307.1_window.fna')
        self.ign_window_dust = Path('tests/Tn7/ign_blast/NZ_CP067307.1_window.dust.fna')
        self.ign_blast = Path('tests/Tn7/ign_blast/NZ_CP067307.1_window.blastn')
        self.meme_background = Path('tests/Tn7/meme/NZ_CP067307.1.background')
        self.meme_fimo = Path('tests/Tn7/meme/fimo.tsv')
        self.output = Path('tests/Tn7/predicted_ends.csv')
        self.motifs = Path('tests/Tn7/tns_motifs.xml')
        
        output_ign_dir = Path('tests/Tn7/ign/')
        if not output_ign_dir.exists():
            output_ign_dir.mkdir(parents = True)
            
        output_blast_dir = Path('tests/Tn7/ign_blast/')
        if not output_blast_dir.exists():
            output_blast_dir.mkdir(parents = True)
            
        output_meme_dir = Path('tests/Tn7/meme/')
        if not output_meme_dir.exists():
            output_meme_dir.mkdir(parents = True)
        
    #def assertIsFile(self, path):
    #    if not Path(path).resolve().is_file():
    #        raise AssertionError("File does not exist: %s" % str(path))
    
    ###
    #The order in which tests are run is defined by the string of the function name (alphabetical)
    #so, I just name them in the order in which I want them to run
    ###
    
    def test1_getinterregions(self):
        """
          
        """
        FU.get_interregions_gbk(self.input, self.ign_fasta, self.ign_csv)
        assert self.ign_fasta.is_file() and self.ign_csv.is_file(), f'Could not extract seqs from {self.input}'
        
    def test2_getwindow(self):
        tn_neighborhood = FU.get_igns_near_tnsB_as_df('NZ_CP067307.1', 5042301, self.ign_csv)
        FU.get_igns_in_window(tn_neighborhood, self.ign_fasta, self.ign_window)
        assert self.ign_window.is_file(), f'couldnt extract seqs up/downstream from {self.ign_fasta}'
        
    def test3_dust(self):
        FU.run_dust(self.ign_window, self.ign_window_dust, 30)
        assert self.ign_window_dust.is_file(), f'couldnt run Dust from the MemeSuite, are you sure its on path?'
        self.ign_window_dust.rename(self.ign_window)
    
    def test4_markov(self):
        M.make_meme_background_file(self.ign_fasta, self.meme_background)
        assert self.meme_background.is_file(), f'running "fasta_get_markov" failed, are you sure its on your path?'
        #self.assertIsFile('tests/background2.meme')
        
    def test5_blast(self):
        B.run_blast(self.ign_window)
        assert self.ign_blast.is_file(), f'problem running BLASTn, are you sure its on path?'
        
    def test6_pipeline(self):
        subprocess.run("python bin/tnacity.py -i tests/NZ_CP067307.1.gb -infmt genbank -o tests/Tn7 -c 5042301 -a NZ_CP067307.1 -k",
                  shell = True)
        assert self.output.is_file(), f"Hmm, a file named {self.output} should be produced, but it isnt. Are you sure Meme is installed?'"
        
    
if __name__ == '__main__':
    unittest.main()