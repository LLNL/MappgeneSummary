"""
mappygene: aux python module for mappgene summary

"""
__all__ = ['__author__', '__date__', '__version__']
__author__ = 'Tyshawn Ferrell and Jose Manuel Martí'
__email__ = 'martimartine1 **AT** llnl.gov'
__status__ = 'Alpha'
__date__ = 'Jan 2022'
__version__ = '0.1.0'

# Python Standard Library
import io
from enum import Enum
import os
from pathlib import Path
from typing import Any, Dict, List, NamedTuple, NewType, Set, Tuple, Union

# 3rd party
import Bio
from Bio import SeqIO, SeqUtils, SeqRecord
from Bio.Seq import Seq
import numpy as np
import pandas as pd

# Here we define some constants and alike
debug: bool = False  # True here will print more stuff useful for debugging
refgen: Path = Path('NC_045512.2.fasta')  # Expected location for R/reticulate
if not refgen.is_file(): # Let's look for the file
    refgen = Path('inst/extdata/NC_045512.2.fasta')
    if not refgen.is_file():  # We are debugging into "/inst/python"
        debug = True
        refgen = Path('../extdata/NC_045512.2.fasta')

# Read the sequence into a SeqIO object
refseq = SeqIO.read(refgen, format='fasta')


# MUTATIONS
# Enumeration to list the possibilities for the target seq (of mutations)
class Target(Enum):
    """Enumeration with options for target element."""
    AA = 0  # AminoAcid
    NT = 1  # Nucleotide
    AMINOACID = AA  # This is an alias for AA
    NUCLEOTIDE = NT  # This is an alias for NT

    def __str__(self):
        return f'{str(self.name)}'


class Mutation(Dict[str, Any]):
    """Class representing a mutation as a dict"""

    def __init__(self, mutation_input: str,
                 mutation_target: Target = Target.AA):
        """Parse a mutation in a string and return a dict with useful info.

           Args:
             mutation_input: A string representing a SARS-CoV2 mutation
             mutation_target: Either MutationTarget.NT or .AA (default)"""

        super().__init__()
        try:
            mutation_gene, mutation_details = mutation_input.split(':')
        except ValueError:
            print(
                f"ERROR in Mutation {mutation_input}! Have you forgotten the gene?")
            raise
        self['gen'] = mutation_gene
        self['ref'] = mutation_details[0]  # ref(erence) is the 1st character
        self['mut'] = mutation_details[-1]  # mut(ated) is the last character
        self['loc'] = int(mutation_details[
                          1:-1])  # loc is all in-between those, and we convert to int
        self['tgt'] = mutation_target

    def __str__(self):
        return f"{self['gen']}:{self['ref']}{self['loc']}{self['mut']}"

# CODONS
aa2codons: Dict[str, List[str]] = {
    'A': ["GCT", "GCC", "GCA", "GCG"],
    'R': ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    'N': ["AAT", "AAC"],
    'D': ["GAT", "GAC"],
    'B': ["AAT", "AAC", "GAT", "GAC"],  # N or D
    'C': ["TGT", "TGC"],
    'Q': ["CAA", "CAG"],
    'E': ["GAA", "GAG"],
    'Z': ["CAA", "CAG", "GAA", "GAG"],  # Q or E
    'G': ["GGT", "GGC", "GGA", "GGG"],
    'H': ["CAT", "CAC"],
    'I': ["ATT", "ATC", "ATA"],
    'L': ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
    'K': ["AAA", "AAG"],
    'M': ["ATG"],
    'F': ["TTT", "TTC"],
    'P': ["CCT", "CCC", "CCA", "CCG"],
    'S': ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    'T': ["ACT", "ACC", "ACA", "ACG"],
    'W': ["TGG"],
    'Y': ["TAT", "TAC"],
    'V': ["GTT", "GTC", "GTA", "GTG"],
    '*': ["TAA", "TGA", "TAG"]
}


# GENE
# Derive a namedtuple to store the information for each gene
class Gene(NamedTuple):
    """Represents a gene."""

    ref_seq: SeqRecord.SeqRecord  # Reference genome
    name: str  # Name of the gene
    start: int  # Start NT position in the genome
    end: int  # End NT position in the genome

    def __len__(self) -> int:
        return self.end - self.start + 1

    def num_codons(self) -> int:
        """Get the number of codons"""
        if self.name == 'ORF1ab':  # Account for frameshift a-1
            return (self.end - 1 - self.start) // 3
        else:
            return (self.end - 2 - self.start) // 3

    def pos_gene2genome(self, pos: int,
                        src: Target = Target.NT,
                        tgt: Target = Target.NT) -> List[int]:
        """Convert a position in the gene to a position in the genome"""
        if src is Target.AA:
            pos = 3 * pos - 2
        if pos >= len(self):
            raise ValueError(f"{self.name}.pos_gene2genome:"
                             f" Position {pos} exceeds the length of the gene")
        pos = pos + self.start - 1
        if tgt is Target.AA:
            pos = (pos + 2) / 3
        elif src is Target.AA:
            return [pos, pos + 1, pos + 2]
        return [pos]

    def pos_genome2gene(self, pos: int,
                        src: Target = Target.NT,
                        tgt: Target = Target.NT) -> List[int]:
        """Convert a position in the genome to a position in the gene"""
        if src is Target.AA:
            pos = 3 * pos - 2
        if pos < self.start or pos > self.end:
            raise ValueError(f"{self.name}.pos_genome2gene:"
                             f" Position {pos} does not belong to the gene")
        pos = pos - self.start + 1
        if tgt is Target.AA:
            pos = (pos + 2) / 3
        elif src is Target.AA:
            return [pos, pos + 1, pos + 2]
        return [pos]

    def translate_(self, start: int = None,
                   end: int = None) -> SeqIO.SeqRecord.seq:
        """Translate all or part to the AA seq"""
        if start is None:
            start: int = self.start - 1
        if end is None:
            end: int = self.end
        return self.ref_seq[start:end].translate().seq

    def valid_mutation(self, mut: Mutation,
                       verb: bool = True) -> bool:
        """Check validity of a mutation (True/False)"""

        def vprint(*arg, **kwarg):
            if verb:
                print(*arg, **kwarg)

        # Check for name of gene
        if mut['gen'] != self.name:
            vprint(f"{self.name}.valid_mutation:"
                   f" Mutation does not belong to this gene!")
            return False
        # Check for nuc/aa in the reference
        ref: str
        if mut['tgt'] is Target.AA:
            ref = self.translate_()[mut['loc'] - 1]
        else:
            ref = self.ref_seq[mut['loc'] - 1]
        if ref != mut['ref']:
            vprint(
                f"Mutation of {mut['ref']} but reference has a {ref}"
                f" (target is {mut['tgt']})!")
            return False
        return True

    def explain_aa_mutation(self, mut: Mutation,
                            verb: bool = True) -> Dict[
        Seq, List[Tuple[int, str, str]]]:
        """Explain an AA mutation in terms of nucleotide mutation(s)

        Returns: Dict with a key for each possible codon and corresponding
          value is a list of possible nucleotide mutations"""

        def vprint(*arg, **kwarg):
            if verb:
                print(*arg, **kwarg)

                # First, check the validity

        if not self.valid_mutation(mut, verb=verb):
            raise ValueError("This mutation is not valid!")
        # Try to explain it
        #ref = self.translate_()[mut['loc'] - 1]
        lst = self.pos_gene2genome(mut['loc'], src=Target.AA, tgt=Target.NT)
        seq_ref = self.ref_seq[lst[0] - 1:lst[-1]].seq
        #        print(seq_ref.translate())
        #        print(aa2codons[mut['mut']])
        nuc_mutations: Dict[Seq, List[Tuple[int, str, str]]] = {}
        for codon in aa2codons[mut['mut']]:
            pos_mut: List[Tuple[int, str, str]] = []
            codon_seq: Seq = Seq(codon)
            vprint(
                f"\nIn {seq_ref} ({seq_ref.translate()}) to "
                f"{codon} ({codon_seq.translate()}): ",
                end='')
            for num, nuc in enumerate(codon):
                if nuc != seq_ref[num]:
                    vprint(f"{seq_ref[num]} could mute to {nuc}; ", end='')
                    pos_mut.append((lst[num - 1] + 1, seq_ref[num], nuc))
                nuc_mutations[codon_seq] = pos_mut
        return nuc_mutations

    def subgene(self, name: str, aa_ini: int, aa_end: int) -> "Gene":
        """Get an inner gene object, such as the coding region of a NSP"""
        nt_ini_abs = self.pos_gene2genome(aa_ini,
                                          src=Target.AA, tgt=Target.NT)[0]
        nt_end_abs = self.pos_gene2genome(aa_end,
                                          src=Target.AA, tgt=Target.NT)[0] + 2
        # Adding 2 NT to account until the last NT of the last codon
        return Gene(self.ref_seq, name, nt_ini_abs, nt_end_abs)

    def __repr__(self) -> str:
        return (f'<Gene {self.name} has {self.num_codons()} codons, '
                f'from {self.start} to {self.end} nt of the genome>')


# Define the genes to store the data
#
# References:
# 0) Wu, F., Zhao, S., Yu, B. et al. A new coronavirus associated with human
# respiratory disease in China. Nature 579, 265–269 (2020).
# 1) Supplementary Note 1 (Resolution of ambiguous gene names) from:
# Jungreis, I., Sealfon, R. & Kellis, M. "SARS-CoV-2 gene content and COVID-19
# mutation impact by comparing 44 Sarbecovirus genomes". Nat Commun 12, 2642 (2021).
# 2) Jungreis, I. et al, "Conflicting and ambiguous names of overlapping ORFs in
# the SARS-CoV-2 genome: A homology-based resolution". Virology 558, 145-151 (2021).
# 3) Khailany, Rozhgar A et al. “Genomic characterization of a novel SARS-CoV-2”.
# Gene reports vol. 19 (2020): 100682. doi:10.1016/j.genrep.2020.100682

NCS5UTR = Gene(refseq, 'NCS5UTR', 1, 265)
ORF1a = Gene(refseq, 'ORF1a', 266, 13483)
ORF1b = Gene(refseq, 'ORF1b', 13468, 21555)
ORF1ab = Gene(refseq, 'ORF1ab', ORF1a.start, ORF1b.end)

S = Gene(refseq, 'S', 21563, 25384)
ORF2b = Gene(refseq, 'ORF2b', 21744, 21860)  # Also named: S.iORF1

ORF3 = Gene(refseq, 'ORF3', 25393, 26220)
ORF3a = Gene(refseq, 'ORF3a', 25393, 26220)  # Assuming same pos than ORF3!
ORF3c = Gene(refseq, 'ORF3c', 25457, 25579)  # Also named: ORF3h, 3a.iORF1, ORF3b
ORF3d = Gene(refseq, 'ORF3d', 25524, 25694)  # Also named: ORF3b
ORF3d_2 = Gene(refseq, 'ORF3d_2', 25596, 25694)  # Also named: 3a.iORF2
ORF3b = Gene(refseq, 'ORF3b', 25814, 25879)  # 5’ end of SARS-CoV ORF3b ortholog

E = Gene(refseq, 'E', 26245, 26472)
M = Gene(refseq, 'M', 26523, 27191)
ORF6 = Gene(refseq, 'ORF6', 27202, 27387)
ORF7a = Gene(refseq, 'ORF7a', 27394, 27759)
ORF7b = Gene(refseq, 'ORF7b', 27756, 27887)
ORF8 = Gene(refseq, 'ORF8', 27894, 28259)

N = Gene(refseq, 'N', 28274, 29533)
ORF9b = Gene(refseq, 'ORF9b', 28284, 28574)  # Also named: ORF9a, N.iORF1
ORF9c = Gene(refseq, 'ORF9c', 28734, 28952)  # Also named: ORF9b, ORF14
ORF10 = Gene(refseq, 'ORF10', 29558, 29674)
NCS3UTR = Gene(refseq, 'NCS3UTR', 29675, 29903)

# Define the proteins (NSP and others) in ORF1a with AA
NSP1 = ORF1a.subgene('NSP1', 1, 180)
NSP2 = ORF1a.subgene('NSP2', 181, 818)
DUF3655 = ORF1a.subgene('DUF3655', 880, 1050)
# NSP3_LIKE = Gene(refseq, 'NSP3_LIKE', 1054, 1177)
# MACRO_SF = Gene(refseq, 'MACRO_SF', 1233, 1358)
# SUD_M = Gene(refseq, 'SUD_M', 1351, 1493)
NSP3 = ORF1a.subgene('NSP3', 819, 2763)
NSP4 = ORF1a.subgene('NSP4', 2764, 3263)
NSP5 = ORF1a.subgene('NSP5', 3264, 3569)
NSP6 = ORF1a.subgene('NSP6', 3570, 3859)
NSP7 = ORF1a.subgene('NSP7', 3860, 3942)
NSP8 = ORF1a.subgene('NSP8', 3943, 4140)
NSP9 = ORF1a.subgene('NSP9', 4141, 4253)
NSP10 = ORF1a.subgene('NSP10', 4254, 4392)
NSP11 = ORF1a.subgene('NSP11', 4393, 4406)

# Define the proteins (NSP and others) in ORF1b
NSP12 = ORF1ab.subgene('NSP12', 4393, 5324)
NSP13 = ORF1ab.subgene('NSP13', 5325, 5925)
NSP14 = ORF1ab.subgene('NSP14', 5926, 6452)
NSP15 = ORF1ab.subgene('NSP15', 6453, 6798)
NSP16 = ORF1ab.subgene('NSP16', 6799, 7095)

genes: List[Gene] = [S, E, M, N]
orfs: List[Gene] = [ORF1a, ORF1b, ORF1ab, ORF2b, ORF3, ORF3a, ORF3b, ORF3c,
                    ORF3d, ORF3d_2, ORF6, ORF7a, ORF7b, ORF8, ORF9b, ORF9b,
                    ORF10]
genes.extend(orfs)
ncs: List[Gene] = [NCS5UTR, NCS3UTR]
genes.extend(ncs)
nsps: List[Gene] = [NSP1, NSP2, DUF3655, NSP3, NSP4, NSP5, NSP6, NSP7, NSP8,
              NSP9, NSP10, NSP12, NSP13, NSP14, NSP15, NSP16]
genes.extend(nsps)  # Extend general gene list with nsp list


def pos2genes(pos:int, gene_name:str=None,
              tgt_genes:List[Gene]=None,
              exclude_gene_name:bool=False,
              out_type:str='dict',
              src: Target = Target.AA,
              tgt: Target = Target.AA,
              ddebug: bool = False,
              ) -> Union[Dict[str,int], List[str], str] :
    """Get dict of genes/rel_pos corresponding to pos in the genome or genes

    pos: position for the search, relative of absolute depending on next param
    gene_name: If not set the position will be absolute for the genome,
      otherwise it will be assumed to be a relative position to that gene.
    tgt_genes: List of genes used for the search, all the defined by default
    exclude_gene_name: if True and gene_name set, exclude it from the output
    out_type: 'dict' for dict of gene:position, 'list' for list of "gene.pos",
      or 'str' for just a string with "gene.pos gene.pos" and so forth.
    src: Target enum indicating kind of position for the source (gene_name)
    tgt: Target enum indicating kind of position for the targets (tgt_genes)
    ddebug: Double debug boolean
    """
    # Initialization
    if tgt_genes is None:
        tgt_genes = genes
    touched:Dict[str,int] = {}
    abs_pos:int  # Absolute position in the genome in NT
    # Get abs_pos
    if gene_name is None:
        if src is Target.AA:
            pos = 3 * pos - 2
        abs_pos = pos
    else:
        try:
            abs_pos = eval(gene_name).pos_gene2genome(
                pos, src=src, tgt=Target.NT)[0]
        except ValueError:  # Assume then the pos is absolute
            if src is Target.AA:
                pos = 3 * pos - 2
            abs_pos = pos
    if ddebug:
        print(f'{gene_name}: {pos} {src} converted to {abs_pos} {Target.NT}')
    # Get the genes and relative position containing the absolute position
    for gene in tgt_genes:
        rel_pos:int
        try:
            rel_pos = gene.pos_genome2gene(
                abs_pos, src=Target.NT, tgt=tgt)[0]
        except ValueError:  # No hit!
            if ddebug:
                print('No hit for', gene, 'for abs_pos =', abs_pos, Target.NT)
            pass
        else:  # Hit!
            if ddebug:
                print('HIT for', gene, 'for abs_pos =', abs_pos, Target.NT)
            touched[gene.name] = rel_pos if tgt is Target.NT else int(rel_pos)
    # Exclude gene_name from the output if it's there
    if exclude_gene_name:
        touched.pop(gene_name, None)
    # Convert the output to the desired format
    if out_type == 'list':
        retouched:List[str] = [
            f'{gen}.{int(pos)}' for (gen, pos) in touched.items()]
        return retouched
    if out_type == 'str':
        retouched:str = ' '.join([
            f'{gen}.{int(pos)}' for (gen, pos) in touched.items()])
        return retouched
    else:
        return touched


def mappgene_summary_populate_df_cols(df:pd.DataFrame,
                                      nsp:bool = True,
                                      out_type:str = 'dict'
                                      ) -> pd.DataFrame:
    """Populate column(s) of (R) dataframe

     * additional col with NSPs or all genes
     * additional col with ORF1ab to/from ORF1a/b

    NOTE: Function interfacing with mappgene_summary via R::reticulate
    """

    def df_apply_pos2genes(tgt_genes: List[Gene], col:str,
                           exclude_gene_name:bool = False) -> None:
        """Aux: use df.apply with pos2genes for each row of the df"""
        # CAUTION: row.POS is absolute (chromosome), no rel to row.GENE
        df[col] = df.apply(
            lambda row: pos2genes(row.POS, # gene_name=row.GENE,
                                  tgt_genes=tgt_genes, out_type=out_type,
                                  exclude_gene_name=exclude_gene_name,
                                  src=Target.NT, tgt=Target.AA),
            axis=1)

    # NSPs or all genes
    if nsp:
        df_apply_pos2genes(nsps, 'NSP')
    else:
        df_apply_pos2genes(genes, 'ALL_GENES')
    # ORF
    df_apply_pos2genes(orfs, 'ORF', exclude_gene_name=False)
    #df_apply_pos2genes(orfs, 'OTHER_ORF', exclude_gene_name=True)
    return df


def do_checks():
    """Perform some basic checks and print results"""
    print('>> Performing basic checks...')

    # Check the reference genome (a SeqIO object)
    print(refseq)

    # Set the input for the example code (a mutation described with AA)
    mut_name:str = "ORF1b:P314L"  # This is an example from outbreak.info to start with
    mut = Mutation(mut_name)
    print(mut)

    print('> Explain S:P681H mutation:', end='')
    S.explain_aa_mutation(Mutation('S:P681H'), verb=True)

    print('\n> Check ORF1ab overlap:')
    for gene in [ORF1a, ORF1b, S]:
        print(gene)
        print(f'{gene.translate_()[:10]}...{gene.translate_()[-10:]}')
    print(
        f'So the ORF1 overlap is from {ORF1b.start} to '
        f'{ORF1a.end} nt of the genome')

    print('\n> Check solution to Jeff\'s example:')  
    # Transform positions
    # -> from: source position is referred to nucleotide number from the start of the genome
    # -> to: target position is referred to aminoacid number from the start of the gene
    # NOTE: You can have multiple solutions for overlapping genes (ORF1ab, ORF1a, ORF1b)
    #
    abs_pos: int = 23367 # 
    src_pos: Target = Target.NT 
    tgt_pos: Target = Target.AA 
        
    for gene in genes:
        rel_pos: int
        try:
            rel_pos = gene.pos_genome2gene(
              abs_pos, src = src_pos, tgt = tgt_pos)[0]
        except ValueError:
            pass
        else:
            print(f'Position {abs_pos} ({src_pos}) in genome corresponds to '
                  f'position {int(rel_pos)} ({tgt_pos}) in {gene.name} gene')

    nsp3_rel_ini_aa: int = int(ORF1ab.pos_genome2gene(
        NSP3.start, src=Target.NT, tgt=Target.AA)[0])
    nsp3_rel_end_aa: int = int(ORF1ab.pos_genome2gene(
        NSP3.end, src=Target.NT, tgt=Target.AA)[0])
    print(
        f'\n> NSP3 has {NSP3.num_codons()} AAs, '
        f'from {nsp3_rel_ini_aa} to {nsp3_rel_end_aa} of ORF1ab.')
    nsp16_rel_ini_aa: int = int(ORF1ab.pos_genome2gene(
        NSP16.start, src=Target.NT, tgt=Target.AA)[0])
    nsp16_rel_end_aa: int = int(ORF1ab.pos_genome2gene(
        NSP16.end, src=Target.NT, tgt=Target.AA)[0])
    print(
        f'> NSP16 has {NSP16.num_codons()} AAs, '
        f'from {nsp16_rel_ini_aa} to {nsp16_rel_end_aa} of ORF1ab.')
    gene_name:str = "ORF1ab"
    mut_aa_pos:int = 4715
    print(
        f'\n> Mut in AA pos {mut_aa_pos} of {gene_name} corresponds to muts:'
        f' {pos2genes(mut_aa_pos, gene_name)}, \n  so to NSPs (if any):'
        f' {pos2genes(mut_aa_pos, gene_name, tgt_genes=nsps, out_type="str")}')
    gene_name:str = "ORF1ab"
    mut_nt_pos:int = 539
    all_muts = pos2genes(mut_nt_pos, gene_name, src=Target.NT, tgt=Target.AA)
    nsp_muts = pos2genes(mut_nt_pos, gene_name,
                         tgt_genes=nsps, out_type="str",
                         src=Target.NT, tgt=Target.AA, ddebug=False)
    print(
        f'\n> Mut in NT pos {mut_nt_pos} of {gene_name} corresponds to muts:'
        f' {all_muts}\n  so to NSPs (if any): {nsp_muts}')
    mut_abs_nt_pos:int = 14408
    all_muts = pos2genes(mut_abs_nt_pos, src=Target.NT, tgt=Target.AA)
    nsp_muts = pos2genes(mut_abs_nt_pos,
                         tgt_genes=nsps, out_type="str",
                         src=Target.NT, tgt=Target.AA, ddebug=False)
    print(
        f'\n> Mut in absolute NT pos {mut_abs_nt_pos} corresponds to muts:'
        f' {all_muts}\n  so to NSPs (if any): {nsp_muts}')
if __name__ == '__main__':
    # Check that we got it right
    if debug:
        do_checks()

