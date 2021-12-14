"""Miscellaneous bioinformatics utility functions."""

from subprocess import PIPE, Popen

from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna


def complement(seq):
    """Returns the complement of the seq.

    Args:
        seq (string): the DNA sequence.
    Returns:
        string: the complement sequence
    """
    dna_alphabet = "ACGT"
    return str(Seq(seq, dna_alphabet).complement())


def reverse_complement(seq):
    """Returns the reverse complement of the given sequence.

    Args:
        seq (string): the DNA sequence.
    Returns:
        string: the reverse complement sequence
    """
<<<<<<< HEAD
    dna_alphabet = "ACGT"
    return str(Seq(seq).reverse_complement())
#return str(Seq(seq, dna_alphabet).reverse_complement())

=======
    return str(Seq(seq).reverse_complement())
>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d

def weblogo(seqs, filename):
    """Generates the sequence logo for the given sequences.

    Uses WebLogo program (http://weblogo.threeplusone.com/).
    """
    # Sequences in FASTA format
    fasta = '\n'.join('>seq%d\n%s' % (i, seq) for i, seq in enumerate(seqs))
    p = Popen(['weblogo',
               '--format', 'png',
               '--fout', filename,
               '--color-scheme', 'classic',
               '--errorbars', 'YES'],
              stdout=PIPE, stderr=PIPE, stdin=PIPE, close_fds=True)
<<<<<<< HEAD
    p.communicate(input=fasta.encode())
    #p.communicate(input=fasta)
=======
    p.communicate(input=fasta.encode())
>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d
