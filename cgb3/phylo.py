"""Methods for multiple sequence alignment and tree construction."""
import copy
import pylab

from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
from Bio import Phylo as BioPhylo
from Bio.Phylo import NewickIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from cached_property import cached_property

from . import visualization
from . import misc

# Structure of a Nexus tree-only file
NEX_TEMPLATE = """\
#NEXUS

BEGIN TREES;
TRANSLATE
%(translate)s
  ;
Tree tree= %(tree)s
END;
"""


class Phylo:
    """Class for phylogeny.

    The class Phylogeny contains the phylogenetic tree initialized with the
    given collection of proteins. For the construction of the tree following
    are used:

    - Clustal Omega command-line tool for multiple sequence alignment
    - Distance matrix to compute the distance between two proteins
      (e.g. BLOSUM62, identity matrix)
    - Tree construction algorithm: UPGMA (Unweighted Pair Group Method with
      Arithmetic Mean) or NJ (Neighbor Joining)

    The class also provides methods for outputting the built phylogenetic tree,
    such as drawing it as a string as well as exporting it to a Newick file.
    """
    def __init__(self, proteins, names, distance_model='identity',
                 tree_algorithm='nj'):
        """Initializes a Phylo object.

        Args:
            proteins ([Protein]): list of Protein objects
            names ([String]): list of names to be used for each species on the
                tree.
            distance_model (string): see DistanceCalculator.protein_models
            tree_algorithm (string): 'nj' or 'upgma'
        """

        self._proteins = proteins
        self._names = names
        self._distance_model = distance_model
        self._tree_algorithm = tree_algorithm

    @property
    def proteins(self):
        """Returns Protein objects."""
        return self._proteins

    @cached_property
    def alignment(self):
        """Returns the multiple sequence alignment."""
        return self._clustalo()

    def proteins_to_fasta_file(self, filename):
        """Writes proteins to a temporary FASTA file."""
        with open(filename, 'w') as f:
            for i, protein in enumerate(self.proteins):
                f.write(protein.to_fasta(self._names[i]))

    def _clustalo(self):
        """Performs Clustal-Omega multiple sequence alignment.

        Args:
            proteins (list): List of Protein objects to be aligned
        Returns:
            MultipleSeqAlignment: A Bio.Align.MultipleSeqAlignment object.
        """
        infile = '/tmp/input.fasta'
        outfile = '/tmp/output.aln'
        self.proteins_to_fasta_file(infile)
        clustalo_cline = ClustalOmegaCommandline(
            'clustalo',             # executable
            infile=infile,          # input file name
            outfile=outfile,        # output file name
            outfmt='clustal',       # output format
            verbose=True,           # verbose output
            auto=True,              # set options automatically
            force=True)             # force file overwriting
        stdout, stderr = clustalo_cline()
        print(stderr)

        align = AlignIO.read(outfile, 'clustal')
        return align

    @cached_property
    def tree(self):
        """Returns a phylogenetic tree constructed from the given alignment."""
        calculator = DistanceCalculator(self._distance_model)
        constructor = DistanceTreeConstructor(calculator, self._tree_algorithm)
        tree = constructor.build_tree(self.alignment)
        # Make the tree rooted.
        tree.root_at_midpoint()
        tree.root.name = 'Root'
        return tree

    @cached_property
    def tree_lookup(self):
        lookup = {}
        for clade in self.tree.get_terminals():
            lookup[clade.name] = clade
        return lookup

<<<<<<< HEAD
    def tree_distance(self, pa, pb):
        """Finds the phylogenetic distance between two proteins."""
        print("dfgshk.jjjjjjjjfdgsssssssssshjk.........")
        
        clade_a = self.tree_lookup[pa.accession_number]
        print("ACESSION NUMBER", self.tree_lookup[pa.accession_number])
        clade_b = self.tree_lookup[pb.accession_number]
        print("HEY!")
        print("Clade A is ", clade_a.name)
        print("Clade B is ", clade_b.name)
=======
    def distance(self, pa, pb):
        """Finds the phylogenetic distance between two proteins."""
        clade_a = self.tree_lookup[pa.accession_number]
        clade_b = self.tree_lookup[pb.accession_number]
>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d
        return self.tree.distance(clade_a, clade_b)

    def draw_ascii(self):
        """Draws the tree in ASCII format."""
        BioPhylo.draw_ascii(self.tree)

    def to_newick(self, filename):
        """Writes the tree to the given file in newick format."""
        BioPhylo.write(self.tree, filename, 'newick')

    def to_nexus(self, filename):
        """Writes the tree to the given file in nexus format.

        This method doesn't call Bio.Phylo.NexusIO as BayesTraitsV2 requires a
        different dialect of Nexus format.
        """
        # Copy the tree before making any changes on it.
        tree = copy.deepcopy(self.tree)
        # BayesTraits requires the Nexus file to have a "Translate" block which
        # declares a number->taxon mapping so that numbers, not long taxa
        # names, are used in the tree descriptions.
        names_to_ints = dict((clade.name, i) for i, clade in enumerate(
            tree.get_terminals(), start=1))
        # Assign numbers to terminal clades
        for node in tree.get_terminals():
            node.name = str(names_to_ints[node.name])
        # Drop names of the inner nodes
        for n in tree.get_nonterminals():
            n.name = None
        # Tree to string
        writer = NewickIO.Writer([tree])
        nexus_tree = NEX_TEMPLATE % {
            'translate': ',\n'.join('%d %s' % (name, id)
                                    for id, name in list(names_to_ints.items())),
            'tree': next(writer.to_strings(plain=False, plain_newick=True))}
        # Write string to file
        with open(filename, 'w') as handle:
            handle.write(nexus_tree)

    def draw(self, filename):
        """Draws tree and saves it into the given file."""
        BioPhylo.draw(self.tree, do_show=False)
        pylab.savefig(filename)

    @property
    def svg_view(self):
        """Draws the tree in SVG format."""
        # Convert the tree from Biopython's Bio.Phylo.Tree to ETE3 Treenode
        t = visualization.biopython_to_ete3(self.tree)
        temp_file = misc.temp_file_name(suffix='.svg')
        visualization.tree_svg_plot(t, temp_file)
        with open(temp_file) as f:
            contents = f.read()
<<<<<<< HEAD
        return contents
=======
        return contents
>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d
