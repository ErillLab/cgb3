"""Module for orthologous groups."""

import csv

from tqdm import tqdm
import networkx as nx

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# from Bio.Alphabet import IUPAC

from . import misc
from . import visualization
from . import bayestraits_wrapper
from .my_logger import my_logger
from .misc import mean
<<<<<<< HEAD
from ete3 import Tree
=======
>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d

from .hmmer import run_COG_hmmscan, process_COG_hmmscan,\
    run_eggNOG_hmmscan, process_eggNOG_hmmscan, \
    run_PFAM_hmmscan, process_PFAM_hmmscan


class OrthologousGroup:
    """Class definition for OrthologousGroup.
    
    The OrthologousGroup class holds a group of genes 
    which are orthologous to each other. Two genes are 
    determined as orthologs if they are best BLAST hits 
    for each other (reciprocal best BLAST hits).

    The genes in the orthologous group are sorted by 
    their posterior probabilities of regulation.

    """
    def __init__(self, genes):
        self._genes = genes
        self._genes.sort(key=lambda g: g.operon.regulation_probability, reverse=True)
<<<<<<< HEAD
        self._weighted_average_prob_regulation = None
=======
>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d
        self._COGs = []
        self._NOGs = []
        self._PFAMs = []
        
    @property
    def genes(self):
        """Returns the list of orthologous genes."""
        return self._genes
<<<<<<< HEAD
    
    @property
    def weighted_average_prob_regulation(self):
        return self._weighted_average_prob_regulation
=======
>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d

    @property
    def description(self):
        """Returns the description for the genes in the group."""
        gene_descs = [gene.product for gene in self.genes if gene.product]
        #default to first description, if any
        ret_desc = gene_descs[0] if gene_descs else ""

        #try to return a meaningful description
        for desc in gene_descs:
            if not(('hypothetical' in desc) or ('hypothetical' in desc)):
                ret_desc = desc
                break
        return ret_desc

    @property
    def COGs(self):
        """Returns the list of NOG IDs assigned to this group."""
        return self._COGs

    @property
    def NOGs(self):
        """Returns the list of NOG IDs assigned to this group."""
        return self._NOGs

    @property
    def PFAMs(self):
        """Returns the list of NOG IDs assigned to this group."""
        return self._PFAMs
    
    def member_from_genome(self, genome_name):
        """Returns the member of the group from the given genome.

        If there are multiple genes in the group from the given genome, the one
        with the maximum posterior probability of regulation is selected.

        Returns None if the specified genome has no genes in the group.

        """
        genes = [g for g in self.genes
                 if g.genome.strain_name == genome_name]
        if genes:
            return genes[0]
        return None

    def all_genes_from_genome(self, genome_name):
        """Returns all genes of the group from the given genome.

        Returns empty list if the given genome has no genes in the group.
        """
        genes = [g for g in self.genes if g.genome.strain_name == genome_name]
        return genes

    def assign_NOGs(self,userin):
        """Calls hmmscan to obtain matches for first protein in the group.
           Invokes process_hmmscan to obtain a list of NOGs according to user
           input preferences (maxNOG number, e-value jump), and assings the
           resulting list to the groups' NOG property.
        """

        #take first proper (protein coding) gene as query
        query=None
        for g in self.genes:
            if g.is_protein_coding_gene:
                query=g
                break
		
		#if the group contains at least a valid protein coding gene
        if query:
			#obtain protein sequence and generate seq object
            query_record=SeqRecord(Seq(query.translate, IUPAC.protein),\
						 id=query.protein_accession_number,name=query.name,\
						 description=query.product)

            #invoke hmmscan
            run_eggNOG_hmmscan(query_record,userin)

            #process result from hmmscan and assign NOGs to orthologous group
            self._NOGs=process_eggNOG_hmmscan(userin)
        else:
            self._NOGs=[]
            
    def assign_PFAMs(self,userin):
        """Calls hmmscan to obtain matches for first protein in the group.
           Invokes process_hmmscan to obtain a list of NOGs according to user
           input preferences (maxNOG number, e-value jump), and assings the
           resulting list to the groups' NOG property.
        """

        #take first proper (protein coding) gene as query
        query=None
        for g in self.genes:
            if g.is_protein_coding_gene:
                query=g
                break
		
		#if the group contains at least a valid protein coding gene
        if query:
			#obtain protein sequence and generate seq object
            query_record=SeqRecord(Seq(query.translate, IUPAC.protein),\
						 id=query.protein_accession_number,name=query.name,\
						 description=query.product)

            #invoke hmmscan
            run_PFAM_hmmscan(query_record,userin)

            #process result from hmmscan and assign NOGs to orthologous group
            self._PFAMs=process_PFAM_hmmscan(userin)
        else:
            self._PFAMs=[]            
            
    def assign_COGs(self,userin):
        """Calls hmmscan to obtain matches for first protein in the group.
           Invokes process_hmmscan to obtain a list of NOGs according to user
           input preferences (maxNOG number, e-value jump), and assings the
           resulting list to the groups' NOG property.
        """

        #take first proper (protein coding) gene as query
        query=None
        for g in self.genes:
            if g.is_protein_coding_gene:
                query=g
                break
		
		#if the group contains at least a valid protein coding gene
        if query:
			#obtain protein sequence and generate seq object
            query_record=SeqRecord(Seq(query.translate, IUPAC.protein),\
						 id=query.protein_accession_number,name=query.name,\
						 description=query.product)

            #invoke hmmscan
            run_COG_hmmscan(query_record,userin)

            #process result from hmmscan and assign NOGs to orthologous group
            self._COGs=process_COG_hmmscan(userin)
        else:
            self._COGs=[]                        
            
    def discretize_regulation_states(self, phylo):
        """Discretizes the trait of regulation for all genes in the orthologous
           group.

        Each gene in the orthologous group has a posterior probability of
        regulation that is computed using the binding model. For each gene in
        the orthologous group, this method chooses one of the two possible
        states of the trait, (1) regulation and (2) not regulation, based on
        the posterior probability of the regulation of the gene.

        Returns: {string: int}: the dictionary containing (key, value) pairs
            where key is the accession number of the regulated gene, and the
            value is 0/1 indicating the regulation trait (or absence A).
        """

        #define the probabilities for each state on terminal nodes
        #(where prob of regulation is known)
        terminal_states = self.get_terminal_states_probs(phylo)
        #initialize discrete trait dictionary
        trait = {}
        #for each terminal node in the phylogeny
        for node in phylo.tree.get_terminals():
            states = ['1', '0', 'A']
        	#grab the node probabilities for each of the 3 states into a vector
            probabilities = [terminal_states[(node.name, state)]
                             for state in states]
        	#sample the vector of states using their associated probablities
        	#and assign the sampled discrete state (1, 0 or A) to the trait
        	#dictionary entry for that node [i.e. species]
            trait[node.name], = misc.weighted_choice(states, probabilities)

        #at the end of the process, we obtain a dictionary in which we have
        #sampled, for each species, one of the three possible states according
        #to the probability of regulation by the TF in that species, for the
        #gene instances belonging to this orthologous group
        return trait

    def bootstrap_traits(self, phylo, sample_size):
        """Sample discrete traits for each gene in the ortho_group, as a list.
        Calls 'sample_size' times the 'discretize_regulation_states(phylo)'
        method, generating 'sample_size' samples of the discretized states of
        each terminal node in the tree.

        This is returned as a list of dictionaries, where each dictionary
        contains entries for each terminal node in the tree and its discrete
        state.

        Each instance of the sample contains the discrete traits of regulation
        associated with each gene.
        """
        return [self.discretize_regulation_states(phylo)
            for _ in range(sample_size)]


    """For a given phylogeny and set of orthologous genes [ortho_group], it
    returns a dictionary, indexed by genome_name and the three possible
    states: "1" for regulation, "0" for non-regulation and "A" for absent,
    encoding the probability of each state in each species.
    
    For each genome, the code looks at the member in the orthologous group.
    If the genome does not have an instance of the gene, it assigns 1 to
    absent ("A"), and 0 to states "1" and "0".
    If the genome does have an instance of the gene, it assigns to "1" the
    probability of regulation in that species, to "0" 1-prob and to "A" 0.
	"""
    def get_terminal_states_probs(self, phylo):
        states = {}
        genome_names = [node.name for node in phylo.tree.get_terminals()]
        #for each genome in the analysis (looking at this [self] ortho_group
        for genome_name in genome_names:
            # Check if the group contains a gene from the current genome
            gene = self.member_from_genome(genome_name)
            if gene:
            	#if so, assign the corresponding probabilities of regulation
                #and zero probability for absence
                p_reg = gene.operon.regulation_probability
                states[(genome_name, '1')] = p_reg
                states[(genome_name, '0')] = 1 - p_reg
                states[(genome_name, 'A')] = 0
            else:
                #No gene in this orthologous group from the genome
                #Assign zero prob to both regulation states, and 1 to Absent
                states[(genome_name, '1')] = 0
                states[(genome_name, '0')] = 0
                states[(genome_name, 'A')] = 1
        return states

    def ancestral_state_reconstruction(self, phylo, user_input):
        """Runs BayesTraits for ancestral state reconstruction.

        It estimates whether the gene is likely to be present in ancestral
        nodes, as well as its regulation, if the gene is present.
        """
        #define possible state list
        states = ['1', '0', 'A']
        #define dictionary with dimensions [states]x[non_terminal_nodes],
        #initialized to zeroes
        bootstrap_inferred_states = {(node.name, state): 0
                                     for state in states
                                     for node in phylo.tree.get_nonterminals()}

        #read sample size from options
        sample_size=user_input.bootstrap_replicates

        #generate 'sample_size' replicates of the bootstrapped discretized
        #states for the terminal nodes in the tree
        for trait in self.bootstrap_traits(phylo, sample_size):
        	#for each of those replicates, call bayestraits and infer ancestral
        	#states
            inferred_states = bayestraits_wrapper.bayes_traits(phylo, trait)
            for node in phylo.tree.get_nonterminals():
                #for each non-terminal node, get the probability assigned by
                #BayesTraits to each state, and add it to that state/node tally
                for state in states:
                    k = (node.name, state)
                    bootstrap_inferred_states[k] += inferred_states.get(k, 0)
        #Normalize the added probabilities for the non-terminal states
        #by dividing by the sample_size (number of boostrap replicates)
        nonterminal_states = {k: v/sample_size
                              for (k, v) in list(bootstrap_inferred_states.items())}
        #join the terminal and non-terminal states into a single dictionary
        #indexed by genome_name and state, and containing the observed (for
        #terminal) and inferred (for non-terminal, coming from boostrap) probab
        #of regulation [and absence]
        all_states = dict(list(self.get_terminal_states_probs(phylo).items()) +
                          list(nonterminal_states.items()))
        # Store these reconstructed+observed ancestral states into the ortholog
        #group '_regulation_states' field, as a dictionary
        self._regulation_states = all_states
<<<<<<< HEAD
        
    def assign_weighted_average_prob_regulation(self, ortho, phylogeny, prior_reg, reg_choice):
        "Obtains the pairwise distances for the members of an orthologous group and computes the average"
        pairwise = [] #To store the pairwise distances
        n = 0 #The number of genes for a species
        ortho_genes = []         
       
        #Get all the genes for each species
        for node in phylogeny.tree.get_terminals():
            matching_genes = [g for g in ortho.genes \
                              if g.genome.strain_name == node.name]
            if len(matching_genes) > 0:
                ortho_gene = {}
                ortho_gene['species'] = node.name
                ortho_gene['gene'] = matching_genes[0]
                ortho_gene['prob'] = ortho_gene['gene'].operon.regulation_probability
                ortho_genes.append(ortho_gene)
                
            elif reg_choice:
                ortho_gene = {}
                ortho_gene['species'] = node.name
                ortho_gene['gene'] = 'None'
                ortho_gene['prob'] = prior_reg
                ortho_genes.append(ortho_gene)
            n = n+1
        
        i = 0
        while i < len(ortho_genes):
            j = i+1
            while j < len(ortho_genes):
                pairwise.append(
                    ((ortho_genes[i]['prob'] + ortho_genes[j]['prob']) *
                     phylogeny.tree.distance(ortho_genes[i]['species'], ortho_genes[j]['species']))
                    /2
                    )
                j+=1
            i+=1

        #Compute the average
        n_choose_two= n * (n - 1)/2
        
        if len(pairwise) == 0:
            print("Len 0:")
            self._weighted_average_prob_regulation = ortho.genes[0].regulation_probability/n_choose_two
        else:
            #self._weighted_average_prob_regulation = (sum(pairwise))/comparisons
            self._weighted_average_prob_regulation = (sum(pairwise)/n_choose_two)
=======
>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d

    @property
    def regulation_states(self):
        """Gets the field _regulation_states"""
        return self._regulation_states

    @property
    def prob_regulation_at_root(self):
        """Returns the probability of regulation at the root of the tree.
           'Root' here standing for the 'genome_name' and '1' for the prob of
           regulation state.
        """
        return self.regulation_states[('Root', '1')]

    def most_likely_state_at(self, node_name):
        """Returns the most likely state at the given node.
           The function calls 'max' using the 'regulation_states' value for each
           state in that node as the 'key' for sorting [and inferring the max]
        """
        return max(['1', '0', 'A'],
                   key=lambda x: self.regulation_states[(node_name, x)])

    def ancestral_state_reconstruction_svg_view(self, phylo):
        """Returns an SVG-parsed view of the reconstructed ancestral states
           using ete3 visualization options to generate a temp SVG file and
           then read it back to return the SVG parse
        """
        temp_file = misc.temp_file_name(suffix='.svg')
        t = visualization.biopython_to_ete3(phylo.tree)
        visualization.view_by_gene(t, self, temp_file)
        with open(temp_file) as f:
            contents = f.read()
        return contents

    def __repr__(self):
        return str(self.genes)


# Class-associated functions
#
# The following functions provide the means to instantiate orthologous groups
# from a pre-defined subset of genes in all genomes under analysis and to
# export them in CSV format.


def construct_orthologous_groups(genes, genomes, cache,h_eval):
    """Constructs orthologous groups starting with the given list of genes.

    For each genome, candidate genes that are identified as likely to be
    regulated are tagged for orthology detection.

    This constructor function receives the genome objects and the list
    of genes from each of these genomes on which reciprocal BLAST will be
    applied to infer orthologs.

    For each gene, it identifies the reciprocal best BLAST hits in other
    genomes and adds the gene and its orthologs to the orthologous group.
	This is done in a serial manner, going through each gene in each
	genome and adding its rbbh's as orthologs (except for those that have
	already been incorporated to another group).
	At the end of this process "each gene" will have spawned an orthologous
	group (or it wil have been incoroporated into a preexisting one).
	These are "primordial" orthologous groups.

	The function merge_orthologous_groups is then called.
	This function takes all these "primordial" orthologous groups and
	generates an interconnection graph with them. If two groups intersect at
	any node, they are collapsed into a bigger group.

    Each orthologous group is hence a list of gene objects that have
	been found to be connected via best-reciprocal BLAST hit relationships.

    The function returns a list of such orthologous groups.
    """
    my_logger.info("Constructing orthologous groups.")
    groups = []
    for gene in tqdm(genes):
        # Check whether gene is already in a group; if it is, it skips the gene
        if any(gene in grp.genes for grp in groups):
            continue
        # If the gene is not in any group, create list of orthologous genes (group)
        # by performing reciprocal BLAST against all other genomes.
        # We will end up with a group for each gene (except those that were
        # picked up by a previous group
        rbhs = [gene.reciprocal_blast_hit(other_genome, cache,h_eval)
                for other_genome in genomes if gene.genome != other_genome]
        # Create the orthologous group
        grp = OrthologousGroup([gene] + [rbh for rbh in rbhs if rbh])
        groups.append(grp)

    # Collapse groups that contains a same gene.
    my_logger.info("Collapsing orthologous groups.")
    return merge_orthologous_groups(groups) #call merge function


def merge_orthologous_groups(groups):
    """Merges intersecting orthologous groups."""
    G = nx.Graph()
    for grp in groups:
        #all genes in group become nodes
        G.add_nodes_from(grp.genes)
        #all genes in group get connected (straight neighbor connectivity)
        G.add_edges_from(list(zip(grp.genes, grp.genes[1:])))

    merged_grps = []
    #for each list of neighbor connections, including any overlaps
    #note that the connectivity might travel backwards to a genome, detecting
    #paralogs through connections in other genomes
    for cc in nx.connected_components(G):
        #create a new orthologous group, based on that connectivity
        merged_grps.append(OrthologousGroup(list(cc)))
    return merged_grps


<<<<<<< HEAD
def orthologous_grps_to_csv(groups, phylogeny, filename, weight_choice):
=======
def orthologous_grps_to_csv(groups, phylogeny, filename):
>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d
    species = phylogeny.tree.find_elements(terminal=True, order='postorder')
    genome_names = [node.name for node in species]
    with open(filename, 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        header_row = (['average_probability',
<<<<<<< HEAD
                       'average_probability_all', 'weighted_probability'
=======
                       'average_probability_all',
>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d
                       'ortholog_group_size', 'description', 'COGs', 'eval', \
                       'NOGs', 'eval', 'PFAMs', 'eval'] +
                      [field for genome_name in genome_names
                       for field in ['probability (%s)' % genome_name,
                                     'locus_tag (%s)' % genome_name,
                                     'protein_id (%s)' % genome_name,
                                     'product (%s)' % genome_name,
                                     'operon id (%s)' % genome_name,
                                     'paralogs (%s)' % genome_name]])
        csv_writer.writerow(header_row)
        csv_rows = []
        for group in groups:
            genes = [group.member_from_genome(genome_name)
                     for genome_name in genome_names]
            # Average regulation probability
            avg_p = mean([g.operon.regulation_probability
                          for g in genes if g])
            # Average regulation probability (p=0 for absent genes in the grp)
            avg_p_all = mean([g.operon.regulation_probability if g else 0
                              for g in genes])
<<<<<<< HEAD
            #Weighted regulation probability
            weighted_p = group.weighted_average_prob_regulation
=======
>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d
            # Orthologous group size
            grp_size = len([g for g in genes if g])
            #COGs and evalues
            COGs=' | '.join([item['ID'] for item in group.COGs])
            COGes=' | '.join([str(item['eval']) for item in group.COGs])            
            #NOGs
            NOGs=' | '.join([item['ID'] for item in group.NOGs])
            NOGes=' | '.join([str(item['eval']) for item in group.NOGs])
            #PFAMs
            PFAMs=' | '.join([item['ID'] for item in group.PFAMs])
            PFAMes=' | '.join([str(item['eval']) for item in group.PFAMs])
            
            #row start
<<<<<<< HEAD
            row = [avg_p, avg_p_all, weighted_p, grp_size,group.description,\
=======
            row = [avg_p, avg_p_all, grp_size,group.description,\
>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d
                   COGs,COGes,NOGs,NOGes,PFAMs,PFAMes]

            for genome_name in genome_names:
                all_genes = group.all_genes_from_genome(genome_name)
                # Write info on the gene of the genome
                if all_genes:
                    gene = all_genes[0]
                    row.extend(['%.3f' % gene.operon.regulation_probability,
                                gene.locus_tag,
                                gene.protein_accession_number \
                                if gene.is_protein_coding_gene else ' ',
                                gene.product,
                                gene.operon.operon_id])
                else:
                    row.extend(['', '', '', '', ''])
                # Write all paralogs into a cell
                paralogs = [':'.join(('%.3f' % g.operon.regulation_probability,
                                      g.locus_tag,
                                      g.protein_accession_number \
                                      if g.is_protein_coding_gene else ' ',
                                      g.product,
                                      str(g.operon.operon_id)))
                            for g in all_genes[1:]]
                row.append('|'.join(paralogs))

            csv_rows.append(row)

        # Sort rows by average probability
<<<<<<< HEAD
        if weight_choice:
            csv_rows.sort(key=lambda row: row[2], reverse=True)
        else:
            csv_writer.writerows(csv_rows)

    """
    def weighted_orthologous_grps_to_csv(groups, phylogeny, filename):
    species = phylogeny.tree.find_elements(terminal=True, order='postorder')
    genome_names = [node.name for node in species]
    with open(filename, 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        header_row = (['average_probability',
                       'average_probability_all', 'weighted_probability'
                       'ortholog_group_size', 'description', 'COGs', 'eval', \
                       'NOGs', 'eval', 'PFAMs', 'eval'] +
                      [field for genome_name in genome_names
                       for field in ['probability (%s)' % genome_name,
                                     'locus_tag (%s)' % genome_name,
                                     'protein_id (%s)' % genome_name,
                                     'product (%s)' % genome_name,
                                     'operon id (%s)' % genome_name,
                                     'paralogs (%s)' % genome_name]])
        csv_writer.writerow(header_row)
        csv_rows = []
        for group in groups:
            genes = [group.member_from_genome(genome_name)
                     for genome_name in genome_names]
            # Average regulation probability
            avg_p = mean([g.operon.regulation_probability
                          for g in genes if g])
            # Average regulation probability (p=0 for absent genes in the grp)
            avg_p_all = mean([g.operon.regulation_probability if g else 0
                              for g in genes])
            weight = group.weighted_average_prob_regulation
            # Orthologous group size
            grp_size = len([g for g in genes if g])
            #COGs and evalues
            COGs=' | '.join([item['ID'] for item in group.COGs])
            COGes=' | '.join([str(item['eval']) for item in group.COGs])            
            #NOGs
            NOGs=' | '.join([item['ID'] for item in group.NOGs])
            NOGes=' | '.join([str(item['eval']) for item in group.NOGs])
            #PFAMs
            PFAMs=' | '.join([item['ID'] for item in group.PFAMs])
            PFAMes=' | '.join([str(item['eval']) for item in group.PFAMs])
            
            #row start
            row = [avg_p, avg_p_all, weight, grp_size,group.description,\
                   COGs,COGes,NOGs,NOGes,PFAMs,PFAMes]

            for genome_name in genome_names:
                all_genes = group.all_genes_from_genome(genome_name)
                # Write info on the gene of the genome
                if all_genes:
                    gene = all_genes[0]
                    row.extend(['%.3f' % gene.operon.regulation_probability,
                                gene.locus_tag,
                                gene.protein_accession_number \
                                if gene.is_protein_coding_gene else ' ',
                                gene.product,
                                gene.operon.operon_id])
                else:
                    row.extend(['', '', '', '', ''])
                # Write all paralogs into a cell
                paralogs = [':'.join(('%.3f' % g.operon.regulation_probability,
                                      g.locus_tag,
                                      g.protein_accession_number \
                                      if g.is_protein_coding_gene else ' ',
                                      g.product,
                                      str(g.operon.operon_id)))
                            for g in all_genes[1:]]
                row.append('|'.join(paralogs))

            csv_rows.append(row)

        # Sort rows by average probability
        csv_rows.sort(key=lambda row: row[2], reverse=True)
        csv_writer.writerows(csv_rows)
        """
=======
        csv_rows.sort(key=lambda row: row[1], reverse=True)
        csv_writer.writerows(csv_rows)

>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d

def ancestral_state_reconstruction(ortho_grps, phylo, user_input):
    """Performs ancestral state reconstruction for all orthologous groups.

    Each orthologous group consists of genes (one from each genome) and
    associated posterior probabilities of regulation. Given each orthologous
    group, this method uses BayesTraits
    (http://www.evolution.rdg.ac.uk/BayesTraits.html) to estimate the state of
    regulation on internal nodes and the root of the phylogenetic tree.

    For each group, it randomly samples trees with discrete states on genes:
    regulation or not regulation. The discretization, setting each gene as
    regulated or not, is done proportional to posterior probability of the gene
    regulation.

    For each sampled tree with discrete states, BayesTraits performs the
    ancestral state reconstruction and computes the probability of regulation
    on each internal node and the root.. The final step is to average all
    probabilities from each run of BayesTraits on each sampled tree.

    Args:
        ortho_grps ([OrthologousGroup]): the list of orthologous groups
        genomes ([Genome]): the list of target genomes.
    Returns:
    """
    my_logger.info("Ancestral state reconstruction")
    #for each orthologous group
    for ortho_grp in tqdm(ortho_grps):
        ortho_grp.ancestral_state_reconstruction(phylo, user_input)
    my_logger.info("Ancestral state reconstruction [DONE]")


def ancestral_states_to_csv(ortho_grps, phylo, filename):
    with open(filename, 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['genes', 'description'] +
                            ['%s P(%s)' % (node.name, state)
                             for node in phylo.tree.find_clades()
                             for state in ['1', '0', 'A']])
        for ortho_grp in ortho_grps:
            states = ortho_grp.regulation_states
            csv_writer.writerow(
                [', '.join(g.locus_tag for g in ortho_grp.genes),
                 ortho_grp.genes[0].product] +
                ['%.2f' % states[(node.name, state)]
                 for node in phylo.tree.find_clades()
                 for state in ['1', '0', 'A']])
