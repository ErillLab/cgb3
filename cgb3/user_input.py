import json
from .my_logger import my_logger
from cached_property import cached_property

class UserInput:
    """Class definition for UserInput.

    UserInput class gives an interface to access the run parameters provided by
    the user. It allows us to make changes in the input format (currently,
    JSON) relatively easily. It also checks the content of the input and sets
    default values for parameters not specified by the user.
    """
    def __init__(self, input_filename):
        """Constructs Input object using the input file."""
        self._input = {}
        with open(input_filename) as f:
            for k, v in list(json.load(f).items()):
                self._input[k] = v

    @property
    def genome_name_and_accessions(self):
        """Returns the list of (genome name, accession_numbers) tuples."""
        return [(g['name'], g['accession_numbers'])
                for g in self._input['genomes']]

    @property
    def genome_names(self):
        """Returns the list of genome names used."""
        return [g['name'] for g in self._input['genomes']]

    @property
    def protein_accessions(self):
        """Returns the protein accession numbers."""
        return [p['protein_accession'] for p in self._input['motifs']]

    @property
    def protein_names(self):
        """Returns the species names for reference proteins."""
        return [p['name'] for p in self._input['motifs']]

    @property
    def sites_list(self):
        """Returns the lists of binding sites."""
        return [m['sites'] for m in self._input['motifs']]

    @property
    def genomes_acc_list(self):
        """Returns the lists of reference genome accessions."""
        return [g['genome_accessions'] for g in self._input['motifs']]

    @property
    def protein_accessions_and_sites(self):
        """Zips protein accessions and binding sites."""
        return list(zip(self.protein_accessions, self.sites_list))

    @property
    def protein_names_and_genome_accessions(self):
        """Zips protein names and genome accessions for motifs."""
        return list(zip(self.protein_names, self.genomes_acc_list))

    @cached_property
    def prior_regulation_probability(self):
        """Returns the prior probability of regulation.

        It is used by the Bayesian estimator of operon regulation as the prior
        probability of an operon being regulated.

        If not provided by the user, it is estimated by predicting operons in
        the genome of the TF for which the motif is provided and dividing the
        number of sites in the motif (assuming that there is one site per
        operon) by the number of operons. For instance, if 30 sites are
        available for LexA in E. coli, then the prior for regulation is
        30/~2300 operons. See "get_prior" in "main.py" for implementation.
        """
        try:
            value = float(self._input['prior_regulation_probability'])
        except KeyError:
            value = None
        except ValueError:
            value = None
        return value

    @property
    def has_prior_probability_set(self):
        """Returns true if the prior probability of regulation is provided."""
        return 'prior_regulation_probability' in self._input

    @cached_property
    def posterior_probability_threshold_for_reporting(self):
        """Returns the threshold for regulation probabilities.

        Only the operons with a regulation probability above the threshold will
        be reported.

        The default value for the probability threshold is 0.5
        """
        try:
            value = self._input['posterior_probability_threshold_for_reporting']
            #limit range
            if value < 0.0:
                my_logger.info("WARNING: "\
                               "posterior_probability_threshold_for_reporting "\
                               "(%d) out"\
                               "of range in input file; will be reset to %d" %
                              (value, 0.0))
                value=0.0
            if value > 1.0:
                my_logger.info("WARNING: "\
                               "posterior_probability_threshold_for_reporting "\
                               "(%d) out"\
                               "of range in input file; will be reset to %d" %
                              (value, 1.0))
                value=1.0
        except KeyError:
            value = 0.5
        return value

    @cached_property
    def phylogenetic_weighting(self):
        """Returns if the phylogenetic-weighting option is on."""
        try:
            value = self._input['phylogenetic_weighting']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "phylogenetic_weighting (%s) not properly "\
                               "defined in input file; will be reset to %d" %
                              (str(value), False))
        except KeyError:
            value = False
        return value

    @cached_property
    def site_count_weighting(self):
        """Returns if the site count-weighting option is on."""
        try:
            value = self._input['site_count_weighting']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "site_count_weighting (%s) not properly defined; "\
                               "will be reset to %d" %
                              (str(value), False))
                value=False
        except KeyError:
            value = False
        return value

    @cached_property
    def operon_prediction_probability_threshold(self):
        """Returns the threshold for posterior probability of regulation used
        for operon prediction. The default value is 0.5.
        """
        try:
            value = self._input['operon_prediction_probability_threshold']
            #limit range
            if value < 0.0:
                my_logger.info("WARNING: "\
                               "operon_prediction_probability_threshold (%d) "\
                               "out of range in input file; will be reset to %d"
                               % (value, 0.0))
                value=0.0
            if value > 1.0:
                my_logger.info("WARNING: "\
                               "operon_prediction_probability_threshold (%d) "\
                               "out of range in input file; will be reset to %d"
                               % (value, 1.0))
                value=1.0
        except KeyError:
            value = 0.5
        return value

    @cached_property
    def operon_prediction_distance_tuning_parameter(self):
        """Returns the sigma value used to calibrate distance threshold used
        for operon prediction. The value should be positive. The default value
        is 1.0. Values greater than 1.0 results in higher intergenic distance
        threshold, therefore, less number of operons.
        """
        try:
            value = self._input['operon_prediction_distance_tuning_parameter']
            #limit range
            if value < 0.5:
                my_logger.info("WARNING: "\
                               "operon_prediction_distance_tuning_parameter "\
                               "(%d) out "\
                               "of range in input file; will be reset to %d" %
                              (value, 0.5))
                value=0.5
            if value > 5.0:
                my_logger.info("WARNING: "\
                               "operon_prediction_distance_tuning_parameter "\
                               "(%d) out "\
                               "of range in input file; will be reset to %d" %
                              (value, 5.0))
                value=5.0
        except KeyError:
            value = 1.0
        return value

    @cached_property
    def ancestral_state_reconstruction(self):

        """Returns True/False which specifies whether ancestral state reconstruction
        will be performed or not.
        """
        try:
            value = self._input['ancestral_state_reconstruction']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "ancestral_state_reconstruction (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), False))
                value=False
        except:
            value = False
        return value

    @cached_property
    def bootstrap_replicates(self):
        """Returns the number of bootstrap replicates to be performed for
           ancestral state reconstruction. Defaults to 100.
        """
        try:
            value = self._input['bootstrap_replicates']
            #limit range
            if value < 0:
                my_logger.info("WARNING: "\
                               "mixing rate alpha (%d) out"\
                               "of range in input file; will be reset to %d" %
                              (value, 0))
                value=0.0
            if value > 10000:
                my_logger.info("WARNING: "\
                               "mixing rate alpha (%d) out"\
                               "of range in input file; will be reset to %d" %
                              (value, 10000))
                value=10000
        except:
            value = 100
        return value

    @cached_property
    def alpha(self):
        """Returns the alpha mixing rate parameter for the mixture distribution
           of scores in regulated promoters. Defaults to 0.03 [1/300].
        """
        try:
            value = self._input['alpha']
            #limit range
            if value < 0.0:
                my_logger.info("WARNING: "\
                               "mixing rate alpha (%d) out"\
                               "of range in input file; will be reset to %d" %
                              (value, 0.0))
                value=0.0
            if value > 1.0:
                my_logger.info("WARNING: "\
                               "mixing rate alpha (%d) out"\
                               "of range in input file; will be reset to %d" %
                              (value, 1.0))
                value=1.0
        except:
            value = 0.03
        return value

    @cached_property
    def promoter_up_distance(self):
        """Returns the maximum distance upstream of a gene predicted TLS to
           search for putative TF-binding sites. Defaults to 300 bp.
        """
        try:
            value = self._input['promoter_up_distance']
            #limit range
            if value < 0:
                my_logger.info("WARNING: "\
                               "promoter_up_distance (%d) out"\
                               "of range in input file; will be reset to %d" %
                              (value, 0))
                value=0
            if value > 1000:
                my_logger.info("WARNING: "\
                               "promoter_up_distance (%d) out"\
                               "of range in input file; will be reset to %d" %
                              (value, 1000))
                value=1000
        except:
            value = 300
        return value

    @cached_property
    def promoter_dw_distance(self):
        """Returns the maximum distance downstream of a gene predicted TLS to
           search for putative TF-binding sites. Defaults to 50 bp.
        """
        try:
            value = self._input['promoter_dw_distance']
            #limit range
            if value < 0:
                my_logger.info("WARNING: "\
                               "promoter_dw_distance (%d) out"\
                               "of range in input file; will be reset to %d" %
                              (value, 0))
                value=0
            if value > 1000:
                my_logger.info("WARNING: "\
                               "promoter_dw_distance (%d) out"\
                               "of range in input file; will be reset to %d" %
                              (value, 1000))
                value=1000
        except:
            value = 50
        return value

    @cached_property
    def heatmap_plot(self):

        """Returns True/False which specifies whether a heatmap plot will be
           generated or not.
        """
        try:
            value = self._input['heatmap_plot']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "heatmap_plot (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), True))
                value=True
        except:
            value = True
        return value

    @cached_property
    def motif_plot(self):

        """Returns True/False which specifies whether a motif plot will be
           generated or not.
        """
        try:
            value = self._input['motif_plot']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "motif_plot (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), True))
                value=True
        except:
            value = True
        return value

    @cached_property
    def gene_regulation_plot(self):

        """Returns True/False which specifies whether a gene regulation plot
           will be generated or not.
        """
        try:
            value = self._input['gene_regulation_plot']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "gene_regulation_plot (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), False))
                value=False
        except:
            value = False
        return value

    @cached_property
    def taxon_regulation_plot(self):

        """Returns True/False which specifies whether a taxon regulation plot
           will be generated or not.
        """
        try:
            value = self._input['taxon_regulation_plot']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "taxon_regulation_plot (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), False))
                value=False
        except:
            value = False
        return value

    @cached_property
    def network_size_plot(self):

        """Returns True/False which specifies whether a network size plot will
           be generated or not.
        """
        try:
            value = self._input['network_size_plot']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "network_size_plot (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), False))
                value=False
        except:
            value = False
        return value

    @cached_property
    def site_printout(self):

        """Returns True/False which specifies whether a a printout of predicted
           sites will be generated or not.
        """
        try:
            value = self._input['site_printout']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "site_printout (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), True))
                value=True
        except:
            value = True
        return value
<<<<<<< HEAD
    
    @cached_property
    def colorblind_compatibility(self):
        "Returns the user's preference for colorblind compatibility in plots"    
        try:
            value = self._input['colorblind_compatibility']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "colorblind_compatibility (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), False))
                value=False
        except:
            value = False
        return value
    
    @cached_property
    def min_size_orthologs(self):
        "Returns the user's number of orthologous groups to be considered"
        try:
            value = self._input['min_size_orthologs']
            #limit range
            if value <= 0:
                my_logger.info("WARNING: "\
                               "min_size_orthologs (%d) out"\
                               "of range in input file; will be reset to %d" %
                              (value, 2))
                value=2
        except:
            value = 2
        return value
            
    @cached_property
    def weighted_average_sorting(self):
        "Returns the user's preference for weighted sorting in plots"    
        try:
            value = self._input['weighted_average_sorting']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "weighted_average_sorting (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), True))
                value=True
        except:
            value = True
        return value
    
    @cached_property
    def use_prior_for_absence(self):
        "Returns the user's prefernce for using the prior probability of regulation\
            for weighted average sorting"
        try:
            value = self._input['use_prior_for_absence']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "use_prior_for_absence (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), True))
                value=True
        except:
            value = True
        return value
=======
>>>>>>> 26423e6e9d1842cc737bcab2b17a6f0e83d1194d

    @cached_property
    def entrez_email(self):
        """Returns the email address used for entrez.

        It is used by the entrez library to fetch genome records, and NCBI to
        warn users of API over-use.

        If not provided by the user, the program halts.
        """
        try:
            value = self._input['entrez_email']
        except KeyError:
            value = None
        except ValueError:
            value = None
        return value

    @cached_property
    def entrez_apikey(self):
        """Returns the API key address used for entrez.

        It is used by the enable more requests per second on NCBI's API

        If not provided by the user, it's set by default to none.
        """
        try:
            value = self._input['entrez_apikey']
        except:
            value = None
        return value

    @cached_property
    def sleep(self):
        """Returns the time set for sleep after NCBI queries
           in case  "HTTP Error 429: Too Many Requests" are received.
           Defaults to 0 s.
        """
        try:
            value = self._input['sleep']
            #limit range
            if value < 0:
                my_logger.info("WARNING: "\
                               "Sleep time (%d) out"\
                               "of range in input file; will be reset to %d" %
                              (value, 0))
                value=0
            if value > 1000:
                my_logger.info("WARNING: "\
                               "Sleep time (%d) out"\
                               "of range in input file; will be reset to %d" %
                              (value, 10))
                value=1000
        except:
            value = 0
        return value


    @cached_property
    def TF_eval(self):
        """Returns the e-value used for identifying TF instances in target
           genomes.
           Defaults to 10**-3.
        """
        try:
            value = self._input['TF_eval']
            #limit range
            if value < 0:
                my_logger.info("WARNING: "\
                               "E-value (%d) must be positive;"\
                               "will be reset to %d" %
                              (value, 0.001))
                value=0.001
            if value > 1:
                my_logger.info("WARNING: "\
                               "e-value (%d) too large for meaningful results"\
                               "; will be reset to %d" %
                              (value, 0.001))
                value=0.001
        except:
            value = 0.001
        return value

    @cached_property
    def homolog_eval(self):
        """Returns the e-value used for identifying homologs among target
           genomes.
           Defaults to 10**-3.
        """
        try:
            value = self._input['homolog_eval']
            #limit range
            if value < 0:
                my_logger.info("WARNING: "\
                               "E-value (%d) must be positive;"\
                               "will be reset to %d" %
                              (value, 0.001))
                value=0.001
            if value > 1:
                my_logger.info("WARNING: "\
                               "e-value (%d) too large for meaningful results"\
                               "; will be reset to %d" %
                              (value, 0.001))
                value=0.001
        except:
            value = 0.001
        return value

    @cached_property
    def hmmer_eval(self):
        """Returns the e-value used for identifying homologs in the eggNOG or
           PFAM databases.
           Defaults to 10**-5.
        """
        try:
            value = self._input['hmmer_eval']
            #limit range
            if value < 0:
                my_logger.info("WARNING: "\
                               "E-value (%d) must be positive;"\
                               "will be reset to %d" %
                              (value, 0.00001))
                value=0.00001
            if value > 1:
                my_logger.info("WARNING: "\
                               "e-value (%d) too large for meaningful results"\
                               "; will be reset to %d" %
                              (value, 0.00001))
                value=0.00001
        except:
            value = 0.00001
        return value

    @cached_property
    def COG_search(self):
        """Returns True/False which specifies whether a HMMER COG search
           will be run or not.
        """
        try:
            value = self._input['COG_search']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "COG_search (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), False))
                value=False
        except:
            value = False
        return value
    
    @cached_property
    def NOG_search(self):
        """Returns True/False which specifies whether a HMMER eggNOG search
           will be run or not.
        """
        try:
            value = self._input['NOG_search']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "NOG_search (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), False))
                value=False
        except:
            value = False
        return value

    @cached_property
    def PFAM_search(self):
        """Returns True/False which specifies whether a HMMER PFAM search
           will be run or not.
        """
        try:
            value = self._input['PFAM_search']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "PFAM_search (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), False))
                value=False
        except:
            value = False
        return value

    @cached_property
    def COG_dbname(self):
        """Returns the path (full or relative, including file name) to
           the COG database.

           If not provided by the user, it's set by default to none.
        """
        try:
            value = self._input['COG_dbname']
        except:
            value = None
        return value
    
    @cached_property
    def eggNOG_dbname(self):
        """Returns the path (full or relative, including file name) to
           the eggNOG database.

           If not provided by the user, it's set by default to none.
        """
        try:
            value = self._input['eggNOG_dbname']
        except:
            value = None
        return value

    @cached_property
    def PFAM_dbname(self):
        """Returns the path (full or relative, including file name) to
           the eggNOG database.

           If not provided by the user, it's set by default to none.
        """
        try:
            value = self._input['PFAM_dbname']
        except:
            value = None
        return value


    @cached_property
    def OGejump(self):
        """Returns the exponent difference in e-value used filter COG/NOG hits.
           Defaults to 5 (e.g. allows jump from 1e-15 to 1e-10).
        """
        try:
            value = self._input['OGejump']
            #limit range
            if value < 0:
                my_logger.info("WARNING: "\
                               "E-value jump (%d) must be positive;"\
                               "will be reset to %d" %
                              (value, 5))
                value=5
        except:
            value = 5
        return value

    @cached_property
    def maxCOG(self):
        """Returns the maximum number of COGs to be reported for any given
           ortholog group.
           Defaults to 1.
        """
        try:
            value = self._input['maxCOG']
            #limit range
            if value < 0:
                my_logger.info("WARNING: "\
                               "Max number of COGs (%d) must be positive;"\
                               "will be reset to %d" %
                              (value, 1))
                value=1
        except:
            value = 1
        return value

    @cached_property
    def maxNOG(self):
        """Returns the maximum number of NOGs to be reported for any given
           ortholog group.
           Defaults to 1.
        """
        try:
            value = self._input['maxNOG']
            #limit range
            if value < 0:
                my_logger.info("WARNING: "\
                               "Max number of NOGs (%d) must be positive;"\
                               "will be reset to %d" %
                              (value, 1))
                value=1
        except:
            value = 1
        return value


    @cached_property
    def maxPFAM(self):
        """Returns the maximum number of PFAMs to be reported for any given
           ortholog group.
           Defaults to 1.
        """
        try:
            value = self._input['maxPFAM']
            #limit range
            if value < 0:
                my_logger.info("WARNING: "\
                               "Max number of PFAMs (%d) must be positive;"\
                               "will be reset to %d" %
                              (value, 1))
                value=1
        except:
            value = 1
        return value

    @cached_property
    def use_up_dist_site_scan(self):
        """Returns True/False which specifies whether the promoter up distance
           will be used during the site scan.
           If the promoter up distance is not used, the site scan will run up
           to the preceding gene end/start. If it is used, it will run up to
           that distance, regardless of any annotated features.
        """
        try:
            value = self._input['use_up_dist_site_scan']
            #test value
            if not(isinstance(value, bool)):
                my_logger.info("WARNING: "\
                               "use_up_dist_site_scan (%s) not "\
                               "properly defined in input file; "\
                               "will be reset to %d" %
                              (str(value), False))
                value=False
        except:
            value = False
        return value