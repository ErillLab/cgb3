
from Bio.motifs import Motif, Instances
from Bio import motifs
from Bio.Seq import Seq
from .lasagna import LASAGNA_alignment
from tqdm import tqdm


def get_alignment_offset(motif, other):
    '''
    Determines the optimal alignment of two motifs by maximizing the information content (ic) in the aligned regions.
    Parameters
    ----------
    motif, other: Motif objects
        The two motifs of interest.
    
    Returns
    -------
    offsets: int
        The offset that results in the maxium ic in the alignment. 
    '''
    max_ic = float('-inf')
    for offset in range(-len(motif) + 1, len(other)):
        if offset < 0:
            ic = ic_at(motif, other, -offset)
        else:
            ic = ic_at(other, motif, offset)

        
        if ic > max_ic:
            max_ic = ic
            max_offset = offset
            
    return max_offset

def ic_at(motif, other, offset):
    '''
    Caculates the information content, ic, for a specific alignment. The approach makes a temporary motif object containing the overlapping sequences in the alignemnt and taking the average of the pssm.
    Parameters
    ----------
    motif, other: Motif objects
        The motifs of interest
    offset: int
        The offset value that results in the alignment of interest. 
    '''

    #Pull the sequences containined in the aligned region of the motifs from each of the motif instances. 
    alignment_len = min(len(motif)-offset, len(other))
    motif_seqs = [site[offset:alignment_len+offset] for site in motif.instances]
    other_seqs = [site[:alignment_len] for site in other.instances]

    # Create the motif and compute the IC
    amotif = Motif(instances=Instances(motif_seqs+other_seqs))
    amotif.pseudocounts = dict(A=0.25, C=0.25, G=0.25, T=0.25)

    #print('Motif Seqs: ' , motif_seqs)
    #print('Other Seqs: ' , other_seqs)
    #print('Offset ', offset)
    #print('IC: ' , amotif.pssm.mean(), '\n\n')

    return amotif.pssm.mean()

def createMotifAligned(motif, other):
    """
    Given 2 motifs with determinated aligned sites, creates a new motif that is the aligned combination of the sites of both initial motifs
    """

    instances1=[]
    instances2=[]
    
    for seq in motif:
        instances1.append(Seq(seq))
        
    for seq in other:
        instances2.append(Seq(seq))
        
    motif=motifs.create(instances1)
    other=motifs.create(instances2)
     
    offset=get_alignment_offset(motif, other)
    non_motif=[]

    if offset<0:
        offset=-offset
    else:
        tmp=motif
        motif=other
        other=tmp
    alignment_len = min(len(motif)-offset, len(other))
    motif_seqs = [site[offset:alignment_len+offset] for site in motif.instances]
    other_seqs = [site[:alignment_len] for site in other.instances]
    amotif = Motif(instances=Instances(motif_seqs+other_seqs))
    for seq in amotif.instances:
        non_motif.append(str(seq))

    return non_motif
    
def alignSites(sites):
    """
    Given the sites of all motifs of the input, it aligns them
    """
    #Count how many sites has every motif
    quantity_of_sequences=[]
    for seq in sites:
        quantity_of_sequences.append(len(seq))
    #Loop inicialization  
    it = iter(sites[0])
    the_len = len(next(it))
    if not all(len(l) == the_len for l in it):
        sites[0]=LASAGNA_alignment(sites[0])
    #loop
    for i in tqdm(range(len(sites)-1)):
        it = iter(sites[i+1])
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            sites[i+1]=LASAGNA_alignment(sites[i+1])
        sites[i+1]=createMotifAligned(sites[i],sites[i+1])
    #now that we have all the sites of all motifs aligned in sites[-1] we have to divide them to storage them in the respective motif
    i=0
    old_quantity=0
    for quantity in quantity_of_sequences:
        sites[i]=sites[-1][old_quantity:old_quantity+quantity]
        i=i+1
        old_quantity=old_quantity+quantity

    return sites
