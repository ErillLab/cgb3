"""Module containing functions using NCBI Entrez utility.

See NCBI Entrez page (http://www.ncbi.nlm.nih.gov/books/NBK3837/) and Biopython
(http://biopython.org/DIST/docs/tutorial/Tutorial.html) tutorial for more
information.
"""

import os
from tqdm import tqdm
from Bio import Entrez
from Bio import SeqIO
from ete3 import NCBITaxa
from .misc import directory
from .my_logger import my_logger
import sys
import time
from collections import Counter
import re


# The directory used to save NCBI records for later use.
ENTREZ_DIRECTORY = directory('entrez_cache')


def set_entrez_email(email_address):
    Entrez.email = email_address

def set_entrez_apikey(api_key):
    Entrez.api_key = api_key

def set_entrez_delay(delay):
    global sleep_time
    sleep_time = delay
    
def set_entrez_retry_number(retries):
    global retry_number
    retry_number = retries

def get_genome_record(accession):
    """Gets the genome record from NCBI RefSeq."""
    i=0
    genbank_file = os.path.join(ENTREZ_DIRECTORY, accession+'.gb')
    if not os.path.isfile(genbank_file):
        while i<retry_number:
            try:
                # Download and save Genbank record
                my_logger.info("Downloading %s" % accession)
                handle = Entrez.efetch(db='nuccore', id=accession,
                                       rettype='gbwithparts', retmode='text')
                record = handle.read()
                # add further delay if NCBI traffic is high and 
                # "HTTP Error 429: Too Many Requests" is received
                with open(genbank_file, 'w') as f:
                    f.write(record)
                break
            except:
                time.sleep(sleep_time+i)
                i=i+1
                continue
    if i==retry_number and i!=0:
        my_logger.warning("Retrieving accessions from NCBI databases is not possible, cannot save on cache")
        sys.exit() #without the genomes the execution has to stop
    handle = open(genbank_file)
    return handle.read()
        


#takes an accession number for the protein, gets the record from NCBI
#(unless it is already stored locally) and saves it to file locally
#returns the object in the local file
def get_protein_record(accession):
    """Fetches the protein record from NCBI Protein database."""
    #global sleep_time
    protein_file = os.path.join(ENTREZ_DIRECTORY, accession+'.gb')
    #if file not locally available, fetch and save locally (in ENTREZ_DIRECTORY cache)
    if not os.path.isfile(protein_file):
        # Download and save file
        try:
            handle = Entrez.efetch(db='protein', id=accession,
                                   rettype='gb', retmode='text')
            record = handle.read()
        except Exception as e:
            my_logger.info("Unable to access %s" % e.geturl())
        # add further delay if NCBI traffic is high and 
        # "HTTP Error 429: Too Many Requests" is received
        time.sleep(sleep_time)
        with open(protein_file, 'w') as f:
            f.write(record)

    #read file and return object
    handle = open(protein_file)
    return handle.read()


def query_assembly(assembly_ID):
    """
    Given a assembly ID returns the list of nucleotide ID linked to that assembly ID.
    First searchs the list of nucleotide ID with linkname="assembly_nuccore_refseq" and if the elink returns nothing tries again without linkname
    """
    nucleotide_ID_list=[]
    nucleotide_ID_list_object=Entrez.read(Entrez.elink(db="nuccore",dbfrom="assembly", id=assembly_ID, linkname="assembly_nuccore_refseq", retmax=9999))  
    if  nucleotide_ID_list_object[0]["LinkSetDb"]==[]:
        nucleotide_ID_list_object=Entrez.read(Entrez.elink(db="nuccore",dbfrom="assembly", id=assembly_ID, retmax=9999))
        
    for seqId in nucleotide_ID_list_object[0]["LinkSetDb"][0]["Link"]:
        nucleotide_ID_list.append(seqId["Id"])
        
    return nucleotide_ID_list

def nucleotide_query_species(species_taxID, minN50):
    """
    Given a species TaxID obtains the best assembly for that specie based on 4 scores:
        -100 --> representative genome
        -10  --> RefSeq complete 
        -1   --> GenBank complete
        -0,contigN50 --> Incomplete assembly
        
    Once the best assembly is retrieved it calls query_assembly() with the assembly ID retrieved as parameter and then returns the result of that call and the score of the assembly 
    The second parameters is minN50 a parameter from the json input
    """
    score_n50=float("-inf")
    my_logger.info("Analyzing species with taxonomic ID = %s"%species_taxID)
   
    #If score is 100 it is representative and finishes 
    representative=False
    genome_ID_object=Entrez.esearch(db="genome", term="txid%s[orgn] AND complete[Status]"%species_taxID)
    genome_ID=Entrez.read(genome_ID_object)['IdList']
    
    if genome_ID!=[]:
        representative=True
        assembly_ID_object=Entrez.read(Entrez.efetch(db="Genome",id=genome_ID,rettype="docsum", retmode="xml"))	
        assembly_ID=assembly_ID_object[0]['AssemblyID']
        if assembly_ID!="0":
            score=100
        else:
            representative=False   
    #The score is not 100 then we it will search for the best score available   
    if not representative:
        assembly_ID_object=Entrez.read(Entrez.esearch(db="assembly", term="txid%s[orgn] AND latest_refseq[filter] AND complete_genome[filter]"%species_taxID))
        assembly_ID=assembly_ID_object["IdList"]
        score=10
        if assembly_ID==[]:
            assembly_ID_object=Entrez.read(Entrez.esearch(db="assembly", term="txid%s[orgn] AND complete_genome[filter]"%species_taxID))
            assembly_ID=assembly_ID_object["IdList"]
            score=1
            if assembly_ID==[]:
                assembly_ID_object=Entrez.read(Entrez.esearch(db="assembly", term="txid%s[orgn]"%species_taxID,retmax=9999))
                assembly_ID_list=assembly_ID_object["IdList"]
                
                for assembly in assembly_ID_list:
                    score=getContigN50(assembly)
                    if score>score_n50 and score>minN50:
                        score_n50=score
                        assembly_ID=assembly 
            
            score=score_n50/pow(10,(len(str(int(score)))))
            
        if assembly_ID==[]:
            return None, None
        
    return query_assembly(assembly_ID), score



def nucleotide_query(parameters):
    """
    The function works in two ways depending if the level of search is species or not.
    If is the level is species it calls the function nucleotide_query_species() for every descendant of the taxID
    If the level is not species it calls getChildren() with parameter level = "species" to get the descendant species for every TaxID in the "descendats" variable and filters the best species as the representative for that TaxID
    """
    output_list=[]
    taxID=parameters[0]['target_taxID']
    level=parameters[1]['target_rank']
    letNoRank=parameters[2]['allow_no_rank']
    letUnClassified=parameters[3]['allow_unclassified']
    letEnviormental=parameters[4]['allow_environmental']
    minN50=parameters[5]['min_N50_contig']

    descendants,names_descendants=getChildren(taxID, level, letUnClassified, letNoRank, letEnviormental)
    if level =="species":
        for species in tqdm(descendants):
            nucleotide_ID_list,score=nucleotide_query_species(species, minN50)
            if score!=None:
                my_logger.info("Species with taxId = %s has an assembly", species)
                output_list.append(nucleotide_ID_list)
            else:
                my_logger.info("SPECIES WITH TAXID = %s HAS NOT AN ASSEMBLY", species)
        return output_list,names_descendants,level
    else:
        return_names_descendants=False
        for descendant in tqdm(descendants):
            score_old=0
            species_list=getChildren(descendant, "species", letUnClassified, letNoRank, letEnviormental, return_names_descendants)
            for species in species_list:
                nucleotide_ID_list,score=nucleotide_query_species(species, minN50)
                if score!=None:
                    if score==100:
                        best=nucleotide_ID_list
                        score_old=100
                        break
                    if score > score_old:
                        score_old=score
                        best=nucleotide_ID_list
            if score_old != 0:
                my_logger.info("The representative species for the taxonomic group has been found")
                output_list.append(best)
            else:
                my_logger.info("THERE ARE NOT SUITABLE REPRESENTATIVE ASSEMBLIES FOR THE TAXONOMIC GROUP")
           
        return output_list,names_descendants,level
    
def getChildren(TaxId,level="species", letUnClassified=True, letNoRank=True, letEnviormental=True,returnNames=True):
    """
    Using ete3 library it returns the descendants of a TaxID filtering by level in the phylogenetic tree
    Also depending on the values of the last 3 parameters of the function it will also filter "no rank", "unclassified" or "environmental samples" descendats
    """
    ncbi = NCBITaxa() 
    descendants=ncbi.get_descendant_taxa(TaxId, intermediate_nodes=True)
    #names_descendant=ncbi.translate_to_names(descendants)
    ranks_of_descendants=ncbi.get_rank(descendants)
    if checkLevel(level, ranks_of_descendants):
        records=[]
        for rank in ranks_of_descendants:
            if ranks_of_descendants[rank]==level:
                number_of_no_ranks_in_lineage=Counter(list(ncbi.get_rank(ncbi.get_lineage(rank)).values()))['no rank']
                if letNoRank!=True and number_of_no_ranks_in_lineage>2: 
                    pass
                elif letUnClassified!=True and "unclassified" in str(list(name for name in ncbi.get_taxid_translator(ncbi.get_lineage(rank)).values() if "unclassified" in name)):
                    pass
                elif letEnviormental!=True and "environmental samples" in str(list(name for name in ncbi.get_taxid_translator(ncbi.get_lineage(rank)).values() if "environmental samples" in name)):
                    pass
                else:
                    records.append(rank)
        if returnNames:
            names_descendant=ncbi.translate_to_names(records)              
            return records,names_descendant
        else:
            return records

    else:
        my_logger.info("Can not search for %s inside %s",level,TaxId)
            
    
    
def getContigN50(assemblyId):
    """
    Using efetch to retrieve the docsum object of a determinated assemblyId It returns the "contig_n50" value
    """
    data_object=Entrez.read(Entrez.efetch(db="assembly", id=assemblyId,rettype="docsum", retmode="xml"))
    data=data_object["DocumentSummarySet"]["DocumentSummary"][0]["Meta"]
    for child in data.split("<"):
        if "contig_n50" in child:
            left, right = child.split(">")
            return float(right)
        
def checkLevel(level, rank):
    """
    Function used to check if the user is searching correctly. That is to say, for example,that it is not searching for orders inside a genus taxId
    """
    if level not in list(rank.values()):
        return False
    else:
        return True
            
def saveGenBankFileCache(gb_obj,accession):
    """
    Saves gb_obj (genbank object) in the cache with SeqIO.write()
    """

    genbank_file = os.path.join(ENTREZ_DIRECTORY, accession+'.gb')
    if not os.path.isfile(genbank_file):
        my_logger.info("Saving %s in cache" %accession)

        with open(genbank_file, 'w') as f:
            SeqIO.write(gb_obj,f,'genbank')

def update_genomes_dictionary(uid_list,group,level):
    """
    Returns the dictionay of user_input updated with the new genomes. It also saves them in cache
    """
    acc_list = convert_uids_to_accessions(uid_list)
    final_acc_list = download_cache_accessions(acc_list)
    sp_name=get_name_species(acc_list[0])
    return({'name': sp_name, str(level):group, 'accession_numbers': final_acc_list})

def get_name_species(accession):
    """
    It fetchs from cache the name of the species given a accession number
    """
    genbank_file = os.path.join(ENTREZ_DIRECTORY, accession+'.gb')
    nucleotide_ID_object=SeqIO.read(genbank_file, "gb")
    return re.sub('[\W]', '_', nucleotide_ID_object.annotations['organism']) #[\W] to remove ALL non-alphanumerics and non-underscore

def convert_uids_to_accessions(uid_list):
    """
    With Entrez.efecth returns the accession numbers of the uid list
    """
    accession_list=[]
    for uid in uid_list:
        with Entrez.efetch(db="nucleotide", rettype="acc", retmode="text", id=uid) as handle:
            accession_list.append(handle.read().rstrip("\n"))
    return accession_list

def processMaterRecordAlt(accession_list,master_record):
    """
    Given a nuccore master record It appends all the range of accession numbers of the contigs submitted as part of that nuccore master record into the "accession_list" variable
    """
    #Start. Here it creates the accession number range
    first=master_record.annotations['wgs'][0] #first accession number
    last=master_record.annotations['wgs'][1] #last accession number
    for i in range(len(first)):#loop to divide the invariable part from the variable part of the accession number
        if first[i] != last[i]:
            index=i
    index=index-1
    prefix_accesion=first[:index]#invariable part of the accession number (prefix)
    range_accesions=range(int(first[index:]),int(last[index:])+1)
    #end
    
    for i in range_accesions: #for every number in the range it creates a accession number string
        accession_tmp=str(prefix_accesion)+str(i).zfill(len(str(last[index:])))+"."+str(master_record.annotations['sequence_version'])
        if accession_tmp not in accession_list:
            accession_list.append(accession_tmp)
    return accession_list   

def download_cache_accessions(accession_list):
    """
    Given a list of accession numbers It downloads and saves in cache the genbank files. It only downloads and saves the ones that are not in cache already
    It checks if the accession number is a master record and It downloads and saves all the accession numbers from the master record. 
    It also updates the accession list every execution. That is done to ensure that the cache has updated and correct information 
    """
    new_accession_list = []
    
    for accession in accession_list:
        accession_cache = os.path.join(ENTREZ_DIRECTORY, accession+'.gb')
        if os.path.isfile(accession_cache):#if already in cache
                new_accession_list.append(accession)
        else:#if not in cache
            i=0
            while i<retry_number:
                try:
                    with Entrez.efetch(db="nucleotide", rettype="gbwithparts", retmode="text", id=accession) as handle:
                        nucleotide_ID_object=SeqIO.read(handle, "gb")
                        if not "wgs" in nucleotide_ID_object.annotations:#if no master record
                            saveGenBankFileCache(nucleotide_ID_object,accession)
                            new_accession_list.append(accession)
                        else:#if master record
                            tmp_accession_list=download_cache_accessions(processMaterRecordAlt(new_accession_list,nucleotide_ID_object))#recursivity with the accession numbers from the master record
                            for accession in tmp_accession_list:
                                if accession not in new_accession_list:
                                    new_accession_list.append(accession)
                    break
                except:
                    time.sleep(sleep_time+i)
                    i=i+1
                    continue
            if i==retry_number and i!=0:
                my_logger.warning("Retrieving accessions from NCBI databases is not possible, cannot save on cache")
                sys.exit() #without the genomes the execution has to stop
    return new_accession_list
