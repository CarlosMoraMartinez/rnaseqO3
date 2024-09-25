import argparse
import os
import glob
import pysam
import gffutils
from typing import List, Set, Dict, Tuple
import pandas as pd
import numpy as np
import json

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class Logger:
    def __init__(self, logname: str = None) -> None:
        self.logname: str = logname
    def set_file(self, logname: str) -> None:
        self.logname = logname
    def log2file(self, msg: str):
        if self.logname is not None:
            with open(self.logname, "a") as fo:
                fo.write(str(msg) + '\n')
    def log(self, msg: str, col: bcolors = bcolors.BOLD) -> None:
        print(col + str(msg) + bcolors.ENDC)
        self.log2file(msg)

logger = Logger()

def run_command(cmd: str):
    logger.log(f"----Running command:\n{cmd}", bcolors.HEADER)
    try:
        exitcode: int = os.system(cmd)
        if exitcode == 0:
            logger.log(f"---Success: {cmd}", bcolors.OKGREEN)
        else:
            logger.log(f"---FAIL {exitcode}: {cmd}", bcolors.FAIL)
    except Exception as e:
        logger.log(e, bcolors.FAIL)
        logger.log(f"---FAIL: {cmd}", bcolors.FAIL)
 
######################

def get_gtf_db(fname, only_genes = False, rewrite_db = False):

    dbname = f"{fname}.genes.db" if only_genes else f"{fname}.trans.db"
    if not os.path.isfile(dbname) or rewrite_db:
        logger.log(f"Creating feature database from {fname}", bcolors.OKBLUE)
        db = gffutils.create_db(fname, dbfn=dbname, disable_infer_transcripts=only_genes)
        logger.log(f"Database created from {fname}", bcolors.OKGREEN)
    else:
        logger.log(f"Database already exists: {fname}", bcolors.WARNING)
    db = gffutils.FeatureDB(dbname, keep_order=True)
    return db


def analyze(bam, db, genelist, add_chr=True):
    results = []
    samfile = pysam.AlignmentFile(bam, "rb")
    genes_found = []
    for gg in db.features_of_type('gene'):
        nn = gg.attributes.get('gene_name', gg.attributes.get('gene_id'))
        if genelist and not any([n in genelist for n in nn]):
            continue
        logger.log(f"Getting reads from gene {nn}", bcolors.OKBLUE)
        genes_found.extend(nn)
        chrom: str = gg.chrom if not add_chr else f"chr{gg.chrom}"

        #Get transcripts in reference from which fastq wewre generated to assign reads to gene
        #query: str= (
        #f"SELECT id, attributes FROM features "
        #f"WHERE featuretype='transcript' AND attributes like '%\"{nn[0]}\"%';"
        #)
        #logger.log(query, bcolors.OKCYAN)
        #cur = db.execute(query)
        #tr_samples_res = cur.fetchall()
        #if not tr_samples_res:
        #    logger.log(f"Gene not found in samples gtf", bcolors.FAIL)
        #    continue
        #transcript_ids = [tt['attributes'][0] for tt in tr_samples_res]
        #print(transcript_ids)
        #logger.log(f"Transcripts in original fastq: {', '.join(transcript_ids)}", bcolors.OKCYAN)
        
        #transcript_ids = [tn for transcr in db.children(gg, featuretype='transcript') for tn in transcr.attributes['transcript_id'] ]
        gene_results = {'gene_name':'|'.join(nn), 
                        'gene_id':'|'.join(gg.attributes.get('gene_id')),
                        'chrom':chrom,
                        'gene_start':gg.start, 
                        'gene_end':gg.end, 
                        'primary_right':0,
                        'primary_wrong':0,
                        'secondary_right':0,
                        'secondary_wrong':0,
                        'confounding_genes':{}
                        }
        
        for read in samfile.fetch(chrom, gg.start, gg.end): 
            kpart1 = 'secondary' if read.is_secondary else 'primary'
            kpart2 = 'right' if any([n == read.query_name.split('.')[0] for n in nn]) else 'wrong'
            gene_results[f"{kpart1}_{kpart2}"] += 1
            if kpart2 == 'wrong':
                original_transcript = f"{read.query_name.split('.')[0]}"
                #original_genes = [n for og in db.parents(original_transcript) for n in og.attributes.get('gene_name', og.attributes.get('gene_id'))]
                wrong_reads = gene_results['confounding_genes'].setdefault(original_transcript, 0)
                gene_results['confounding_genes'][original_transcript] = wrong_reads + 1
        print(gene_results)
        results.append(gene_results)
    if genelist:
        logger.log(f"Genes found: {', '.join(genes_found)}", bcolors.OKGREEN)
        logger.log(f"Genes not found: {', '.join([i for i in genelist if i not in genes_found])}", bcolors.FAIL)
    return results


parser: argparse.ArgumentParser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
                                                          description='Count read origin in bam calculated from simulated reads.')
parser.add_argument('-o', '--outdir', type = str, help='Output directory', default="")
parser.add_argument('-u', '--outfname', type = str, help='Output file name', default="")
parser.add_argument('-b', '--bam', type = str, help='Input bam or directory with bam files.')
parser.add_argument('-l', '--glob', type = str, help='Glob to find bam files in input dir.', default="*.bam")
parser.add_argument('-g', '--genelist', type = str, help='List of genes of interest.', default="")
parser.add_argument('-a', '--annot', type = str, help = 'Annotation file of reference')
parser.add_argument('-x', '--only_genes',  action="store_true" , help = 'Argument disable_infer_transcripts for gffutils.create_db')
parser.add_argument('-r', '--rewrite_db',  action="store_true", help = 'Rewrite database')
parser.add_argument('-c', '--add_chr',  action="store_true", help = 'Add \'chr\' to GTF chromosome when searching bam.')


# python assess_mapping.py -o test -b /home/carmoma/Documents/RNAseq_Rebeca2024/results_simulations/sim1_error20/mg02_star_align/allcDNA_Aligned.sortedByCoord.out.bam -a /home/carmoma/Documents/RNAseq_Rebeca2024/annotations/mm39.ncbiRefSeq.gtf -g genelist.txt

def main():
    args = parser.parse_args()
    outdir: str = args.outdir
    outfname: str = args.outfname
    bam: str = args.bam
    globstr: str = args.glob
    genelistf: str = args.genelist
    gtf_ref: str = args.annot
    only_genes: bool = args.only_genes
    rewrite_db: bool = args.rewrite_db
    add_chr: bool = args.add_chr

    run_command(f"mkdir {outdir}")
    
    if os.path.isfile(genelistf):
        genelist = [i.strip() for i in open(genelistf)]
    else:
        logger.log("No gene list.", bcolors.WARNING)
        genelist = []
    db_ref = get_gtf_db(gtf_ref, only_genes, rewrite_db)
    if bam.endswith(".bam"):
        logger.log(f"Analyzing single bam: {bam}", bcolors.OKBLUE)
        results = analyze(bam, db_ref, genelist, add_chr=add_chr)
        results = pd.DataFrame(results)
        results.to_csv(f"{outdir}/{outfname}", index=False, sep="\t")
    elif os.path.isdir(bam):
        logger.log(f"Analyzing bams in directory: {bam}", bcolors.OKBLUE)
        for bamf in glob.glob(f"{bam}/{globstr}"):
            logger.log(f"****\nAnalyzing bam: {bamf}\n****", bcolors.OKBLUE)
            results = analyze(bamf, db_ref, genelist, add_chr=add_chr)
            results = pd.DataFrame(results)
            results.to_csv(f"{outdir}/{os.path.basename(bamf)}.tsv", index=False, sep="\t")
    
                                
if __name__ == "__main__":
    main()
