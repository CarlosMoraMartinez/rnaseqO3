import os
import glob
from typing import List, Set, Dict, Tuple

DEFAULT_N_READS: int = 1000000
DEFAULT_READ_LENGTH: int = 150
DEFAULT_FRAGMENT_LENGTH: int = 300
DEFAULT_FRAGMENT_SD: int = 50
DEFAULT_BASE_ERROR_RATE: float = 0.0
DEFAULT_MUTATION_RATE: float = 0.001
DEFAULT_INDEL_FRACTION: float = 0
DEFAULT_PROB_INDEL_EXT: float = 0
DEFAULT_SEED: int = 123

indir = "fasta_by_gene_geneNames/tanda2"
outdir = "simulated_fastq_bygene_error0"


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
 
 
def wgsim_get_sample(input_fasta: str, output_fastq: str,
                     num_reads: int, read_length_r1: int, read_length_r2: int, 
                     fragment_length: int, fragment_sdesv: int,
                     base_error_rate: float, mutation_rate: float, 
                     indel_fraction: float, prob_indel_ext: float,
                     rewrite_sim_fastq: bool,
                     seed: int) -> Tuple[str, str]:
    """Creates and runs wgsim command to simulate a given number of reads from an input fasta.

    Args:
        input_fasta (str): Optionally compressed fasta file to simulate reads from. 
        output_fastq (str): Name of the output fastq files.
        num_reads (int): Number of reads to generate. Read wgsim docs.
        read_length_r1 (int): Length of read 1 in each read pair. Read wgsim docs.
        read_length_r2 (int): Length of read 2 in each read pair. Read wgsim docs.
        fragment_length (int): Read wgsim docs.
        fragment_sdesv (int): Read wgsim docs.
        base_error_rate (float): Read wgsim docs.
        mutation_rate (float): Read wgsim docs.
        indel_fraction (float): Read wgsim docs
        prob_indel_ext (float): Read wgsim docs.
        rewrite_sim_fastq (bool): Rewrite output file if it exists.
        seed (int): Random seed for wgsim.

    Returns:
        Tuple[str, str]: Names of the two compressed fastq files generated. 
    """
    r1: str = f"{output_fastq}_R1.fastq"
    r2: str = f"{output_fastq}_R2.fastq"
    if os.path.isfile(r1 + '.gz') or os.path.isfile(r2 + '.gz'):
        if rewrite_sim_fastq:
            cmd_rm: str = f"rm {r1}.gz {r2}.gz"
            logger.log(f"----Species files already exist. Removing", bcolors.WARNING)
            run_command(cmd_rm)
        else:
            logger.log(f"----Species files already exist. Skipping WGSIM", bcolors.WARNING)
            return(r1 + '.gz', r2 + '.gz')
    
    logfile: str = f"{output_fastq}_WGSIM.log"
    errfile: str = f"{output_fastq}_WGSIM.err"
    cmd: str = (f"wgsim -e {base_error_rate}"
                f" -d {fragment_length}"
                f" -s {fragment_sdesv}"
                f" -N {num_reads}"
                f" -1 {read_length_r1}"
                f" -2 {read_length_r2}"
                f" -r {mutation_rate}"
                f" -R {indel_fraction}"
                f" -X {prob_indel_ext}"
                f" -S {seed}"
                f" {input_fasta} {r1} {r2} >{logfile} 2>{errfile};" 
                f" pigz {r1}; pigz {r2}")
    run_command(cmd)
    return(r1 + '.gz', r2 + '.gz')
    
    
    
def main():
    
    try:
        run_command(f"mkdir {outdir}")
    except:
        logger.log("Output dir already exists", bcolors.WARNING)
        
    files = glob.glob(f"{indir}/*.fa")

    for f in files:
        logger.log(f"--Simulating: {f}", bcolors.HEADER)
        output_fastq = os.path.join(outdir, os.path.basename(f).replace('.fa', ''))
        logger.log(f"----Output: {output_fastq}", bcolors.OKCYAN)
        sim_files = wgsim_get_sample(input_fasta=f, 
                                     output_fastq=output_fastq,
                                     num_reads=DEFAULT_N_READS, 
                                     read_length_r1=DEFAULT_READ_LENGTH,
                                     read_length_r2=DEFAULT_READ_LENGTH, 
                                     fragment_length=DEFAULT_FRAGMENT_LENGTH, 
                                     fragment_sdesv=DEFAULT_FRAGMENT_SD,
                                     base_error_rate=DEFAULT_BASE_ERROR_RATE, 
                                     mutation_rate=DEFAULT_MUTATION_RATE, 
                                     indel_fraction=DEFAULT_INDEL_FRACTION, 
                                     prob_indel_ext=DEFAULT_PROB_INDEL_EXT,
                                     rewrite_sim_fastq=True,
                                     seed=DEFAULT_SEED)
                                     
if __name__ == "__main__":
    main()

