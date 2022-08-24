from __future__ import annotations
import logging
from pathlib import Path
from haptools.data.haplotypes import Haplotypes
from pysam import tabix_index
from haptools import data
import tempfile

""""
doc string
ouput if not specified can override the input -- if no output specified then write over where th e
original input was  if the file is already compressed  but if not then we should write a temp file
for temp sort, then write to temp, the bgzip it, then move it to output path -- if output not specific then move it to directory with same input path

also write doc string for less than and sort for haplotypes and haplotype and variant class
"""

def index_haps(
    haplotypes: Path,
    log: Logger = None,
    output: Path = None,
):

    if log is None:
        log = logging.getLogger("run")
        logging.basicConfig(
            format="[%(levelname)8s] %(message)s (%(filename)s:%(lineno)s)",
            level="ERROR",
        )

    log.info("Loading haplotypes")
    #creates instance of haplotypes class
    hp = data.Haplotypes(haplotypes, log=log)
    #load the file into memory
    hp.read()
    hp.sort()
    #before tabix we need to write it as a temp file
    with tempfile.NamedTemporaryFile() as tmp:
        #tmp.name creates unique file when ran
        hp.fname = Path(tmp.name)
        hp.write()

     #compresses and indexes file   
    tabix_index(str(hp.fname), seq_col=2, start_col=3, end_col=4)

    #move temp path to output path
    #check if it is none
    if output is None:
        if haplotypes.suffix.endswith(".gz"):
            output = haplotypes
        #if none we want to change output to be same as input
        else:
            output = haplotypes.append_suffix(".gz")
    #rename/move as a path        
    hp.fname.rename(output)

    #repeat process for .tbi
    hp.fname.append_suffix(".tbi").rename(output.append_suffix(".tbi"))
    


    
    



