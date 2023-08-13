from __future__ import annotations
import shutil
import logging
import tempfile
from pathlib import Path

from pysam import tabix_index

from . import data
from .logging import getLogger


def append_suffix(
    path: Path,
    suffix: str,
):
    """
    Used as a helper method for index_haps. Appends a given suffix to a Path instance.

    Parameters
    ----------

    path : Path
        The path to a file
    suffix : str
        A string to append to the end of the given Path. For example, ".gz" or ".gz.tbi"
    """
    return path.with_suffix(path.suffix + suffix)


def index_haps(
    haplotypes: Path,
    sort: bool = False,
    output: Path = None,
    log: logging.Logger = None,
):
    """
    Takes in an unsorted .hap file and outputs it as a .gz and a .tbi file

    Parameters
    ----------

    haplotypes : Path
        The path to the haplotypes in a .hap file
    output : Path, optional
        The location to which to write output. If an output location is not specified,
        the output will have the same name as the input file.
    log : Logger, optional
        A logging module to which to write messages about progress and any errors
    """
    if log is None:
        log = getLogger(name="index", level="ERROR")

    log.info("Loading haplotypes")
    hp = data.Haplotypes(haplotypes, log=log)

    if sort:
        hp.read()
        hp.sort()
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            hp.fname = Path(tmp.name)
            log.debug(f"writing haplotypes to {hp.fname}")
            hp.write()
    else:
        # copy the file to a tmp location in case the input is /dev/stdin
        # or a file that might otherwise be deleted by tabix_index afterward
        with tempfile.NamedTemporaryFile(delete=False, mode="wt") as tmp:
            with data.Data.hook_compressed(str(hp.fname), mode="r") as haps:
                hp.fname = Path(tmp.name)
                tmp.write(haps.read())

    try:
        tabix_index(str(hp.fname), seq_col=1, start_col=2, end_col=3)
    except OSError as e:
        # check if the error message matches what we expect if the file is unsorted
        if str(e).startswith("building of index for "):
            log.error("Indexing failed. Is your file properly sorted?")
        else:
            # otherwise, re-raise it
            raise

    hp.fname = append_suffix(hp.fname, ".gz")

    if output is None:
        if haplotypes.suffix.endswith(".gz"):
            output = haplotypes
        else:
            output = append_suffix(haplotypes, ".gz")
    # use shutil instead of reanme b/c it won't error out if /tmp is mounted elsewhere
    shutil.copy(str(hp.fname), str(output))
    hp.fname.unlink()
    tbi_file = append_suffix(hp.fname, ".tbi")
    shutil.copy(str(tbi_file), str(append_suffix(output, ".tbi")))
    tbi_file.unlink()
