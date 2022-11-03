from __future__ import annotations
import logging


def getLogger(name: str = None):
    """
    Retrieve a Logger object

    Parameters
    ----------
    name : str, optional
        The name of the logging object
    """
    if name is None:
        pass

    log = logging.getLogger("haptools " + name)
    db_time = "|%(asctime)s" if verbosity == "DEBUG" else ""
    logging.basicConfig(
        format="[%(levelname)8s" + db_time + "] %(message)s (%(filename)s:%(lineno)s)",
        level=verbosity,
        datefmt="%H:%M:%S",
    )
