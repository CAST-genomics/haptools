from __future__ import annotations
import logging


def getLogger(name: str = None, level: str = "ERROR", exact_time=False):
    """
    Retrieve a Logger object

    Parameters
    ----------
    name : str, optional
        The name of the logging object
    level : str, optional
        The level of verbosity for the logger
    """
    if name is None:
        name = ""
    else:
        name = "." + name

    # create logger
    logger = logging.getLogger("haptools" + name)
    logger.setLevel(level)

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(level)

    # create formatter
    db_time = (
        ("|%(asctime)s" + (".%(msecs)03d" if exact_time else ""))
        if level == "DEBUG"
        else ""
    )
    formatter = logging.Formatter(
        fmt="[%(levelname)8s" + db_time + "] %(message)s (%(filename)s:%(lineno)s)",
        datefmt="%H:%M:%S",
    )

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)

    return logger
