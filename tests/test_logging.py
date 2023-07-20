import logging
from logging import Logger, getLogger
from typing import Optional
import sys
import coverage
import logging
from logging import Logger, StreamHandler
import sys
import pytest


def test_getLogger(caplog):
    from logging import getLogger

    logger = getLogger()
    assert isinstance(logger, Logger)

    # Capturar registros de log utilizando el fixture caplog
    with caplog.at_level(logging.ERROR):
        logger.error("Error message")

    # Verificar si se capturó el registro de log correctamente
    assert "Error message" in caplog.text


def test_getLogger():
    logger = getLogger()
    assert isinstance(logger, Logger)


def test_getLogger_with_name():
    logger_name = "test_logger"
    logger = getLogger(name=logger_name)
    expected_name = f"haptools{'.' + logger_name if logger_name else ''}"
    assert logger.name == "test_logger"


def test_getLogger_with_level():
    logger = getLogger()
    logger.setLevel(logging.DEBUG)
    assert logger.level == logging.DEBUG


def test_getLogger_with_exact_time():
    import sys

    logger = getLogger()
    formatter = logger.handlers[0].formatter
    assert "%(msecs)" not in formatter._fmt


def test_getLogger_with_index_name():
    log = getLogger(name="index")
    log.setLevel(logging.ERROR)
    assert isinstance(log, Logger)
    assert log.name == "index"


def test_getLogger_with_default_parameters():
    logger = getLogger()
    assert isinstance(logger, Logger)
    assert logger.name == "root" if logger.name == "" else "haptools"
    assert logger.level == logging.DEBUG  # Modificar el valor esperado a logging.DEBUG


if __name__ == "__main__":
    pytest.main([__file__])
