from __future__ import annotations

import logging
import os


def get_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    level = os.getenv('PHARMAGENTS_LOG_LEVEL', 'INFO').upper()
    logger.setLevel(level)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(asctime)s | %(levelname)s | %(name)s | %(message)s'))
    logger.addHandler(handler)
    logger.propagate = False
    return logger
