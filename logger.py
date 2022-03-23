#! /usr/bin/env python
import logging
import sys
import datetime

class Logger(object):
    def __init__(self, logfile=None, verbose=True):
        """
        A simple logger.
        """
        self.console = sys.stdout
        self.verbose = verbose
        if logfile is not None:
            self.log = open(logfile, 'w')
        else:
            self.log = None

        time = datetime.datetime.now()
        self.write((f'{time.strftime("%A")},'
                    f' {time.strftime("%b")}' 
                    f' {time.strftime("%d")},'
                    f' {time.strftime("%Y")}' 
                    f' {time.strftime("%I")}:'
                    f'{time.strftime("%M")}'
                    f' {time.strftime("%p")}'
        ))

    def write(self, message):
        if self.verbose:
            self.console.write(message+'\n')
        if self.log is not None:
            self.log.write(message+'\n')
            self.log.flush()

    def verbose(self, verbose=True):
        self.verbose = verbose
        
def create_logger(logfile):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.ERROR)
    console_handler_format = '%(asctime)s | %(levelname)s: %(message)s'
    console_handler.setFormatter(logging.Formatter(console_handler_format))
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(logfile)
    file_handler.setLevel(logging.INFO)
    file_handler_format = '%(asctime)s | %(levelname)s | %(lineno)d: %(message)s'
    file_handler.setFormatter(logging.Formatter(file_handler_format))
    logger.addHandler(file_handler)
    return logger
