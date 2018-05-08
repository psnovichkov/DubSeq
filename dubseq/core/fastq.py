#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Clases to read FASTQ files
"""

import os
import sys
import subprocess

try:
    import io
    io_method = io.BytesIO
except:
    import cStringIO
    io_method = cStringIO.StringIO


class FastqRecord:
    """ Represents a record of data in FASTQ format

    Attributes:
        id: id fo the FASTQ record.  The first line of the FASTQ record (should start from @ symbol)
        sequence: nucelcotide sequence. The second line of the FASTQ record
        description: description of the record. The thrid line of the FASTQ record (should start from + symbol)
        quality: the quality of the sequence (from the second line). The forth line of the FASTQ record
    """

    __slots__ = ('__id', '__sequence', '__description', '__quality')

    def __init__(self, id='', sequence='', description='', quality=''):
        """ Inits  FastqRecord with values"""
        self.__id = id
        self.__sequence = sequence
        self.__description = description
        self.__quality = quality

    def validate(self):
        """ Checks if first and third line starts with expected symbols
        Args:
        Returns:
        Raises:
            ValueError: occurs when the format is wrong 
        """
        if not self.__id or self.__id[0] != '@':
            raise ValueError(
                'FASTQ format: the first line of a record should start with @ symbol.'
                ' Line value = ' + str(self.__id))
        if not self.__description or self.__description[0] != '+':
            raise ValueError(
                'FASTQ format: the thrid line of a record should start with + symbol.'
                ' Line value= '
                + str(self.__description))

    def __call__(self, id='', sequence='', description='', quality=''):
        """ Set values of the record
        Args:
            id: id of the record
            sequence: nucleotide sequence
            description: description of the record
            quality: quality of the nucleotide record
        Returns:
        Raises:
        """
        self.__id = id
        self.__sequence = sequence
        self.__description = description
        self.__quality = quality

    @property
    def id(self):
        return self.__id

    @property
    def sequence(self):
        return self.__sequence

    @property
    def description(self):
        return self.__description

    @property
    def quality(self):
        return self.__quality


class FastqReader:
    __FILE_TYPE_GZ = 'gz'
    __FILE_TYPE_TXT = 'txt'

    def __init__(self, file_name):
        self.__file_name = file_name
        self.__file_type = None
        self.__file = None

        if file_name.endswith('.gz'):
            self.__file_type = FastqReader.__FILE_TYPE_GZ
            p = subprocess.Popen(['gunzip', '-c', file_name],
                                 stdout=subprocess.PIPE)
            self.__file = io_method(p.communicate()[0])
        else:
            self.__file_type = FastqReader.__FILE_TYPE_TXT
            self.__file = open(file_name, 'r')

    def __nextline(self):
        line = self.__file.readline()
        if not line:
            raise EOFError()

        if self.__file_type == FastqReader.__FILE_TYPE_GZ:
            line = line.decode("utf-8")

        return line.strip()

    @property
    def file_type(self):
        return self.__file_type

    @property
    def file_name(self):
        return self.__file_name

    def next_record(self, record):
        try:
            record(
                self.__nextline(),
                self.__nextline(),
                self.__nextline(),
                self.__nextline())
            record.validate()
        except EOFError:
            return False

        return True

    def next(self):
        try:
            record = FastqRecord(
                self.__nextline(),
                self.__nextline(),
                self.__nextline(),
                self.__nextline())
            record.validate()
        except EOFError:
            return None

        return record

    def close(self):
        self.__file.close()


class FastqFileStat:

    @staticmethod
    def header(sep='\t'):
        return sep.join([
            'total_reads_count',
            'barcode_extracted_reads_count',
            'barcode_extracted_reads_percent'])

    def __init__(self):
        self.__total_reads_count = 0
        self.__barcode_extracted_reads_count = 0

    def total_reads_inc(self):
        self.__total_reads_count += 1

    def barcode_extracted_reads_inc(self):
        self.__barcode_extracted_reads_count += 1

    def __str__(self):
        return "\t".join(str(x) for x in [
            self.__total_reads_count,
            self.__barcode_extracted_reads_count,
            round(self.__barcode_extracted_reads_count * 100.0 / self.__total_reads_count)])
