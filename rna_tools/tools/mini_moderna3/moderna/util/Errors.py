#!/usr/bin/env python
#
# Errors.py
#
# Contains exception classes.
# 
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"

"""
The exception model of Moderna contains
one separate exception class for each class
in the Moderna program. This allows to trace 
back easily where something went wrong.
There are a few additional classes for important
functionalities of Moderna, in particular all things
written on the SocBin2008 poster.

In general, each time an exception derived from
ModernaError occurs, this could mean that something is
wrong with input. If some other Python exception occurs,
it is always the fault of the developers.
"""

class ModernaError(Exception): pass
# KR: we could __init__ have write a message to log.
# this would make the logfile calls in moderna.py obsolete.

class AlignmentError(ModernaError): pass

class AlphabetError(ModernaError): pass

class IsostericityError(ModernaError): pass

class LirError(ModernaError): pass

class LirRecordError(ModernaError): pass

class SearchLirError(ModernaError): pass

class LirCandidatesError(ModernaError): pass

class RNAResidueError(ModernaError): pass

class ModernaResidueError(RNAResidueError): pass

class RNAChainError(ModernaError): pass

class ModernaStructureError(RNAChainError): pass

class ModernaSuperimposerError(ModernaError): pass

class ProcessPDBError(ModernaError): pass

class SequenceError(ModernaError): pass

class RenumeratorError(ModernaError): pass


class AddModificationError(ModernaResidueError): pass

class ExchangeBaseError(ModernaResidueError): pass

class RemoveModificationError(ModernaResidueError): pass


class ModernaFragmentError(ModernaStructureError): pass

class RnaModelError(ModernaStructureError): pass

class TemplateError(ModernaStructureError): pass

class ModernaAdunModelError(ModernaStructureError): pass


class CopyResidueError(RnaModelError): pass

class InsertLoopError(RnaModelError): pass

class BaseRecognitionError(RnaModelError): pass

class ParameterError(ModernaError): pass


