
# toplevel functions for working with RNA sequences

from rna_tools.tools.mini_moderna3.moderna.util.decorators import toplevel_function
from rna_tools.tools.mini_moderna3.moderna.util.validators import validate_alignment, validate_seq, \
    validate_filename, validate_path, \
    validate_alphabet, validate_alphabet_list
from rna_tools.tools.mini_moderna3.moderna.lphabet import alphabet
from rna_tools.tools.mini_moderna3.moderna.equence import Sequence


@toplevel_function
def load_alignment(file_path):
    """*load_alignment(file_path)*

    Loads a sequence alignment from a FASTA file.
    Produces an Alignment object that can be saved in a variable.

    ModeRNA expects, that the alignment file contains exactly two sequences.
    The first is for the target for which the model is to be built,
    the second for the structural template.

    Standard RNA bases should be written in upper case (ACGU).
    Standard DNA bases should be written in lower case (acgt).
    For modified bases, see the 'concepts' section of the manual.

    :Arguments:
        * path+filename of a FASTA file
    """
    file_path = validate_filename(file_path)
    return read_alignment(file_path)


@toplevel_function
def match_target_with_alignment(alignment, model):
    """*match_alignment_with_model(alignment, model)*

Checks, if the sequence of a model structure is equal to the first sequence in the alignment.
Writes an according message and returns True or False.

Both sequences also count as equal if one has modified nucleotides, and
the other the corresponding unmodified nucleotides in the same position,
or if one of the sequences contains the wildcard symbol '.'.
Thus, the sequence "AGU" is equal to both "A7Y" and "A.U".

:Arguments:
    * Alignment object
    * RnaModel object
    """
    alignment = validate_alignment(alignment)
    model = validate_model(model)

    log.write_message('\nChecking whether alignment matches with model.')
    am = AlignmentMatcher(alignment)
    seq = model.get_sequence()
    result = am.is_target_identical(seq)
    if not result:
        log.write_message("alignment and model match.\n")
    else:
        log.write_message("ALIGNMENT AND MODEL SEQUENCES DO NOT MATCH YET !!!\n")
    return result


@toplevel_function
def match_template_with_alignment(template, alignment):
    """*match_template_with_alignment(template, alignment)*

Checks, if the sequence of the template structure is equal
to the second sequence in the alignment. Writes an according message and returns True or False.
Small inconsistencies between both sequences, e.g. backbone breaks, or missing modification symbols
are corrected automatically in the alignment, and changes are reported in the logfile.

Both sequences also count as equal if one has modified nucleotides, and
the other the corresponding unmodified nucleotides in the same position,
or if one of the sequences contains the wildcard symbol '.'.
Thus, the sequence "AGU" is equal to both "A7Y" and "A.U".

:Arguments:
    * Template object
    * Alignment object
    """
    template = validate_template(template)
    alignment = validate_alignment(alignment)

    log.write_message('Checking whether template matches with alignment.')
    am = AlignmentMatcher(alignment)
    seq = template.get_sequence()
    am.fix_template_seq(seq)
    result = am.is_template_identical(seq)
    if result:
        log.write_message("template and alignment match.\n")
    else:
        log.write_message("TEMPLATE AND ALIGNMENT DO NOT MATCH!\n")
    return result

