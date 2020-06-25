#!/usr/bin/env python
#
# commands.py
#
# commands for RNA 3D structure prediction.
# 
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "1.7.1"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"

"""
 ModeRNA

 A program for comparative RNA structure modeling

 Authors: 
    Magdalena Rother
    Kristian Rother
    Tomasz Puton
    Janusz Bujnicki

 (c) 2008
 
 www.genesilico.pl 

"""

from rna_tools.tools.mini_moderna3.moderna.util.LogFile import log
from rna_tools.tools.mini_moderna3.moderna.RNAModel import RnaModel
from rna_tools.tools.mini_moderna3.moderna.Template import Template
from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaAlphabet import Alphabet
from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaSequence import Sequence
from rna_tools.tools.mini_moderna3.moderna.sequence.RNAAlignment import read_alignment
from rna_tools.tools.mini_moderna3.moderna.ModelingRecipe import RecipeMaker
from rna_tools.tools.mini_moderna3.moderna.ModernaStructure import ModernaStructure
from rna_tools.tools.mini_moderna3.moderna.ModernaFragment import ModernaFragment53, ModernaFragment5, \
    ModernaFragment3, ModernaFragmentStrand, ModernaFragment2D
from rna_tools.tools.mini_moderna3.moderna.Helix import HelixBuilder, HelixFragmentBuilder
from rna_tools.tools.mini_moderna3.moderna.fragment_library.SearchLIR import *
from rna_tools.tools.mini_moderna3.moderna.analyze.ClashRecognizer import ClashRecognizer
from rna_tools.tools.mini_moderna3.moderna.analyze.GeometryAnalyzer import GeometryAnalyzer
from rna_tools.tools.mini_moderna3.moderna.analyze.StackingCalculator import StackingCalculator
from rna_tools.tools.mini_moderna3.moderna.builder.BackboneBuilder import BackboneBuilder
from rna_tools.tools.mini_moderna3.moderna.builder.ChiRotator import rotate_chi as rc
from rna_tools.tools.mini_moderna3.moderna import modifications
from rna_tools.tools.mini_moderna3.moderna.CheckPdb import PdbController
from rna_tools.tools.mini_moderna3.moderna.sequence.AlignmentMatcher import AlignmentMatcher

from rna_tools.tools.mini_moderna3.moderna.util.Errors import ModernaError

from rna_tools.tools.mini_moderna3.moderna.Constants import NUMBER_OF_FRAGMENT_CANDIDATES

from rna_tools.tools.mini_moderna3.moderna.util.decorators import toplevel_function
from rna_tools.tools.mini_moderna3.moderna.util.validators import validate_structure, validate_template, validate_model, \
    validate_alignment, validate_seq, validate_resnum, validate_resnum_list, \
    validate_resi, validate_fragment, validate_frag_candidate_list,  \
    validate_resi_list, validate_filename, validate_path, \
    validate_secstruc, validate_alphabet, validate_alphabet_list


##################################### SETTINGS #####################################

@toplevel_function
def add_all_modifications(template, alignment, model):
    """*add_all_modifications(template, alignment, model)*

Adds all modifications that occur in the target sequence, and are
not present in the template.
    
:Arguments:
    * Template object
    * Alignment object
    * RnaModel object
    """
    template = validate_template(template)
    alignment = validate_alignment(alignment)
    model = validate_model(model)
    
    if template: model.template = template
    if alignment:
        model.alignment = alignment
        model.recipe = RecipeMaker(alignment).recipe
    model.add_all_modifications_copy()    


@toplevel_function
def add_modification(residue, modification_name, model=None, residue_number=None):
    """*add_modification(residue, modification_name, model=None, residue_number=None)*

Adds a modification to a single residue. An RNAModel can be given optionally. If this is given, the modified residue is subsequently copied to it.

Note that deoxynucleotides count as modified bases as well.

:Arguments:
    * residue from a Template or RNAModel object
    * modification name (abbreviated as e.g. 'm6Am', 'm1G', 'm3U', 'mnm5s2U')
    * model (optional; if this is given, the residue is copied)
    * residue position in model (optional; by default the residue number is kept)
    """
    residue = validate_resi(residue)
    modification_name = validate_alphabet(modification_name)
    if model: model = validate_model(model)
    if residue_number: validate_resnum(residue_number)
    
    if model != None:
        model.add_one_modification_copy(residue, modification_name, residue_number)
    else:
        modifications.add_modification(residue, modification_name)


@toplevel_function
def add_pair_to_base(model,  anchor_id,  new_id=None,  new_sequence=None):
    """*add_pair_to_base(model,  anchor_id,  new_id=None,  new_sequence=None)*

Adds Watson-Crick paired second strand to single nucleotide.

:Arguments:
    * RNAModel object
    * Identifier of residue to which a new base pair will be added
    * New identifier for second strand
    * New sequence for second strand
    """
    model = validate_model(model)
    anchor_id = validate_resnum(anchor_id)
    if new_id: new_id = validate_resnum(new_id)
    if new_sequence: new_sequence = validate_seq(new_sequence)
    
    fr = ModernaFragmentStrand(anchor=model[anchor_id], identifier=new_id, new_sequence=new_sequence)  
    model.insert_fragment(fr)
    

@toplevel_function
def analyze_geometry(struc, file_name=None):
    """*analyze_geometry(struc, file_name=None)*

Analyses geometry of a structure. Checks whether all angles and bonds length values
are close to reference values determined from the PDB. 
Writes the result of the analysis to log file or to an output file when specified.

:Arguments:
    * RNAModel object
    * Name of an output file
    """
    struc = validate_structure(struc)
    if file_name: file_name = validate_filename(file_name)
    
    analyzer = GeometryAnalyzer(struc)
    analyzer.analyze()
    raport = str(analyzer)
    if file_name: 
        f = open(file_name, 'w')
        f.write(raport)
    log.write_message(raport) 
    return analyzer
    

@toplevel_function
def apply_alignment(template, alignment, model):
    """*apply_alignment(template, alignment, model)*

Applies all operations on single bases that are in an alignment:
Copying identical residues, exchanging mismatches, adding modifications,
and removing modifications. As a result, all these residues are copied
into a model (which may be empty).
For gaps in the alignment, nothing is done.

:Arguments:
    * Template object
    * Alignment object
    * RNAModel object
    """
    template = validate_template(template)
    alignment = validate_alignment(alignment)
    model = validate_model(model)
    
    model.alignment = alignment
    model.recipe = RecipeMaker(alignment).recipe
    model.template = template
    model.apply_alignment()


@toplevel_function
def apply_all_indels(alignment, model): 
    """*apply_all_indels(alignment, model)*
    
Finds and inserts the best fragment candidate for each gap region in the alignment. 
This command can be used to fill in gaps automatically. 
The procedure automatically retrieves the 20 best fitting candidates 
from the fragment database, checks whether they produce any clashes upon insertion, 
and inserts the best scoring candidate.

It depends, however on
    * gaps being shorter than 17 bases.
    (if longer a database of fragments up to 100 nt is available)

:Arguments:
    * Alignment object
    * RnaModel object
    """
    alignment = validate_alignment(alignment)
    model = validate_model(model)
    
    if alignment:
        model.alignment = alignment
        model.recipe = RecipeMaker(alignment).recipe
    model.insert_all_fragments()


@toplevel_function
def apply_indel(model, res5_number, res3_number, sequence):
    """*apply_indel(model, res5_number, res3_number, sequence)*

Finds the best fitting fragment between two residues, and inserts it into the model. 
The procedure automatically retrieves the 20 best fitting candidates from the fragment database, 
checks whether they produce any clashes upon insertion, and inserts the best scoring candidate.

:Arguments:
    * RNAmodel object
    * number of the residue on the 5' end of the fragment
    * number of the residue on the 3' end of the fragment
    * sequence of the fragment (not including the two residues specified above).
    """
    model = validate_model(model)
    res5_number = validate_resnum(res5_number)
    res3_number = validate_resnum(res3_number)
    sequence = validate_seq(sequence)
    
    model.insert_best_fragment(res5_number, res3_number, sequence) 


@toplevel_function
def apply_missing_ends(alignment,  model):
    """*apply_missing_ends(alignment,  model)*
    
Attaches an extended fragment of structure at both
the 5' and 3' ends of the target, where there is no corresponding
part in the template.

:Arguments:
    * Alignment object
    * model
    """
    alignment = validate_alignment(alignment)
    model = validate_model(model)
    
    if alignment:
        model.alignment = alignment
        model.recipe = RecipeMaker(alignment).recipe
    model.add_missing_5p()
    model.add_missing_3p()

@toplevel_function
def change_sequence(struc, sequence):
    """*change_sequence(struc, sequence)*
    
Changes the entire sequence of the given structure object, exchanging all nucleotides that differ. 

:Arguments:
    * Structure object (RnaModel or Template)
    * sequence
    """
    struc = validate_structure(struc)
    sequence = validate_seq(sequence)
    struc.change_sequence(sequence)

@toplevel_function
def clean_structure(st,  write_structure=False):
    """*clean_structure(structure, write_structure=False)*

Eliminates features that may cause problems during modeling from a template or model structure:
    * water molecules, ions, amino acids, and unidentified residues are deleted.
    * old atom names are replaced (C1* becomes C1', O1P becomes OP1). 
    * missing phosphate groups are added (also adds the OP3 atom for these)
    * reports whether the chain is continuous.

In case some feature cannot be fixed (e.g. chain discontinuity) this is written to the logfile. 
It is  recommended to include such features in the alignment
('.' characters for strange residues, and '_' for backbone breaks).

:Arguments:
    * Stucture object (RnaModel or Template)
    * True/False - whether structure should be written to a PDB file (optional)
    """
    struc = validate_structure(st)
    
    pc = PdbController(st)
    pc.clean_structure()
    log.write_message(str(pc))
    if write_structure: pc.write_structure('fixed_structure.pdb')
    return pc


##################################### COPYING #####################################

@toplevel_function
def copy_identical_residues(template, alignment, model, strict=True, modifications=True):
    """*copy_identical_residues(template, alignment, model, strict=True, modifications=True)*

Copies all bases identical identical in the alignment
from the template structure to a model. The residue numbers
do not change.

:Arguments:
    * Template object
    * Alignment object
    * RnaModel object
    * strict=1 complains and stops on any problem encountered (default); strict=0 goes on if single residues fail to be copied
    * modifications=0 does not copy modified bases; modifications=1 is default.
    """
    template = validate_template(template)
    alignment = validate_alignment(alignment)
    model = validate_model(model)
    
    if template: model.template = template
    if alignment:
        model.alignment = alignment
        model.recipe = RecipeMaker(alignment).recipe
    model.copy_all_residues(strict=strict, modifications=modifications)
   
  
@toplevel_function
def copy_single_residue(residue, model, new_number=None, strict=True):
    """*copy_single_residue(residue, model, new_number=None, strict=True)*
    
Copies a single residue into the model.
The residue can be taken from a Template, or another RnaModel.

Residues are specified by their PDB residue number in quotes and square brackets, e.g. template['5'].
eventually accompanied by the insertion letter (abbreviated as number).
By default, the number of the residue is kept, but it can be given optionally. 

:Arguments:
    * Residue from a Template or RnaModel
    * RnaModel object
    * new residue number in the model (optional)
    * strict=1 complains when a residue with the given number exists already (default); strict=0 copies anyway
    """
    residue = validate_resi(residue)
    model = validate_model(model)
    if new_number: new_number = validate_resnum(new_number)
    
    model.copy_residue(residue, new_number, strict=strict)


@toplevel_function
def copy_some_residues(residue_list, model, new_numbers=None, strict=True):
    """*copy_some_residues(residue_list, model, new_numbers=None, strict=True)*

Copies a set of residues into the model.
The residues can be taken from a Template, or another RnaModel.

Residues are specified by enumerating their numbers in the PDB files,
or specifying a range of numbers (by 'num1':'num2').
By default, the numbers of the residues are kept,
but they can be given optionally.
Note that all residue numbers in ModeRNA have to be quoted.

:Arguments:
    * List of residues from a Template or RnaModel
    * RnaModel object
    * List of new residue numbers in the model (optional)
    * strict=1 complains when a residue with the given number exists already (default); strict=0 copies anyway
    """
    residue_list = validate_resi_list(residue_list)
    model = validate_model(model)
    if new_numbers: new_numbers = validate_resnum_list(new_numbers)
    
    model.copy_list_of_residues(residue_list, new_numbers, strict=strict)


@toplevel_function
def create_fragment(pdb_path, anchor5=None, anchor3=None, chain_name='A', sequence=None):
    """*create_fragment(pdb_path, anchor5=None, anchor3=None, chain_name='A', sequence=None)*

To insert small pieces of custom PDB structures into models, fragments can be used. 
This command loads a fragment, and defines one or two connection points, to which 
the fragment will be inserted by superposition. It returns a ModernaFragment object.

To add a ModernaFragment object to a model, the insert_fragment command should be used.
The scenarios for adding either involve superposition of a single residue with the model
on the 5' end of the fragment, or on the 3' end, or both.
   
:Arguments:
    * path to pdb file with the fragment
    * residue to which the 5' end of the fragment will be superimposed
    * residue to which the 3' end of the fragment will be superimposed
    * chain name of the fragment in the file (optional; A by default)
    * sequence that should be modeled onto the fragment upon insertion (optional; the sequence should not include the 1-2 anchor residues, it therefore has to be shorter than the fragment)
    """
    pdb_path = validate_path(pdb_path)
    if anchor5: anchor5 = validate_resi(anchor5)
    if anchor3: anchor3 = validate_resi(anchor3)
    if sequence: sequence = validate_seq(sequence)
    struc = ModernaStructure(data_type='file', data=pdb_path, chain_name=chain_name)
    
    if anchor5 and anchor3:
        fr = ModernaFragment53(struc, anchor5=anchor5, anchor3=anchor3, new_sequence=sequence, strict=False)
    elif anchor5:
        fr = ModernaFragment5(struc, anchor5=anchor5, new_sequence=sequence, strict=False)
    elif anchor3:
        fr = ModernaFragment3(struc, anchor3=anchor3, new_sequence=sequence, strict=False)
    else:
        raise ModernaError("Anchor residues need to be specified.")
    return fr


@toplevel_function
def create_model(template=None, alignment=None, model_chain_name=None):
    """*create_model(template=None, alignment=None, model_chain_name=None)*
  
Creates a RNA structure model.
Produces a RnaModel object that can be saved
in a variable (see example).

If no arguments are given, an empty RNAModel is created.
If both a Template and an Alignment are given as arguments,
the complete model is built automatically.
The chain id that the model should have can be specified optionally.

:Arguments:
    * Template object (optional)
    * Alignment object (optional)
    * chain id of the model to be built
    """
    if template: template = validate_template(template)
    if alignment: alignment = validate_alignment(alignment)
    
    model = RnaModel(template=template, alignment=alignment, model_chain_name=model_chain_name, data_type=None, data=None)
    if template and alignment:
        if not match_template_with_alignment(template,alignment):
            raise ModernaError("""Template and alignment sequences do not match!
The model cannot be built automatically.
Template : %s
Alignment: %s
"""%(template.get_sequence(),alignment.template_seq))
        model.create_model()
        match_alignment_with_model(alignment,model)
    return model



@toplevel_function
def delete_residue(residue_number, model):
    """*delete_residue(residue_number, model)*
    
Removes a residue from a RNAModel object.
   
:Arguments:
    * residue identifier e.g. '1', '6', '6A'
    * RnaModel object
    """
    model = validate_model(model)
    residue_number = validate_resnum(residue_number)
    
    model.remove_residue(residue_number)


@toplevel_function
def examine_structure(st, ex_log=None,  verbose=True):
    """*examine_structure(structure)*

Checks whether the given structure has any features that may cause any problems during the modeling process.
The user needs to give a structure object, and optionally a name of the file the report is written to.

:Arguments:
    * Stucture object
    * name of logfile (optional)
    """
    struc = validate_structure(st)

    pc = PdbController(st)
    if ex_log: pc.write_log(ex_log)
    else: log.write_message(str(pc))
    if verbose: print(pc)
    return pc


##################################### EXCHANGING #####################################

@toplevel_function
def exchange_mismatches(template, alignment, model):
    """*exchange_mismatches(template, alignment, model)*
    
Exchanges the bases for all standard base mismatches in an alignment.
This applies all exchanges of one standard base by another.
Modified bases are not changed (use the apply_alignment command for that).
 
:Arguments:
    * Template object
    * Alignment object
    * RnaModel object
    """
    template = validate_template(template)
    alignment = validate_alignment(alignment)
    model = validate_model(model)
    
    if template: model.template = template
    if alignment:
        model.alignment = alignment
        model.recipe = RecipeMaker(alignment).recipe
    model.exchange_all_bases()


@toplevel_function
def exchange_single_base(residue, new_name, model=None, new_number=None):
    """*exchange_single_base(residue, new_name, model=None, new_number=None)*
    
Exchanges standard RNA bases in a residue.
This can be used to modify residues already in a model,
or to exchange and copy them in one step, if an RnaModel is given.
It is also possible to modify the template this way!

The residue number does not change by default, but a new one can be given.

:Arguments:
    * Residues from a Template or RnaModel
    * Abbreviation of the new base ('A', 'G', 'C' or 'U')
    * RnaModel object (optional)
    * New residue number in the model (optional)
    """
    residue = validate_resi(residue)
    new_name = validate_alphabet(new_name)
    if model: model = validate_model(model)
    if new_number: new_number = validate_resnum(new_number)

    if model:
        model.copy_residue(residue, new_number)
        number = new_number or residue.identifier
        modifications.exchange_base(model[number], new_name)
    else:
        modifications.exchange_base(residue, new_name)

@toplevel_function
def exchange_some_bases(residue_list, new_names, model=None, new_numbers=None):
    """*exchange_some_bases(residue_list, new_names, model=None, new_numbers=None)*

Exchanges standard RNA bases in a set of residues.
This can be used to modify residues already in a model,
or to exchange and copy them in one step, if a RNAModel is given.
It is also possible to modify the template this way!

The residue numbers are not changed by default, but a list of
new one can be given.

:Arguments:
    * Residues from a Template or RnaModel
    * List of names of the new bases ('A', 'G', 'C' or 'U')
    * RnaModel object (optional)
    * New residue number in the model
    """
    residue_list = validate_resi_list(residue_list)
    new_names = validate_alphabet_list(new_names, len(residue_list))
    if model: model = validate_model(model)
    if new_numbers: new_numbers = validate_resnum_list(new_numbers, len(residue_list))
    
    if model: 
        model.exchange_list_of_bases(residue_list, new_names, new_numbers)
    else:
        for resi, newname in zip(residue_list, new_names):
            modifications.exchange_base(resi, newname)


@toplevel_function
def extend_helix(model,  anchor5_id,  anchor3_id,  new_sequence=None):
    """*extend_helix(model, anchor5_id, anchor3_id, new_sequence=None)*

Makes a helix longer by adding A-type RNA Watson-Crick base pairs. 
A helix of any size may be added to the model.
A helix may be added at the end or inside an already existing helix.

:Arguments:
    * Template or RnaModel
    * Identifier of the residue on the 5' side to which helix will be added
    * Identifier of the residue on the 3' side to which helix will be added
    * Helix sequence as Sequence obiect e.g. Sequence('AAAA_UUUU')
    """
    model = validate_model(model)
    anchor5_id = validate_resnum(anchor5_id)
    anchor3_id = validate_resnum(anchor3_id)
    if new_sequence: new_sequence = validate_seq(new_sequence)

    hb = HelixFragmentBuilder()
    frag = hb.build_fragment(anchor5= model[anchor5_id], \
                    anchor3=model[anchor3_id], sequence=new_sequence, \
                    model=model)
    model.insert_fragment(frag)


@toplevel_function
def find_clashes(residues):
    """*find_clashes(residues)*
    
Detects interatomic clashes (overlapping atom VdW radii)
between atoms from two distinct residues. As input,
a template, a model or a list of residues is taken.

Returns a list of clashing residue pairs. If no clashes
are detected, the list is empty.

:Arguments:
    * list of residues or ModernaStructure object
    """
    residues = validate_resi_list(residues)
    
    # KR: result should be written to logfile
    cr = ClashRecognizer()
    clashes=cr.find_clashes_in_residues(residues)
    if clashes: log.write_message('Checking clashes: some residues are clashing %s' %str(clashes))
    else: log.write_message('Checking clashes: there are no clashes')
    return clashes



@toplevel_function
def find_fragment(model, res5_number, res3_number, sequence, number_of_candidates=NUMBER_OF_FRAGMENT_CANDIDATES, secstruc=None): 
    """*find_fragment(model, res5_number, res3_number, sequence, number_of_candidates=20, secstruc=None)*

Searches fragment candidate structures that fit between two given residues. 
Produces a list of fragments candidate objects. The fragment search sorts candidates 
according to the geometrical fit of their backbones.

As an optional argument, a secondary structure can be specified in dot-bracket format. 
This restricts the found fragment candidates to those having exactly the same Watson-Crick pairs.

:Arguments:
    * RnaModel object
    * number of the residue at the 5' end of the fragment
    * number of the residue at the 3' end of the fragment
    * sequence of the fragment to be inserted (not including the two residues above)
    * number of candidate fragments to be returned (default=20)
    * secondary structure of the fragment in dot-bracket format (optional).
    """
    model = validate_model(model)
    res5_number = validate_resnum(res5_number)
    res3_number = validate_resnum(res3_number)
    sequence = validate_seq(sequence)
    if secstruc: secstruc = validate_secstruc(secstruc)
    
    r5 = model[res5_number]
    r3 = model[res3_number]
    return model.find_fragment_candidates(r5, r3, sequence, number_of_candidates, secstruc=secstruc)


@toplevel_function
def find_modifications(structure):
    """*find_modifications(structure)*
    
Reports all modified nucleotides that occur in a structure
loaded or created with ModeRNA (models, templates).

This command returns all modified residues as a Python dictionary.
    {'1':<res>, '6': <res>, '6a': <res>, ...}

:Arguments:
    * RnaModel or Template object
    """
    structure = validate_structure(structure)
    return structure.get_modified_residues()

######################################  FIX BACKBONE ######################################

@toplevel_function
def fix_backbone(model, resi1_number=None, resi2_number=None):
    """*fix_backbone(structure)*

Fixes interrupted backbones between adjacent residues.

This function either examines the backbone of all nt in a structure (and repairs them if broken), 
or only does this for the two residues specified. For repairing, first only the phosphate and adjacent oxygens are remodeled. 
In case this fails, the FCCD Loop Closing Algorithm (Boomsma/Hamelryck 2005) is used for the entire O3'..C4' segment.
If the two residues are far from each other or oriented in a weird angle, the result will be awkward.
    
:Arguments:
    * RnaModel object
    * residue number on the 5' side of the broken backbone (optional).
    * residue number on the 3' side of the broken backbone (optional).
    """
    model = validate_model(model)
    if resi1_number: resi1_number = validate_resnum(resi1_number)
    if resi2_number: resi2_number = validate_resnum(resi2_number)
    
    if resi1_number and resi2_number:
        model.fix_backbone_between_resis(model[resi1_number], model[resi2_number])
    else:
        model.fix_backbone()


##################################### BASE PAIRS #######################################

@toplevel_function
def get_base_pairs(structure):
    """*get_base_pairs(structure)*
    
Returns a dictionary where keys are residues numbers and values are tuples which each one contains:
- residue number interacting with the key-residue
- interaction type

:Arguments:
    * ModernaStructure object (template or model).
    """
    structure = validate_structure(structure)
    return structure.get_base_pairs()


@toplevel_function
def get_sequence(structure):
    """*get_sequence(structure)*

Retrieves the one-letter-sequence from the coordinates of a structure (template or model). In the sequence, standard RNA bases are denoted by 
upper case letters, DNA bases by lowercase. For many modifications, one-letter ASCII abbreviations exist (according to McCloskey). 
All nucleotides that do not have a one-letter abbreviation are represented by the x letter.
Nucleotides that cannot be recognized (e.g. base missing) are represented by the '.' character.

For determining the sequence, the residues are processed according to their numbers. 
If an unusually long bond is found anywhere in the backbone between two bases, 
an additional "_" symbol is inserted in the sequence to mark this discontinuity.

Also see http://www.genesilico.pl/modomics

:Arguments:
    * structure - a Template or RnaModel object
    """
    # KR: sequence should be written to logfile.
    # KR: discontinuities should get an extra logfile message.
    # MM: is it good place for that or should be sone in get_sequence in ModernaStructure?
    structure = validate_structure(structure)
    seq = structure.get_sequence()
    log.write_message('Checking sequence: \n%s\n' %seq.seq_with_modifications)
    if '_' in seq.seq_with_modifications:
        log.write_message('SEQUENCE IS DISCONTINUOUS !!!\n')
    return seq


@toplevel_function
def get_secstruc(structure):
    """*get_secstruc(structure)*

Retrieves the dot-bracket secondary structure from the coordinates of a structure (template or model). 
Base pairs (only Watson-Crick pairings AU and GC, and GU Wobble pairs) \ are indicated by round brackets,
all other residues by dots.

:Arguments:
    * structure - a Template or RNAModel object
"""
    secstruc = structure.get_secstruc()
    log.write_message('Secondary structure: \n%s\n' %secstruc)
    return secstruc


@toplevel_function
def get_stacking(structure):
    """*get_stacking(structure)*
    
Returns a list of stacking interactions found in the structure, 
using the types <<, >>, <>, >< as defined in
MC-Annotate paper (JMB 2001, 308, p.919ff).

:Arguments:
    * ModernaStructure object (template or model).
    """
    structure = validate_structure(structure)
    scalc = StackingCalculator()
    return [(s.resi1.identifier, s.resi2.identifier, s.type) for s in scalc.get_stacking(structure)]


#@toplevel_function
def insert_two_strand_fragment(model, anchor5_id, anchor3_id, anchor5_upper_id, anchor3_upper_id, \
                 frag_anchor5_id,  frag_anchor3_id, frag_anchor5_upper_id, frag_anchor3_upper_id,  fragment_file,  chain_name='A'):
    """*insert_2D_fragment(model, anchor5_id, anchor3_id, anchor5_upper_id, anchor3_upper_id, \
                 frag_anchor5_id,  frag_anchor3_id, frag_anchor5_upper_id, frag_anchor3_upper_id,  fragment_file,  chain_name='A')*

Inserts a given 2D fragment between two pairs of indicated residues.
Requires to indicate connecting residues from both model (4) and fragment (4).

:Arguments:
    * model - RNAModel object
    * anchor5_id - anchor residue in the model on the 5' side of the fragment.
    * anchor3_id - anchor residue in the model on the 3' side of the fragment.
    * anchor5_upper_id - anchor residue in the model on the 5' side of the middle part of the fragment.
    * anchor3_upper_id - anchor residue in the model on the 3' side of the middle part of the fragment.
    * frag_anchor5_id - anchor residue in the fragment on the 5' side of the fragment.
    * frag_anchor3_id - anchor residue in the fragment on the 3' side of the fragment.
    * frag_anchor5_upper_id - anchor residue in the fragment on the 5' side of the middle part of the fragment.
    * frag_anchor3_upper_id - anchor residue in the fragment on the 3' side of the middle part of the fragment.
    * fragment_file - name of the fragment PDB file
    * chain_name - chain ID (default: 'A')
    """   
    model = validate_model(model)
    anchor5_id = validate_resnum(anchor5_id)
    anchor3_id = validate_resnum(anchor3_id)
    anchor5_upper_id = validate_resnum(anchor5_upper_id)
    anchor3_upper_id = validate_resnum(anchor3_upper_id)
    frag_anchor5_id = validate_resnum(frag_anchor5_id)
    frag_anchor3_id = validate_resnum(frag_anchor3_id)
    frag_anchor5_upper_id = validate_resnum(frag_anchor5_upper_id)
    frag_anchor3_upper_id = validate_resnum(frag_anchor3_upper_id)
    fragment_file = validate_filename(fragment_file)
    
    struc = ModernaStructure(data_type='file', data=fragment_file, chain_name=chain_name)
    mf = ModernaFragment2D(struc, \
                anchor5=model[anchor5_id], anchor3=model[anchor3_id], \
                anchor5_upper=model[anchor5_upper_id], \
                anchor3_upper=model[anchor3_upper_id], \
                frag5_upper=struc[frag_anchor5_upper_id], \
                frag3_upper=struc[frag_anchor3_upper_id], \
                model=model)
    model.insert_fragment(mf)


@toplevel_function
def insert_fragment(model, fragment): 
    """*insert_fragment(model, fragment)*

Inserts a fragment object into the model. 
The model should be the same object, from which the 
two anchor residues were defined when the fragment was created. 
(See also documentation of the create_fragment/find_fragment functions.)

:Arguments:
    * RnaModel object
    * ModernaFragment or LirHit object
    """
    model = validate_model(model)
    fragment = validate_fragment(fragment)
    
    if type(fragment) == LirHit: model.insert_lir_candidate(fragment)
    else: model.insert_fragment(fragment)


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
def load_model(data, chain_name='A',data_type='file',template=None,alignment=None):
    """*load_model(file_path, chain_name='A', data_type='file', template=None, alignment=None)*
    
Loads a structure model that has been built previously,
or any PDB struture that is then to be modified.
Produces a RnaModel object that can be saved in a variable.
Each model in ModeRNA contains only one chain.
Multi-chain structures can be modeled by using more than one
Template/Alignment/RNAModel at a time.

By default, RNAModels are created by reading files.
They can also be created from BioPython PDB objects
(precisely Bio.PDB.Structure.Structure and Bio.PDB.Chain.Chain objects)

:Arguments:
    * path+filename of a PDB structure file; or a Structure or Chain object from BioPython.
    * chain id (by default 'A')
    * data type ('file' by default; if set to 'structure' or 'chain', Bio.PDB objects are read)
    * Template object to be used for this model (optional)
    * Alignment object to be used for this model (optional)
    """
    if data_type == 'file': data = validate_filename(data)
    elif data_type == 'residues': data = validate_resi_list(data)
    if template: template = validate_template(template)
    if alignment: alignment = validate_alignment(alignment)
    return RnaModel(template=template, alignment=alignment, model_chain_name=chain_name, data_type=data_type, data=data)


@toplevel_function
def load_template(file_path, chain_name='A'):
    """*load_template(file_path, chain_name='A')*
    
Loads a template structure from a PDB file.
Produces a Template object that can be saved in a variable.

Each template in ModeRNA has only one chain. By default, the chain
with id 'A' is loaded. Another chain id can be specified optionally

:Arguments:    
    * path+filename of a PDB structure file
    * chain id (by default 'A')
    """
    file_path = validate_filename(file_path)
    
    t = Template(file_path, 'file', chain_name)
    log.write_message('Template loaded from %s.\nTemplate sequence:\n%s\n' %(file_path, t.get_sequence().seq_with_modifications))
    return t 
    # KR: Why is the order of arguments in Template different from that
    #     in ModernaStructure ('file',name)?


@toplevel_function
def match_alignment_with_model(alignment, model):
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


##################################### REMOVING MODIFICATION #####################################

@toplevel_function
def remove_all_modifications(model):
    """*remove_all_modifications(model)*

Removes all base modifications from a given model. The nucleotides are
transformed into standard bases from which the modifications originated.
     
:Arguments:
    * RnaModel object
    """
    model = validate_model(model)
    model.remove_all_modifications()


@toplevel_function
def remove_mismatching_modifications(template, alignment, model):
    """*remove_mismatching_modifications(template, alignment, model)*

Removes all nucleotide modifications that occur in the template, 
and are not present in the target sequence. 
The nucleotides are transformed into standard bases from which the modifications originated.
   
:Arguments:
    * Template object
    * Alignment object
    * RnaModel object
    """
    template = validate_template(template)
    alignment = validate_alignment(alignment)
    model = validate_model(model)
    
    model.template = template
    model.alignment = alignment
    model.recipe = RecipeMaker(alignment).recipe
    model.remove_all_modifications_copy()


@toplevel_function
def remove_modification(residue, model=None, new_number=None):
    """*remove_modification(residue, model=None, new_number=None)*
    
Removes base modifications from a single residue. The nucleotide is
transformed into the standard base from which the modification originated.
A RnaModel can be given optionally. If this is given, the
modified residue is subsequently copied to it.

Note that desoxynucleotides count as modified bases as well.
  
:Arguments:
    * residue from a Template or RnaModel object
    * RnaModel object (optional)
    * new residue number after copying (optional)
    """
    residue = validate_resi(residue)
    if model: model = validate_model(model)
    if new_number: new_number = validate_resnum(new_number)
    
    if model: 
        model.remove_one_modification_copy(residue, new_number)
    else:
        modifications.remove_modification(residue)


@toplevel_function
def renumber_chain(struc,  start_identifier='1'):
    """*renumber_chain(struc, start_id)*

Changes numeration of residues in the chain.
Starts with the given identifier. The new numeration is continuous.

:Arguments:
    * RnaModel or ModernaStructure object
    * identifier for first residue
    """
    struc = validate_structure(struc)
    start_identifier = validate_resnum(start_identifier)
    struc.renumber_chain(start_identifier)
    

@toplevel_function
def rotate_chi(residue, angle=90):
    """*rotate_chi(residue, angle=90)*
    
Rotates a base in a given nucleotide around the glycosidic bond
(the chi torsion angle, between C1' and N1/N9 atoms).

:Arguments:
    * residue (from a template or model)
    * angle in degrees
    """
    residue = validate_resi(residue)
    rc(residue, angle)
    
    
@toplevel_function
def shrink_helix(model,  anchor5_id,  anchor3_id,  anchor5_upper_id,  anchor3_upper_id):
    """*shrink_helix(model,  anchor5_id,  anchor3_id,  anchor5_upper_id,  anchor3_upper_id)*

Makes helix shorter - cuts out a helix fragment between given residues from the model.

:Arguments:
    * Template or RnaModel
    * anchor5 - residue id on the 5' side belonging to the first base pair that will remain
    * anchor3 - residue id on the 3' side belonging to the first base pair that will remain
    * anchor5_upper - residue id on the 5' side belonging to the pair that will end up next to the first.
    * anchor3_upper - residue id on the 3' side belonging to the pair that will end up next to the first.
    """
    model = validate_model(model)
    anchor5_id = validate_resnum(anchor5_id)
    anchor3_id = validate_resnum(anchor3_id)
    anchor5_upper_id = validate_resnum(anchor5_upper_id)
    anchor3_upper_id = validate_resnum(anchor3_upper_id)
    
    hbuilder = HelixBuilder()
    helix = hbuilder.build_helix(Sequence('AA_UU'))
    fr = ModernaFragment2D(helix, anchor5=model[anchor5_id], \
                           anchor3=model[anchor3_id], \
                           anchor5_upper=model[anchor5_upper_id], \
                           anchor3_upper=model[anchor3_upper_id], \
                           frag5_upper = helix['2'], frag3_upper = helix['401'], \
                           model=model)
    model.insert_fragment(fr)

##################################### WRITING RESULTS ##################################

def write_logfile(logfile_name='moderna.log'):
    """*write_logfile(logfile_name='moderna.log')*
    
Writes a text file with all log messages that ModeRNA generated
up to that moment. The log messages in the memory are cleared afterwards.
When ModeRNA terminates, any remaining messages will be written to 'moderna.log'

:Arguments:
    * name of the log file (optional; by default moderna.log)
     """
    logfile_name = validate_filename(logfile_name)

    try:
        log.set_filename(logfile_name)
        log.write_logfile()
    except IOError:
        print(('problem writing log file to "%s"'%logfile_name))
    

@toplevel_function
def write_fragment_candidates(fragment_candidates, output_directory='fragment_candidates', with_model=True, fragment_log_file=True): 
    """*write_fragment_candidates(fragment_candidates, output_directory='fragment_candidates')*

Writes a list of fragment candidates to a set of PDB files. 
The candidates are numbered according to the geometrical fit of their backbones.

:Arguments:
    * fragment candidates list (obtained by the find_fragment command)
    * output directory name
    """
    fragment_candidates = validate_frag_candidate_list(fragment_candidates)
    output_directory = validate_path(output_directory)
    
    fragment_candidates.write_fragment_candidates(output_directory, True, with_model, False,  fragment_log_file)
    log.write_message('Fragment candidates written to %s.' %output_directory)


@toplevel_function
def write_model(model, pdb_file_name='moderna_model.pdb'):
    """*write_model(model, pdb_file_name='moderna_model.pdb')*
    
Writes a model to a PDB file. The residues in the file are sorted.
All residues keep their numbers as last assigned.

:Arguments:    
    * Structure object (model or template)
    * name of the PDB file (optional; by default moderna_model.pdb)
    """
    model = validate_structure(model)
    pdb_file_name = validate_filename(pdb_file_name)

    if model.__class__ == RnaModel: model.refine_model()
    model.write_pdb_file(pdb_file_name)
    log.write_message('Model written to %s'%pdb_file_name)
    
    
@toplevel_function
def write_secstruc(struct, file_name='secstruc.vienna'):
    """*write_secstruc(model, file_name='secstruc.vienna')*
    
Writes secondary struture to a vienna file. 

:Arguments:    
    * Structure object (model or template)
    * name of the vienna file (optional; by default secstruc.vienna)
    """
    struct = validate_structure(struct)
    file_name = validate_filename(file_name)

    struct.write_secstruc(file_name)
    log.write_message('Secondary structure written to %s'%file_name)
    

