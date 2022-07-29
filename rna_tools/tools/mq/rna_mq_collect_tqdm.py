#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""mqaprna.py - a script for running all wrapers on each PDB file in a specified directory
saves results to a CSV file.

ss_agreement is ...

The code is full of # hack and tricks.

.. warning:: Uses global variables

Install: 

    csvsort

Cmd::

     # find . -iname 'FARFAR2*.csv' -exec cat {} + > FARFAR2_hires.csv
     $ rna_mq_collect.py -t FARFAR2_hires -m 4 -f -o FARFAR2_hires.csv -l all.txt x.pdb
     # fake x.pdb when -l is used, -l gets a list of files
     x.pdb
     y.pdb
     z.pdb

88% (49329 of 55689) |############### | Elapsed Time: 0:45:23 ETA:  2 days, 18:42:16

"""
MP_VERBOSE = False
DEBUG_MODE = False

################################################################################
import sys
#sys.path.insert(0, "/Users/magnus/work/src/rna-tools/rna_tools/tools/mq/")  # ugly!
import progressbar
# import mqaprna_score as mqs
import time
import os
import copy
from csvsort import csvsort

import rna_tools.tools.mq.lib.shellgraphics.shellgraphics as sg
sg.color_mode = False
from rna_tools.tools.mq.lib.timex import timex
#import rna_tools.tools.mq.mqaprna_config as Config
import rna_tools.rna_tools_config as Config
################################################################################

import rna_tools
__version__ = rna_tools.__version__

import os
import sys
DIRNAME = os.path.dirname(__file__)
PARENT_DIRNAME = os.path.abspath(os.path.join(DIRNAME, os.path.pardir))
sys.path.append(DIRNAME)
sys.path.append(PARENT_DIRNAME)
import csv
import imp

from optparse import OptionParser, OptionGroup
from ctypes import c_int

from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')


#import lib.rmsd_calc.rmsd_calc as rmsd_calc
from multiprocessing import Pool, Lock, Value

try:
    from wrappers.mqap_score.mqap_score import MqapScore
except ImportError:
    pass

# super-verbose logging
MP_VERBOSE = 0
if MP_VERBOSE:
    import multiprocessing
    logger = multiprocessing.log_to_stderr()
    logger.setLevel(multiprocessing.SUBDEBUG)

# create wrappers for all the methods
MODULES = {}
for m in Config.METHOD_LIST:
    if m.find('_') > -1:
        m,n = m.split('_')
    wrapper_path = os.path.join(Config.WRAPPERS_PATH, m, m + '.py')
    module = imp.load_source(m, wrapper_path)
    MODULES[m] = module

# global variable
c = 0
methods = Config.METHOD_LIST
cleanup = True

counter = Value(c_int)
counter_lock = Lock()

# ['farna_rna_base_axis', 'farna_rna_backbone_backbone', 'farna_rna_base_stack_axis', 'farna_rna_base_stagger', 'farna_rna_base_stack', 'farna_rna_base_pair', 'farna_rna_repulsive', 'farna_rna_vdw', 'farna_rna_base_backbone', 'farna_score_lowres', 'farna_rna_data_backbone', 'farna_linear_chainbreak', 'farna_rna_rg', 'farna_atom_pair_constraint'],

steps = '0' #
attributes = {
    'QRNA' : [ 'qrna_' + steps + '_electro', 'qrna_' + steps ],
    #'RASP' : [ 'rasp_all_pdb_energy', 'rasp_all_no_contacts', 'rasp_all_norm_energy', 'rasp_all_mean_energy', 'rasp_all_sd_energy', 'rasp_all_zscore']
    'RASP' : ['rasp_c3_pdb_energy', 'rasp_c3_no_contacts', 'rasp_c3_norm_energy', 'rasp_c3_mean_energy', 'rasp_c3_sd_energy', 'rasp_c3_zscore', 'rasp_bb_pdb_energy', 'rasp_bb_no_contacts', 'rasp_bb_norm_energy', 'rasp_bb_mean_energy', 'rasp_bb_sd_energy', 'rasp_bb_zscore', 'rasp_bbr_pdb_energy', 'rasp_bbr_no_contacts', 'rasp_bbr_norm_energy', 'rasp_bbr_mean_energy', 'rasp_bbr_sd_energy', 'rasp_bbr_zscore', 'rasp_all_pdb_energy', 'rasp_all_no_contacts', 'rasp_all_norm_energy', 'rasp_all_mean_energy', 'rasp_all_sd_energy', 'rasp_all_zscore'],
    'FARNA_hires' : ['farna_score_hires', 'farna_fa_atr', 'farna_fa_rep', 'farna_fa_intra_rep',
                                       'farna_lk_nonpolar',
                                       'farna_fa_elec_rna_phos_phos',
                                       'farna_ch_bond',
                                       'farna_rna_torsion',
                                       'farna_rna_sugar_close',
                                       'farna_hbond_sr_bb_sc',
                                       'farna_hbond_lr_bb_sc',
                                       'farna_hbond_sc',
                                       'farna_geom_sol',
                                       'farna_atom_pair_constraint_hires',
                                       'farna_linear_chainbreak_hires'],

    'SimRNA_0' : ['simrna_steps', 'simrna_total_energy', 'simrna_base_base', 'simrna_short_stacking', 'simrna_base_backbone',  'simrna_local_geometry', 'simrna_bonds_dist_cp', 'simrna_bonds_dist_pc', 'simrna_flat_angles_cpc', 'simrna_flat_angles_pcp', 'simrna_tors_eta_theta', 'simrna_sphere_penalty', 'simrna_chain_energy'],
    'RNAkb' : ['rnakb_bond', 'rnakb_angle', 'rnakb_proper_dih', 'rnakb_improper_dih', 'rnakb_lj14', 'rnakb_coulomb14', 'rnakb_lj_sr', 'rnakb_coulomb_sr',
               'rnakb_potential', 'rnakb_kinetic_en', 'rnakb_total_energy'],
    'RNAkb_all' : ['rnakb_bond_all', 'rnakb_angle_all', 'rnakb_proper_dih_all', 'rnakb_improper_dih_all', 'rnakb_lj14_all', 'rnakb_coulomb14_all', 'rnakb_lj_sr_all', 'rnakb_coulomb_sr_all',
               'rnakb_potential_all', 'rnakb_kinetic_en_all', 'rnakb_total_energy_all'],

    'RNAscore' : ['x3rnascore'],
    'AnalyzeGeometry' : ['analyze_geometry'],
    'SSAgreement' : ['ss_disagreement'],
    'ClashScore' : ['clash_score'],
    'Ernwin_1' : [ 'ernwin_1' ],
    'Ernwin_1k' : [ 'ernwin_1k' ],
    'eSCORE' : ['escore'],
    'RNA3DCNN' : ['rna3dcnn'],
    'Dfire' : ['dfire'],

    'FARNA':                           ['farna_score_lowres',
                                       'farna_rna_data_backbone',
                                       'farna_rna_vdw',
                                       'farna_rna_base_backbone',
                                       'farna_rna_backbone_backbone',
                                       'farna_rna_repulsive',
                                       'farna_rna_base_pair',
                                       'farna_rna_base_axis',
                                       'farna_rna_base_stagger',
                                       'farna_rna_base_stack',
                                       'farna_rna_base_stack_axis',
                                       'farna_rna_rg',
                                       'farna_atom_pair_constraint',
                                       'farna_linear_chainbreak'],
     
    'FARFAR2_hires': 'ff2_score_hires,ff2_fa_atr,ff2_fa_rep,ff2_fa_intra_rep,ff2_lk_nonpolar,ff2_fa_elec_rna_phos_phos,ff2_rna_torsion,ff2_suiteness_bonus,ff2_rna_sugar_close,ff2_fa_stack,ff2_stack_elec,ff2_geom_sol_fast,ff2_bond_sr_bb_sc,ff2_hbond_lr_bb_sc,ff2_hbond_sc,ff2_ref,ff2_free_suite,ff2_free_2HOprime,ff2_intermol,ff2_other_pose,ff2_loop_close,ff2_linear_chainbreak_hires'.split(','),
    
    #'SimRNA_0' : ['', 'simrna', '', '', '',  '', '', '', '', '', '', '', ''],
    'rmsd_all': ['rmsd_all'],
}


def single_run(lst):
    """Start a mqaprna run for a given file
    with all methods (according to config file).

    [!] Use global cleanup = False to block cleaning up

    .. warning:: The function uses global variable.
    """
    filename, c, verbose, methods, opt, ref_seq = lst
    all_results = {}

    for m in methods:

            arguments = ''
            #if DEBUG_MODE: print 'method', m, arguments
            mfull = m

            if verbose: print(m + '...') # show method 'eSCORE...'

            if m == 'FARNA':
                mfull = m
                arguments = [filename] + [False]

            if m == 'FARNA_hires':
                m = 'FARNA'
                mfull = 'FARNA_hires'
                arguments = [filename] + [True]

            if m == 'FARFAR2':
                m = 'FARFAR2'
                mfull = 'FARFAR2'
                arguments = [filename] + [False]

            if m == 'FARFAR2_hires':
                m = 'FARFAR2'
                mfull = 'FARFAR2_hires'
                arguments = [filename] + [True]
                
            if m == 'RNAkb_all':
                m = 'RNAkb'
                mfull = 'RNAkb_all'
                arguments = [filename] + ['aa']

            if m.find('_') > -1:
                m, n = m.split('_')
                n = n.replace('n', '') # n_XXX
                n = n.replace('k', '000')
                n = n.replace('m', '000000')
                arguments = [filename] + [n]

            if not arguments:
                arguments = [filename] + Config.WRAPPER_OPTIONS[m]

            if m == 'escore':
                m = 'eSCORE'
            wrapper = getattr(MODULES[m], m)()#verbose) # ref_seq, ref_ss, verbose)  # for all wrappers but SSAgrement '','' is OK

            if m == 'NAST_pyro':
                lock.acquire()

            if DEBUG_MODE:
                result = wrapper.run(*arguments)
                if verbose: print(m, result) # ClashScore 12.256669
                all_results[mfull] = result
                if cleanup: wrapper.cleanup()
            else:
                try:
                    result = wrapper.run(*arguments)
                    all_results[mfull] = result
                    if cleanup: wrapper.cleanup()
                except:
                    all_results[mfull] = 'error'
                    if cleanup: wrapper.cleanup()

            # {'ClashScore': 12.256669}
            # {'ClashScore': 12.256669, 'AnalyzeGeometry': 32.5581}
            # {'ClashScore': 12.256669, 'AnalyzeGeometry': 32.5581, 'FARNA': '-20.008,-2.739,-13.175,-77.67,-10.652,-158.51,9.547,8.39,-16.246,-263.281,0.0,0.0,17.782,0.0'}
            #if verbose: print 'all_results:', all_results # this every each method showed

            if m == 'NAST_pyro':
                lock.release()

    # get rmsd
    if opt.native_pdb_filename:
        rmsd = rmsd_calc.get_rmsd(opt.native_pdb_filename,
                                     filename)
        all_results['rmsd'] = rmsd
        methods = methods + ['rmsd']
    else:
        methods = methods

    # length
    length = len(ref_seq)
    all_results['length'] = length

    if opt.mqapscore:
        # meta-score
        ms =  MqapScore(all_results)
        mqap_score = ms.get_score()
        methods = methods + ['SCORE']
        all_results['SCORE'] = mqap_score

    if True:
        # lock.acquire()

        global counter_lock
        #with counter_lock:
        counter.value += 1

        if counter.value != 1:
            # @todo does not work
            #sys.stdout.write('\033[F')
            #sys.stdout.write('\033[F')
            pass

        #results = [str(round(all_results[mfull],2)).strip().rjust(9) for m in methods]

        results_str = str(all_results) # "{'AnalyzeGeometry': 0.0, 'eSCORE': 0.10661, 'FARNA': ['-2.411', '0.0', '0.0', '-9.672', '0.0', '-25.678', '0.0', '1.061', '0.0', '-32.098', '0.0', '0.0', '4.601', '0.0'], 'ClashScore': 36.458333, 'length': 0, 'SimRNA_0': ['0', '67.345305', '-37.428', '-23.073', '0.248', '104.524975', '87.955', '9.938', '5.669', '1.089', '-0.126', '', '67.345305'], 'FARNA_hires': ['0.0', '-13.107', '-0.711', '0.0', '5.22', '2.734', '-30.044', '0.223', '-10.511', '-0.173', '-4.719', '1.143', '0.0', '14.371', '9.358'], 'RNAscore': 8.11007, 'RASP': ['-0.1382', '15', '-0.00921333', '-0.0845115', '0.454033', '-0.118248', '-277.666', '949', '-0.292588', '-273.37', '2.51163', '-1.71042', '-584.451', '2144', '-0.272598', '-564.143', '5.77609', '-3.51588', '-1616.08', '6700', '-0.241206', '0', '0', '0'], 'RNAkb': -1}"

        results = [all_results[mfull] for m in methods]
        # progress bar
        #sys.stdout.write('\r')
        #sys.stdout.flush()
        #sys.stdout.write('\r' + ' ' * 110 + '\r' + filename.split(os.sep)[-1].ljust(50) + ' ' + ' '.join(results))

        ########### line with resluts ######################
        #bar.update(counter.value)
        ## my old progress bar here:
        # print(sg.pprogress_line(counter.value, filename_length, ''))# ,
        ## print results, use --verbose now
        if verbose:
            print(filename.split(os.sep)[-1].ljust(20) + ' ' + results_str)
        
        ## [          ]   1 7.14 % 14 3_solution_1.pdb     {'AnalyzeGeometry': 0.0, 'eSCORE': 1.70264, 'FARNA': ['-31.498', '-11.589', '-32.7', '-123.708', '-25.514', '-271.337', '33.563', '2.957', '-36.699', '-471.864', '0.0', '0.0', '24.659', '0.0'], 'ClashScore': 2.201835, 'length': 0, 'SimRNA_0': ['0', '-1016.539381', '-599.475', '-223.162', '-3.935', '-413.129576', '-65.066', '-71.505', '-68.947', '-45.989', '-161.622', '', '-1016.539381'], 'FARNA_hires': ['0.0', '-541.374', '-0.59', '0.0', '1.85', '8.12', '-433.113', '17.811', '-229.203', '3.074', '-140.106', '13.875', '-17.245', '226.762', '7.39'], 'RNAscore': 26.7066, 'RASP': ['-9.3599', '987', '-0.00948318', '8.16333', '3.95157', '-4.4345', '-7976.88', '60547', '-0.131747', '-7274.73', '52.7448', '-13.3123', '-17537.5', '138719', '-0.126424', '-15578.4', '106.602', '-18.3777', '-34270.8', '483436', '-0.07089', '0', '0', '0'], 'RNAkb': -0.019507621989000006}

        #sys.stdout.flush()
        #print
        #sys.stdout.write(sg.pprogress_line(counter.value, filename_length))
        #print sg.pprogress_line(counter.value, filename_length)
        #sys.stdout.flush()

        ## for graphics debugging
        #import time
        #time.sleep(1)

        #format_line([filename.split(os.sep)[-1] + [all_results[m] for m in methods]])  # @todo Nice print with ShellGraphics
        cells = [c, filename.split(os.sep)[-1]] # add id 
        for m in methods:
            if type(all_results[m]) == list:
                cells.extend(all_results[m])
            else:
                cells.append(all_results[m])
        #csv_writer.writerow(cells)
        return cells
        #print 'mqaprna::filename: %i %s' % (counter.value, filename)
        #csv_file.flush()
        #lock.release()

    # hack
    try:
        methods.remove('SCORE')
    except ValueError:
        pass

    try:
        methods.remove('rmsd')
    except ValueError:
        pass


def option_parser():
    """Get options or show usage msg.
    """
    description = ''
    version = __version__
    usage = '\t%prog [-m <number_processes>] [-n <native_pdb_filename>] [-s <seq_ss_filename>] [-g <ignore_pdb_filename>] \ \n\t -o <output csv> <dir/*> # [!] no .csv! the file will get version of mqaprna \n\t' + __version__
    parser = OptionParser(description=__doc__,
                              version=version,
                              usage=usage)

    parser.add_option("-q", "--mQapscore",
                     action="store_true", default=False, dest="mqapscore", help="calculate mqapscore")

    parser.add_option("-v", "--verbose",
                     action="store_true", default=False, dest="verbose", help="verbose")

    parser.add_option("-f", "--no-filename-version",
                     action="store_true", default=False, dest="no_filename_version", help="don't add version of tool to csv filename")


    parser.add_option("-n", "--native_pdb_filename",
                     action="store", type="string", dest="native_pdb_filename", help="native structure in PDB format to calculate RMSD")

    parser.add_option("-m", "--multiprocessing",
                     action="store", type="int", dest="number_processes", default=1,
                      help="set a number of processes, default=8, 0 is no multiprocessing")

    group2 = OptionGroup(parser, "Ignore pdbs, don't have empty lines here! Example",
                        """1xjrA_output3-000142_AA.pdb
                         1xjrA_output3-000208_AA.pdb
                         1xjrA_output3-000166_AA.pdb""")

    group2.add_option("-g", "--ignore-pdbs",
                  action="store", type="string", dest="ignore_pdb_filename")

    group = OptionGroup(parser, "Seq-SS. Example",
                        """>1xjrA
                        GAGUUCACCGAGGCCACGCGGAGUACGAUCGAGGGUACAGUGAAUU
                        .(((((((...((((.((((.....))..))..))).).)))))))""")

    group.add_option("-t", "--methods",
                  action="store", type="string", dest="methods", help=', '.join(['RASP',  'SimRNA', 'AnalyzeGeometry','FARNA', 'QRNA', 'NAST_pyro',
                                                                                 'radius_of_gyration', 'SSAgreement', 'ClashScore', 'RNAkb',
                                                                                 'RNAkb_all', 'FARNA_hires', 'FARNA', 'FARFAR2',
                                                                                 'FARFAR2_hires', 'Dfire', 'RNA3DCNN', 'eSCORE']))

    group.add_option("-s", "--seq-ss",
                  action="store", type="string", dest="seq_ss_filename", help="")

    group.add_option("-o", "--output",
                  action="store", type="string", dest="output", help="output csv file")

    group.add_option("-l", "--list-of-files",
                  action="store", type="string", dest="list_of_files", help="list of files")


    parser.add_option_group(group)
    parser.add_option_group(group2)

    (opt, arguments) = parser.parse_args()

    arguments = [f for f in arguments if f.endswith('.pdb')]

    if len(arguments) == 0:
        parser.print_help()
        print('\n   Curr methods: ', ','.join(methods), end=' ')
        sys.exit(1)

    return arguments, opt


class RunAllDirectory():
    """Class for running wrappers for all files in a directory
    """
    def __init__(self):
        pass

    def run(self, filenames, csv_path, opt):
        """Open csv (with appropriate headers), run methods, print & save csv

        There are two modes of execution:
         * multiprocessing
         * single

        .. warning:: Works on global variables: ref_seq, ref_ss, methods, lock, c
        """
        global ref_seq, ref_ss, verbose, methods, lock, c

        if opt.seq_ss_filename:
            pdb_id, ref_seq, ref_ss = [x.strip() for x in open(opt.seq_ss_filename).read().strip().split('\n')]
            #sg.phr_text('FASTA SEQ/SS')
            sg.poptions({'AnalyzeGeometry': True, 'SSAgreement' : True})
            sg.poption('pdb_id', pdb_id)
            sg.poption('ref_seq', ref_seq)
            sg.poption('ref_ss', ref_ss)
        else:
            pdb_id, ref_seq, ref_ss = ['', '', '']
            sg.poptions({'SSAgreement' : True})
            # hack
            try: # if it's not on the list
                methods.remove('SSAgreement')
            except ValueError:
                pass

        verbose = opt.verbose

        global csv_file, csv_writer  # hack
        # csv open & add header
        csv_file = open(csv_path, 'a')
        csv_writer = csv.writer(csv_file,  delimiter=',')
        # make header
        headers = ['id', 'fn']
        for m in methods:
            headers += attributes[m]

        if opt.native_pdb_filename:
            headers += ['RMSDALL']
        if opt.mqapscore:
            headers += ['SCORE']
        csv_writer.writerow(headers)
        csv_file.flush()

        # remove ~ and remove .out
        for f in copy.copy(filenames):
            if f.endswith('~'):
                filenames.remove(f)
            if f.endswith('.out'):
                filenames.remove(f)
            if f.find('._')>-1:
                filenames.remove(f)

        files_to_ignore = []
        # or if not provided
        import glob
        opt.ignore_pdb_filename = glob.glob('*' + opt.methods + '*.csv')
        for f in opt.ignore_pdb_filename:  # do it for the list, that's nice!
            fn = open(f)
            for f in fn.read().strip().split('\n'):
                if 'error' in f:
                    continue  # don't add files with errors, so the program will be re-run for them
                # if there is an error, this will give error again quickly
                # but this solves when you kill the job, you get erros, but it's not rally errors
                # but stopped jobs
                if f.find('\t') > -1:
                    f = f.split('\t')[1] # id, fn
                if f.find(',') > -1:
                    f = f.split(',')[1] # id, fn
                files_to_ignore.append(os.path.basename(f))

        ## files to ignore
        print(' to ignore', len(files_to_ignore), files_to_ignore[:4])

        filenames = []
        for i, f in enumerate(input_files):
            # print(i, f)
            if '/_' in f:  # skip
                continue
            if os.path.basename(f) not in files_to_ignore:
                filenames.append(f)

        with open('_mq_to_run_.txt', 'w') as f:
            f.write('\n'.join(filenames))
        print(' save filenames to run to _mq_to_run_.txt')

        ## for fi in files_to_ignore:
        ##     for fn in copy.copy(filenames):
        ##         if os.path.basename(fn).startswith('._'):
        ##             filenames.remove(fn)
        ##         if os.path.basename(fn).startswith(fi.split('\t')[0]): # # hack,  @todo <- re could be used here!  to ignore ['fn,RASP,SimRNA,FARNA,NAST_pyro\r', '1ykv_1_ba_c.pdb,-0.104705,-504.468933,-306.245,122.7\r', '2esj_1_ba_c.pdb,-0.1522,-1,-266.217,46.7\r', '2quw_1_ba_c.pdb,-0.103789,-729.386726,-419.047,984.0\r
        ##             filenames.remove(fn)
        print(' files to analyze: %s' % len(filenames), filenames[:5])
        ## headers
        methods_to_print = copy.copy(methods)
        if opt.native_pdb_filename:
            methods_to_print += ['RMSDALL']
        if opt.mqapscore:
            methods_to_print += ['SCORE']

        ## if verbose: print ''.ljust(80), ''.join([m[:9].ljust(10) for m in methods_to_print]) ## print headers

        sg.phr()

        lock = Lock()
        
        counter.value = len(files_to_ignore)

        flist = []
        c  = 1
        # two running modes
        global filename_length
        filenames_length = len(filenames) + len(files_to_ignore)

        global bar
        bar = progressbar.ProgressBar(max_value=filenames_length)
        bar.update(len(files_to_ignore))

        fl = []
        for f in filenames:
            fl.append([f,filenames_length])

        lst = []
        for f in fl:
            # ['test/1xjrA_M1.pdb', 1, True, ['RASP']]
            lst.append([f[0], f[1], verbose, methods, opt, ref_seq])

        if int(opt.number_processes) > 1:
            pool = Pool(opt.number_processes)
            from tqdm.contrib.concurrent import process_map
            #pool.map(single_run, lst)
            outputs = process_map(single_run, lst, max_workers=2)
            pool.close()

            for cells in outputs:
                csv_writer.writerow(cells)
        else:
            for l in lst:
                single_run(l)

#main
if __name__ == '__main__':
    from icecream import ic
    import sys
    ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
    ic.configureOutput(prefix='> ')

    t = timex.Timex()
    t.start()

    arguments, opt = option_parser()

    # files
    input_files = arguments[:]
    if opt.list_of_files:
       for l in open(opt.list_of_files):
           input_files.append(l.strip())
    #ic(input_files)
    
    if not opt.methods:
       opt.methods = ','.join(Config.METHOD_LIST)

    if opt.no_filename_version:
        output_csv = opt.output
    else:
        import platform
        platform = platform.node()
        if opt.output:
            output_csv = opt.output.replace('.csv','') + '-' + __version__ + '-' + platform + '.csv'
        else:
            output_csv = opt.methods + '-' + __version__  + '-' + platform + '.csv'

    sg.pbanner_simply(os.path.basename(sys.argv[0]))

    try:
        rnakb_option = Config.WRAPPER_OPTIONS['RNAkb'][0]
    except KeyError:
        rnakb_option = None
    try:
        rasp_option = Config.WRAPPER_OPTIONS['RASP'][0]
    except KeyError:
        rasp_option = None

    if opt.methods:
        methods = [x.strip() for x in opt.methods.split(',')]

    print('ver:',  __version__ + '\n')
    print('start ', time.strftime("%Y-%m-%d %H:%M:%S"))
    
    opts = {
        'Input files': '#' + str(len(input_files)) + ' ' + str(input_files[:3]),
        'Multiprocessing': True if opt.number_processes > 1 else False,
        'Output csv': output_csv,
        'Seq ss fn': opt.seq_ss_filename,
        'Ignore pdb fn': opt.ignore_pdb_filename,
        'Native pdb': opt.native_pdb_filename,
        'RNAkb' : rnakb_option,
        'RASP' : rasp_option,
     #   'rmsd' : rmsd_calc.RMSD_DEFAULT_METHOD,
        'Model path' : Config.ML_MODEL_PATH,
        'Methods' : ','.join(methods),
        'Verbose' : opt.verbose,
    }
    sg.poptions(opts)

    runner = RunAllDirectory()
    runner.run(input_files, output_csv, opt)
    # meta-scoring
    #output_csv = "test_data/1xjr_m500_m1.csv"
    #mqs.do_scoring(output_csv)

    log = t.end('process: %i' % opt.number_processes)
    print('\n', log)
    print('Output: %s \n' % output_csv)
    ## log
    log_fn = output_csv.replace('.csv', '.log')
    f = open(log_fn, 'w')
    f.write(log + '\n')
    f.write(str(opts) + '\n')
    f.write('Output: %s\n' % output_csv)
    f.close()
    print('logging: %s' % log_fn)
    print('logging wrappers %s' % Config.LOG_DIRECTORY + os.sep)

    #with open(output_csv) as f:
    #    print(f.read())
