/*
 * ClaRNAwd_output_parser.h
 *
 *  Created on: 2015-11-14
 *      Author: Michal Boniecki
 *      based on: RNAView_output_parser.h //but separate, recoded 2016-10-25
 */

#ifndef _CLARNA_OUTPUT_PARSER_
#define _CLARNA_OUTPUT_PARSER_

#define MAX_N_PAIRS 10000
#define MAX_LINE_LENGTH 1024
#define MAX_NAME_LENGTH 32
#define MAX_N_CHAINS 52
#define MAX_N_ITEMS 64

#define CONTACT_EDGE_EDGE 0
#define CONTACT_STACKING  1
#define CONTACT_BASE_BR 2
#define CONTACT_BASE_PH 3
#define CONTACT_OTHER   4

#define MAX_SS_LINES 25

typedef struct
{
    int nucl_number_in_pdb;
    char chain_id, nucl_name;
}  T_ClaRNAwd_nucleotide;

typedef struct
{
    T_ClaRNAwd_nucleotide nucl_1, nucl_2;
    char contact_type[MAX_NAME_LENGTH];
    char cis_or_trans[MAX_NAME_LENGTH];
    int contact_index;
    double confidence;
}  T_ClaRNAwd_pair;

typedef struct
{
    char chain_id;
    int begin_pdb, term_pdb;
    int begin_renum, term_renum;
} T_chain_info;

class T_ClaRNAwd_data
{

protected:
    T_ClaRNAwd_pair v_nucleotide_pairs[MAX_N_PAIRS];
    T_chain_info v_chains[MAX_N_CHAINS];
    int n_nucleotide_contacts, n_nucleotides, n_chains;
    
    int *v_ss_lines[MAX_SS_LINES]; //<--- in case of psedoknots next line will be allocated, that's why is a table of tables here
    int n_ss_lines;
    
    void add_pair(char **in_parsed_line, int in_n_items);
    int parse_chains_info(char *curr_line);

public:
    T_ClaRNAwd_data(const char *filename);
    ~T_ClaRNAwd_data();
    int is_pair_canonical(T_ClaRNAwd_pair *curr_pair);
    int is_chain_break(int inp_i);
    int renum_nucl_number(T_ClaRNAwd_nucleotide *curr_nucleotide);
    int get_n_nucleotide_contacts();
    int get_n_nucleotides();
    int print_ss();
    int print_ss(const char *filename);
    int print_ss(FILE *outfile);
    int calc_ss();
    T_ClaRNAwd_pair *get_ClaRNAwd_pair_ptr(int i_pair);
};

#endif
