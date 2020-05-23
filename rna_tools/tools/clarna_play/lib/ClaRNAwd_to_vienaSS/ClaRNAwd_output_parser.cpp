/*
 * ClaRNAwd_output_parser.c
 *
 *  Created on: 2015-12-10
 *      Author: Michal Boniecki
 *      based on: RNAView_output_parser.c //but separate, recoded 2016-10-25
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ClaRNAwd_output_parser.h"

#define TREAT_GU_AS_CANON  //extended ss means WC pairs considering not only A-U/U-A and C-G/G-C pairings but also G-U/U-G pairings

#define DATA_LINE_LENGTH 65
#define MIN_N_ITEMS 9

//#define VERBOSE

int split_string(char *string2parse, const char *delim_string, char **output_list, int max_n_items)
{
    char *curr_token,*saveptr;
    int n_items=0;
  
    curr_token=strtok_r(string2parse,delim_string,&saveptr);
    while(curr_token!=NULL)
    {
	output_list[n_items]=curr_token;
	n_items++;
	curr_token=strtok_r(NULL,delim_string,&saveptr);
	if(n_items>=MAX_N_ITEMS)
	{
	    fprintf(stderr,"too many items in output list from parsing in function split_string(char*,const char*,char**,int)\n");
	    fprintf(stderr,"increase constant MAX_N_ITEMS and recompile\n");
	}
    }
    return n_items;
}

T_ClaRNAwd_data::T_ClaRNAwd_data(const char *filename){

    FILE *inpfile;
    char curr_line[MAX_LINE_LENGTH];
    char *splitted_line[MAX_N_ITEMS];
    char chains_info_tmp[8]; // it will contain 7 chars + terminator, for checking "chains:" tag
    int i, n_items, chains_info_found=0;
    
    n_chains = 0;
    n_nucleotides = 0;
    n_nucleotide_contacts = 0;
    
    for(i=0; i<MAX_SS_LINES; i++)
	v_ss_lines[i] = NULL;
    n_ss_lines = 0;

    if(!(inpfile = fopen(filename,"r"))){
	fprintf(stderr,"error opening file %s in T_ClaRNAwd_data::T_ClaRNAwd_data(const char*)\n",filename);
	exit(EXIT_FAILURE);
    }
    
    while(fgets(curr_line,MAX_LINE_LENGTH,inpfile) != NULL){

	int curr_length = strlen(curr_line);

	if(curr_length>0){
	    curr_line[curr_length-1] = 0;
	    if(curr_length >= 7){
		strncpy(chains_info_tmp, curr_line, 7);
//		fprintf(stderr, "chains_info_tmp: _%s_", chains_info_tmp);
		chains_info_tmp[7] = 0;
		if(strcmp(chains_info_tmp, "chains:") == 0){
		    if(chains_info_found == 0){
			chains_info_found = 1;
			parse_chains_info(curr_line+7);
		    }
		    else{
			fprintf(stderr, "chains: ... exists at least two times in file: %s, it must exist only once, check file format, exiting program ...\n", filename);
			exit(EXIT_FAILURE);
		    }
		}
	    }
	}

	if(curr_length >= DATA_LINE_LENGTH){
//	    fprintf(stderr,"%s\n",curr_line);
//	fprintf(stderr,"length: %3d\n",curr_length);

	    n_items = split_string(curr_line, " ", splitted_line, MAX_N_ITEMS);
	    if(n_items < MIN_N_ITEMS){
		fprintf(stderr, "too few items in line:\n%s\n", curr_line);
		exit(EXIT_FAILURE);
	    }
	
	    add_pair(splitted_line, n_items);
	    if(n_items > MIN_N_ITEMS){
		splitted_line[7] = splitted_line[9];
		splitted_line[8] = splitted_line[10];
		add_pair(splitted_line, n_items);
	    }
	}
    }

    if(chains_info_found == 0){
	fprintf(stderr, "there is no chains information in file: %s\n", filename);
	fprintf(stderr, "chains field should look like:\n");
	fprintf(stderr, "chains:  A 1 47  B 52 62\n");
	fprintf(stderr, "where A and B are chain ids, 1 47 and 52 62 are first and last number of nucleotide (as it is in pdb) for chain A and B, respectively\n");
	exit(EXIT_FAILURE);
    }

    fclose(inpfile);
}

T_ClaRNAwd_data::~T_ClaRNAwd_data(){

    for(int i=0; i<MAX_SS_LINES; i++)
	if(v_ss_lines[i] != NULL)
	    delete [] v_ss_lines[i];
}

int  T_ClaRNAwd_data::parse_chains_info(char *curr_line){

//    fprintf(stderr, "calling function: int  T_ClaRNAwd_data::parse_chains_info(const char *curr_line)\n");
//    fprintf(stderr, "line to parse: _%s_\n", curr_line);

    int i, n_items;
    char *splitted_line[MAX_N_ITEMS];
    char curr_line_backup[MAX_LINE_LENGTH];

    strcpy(curr_line_backup, curr_line); // in order to have it intact for displaying error messages

    n_items = split_string(curr_line, " ", splitted_line, MAX_N_ITEMS);
    if(n_items < 3){
	fprintf(stderr, "there is no chains information\n");
	fprintf(stderr, "chains field should look like:\n");
	fprintf(stderr, "chains: A 1 47  B 52 62\n");
	fprintf(stderr, "where A and B are chain ids, 1 47 and 52 62 are first and last number of nucleotide (as it is in pdb) for chain A and B, respectively\n");
	exit(EXIT_FAILURE);
    }

    if(n_items%3 != 0){
	fprintf(stderr, "parsing error of line \"chains:\", number of items should be dividable by 3, because of the line format which is:\n");
	fprintf(stderr, "chains: chain_id_1 first_nucl_num_in_chain_1 last_nucl_num_in_chain_1 ...\n");
	fprintf(stderr, "i.e.\n");
	fprintf(stderr, "chains:  A 1 47  B 52 62\n");
	fprintf(stderr, "where A and B are chain ids, 1 47 and 52 62 are first and last number of nucleotide (as it is in pdb) for chain A and B, respectively. They should be separated by space.\n\n");
	fprintf(stderr, "in fact chains line is:\n%s\n", curr_line);
	fprintf(stderr, "n_items after \"chain:\" is %d\n", n_items);
    }

    n_chains = n_items/3; //n_chains is actually number of triples
    if(n_chains >= MAX_N_CHAINS){
	fprintf(stderr, "field: chains: indicates too many chains ... max number of chains is: %d, but actually is: %d, probably data error\n", MAX_N_CHAINS, n_chains);
	exit(EXIT_FAILURE);
    }
    int curr_begin_renum = 0;
    for(i=0; i<n_chains; i++){ // parsing consecutive chain fields in chains: ...

	int curr_begin_pdb, curr_term_pdb, curr_chain_length;

	if(strlen(splitted_line[i*3]) != 1){
	    fprintf(stderr, "error parsing of chain line: %s,\nchain_id should span ONE character, example: A 1 47\n", curr_line_backup);
	    exit(EXIT_FAILURE);
	}
	v_chains[i].chain_id = splitted_line[i*3][0];
//fprintf(stderr, "chain_id: %c\n", v_chains[i].chain_id);

	if(sscanf(splitted_line[i*3+1], "%d", &curr_begin_pdb) != 1){
	    fprintf(stderr, "error parsing of chain line: %s,\nerror converting: %s into int\n", curr_line_backup, splitted_line[i*3+1]);
	    exit(EXIT_FAILURE);
	}
	if(sscanf(splitted_line[i*3+2], "%d", &curr_term_pdb) != 1){
	    fprintf(stderr, "error parsing of chain line: %s,\nerror converting: %s into int\n", curr_line_backup, splitted_line[i*3+2]);
	    exit(EXIT_FAILURE);
	}
	if(curr_term_pdb <= curr_begin_pdb){
	    fprintf(stderr, "data inconsistency in chain line: %s, last_nucl_number is lower than first_nucl_number: chain_id: %c first: %d, last: %d\n", curr_line_backup, v_chains[i].chain_id, curr_begin_pdb, curr_term_pdb);
	    exit(EXIT_FAILURE);
	}

	curr_chain_length = curr_term_pdb - curr_begin_pdb + 1;

	v_chains[i].begin_pdb = curr_begin_pdb;
	v_chains[i].term_pdb = curr_term_pdb;
	v_chains[i].begin_renum = curr_begin_renum;
	v_chains[i].term_renum = curr_begin_renum + curr_chain_length - 1;
	curr_begin_renum += curr_chain_length;
	
//fprintf(stderr, "begin_term_pdb: %d %d\n", curr_begin_pdb, curr_term_pdb);
//fprintf(stderr, "begin_term_ren: %d %d\n", v_chains[i].begin_renum, v_chains[i].term_renum);

//	fprintf(stderr, "%s\n", splitted_line[i]);
    }
    n_nucleotides = v_chains[i-1].term_renum+1;

    return n_chains;
}

void T_ClaRNAwd_data::add_pair(char **in_parsed_line, int in_n_items){

    int  n_items;
    char *splitted_field[MAX_N_ITEMS];
/*
    for(int i=0; i<in_n_items; i++){
	fprintf(stderr, " _%s_ ", in_parsed_line[i]);
    }
*/
//    fprintf(stderr, "\n");

    v_nucleotide_pairs[n_nucleotide_contacts].contact_index = CONTACT_OTHER; //assigning contact type: other as a default, to be changes (or leave like this)

    v_nucleotide_pairs[n_nucleotide_contacts].nucl_1.chain_id = in_parsed_line[0][0];
    v_nucleotide_pairs[n_nucleotide_contacts].nucl_2.chain_id = in_parsed_line[2][0];

    v_nucleotide_pairs[n_nucleotide_contacts].nucl_1.nucl_name = in_parsed_line[5][0];
    v_nucleotide_pairs[n_nucleotide_contacts].nucl_2.nucl_name = in_parsed_line[6][0];

    if( sscanf(in_parsed_line[1],"%d",&(v_nucleotide_pairs[n_nucleotide_contacts].nucl_1.nucl_number_in_pdb)) != 1 ){
	fprintf(stderr,"error processing data outCR file, error converting %s into int variable\n", in_parsed_line[1]);
	exit(EXIT_FAILURE);
    }

    if( sscanf(in_parsed_line[3],"%d",&(v_nucleotide_pairs[n_nucleotide_contacts].nucl_2.nucl_number_in_pdb)) != 1 ){
	fprintf(stderr,"error processing data outCR file, error converting %s into int variable\n", in_parsed_line[3]);
	exit(EXIT_FAILURE);
    }
#ifdef VERBOSE
    fprintf(stderr, "adding pair: %c %c %d <---> %c %c %d --- type: %s %s\n",
    v_nucleotide_pairs[n_nucleotide_contacts].nucl_1.chain_id, v_nucleotide_pairs[n_nucleotide_contacts].nucl_1.nucl_name, v_nucleotide_pairs[n_nucleotide_contacts].nucl_1.nucl_number_in_pdb,
    v_nucleotide_pairs[n_nucleotide_contacts].nucl_2.chain_id, v_nucleotide_pairs[n_nucleotide_contacts].nucl_2.nucl_name, v_nucleotide_pairs[n_nucleotide_contacts].nucl_2.nucl_number_in_pdb,
    in_parsed_line[7], in_parsed_line[8]);
#endif
    n_items = split_string(in_parsed_line[7], "_", splitted_field, MAX_N_ITEMS);
//    fprintf(stderr, "n_items: %d\n", n_items);

    if(n_items == 2){

	    if( (strcmp(splitted_field[1], "cis") == 0) || (strcmp(splitted_field[1], "tran") == 0) ){
		if(strlen(splitted_field[0]) == 2){
		    strcpy(v_nucleotide_pairs[n_nucleotide_contacts].contact_type, splitted_field[0]);
		    strcpy(v_nucleotide_pairs[n_nucleotide_contacts].cis_or_trans, splitted_field[1]);
		    v_nucleotide_pairs[n_nucleotide_contacts].contact_index = CONTACT_EDGE_EDGE;
		}
		else fprintf(stderr, "warning: incorrect type of base contact: _%s_\n", in_parsed_line[7]);
//fprintf(stderr, "--- cis or trans detected ---\n");
	    }
	    else{
//fprintf(stderr, "--------------------------------\n");
		char tmp_str[DATA_LINE_LENGTH];
		strncpy(tmp_str, splitted_field[1]+strlen(splitted_field[1])-2, 2);
		tmp_str[2]=0;
		if(strcmp(tmp_str, "BR") == 0){
		    strcpy(v_nucleotide_pairs[n_nucleotide_contacts].contact_type, tmp_str);
		    strcpy(v_nucleotide_pairs[n_nucleotide_contacts].cis_or_trans, "");
		    v_nucleotide_pairs[n_nucleotide_contacts].contact_index = CONTACT_BASE_BR;
//fprintf(stderr, "--- contact BR detected ---\n");
		}
		else if(strcmp(tmp_str, "Ph") == 0){
		    strcpy(v_nucleotide_pairs[n_nucleotide_contacts].contact_type, tmp_str);
		    strcpy(v_nucleotide_pairs[n_nucleotide_contacts].cis_or_trans, "");
		    v_nucleotide_pairs[n_nucleotide_contacts].contact_index = CONTACT_BASE_PH;
//fprintf(stderr, "--- contact PH detected ---\n");
		}
//fprintf(stderr, "_%s_\n", tmp_str);
	    }

    }
    else if(n_items == 1){

	char tmp_str[DATA_LINE_LENGTH];
	strcpy(tmp_str, in_parsed_line[7]);

	if(strlen(tmp_str) == 2){
	    if( (strcmp(tmp_str,"<<") == 0) || (strcmp(tmp_str,">>") == 0) || (strcmp(tmp_str,"<>") == 0) || (strcmp(tmp_str,"><") == 0) ){
		strcpy(v_nucleotide_pairs[n_nucleotide_contacts].contact_type, tmp_str);
		strcpy(v_nucleotide_pairs[n_nucleotide_contacts].cis_or_trans, "");
		v_nucleotide_pairs[n_nucleotide_contacts].contact_index = CONTACT_STACKING;
//fprintf(stderr, "--- contact stacking detected ---\n");
	    }
	    else fprintf(stderr, "warning: incorrect base-base contact detected: %s\n", tmp_str);
	}
    }

    n_nucleotide_contacts++;

}

int T_ClaRNAwd_data::get_n_nucleotides(){

    return n_nucleotides;
}

int T_ClaRNAwd_data::get_n_nucleotide_contacts(){

    return n_nucleotide_contacts;
}

T_ClaRNAwd_pair* T_ClaRNAwd_data::get_ClaRNAwd_pair_ptr(int i_pair){
    
    if(i_pair < 0 || i_pair >= n_nucleotide_contacts){
	fprintf(stderr,"i_pair exceeds the data table in T_ClaRNAwd_pair* T_ClaRNAwd_data::get_ClaRNAwd_pair_ptr(int), should be <0,%d>, but is %d\n", n_nucleotide_contacts-1, i_pair);
	exit(EXIT_FAILURE);
    }
    
    return v_nucleotide_pairs+i_pair;
}

int T_ClaRNAwd_data::print_ss(){

    return print_ss(stdout);
}

int T_ClaRNAwd_data::print_ss(const char *filename){

    FILE *outfile;
    int result;
    
    if((outfile=fopen(filename,"w")) == NULL){
	fprintf(stderr,"error opening file %s in int T_ClaRNAwd_data::print_ss(const char*)\n",filename);
	exit(EXIT_FAILURE);
    }
    
    result = print_ss(outfile);
    
    fclose(outfile);
    
    return result;
}

int T_ClaRNAwd_data::print_ss(FILE *outfile){

    int result, i, i2;
    
    result = calc_ss();
    
    for(i=0; i<n_ss_lines; i++){
	for(i2=0; i2<n_nucleotides; i2++){
	    if(v_ss_lines[i][i2] == 0){
		fprintf(outfile,".");
		if(is_chain_break(i2))
		    fprintf(outfile," ");
	    }
	    else if(v_ss_lines[i][i2] == -1)
		fprintf(outfile,"(");
	    else if(v_ss_lines[i][i2] == 1)
		fprintf(outfile,")");
	    else{
		fprintf(stderr,"error processing data in int T_ClaRNAwd_data::print_ss(FILE*)\nv_ss_lines[][] can contain only values 0, -1 or 1, but in v_ss_lines[%d][%d] is %d\n",i,i2,v_ss_lines[i][i2]);
		exit(EXIT_FAILURE);
	    }
	}
	fprintf(outfile,"\n");
    }

    return result;
}


int T_ClaRNAwd_data::is_pair_canonical(T_ClaRNAwd_pair *curr_pair){
// it could be coded better, but in this case is not critical
    if(strcmp(curr_pair->contact_type, "WW") != 0)
	return 0;
    if(strcmp(curr_pair->cis_or_trans, "cis") != 0)
	return 0;

    if(curr_pair->nucl_1.nucl_name == 'C' && curr_pair->nucl_2.nucl_name == 'G')
	return 1;
    if(curr_pair->nucl_1.nucl_name == 'G' && curr_pair->nucl_2.nucl_name == 'C')
	return 1;
    if(curr_pair->nucl_1.nucl_name == 'A' && curr_pair->nucl_2.nucl_name == 'U')
	return 1;
    if(curr_pair->nucl_1.nucl_name == 'U' && curr_pair->nucl_2.nucl_name == 'A')
	return 1;
#ifdef TREAT_GU_AS_CANON
    if(curr_pair->nucl_1.nucl_name == 'G' && curr_pair->nucl_2.nucl_name == 'U')
	return 1;
    if(curr_pair->nucl_1.nucl_name == 'U' && curr_pair->nucl_2.nucl_name == 'G')
	return 1;
#endif
    return 0;
}

int T_ClaRNAwd_data::renum_nucl_number(T_ClaRNAwd_nucleotide *curr_nucleotide){

    int i, nucl_renum;

    for(i=0; i<n_chains; i++)
	if(curr_nucleotide->chain_id == v_chains[i].chain_id){
	
	    nucl_renum = curr_nucleotide->nucl_number_in_pdb - v_chains[i].begin_pdb + v_chains[i].begin_renum;

	    if(nucl_renum < v_chains[i].begin_renum){
		fprintf(stderr, "data inconsistensy: nucleotide in contact list has number lower than declared in chains line, for this chain.\n");
		fprintf(stderr, "chain %c spans: %d %d, but there is nucleotide: chain: %c number: %d\n",\
		v_chains[i].chain_id, v_chains[i].begin_pdb, v_chains[i].term_pdb, curr_nucleotide->chain_id, curr_nucleotide->nucl_number_in_pdb);
		exit(EXIT_FAILURE);
	    }

	    if(nucl_renum > v_chains[i].term_renum){
		fprintf(stderr, "data inconsistensy: nucleotide in contact list has number greater than declared in chains line, for this chain.\n");
		fprintf(stderr, "chain %c spans: %d %d, but there is nucleotide: chain: %c number: %d\n",\
		v_chains[i].chain_id, v_chains[i].begin_pdb, v_chains[i].term_pdb, curr_nucleotide->chain_id, curr_nucleotide->nucl_number_in_pdb);
		exit(EXIT_FAILURE);
	    }
	    return nucl_renum;
	}

    fprintf(stderr, "data inconsistency: nucleotide present in base-pair list has pdb_id: %c, which hasn't been found among chains ids defined in field \"chain:\"\n", curr_nucleotide->chain_id);
    fprintf(stderr, "chain ids defined in field \"chains:\" are: ");
    for(i=0; i<n_chains; i++)
	fprintf(stderr, " %c", v_chains[i].chain_id);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
    return -1; //just to shut up compiler ;-)
}

int T_ClaRNAwd_data::is_chain_break(int inp_i){

    int i;

    for(i=0; i<n_chains-1; i++) //n_chains-1 because we don't want to have break after last chain
	if(inp_i == v_chains[i].term_renum)
	    return 1;

    return 0;
}

int T_ClaRNAwd_data::calc_ss(){
    
    int i, i2, j, k;
    int open_bracket_index, close_bracket_index;
    int partial_sum;

    if(n_nucleotides && n_nucleotide_contacts){
	
	for(i=0; i<MAX_SS_LINES; i++){
	    v_ss_lines[i] = new int[n_nucleotides];
	    for(i2=0; i2<n_nucleotides; i2++)
		v_ss_lines[i][i2]=0;
	}
	n_ss_lines=1;

	for(i=0; i<n_nucleotide_contacts; i++){

	    if(is_pair_canonical(v_nucleotide_pairs+i) == 1){
//fprintf(stderr,"%s %s %c %c\n",v_nucleotide_pairs[i].contact_type, v_nucleotide_pairs[i].cis_or_trans, v_nucleotide_pairs[i].nucl_1.nucl_name, v_nucleotide_pairs[i].nucl_2.nucl_name);

		open_bracket_index = renum_nucl_number(&v_nucleotide_pairs[i].nucl_1);
		close_bracket_index = renum_nucl_number(&v_nucleotide_pairs[i].nucl_2);

//fprintf(stderr,"open_bracket_index = %d\n",open_bracket_index);
//fprintf(stderr,"close_bracket_index = %d\n\n",close_bracket_index);
		
		for(j=0; j<MAX_SS_LINES; j++){
		    
		    partial_sum = 0; //partial_sum indicates what's going on between opening and closing braket, during summation partial_sum cannot be positive value --- it means that there are closing brackets first
		    if(v_ss_lines[j][open_bracket_index]==0 && v_ss_lines[j][close_bracket_index]==0){
			for(k=open_bracket_index; k<=close_bracket_index; k++){ //also when partial sum is not equal 0 it means that there is a pseudoknot and program should go farther to the next lines
			    partial_sum += v_ss_lines[j][k];
			    if(partial_sum > 0){
				partial_sum += 100000;
			        break;
			    }
			}
		    }
		    else partial_sum = -1;
//fprintf(stderr,"j = %d, partial_sum = %d\n",j,partial_sum);
			
		    if(partial_sum == 0){
			v_ss_lines[j][open_bracket_index]= -1;
			v_ss_lines[j][close_bracket_index]= 1;
			break;
		    }
		    
		    if(n_ss_lines==j+1) n_ss_lines++;
//fprintf(stderr,"n_ss_lines = %d\n",n_ss_lines);
		}

		if(j == MAX_SS_LINES){
		    fprintf(stderr,"too many pseudoknots detected, increase MAX_SS_LINES and recompile, or error in structure or ClaRNAwd classication is possible - check data\n");
		    exit(EXIT_FAILURE);
		}
	    }
	}
    }
    else
	return 0;

    return n_ss_lines;
}


int main(int argc, char *argv[]){

    if(argc<2){
	fprintf(stderr, "usage: ClaRNAwd_output_parser_get_SS file.pdb.outCR\n");
	exit(EXIT_FAILURE);
    }

    T_ClaRNAwd_data *curr_ClaRNAwd_data;
    int n_nucleotide_contacts, n_nucleotides;
    
    curr_ClaRNAwd_data = new T_ClaRNAwd_data(argv[1]);
#ifdef VERBOSE
fprintf(stderr, "\nDATA HAVE BEEN LOADED\n\n");
#endif

    n_nucleotides = curr_ClaRNAwd_data->get_n_nucleotides();
#ifdef VERBOSE
fprintf(stderr,"n_nucleotides = %d\n",n_nucleotides);
#endif

    n_nucleotide_contacts = curr_ClaRNAwd_data->get_n_nucleotide_contacts();
#ifdef VERBOSE
fprintf(stderr,"n_nucleotide_contacts = %d\n",n_nucleotide_contacts);
#endif

#ifdef VERBOSE
    for(int i=0;i<n_nucleotide_contacts;i++){
	
	curr_ClaRNAwd_pair = curr_ClaRNAwd_data->get_ClaRNAwd_pair_ptr(i);
//	printf("%3d\n",curr_ClaRNAwd_pair->nucl_1.nucl_number_from_1);
	printf("%3d %c %c %s %s %3d %c %c\n", curr_ClaRNAwd_pair->nucl_1.nucl_number_in_pdb, curr_ClaRNAwd_pair->nucl_1.chain_id, curr_ClaRNAwd_pair->nucl_1.nucl_name, \
curr_ClaRNAwd_pair->contact_type, curr_ClaRNAwd_pair->cis_or_trans, curr_ClaRNAwd_pair->nucl_2.nucl_number_in_pdb, curr_ClaRNAwd_pair->nucl_2.chain_id, curr_ClaRNAwd_pair->nucl_2.nucl_name);
    }
#endif

    curr_ClaRNAwd_data->print_ss();

    delete curr_ClaRNAwd_data;

    return 0;
}
