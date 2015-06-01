/****************************************************/
/***   ACS au complet et en detail, C++ version  ****/
/****************************************************/

#include "cdawg_utils.h"
#include "dynamic.h"


/* balenced parenthesis of sdsl */
#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>
#include <sdsl/bp_support.hpp>


#include <cstdlib>
#include <iostream>
#include <vector>

#include <sdsl/wavelet_trees.hpp>
#include <rlcsa.h>


using namespace std;
using namespace sdsl;
using namespace CSA;

#include <getRSS.h>



typedef sdsl::bit_vector::size_type size_type;



/* decompresses a node in the succint structure */
/* WARNING: res must point to an existing structure */
void decompress_node_succ( sdsl::bit_vector * balanced, sdsl::sd_vector<> * node_to_pos, sdsl::bp_support_sada<> * bps, 
			   sdsl::sd_vector<>::select_1_type * sdb_sel,
			   sacs * res, int number_node, int node_final_node, uchar * compressed) {
  int j;
  uint32_t code_trans,ltrans;
  int toadd;
  int nb_succinct_trans=0;
  int position,number_target_succ=0;
  /* decoding the byte coding transitions */

  if (number_node>0)
    position=(*sdb_sel)(number_node+1);  /* the position of the current node in compress */
  else position=0;
     
  code_trans=compressed[position];
  if (code_trans==255) { /* terminal node */
    res->s=1;
    position++;
    code_trans=compressed[position];
  } else {
     res->s=0;
  }

  for (j=0; j<MAXALPHA; j++) {
    res->arc[j]=code_trans%3;
    code_trans/=3;
  }

#ifdef DEBUGSEARCHSUCC
  fprintf(stdout,"\nDecoding node numner %d at pos %d, final %d :",number_node,position,res->s);
  
  for (j=0; j<MAXALPHA; j++) {
    fprintf(stdout,"[%c - %d]",indtolet(j),res->arc[j]);
  }
  fprintf(stdout,"\n");
#endif

  position++;

  for (j=0; j<MAXALPHA; j++) {
    switch (res->arc[j]) {
    case 1:  /* succint encoding in the trie */
     /* the number of the transition in the list of succint transitions from the same node */
      res->larc[j]=length_decode_int(&position,compressed+position); /* length of the transition */
      /* question: target node -> back to the tree */

 #ifdef DEBUGSEARCHSUCC
      fprintf(stdout,"CASZ 1-bis: number_target_succ %d\n",nb_succinct_trans);
#endif

      if (nb_succinct_trans==0) {    /* first child in the tree */
      //number_target_succ=bps->find_open(number_node+1)+1;
      //      number_target_succ=bps->rank(bps->select(number_node+1)+1)-1;
	number_target_succ=number_node+1;

#ifdef DEBUGSEARCHSUCC
	fprintf(stdout," -> number node %d CASZ 1-bis: first child %d\n",number_node,number_target_succ);
#endif


      }
      else {

#ifdef DEBUGSEARCHSUCC
	fprintf(stdout,"CASZ 1-bis: before next child %d\n",number_target_succ);
#endif


	
#ifdef DEBUGSEARCHSUCC
	int sel = bps->select(number_target_succ);
	int sel2 = bps->select(number_target_succ+1);
	int tmp=bps->find_close(number_target_succ);
	int tmp2=bps->find_close(number_target_succ+1);
        int tmp3 = bps->select(bps->find_close(number_target_succ));
	
	fprintf(stdout,"CASZ 1-BIS: %d, select %d, select+1 %d, find_close_num %d, find_close+1 %d\n",number_target_succ, sel, sel2, tmp,tmp2);
	fprintf(stdout,"CASZ 1-YES: select + close %d\n",  tmp3);
        fprintf(stdout,"CASZ 1-YES: select 1+1 %d, find_close de select 1+1 %d\n",bps->select(1+1),  bps->find_close(bps->select(1+1)));
	fprintf(stdout,"CASZ 1-YES: rank( find_close de select 1+1 + 1) %d\n",bps->rank(bps->find_close(bps->select(1+1))+1)-1);

#endif
	//	number_target_succ=bps->rank(bps->find_close(bps->select(number_target_succ))+1)-1;
	//	number_target_succ=bps->rank(bps->find_close(number_target_succ+1)-1;
	number_target_succ=bps->rank(bps->find_close(bps->select(number_target_succ+1))+1)-1;

#ifdef DEBUGSEARCHSUCC
	fprintf(stdout,"CASZ 1-bis: next child %d\n",number_target_succ);
#endif
      

      }
      res->arc[j]=number_target_succ;
      nb_succinct_trans++; 

#ifdef DEBUGSEARCHSUCC
      fprintf(stdout,"CASZ 1: trans %c -> num target node %d, length %d\n",indtolet(j),res->arc[j],res->larc[j]);
#endif

      break;

    case 2:  /* the number target of the node is encoded with the transition */

      ltrans=length_decode_int(&position,compressed+position);

#ifdef DEBUGSEARCHSUCC
      //fprintf(stdout,"after trans %c -> ltrans %d\n",indtolet(j),ltrans);
#endif

      toadd = ltrans%2;
      res->larc[j]=(ltrans>>1);   /* length of the transition */
      if (toadd==1)  { 
	/* res->arc[j]=length_decode_int(&position,compressed+position);  */
	res->arc[j]=length_decode_int(&position,compressed+position);  /* the number of the node, to search in sd_structure */
	      }
     else res->arc[j]=node_final_node; /* number of target node */ 
      
#ifdef DEBUGSEARCHSUCC
      fprintf(stdout,"trans %c -> toadd %d, num target node %d, length %d\n",indtolet(j),toadd,res->arc[j],res->larc[j]);
#endif
      
      break;
    default: res->arc[j]=INDEFINI;
    }
  }


return;
}



/* search in the succint compressed file, returning the number of ocurrences*/
/* warning: sastmp must point to an existing sas structure */
/* node_to_pos is the vector of node/position */
/* bps is the sasakane structure to index the parenthesis */
/* balanced is the bit vector of parenthesis */
/* cdawg_pos: position of the result */
int search_compressed_succ( sdsl::bit_vector * balanced, sdsl::sd_vector<> * node_to_pos, sdsl::bp_support_sada<> * bps,  
			    sdsl::sd_vector<>::select_1_type * sdb_sel,
			    sacs * sastmp,char * pattern, uchar * compressed, dynamic_array<int> * theocc, dynamic_array<locStack> * stack_acs,int seq_length, int number_final_node) {
  int lpat,pos_test,j,node_target,length_target,current_length,number_initial_state;
  locStack * Lacs;

#ifdef DEBUGSEARCHTEST
    fprintf(stdout, "searching %s\n",pattern);
#endif


     number_initial_state=0;   /* begin of the buffer = first node */ 

     /* succint searching of the position of the final node */
     /* pos_final_node=sdb_sel(number_final_node); */

#ifdef DEBUGSEARCHSUCC
    fprintf(stdout, "number of the final node : %d\n",number_final_node);
#endif

     lpat=strlen(pattern); /* length of the pattern. WARNING: the pattern must be in the sequence !! */

     /* first node */
     decompress_node_succ(balanced,node_to_pos,bps,sdb_sel,sastmp,number_initial_state, number_final_node, compressed);
     pos_test=0;
     
     while (pos_test<lpat) {    /* fastfind */
       j=lettoind(pattern[pos_test]);
       node_target=sastmp->arc[j];
       pos_test+=sastmp->larc[j];
       decompress_node_succ(balanced,node_to_pos,bps,sdb_sel,sastmp,node_target, number_final_node, compressed);
     }

#ifdef DEBUGSEARCHSUCC
    fprintf(stdout, "HERE found node at the end of the search %d\n",pos_test);
#endif

  
    // empty the stack
    stack_acs->makeempty();
    Lacs=stack_acs->push_empty();  //     stackplusplus(stack_acs);

    for (j=0;j<MAXALPHA;j++) {

#ifdef DEBUGSEARCHSUCC
      fprintf(stdout, "[%d - %d]\n",sastmp->arc[j],sastmp->larc[j]);
#endif

      Lacs->acs.arc[j]=sastmp->arc[j];
      Lacs->acs.larc[j]=sastmp->larc[j];

    }
    Lacs->acs.s=sastmp->s;
    Lacs->children=-1;
    Lacs->lengths=pos_test;
    Lacs->acs.lg=0;
    Lacs->acs.pos=INDEFINI;

    // copy of s, final or not !!
  
#ifdef DEBUGSEARCHSUCC
    fprintf(stdout, "length of the node at the end of the search %d\n\n",pos_test);
#endif
   
    while (!stack_acs->isempty()) {  // while there is a remaining state in the stack
     /* */
      Lacs=stack_acs->point_top();

      if (Lacs->children<MAXALPHA-1) {

	current_length=Lacs->lengths; // current length. Warning, should be here because nb_state_max can move
	Lacs->children+=1;
	
	node_target= Lacs->acs.arc[Lacs->children];


	if (node_target!=INDEFINI) { /* there is a transition */
	  
	  length_target=Lacs->acs.larc[Lacs->children];
	  
	  if (node_target==number_final_node) { /* WARNING transition to the final node !*/

	    theocc->push(seq_length-(current_length+length_target)+1);

#ifdef DEBUGSEARCHSUCC
	    fprintf(stdout, "[occurence at %d]\n",seq_length-(current_length+length_target)+1);
	     fprintf(stdout, "nb %d\n",stack_acs->nb_in);
#endif


	  } else { /* not to the last node ! */
	    
	    Lacs=stack_acs->push_empty();  //     stackplusplus(stack_acs);
	 
	    decompress_node_succ(balanced,node_to_pos,bps,sdb_sel, &(Lacs->acs) ,node_target, number_final_node, compressed);

	    Lacs->lengths=current_length+length_target;
	    Lacs->children=-1;

	  }
	}
      } else {

#ifdef DEBUGSEARCHSUCC
	// fprintf(stdout, "node on the top of the stak final ? %d\n",stack_acs->acs[nb_state_stack-1].s);
#endif	

	/* managing final nodes */
	if (Lacs->acs.s==1)  { /* the node is final ! */
	   theocc->push(seq_length-Lacs->lengths+1);

#ifdef DEBUGSEARCHSUCC
	    fprintf(stdout, "[occurence final at %d]\n",seq_length-stack_acs->lengths[nb_state_stack-1]+1);
#endif
	}

 #ifdef DEBUGSEARCHSUCC
	     fprintf(stdout, "back on stack\n");
#endif

	stack_acs->pop();
	Lacs=stack_acs->point_top();	
      }
    }

#ifdef DEBUGSEARCHTEST
    fprintf(stdout, "occurrences: %d\n",theocc->nbelem());
#endif
    return theocc->nbelem();
}




/***************************/
/*** FONCTION PRINCIPALE ***/
int main( int c, char* v[])
{

  int n,m;
 uchar * buf_coding;
 FILE *binFile;
 struct stat sdata; 
 clock_t start, end;
 double cpu_time_used;
 dynamic_array<int> * theocc_succ;
 dynamic_array<locStack> * stack_acs_succ;
 int seq_length2,number_final_node, size_comp;
 sacs * sastmp_succ;
 sdsl::bit_vector bp_succ;
 sdsl::sd_vector<> pos_succ;
 int pos_code=0;
 float statnbocc=0;

 if (c<2) { 
   fprintf(stderr,"program %s need 2 arguments: <input file><patterns file> \n", v[0]); 
   exit(-1);
 }



 

 srand(time(NULL));

 start = clock();

 
 std::string idx_basename(v[1]);
 std::string patterns(v[2]);
 

 /************************************************************************************************************/
 /*************************************  LOADING RLCSA  ****************************************************/
 /***********************************************************************************************************/

 std::cout << "Loading rlcsa file :" << idx_basename << std::endl;
 RLCSA * rlcsa_idx = new RLCSA(idx_basename);
 std::cout << "done." << std::endl;

 /************************************************************************************************************/
 /*************************************  LOADING SUCC FILE **************************************************/
 /***********************************************************************************************************/


 std::string name_comp(v[1]);
 name_comp+=".cdawg-succ";

 fprintf(stdout,"\n Loading comp file %s\n", name_comp.c_str());

 /* read buffer */ 
 if (stat((char *)name_comp.c_str(),&sdata) != 0) 
   fatal_error ("Could not open compressed (comp) file\n"); 
 size_comp = sdata.st_size; 
 buf_coding = (uchar*) malloc (size_comp+8); 
 if (buf_coding == NULL) 
   fatal_error ("Cannot allocate memory for the compressed file\n"); 
 
 binFile = fopen ((char *)name_comp.c_str(),"r"); 
  if (fread (buf_coding,size_comp,1, binFile) != 1) 
    fatal_error ("Could not read file\n"); 
  fclose(binFile); 

  fprintf(stdout,"done.\n");
 

 /************************************************************************************************************/
 /*************************************  SEARCHING PATTERNS *************************************************/
 /***********************************************************************************************************/
  
 theocc_succ = new dynamic_array<int>(NBOCC);
 stack_acs_succ = new  dynamic_array<locStack>(NBITEM);

 sastmp_succ=(sacs *) malloc(sizeof(sacs));  


 /* loading the bitbector of balanced parenthesis */
 std::string test5(v[1]);
 test5+=".cdawg-scc-bp";

#ifdef DEBUGOCCSUCC
 std::cout << "loading " << test5 << std::endl;
#endif

 sdsl::bit_vector b;
 sdsl::load_from_file(b, test5);



 /* build sadakane structure on top of b */
 sdsl::bp_support_sada<> bps(&b);

/* loading the bitbector of positions */
#ifdef DEBUGOCCSUCC
  fprintf(stdout,"\n sadakane done\n");
#endif

    // loading the elias-fano code of positions
 std::string test3(v[1]);
 //replaceExt(test3, "scc-ps");
 test3+=".cdawg-scc-ps";

 //std::ifstream in2(test3);
#ifdef DEBUGOCCSUCC
 std::cout << "loading " << test3 << std::endl;
#endif

 sdsl::sd_vector<> b_sd;
 sdsl::load_from_file(b_sd, test3);

#ifdef DEBUGOCCSUCC
 std::cout << " loaded !! " << b_sd.size() << std::endl;
 std::cout << " loaded !! " << b_sd << std::endl;
#endif
 
 sdsl::sd_vector<>::select_1_type sdb_sel(&b_sd);


#ifdef DEBUGOCCSUCC
 
   fprintf(stdout,"\n------------------------------\n");
   fprintf(stdout,"Structures loaded, Sadakane built\n");

#endif

   end = clock();
 cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

 fprintf(stdout,"------------------------------\n");
 fprintf(stdout,"[Succint] Loading in time %f\n",cpu_time_used);
   

 /************************************************************************************************************/
 /*************************************  SEARCHING           **************************************************/
 /***********************************************************************************************************/

 statnbocc=0;
 start = clock();

 
 std::cout << "searching patterns ... " << std::endl;
 std::ifstream ifs(patterns);

 std::string header;
 std::getline(ifs, header);

 n = get_number_of_patterns_bis(header);
 m = get_patterns_length_bis(header);
 
 int last_perc = 0;
 

   /* length of the sequence, nedeed for calculating the poistions of the occurrences */
 seq_length2=length_decode_int(&pos_code,buf_coding);
 /* position of the final node */
 number_final_node=length_decode_int(&pos_code,buf_coding+pos_code);

 for(int i=0;i<n;++i){
 
   int perc = (100*i)/n;
   if(perc>last_perc){
     std::cout << perc << "% done ..." << std::endl;
     last_perc=perc;
   }

   std::string p = std::string();

   for(int j=0;j<m;++j){
     char c;
     ifs.get(c);
     p+=c;
   }
     //locate with rlcsa
   auto rlcsa_range = rlcsa_idx->count(p);
   
   if(rlcsa_range.second >= rlcsa_range.first){
     
     // the search !!
     theocc_succ->makeempty();
  
    statnbocc+=search_compressed_succ(&b, &b_sd, &bps, &sdb_sel, sastmp_succ,(char *)p.c_str(),buf_coding+pos_code,theocc_succ,stack_acs_succ,seq_length2,number_final_node)/(float)n;

    end = clock();
   }
 }
  
 cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
 fprintf(stdout,"Ave. nb occur: %f; total time: %f; ave. time: %f\n",statnbocc,cpu_time_used,cpu_time_used/(float)n);

 /* closing */
  ifs.close();		
  delete rlcsa_idx;
  delete stack_acs_succ;
  delete theocc_succ;
 
}
