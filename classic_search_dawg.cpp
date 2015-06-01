/****************************************************/
/***   ACS au complet et en detail, C++ version  ****/
/****************************************************/

#include "cdawg_utils.h"
#include "dynamic.h"

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



/* decompresseses a node */
/* WARNING: res must point to an existing structure */
void decompress_node(sacs * res, int position, uchar * compressed) {
  int toadd,j;
  uint32_t code_trans,ltrans;

  /* decoding the byte coding transitions */
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

#ifdef DEBUGSEARCH
  fprintf(stdout,"Decoding node at pos %d, final %d :",position,res->s);
  
  for (j=0; j<MAXALPHA; j++) {
    fprintf(stdout,"[%c - %d]",indtolet(j),res->arc[j]);
  }
  fprintf(stdout,"\n");
#endif

  position++;

  for (j=0; j<MAXALPHA; j++) {
    switch (res->arc[j]) { /* right now, no succint encoding */
    case 1: 
    case 2:


      ltrans=length_decode_int(&position,compressed+position);

#ifdef DEBUGSEARCH
      //fprintf(stdout,"after trans %c -> ltrans %d\n",indtolet(j),ltrans);
#endif

      toadd = ltrans%2;
      res->larc[j]=(ltrans>>1);   /* length of the transition */
      if (toadd==1)  res->arc[j]=length_decode_int(&position,compressed+position); else res->arc[j]=POSFINALNODE; /* position of target node */ 
      
#ifdef DEBUGSEARCH
      fprintf(stdout,"trans %c -> toadd %d, pos target node %d, length %d\n",indtolet(j),toadd,res->arc[j],res->larc[j]);
#endif
      
      break;
    default: res->arc[j]=INDEFINI;
    }
  }



return;
}


   
/* search in the compressed file, returning the number of ocurrences*/
/* warning: sastmp must point to an existing sas structure */
/* cdawg_pos: position of the result */
/* dynamic_array<acs> * stack_acs already allocated */
int search_compressed(sacs * sastmp,char * pattern, uchar * compressed, dynamic_array<int> * theocc, 
		      dynamic_array<locStack> * stack_acs,int seq_length) {

     uint32_t pos_initial_state;   /* position of the initial node of the dawg */
     int lpat,pos_test,j,pos_target=0,length_target,current_length;
     locStack * Lacs;

#ifdef DEBUGSEARCHTEST
    fprintf(stdout, "searching %s\n",pattern);
#endif

     pos_initial_state=decode_int(compressed); 

     lpat=strlen(pattern); /* length of the pattern. WARNING: the pattern must be in the sequence !! */

     /* first node */
     decompress_node(sastmp,pos_initial_state, compressed);
     pos_test=0;
     
     while (pos_test<lpat) {    /* fastfind */
       j=lettoind(pattern[pos_test]);
       pos_target=sastmp->arc[j];
       pos_test+=sastmp->larc[j];
       decompress_node(sastmp,pos_target, compressed);
     }

#ifdef DEBUGSEARCH
    fprintf(stdout, "found node at the end of the search %d\n",pos_test);
#endif

     /* sastmp is the node from which we must explore all the subtree, and pos_test the length until the node */
     /* need to manage the stack. basic version: stack of sas, no saving on decompression !! */
    
    if (pos_target==POSFINALNODE) {
      theocc->push(seq_length-pos_test+1);
      return 1;
    }
 
    
    // empty the stack
    stack_acs->makeempty();
    

    Lacs=stack_acs->push_empty();  //     stackplusplus(stack_acs);

   
   

    for (j=0;j<MAXALPHA;j++) {

#ifdef DEBUGSEARCH
      fprintf(stdout, "[%d - %d]\n",sastmp->arc[j],sastmp->larc[j]);
#endif
      
      Lacs->acs.arc[j]=sastmp->arc[j];
      Lacs->acs.larc[j]=sastmp->larc[j];

      //      stack_acs->acs[0].arc[j]=sastmp->arc[j];
      //      stack_acs->acs[0].larc[j]=sastmp->larc[j];
    }
    // copy of s, final or not !!
    // stack_acs->acs[0].s=sastmp->s;
    Lacs->acs.s=sastmp->s;
    Lacs->children=-1;
    Lacs->lengths=pos_test;
    Lacs->acs.lg=0;
    Lacs->acs.pos=INDEFINI;
    
    // stack_acs->children[0]=-1;  /* no child taken rignt now first child ! */
    //stack_acs->lengths[0]=pos_test;  /* length to the first node in the stack */
    //stack_acs->nb_in+=1;

#ifdef DEBUGSEARCH
    fprintf(stdout, "length of the node at the end of the search %d\n\n",pos_test);
#endif
   
    while (!stack_acs->isempty()) {  // while there is a remaining state in the stack
     /* */
     
      Lacs=stack_acs->point_top();

      if (Lacs->children<MAXALPHA-1) {

	current_length=Lacs->lengths; // current length. Warning, should be here because nb_state_max can move
	Lacs->children+=1;
	
	pos_target= Lacs->acs.arc[Lacs->children];
	if (pos_target!=INDEFINI) { /* there is a transition */
	  	  	  
	  length_target=Lacs->acs.larc[Lacs->children];

	  if (pos_target==POSFINALNODE) { /* WARNING transition to the final node !*/
	    /* new occurrence */
	    theocc->push(seq_length-(current_length+length_target)+1);

	    

#ifdef DEBUGSEARCH
	    fprintf(stdout, "[occurence at %d]\n",seq_length-(current_length+length_target)+1);
#endif


	  } else { /* not to the last node ! */

	    
	    Lacs=stack_acs->push_empty();  //     stackplusplus(stack_acs);
	   
	    decompress_node(&(Lacs->acs),pos_target,compressed); /* decompression of the target node in the stack */
	    Lacs->lengths=current_length+length_target;
	    Lacs->children=-1;

	  }
	}
      } else {

#ifdef DEBUGSEARCH
	// fprintf(stdout, "node on the top of the stak final ? %d\n",stack_acs->acs[nb_state_stack-1].s);
#endif	

	/* managing final nodes */
	if (Lacs->acs.s==1)  { /* the node is final ! */
	  theocc->push(seq_length-Lacs->lengths+1);

#ifdef DEBUGSEARCH
	    fprintf(stdout, "[occurence final at %d]\n",seq_length-Lacs->lengths+1);
#endif
	}

 #ifdef DEBUGSEARCH
	     fprintf(stdout, "back on stack \n");
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



// main function
int main( int c, char* v[])
{
 uchar * buf_coding;
 FILE *binFile;
 struct stat sdata; 
 clock_t start, end;
 double cpu_time_used;
 dynamic_array<int> * theocc;
 dynamic_array<locStack> * stack_acs;
 int size_comp,n,m;
 sacs * sastmp;
 sdsl::bit_vector bp_succ;
 sdsl::sd_vector<> pos_succ;
 float statnbocc=0;

if (c<2) { 
   fprintf(stderr,"program %s need 2 arguments: <input file> <patterns file>\n", v[0]); 
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
 /*************************************  LOADING COMP FILE **************************************************/
 /***********************************************************************************************************/


 std::string name_comp(v[1]);
 name_comp+=".cdawg-comp";

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


  theocc = new dynamic_array<int>(NBOCC);
  stack_acs = new  dynamic_array<locStack>(NBITEM);


  sastmp=(sacs *) malloc(sizeof(sacs));  


  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  fprintf(stdout,"------------------------------\n");
  fprintf(stdout,"[Classic] Loading in time (with RCCSA) %f\n",cpu_time_used);
  
  statnbocc=0;
  start = clock();
 
  std::cout << "searching patterns ... " << std::endl;
  std::ifstream ifs(patterns);
  
  std::string header;
  std::getline(ifs, header);
  
  n = get_number_of_patterns_bis(header);
  m = get_patterns_length_bis(header);
  
  int last_perc = 0;
  
  //extract patterns from file and search them in the index
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
      //  p+='\0';
      
      theocc->makeempty();
      
      statnbocc+=search_compressed(sastmp,(char *)p.c_str(),buf_coding,theocc,stack_acs,size_comp)/(float)n;
          
    }
    
 }
  
 end = clock();
 
 cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
 fprintf(stdout,"Ave. nb occur: %f; total time: %f; ave. time: %f\n",statnbocc,cpu_time_used,cpu_time_used/(float)n);


 /* closing */
 ifs.close();		
 delete rlcsa_idx;
 delete stack_acs;
 delete theocc;
  
}
