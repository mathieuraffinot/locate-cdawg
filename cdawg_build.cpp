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


#define INDEFINI -1
#define INITIAL 0
#define FINAL 1
#define MAXALPHA 5 /* alphabet total */
#define TRUE 1
#define FALSE 0

// #define DEBUG 1 
// #define DEBUGSEARCH 1
// #define DEBUGOCC 1
// #define DEBUGPARENTHESIS 1
// #define DEBUGSUCCINT 1
// #define DEBUGSEARCHSUCC 1
// #define DEBUGOCCSUCC 1

#define NBITEM 10000
#define NBOCC 1000     /* initial number of occurrences for the dynamical list */

#define INITCODING    10000    /* first value of the number of char of the coding buffer (dynamic buffer) */
#define SIZEBUFFUCHAR  200    

#define POSFINALNODE 6        /* beginning of the coding of the final node of the automation in BYTES */


typedef unsigned char uchar;
typedef sdsl::bit_vector::size_type size_type;


extern int tobasethree[];



// dynamical increment of the number of states in the dynamic array
void nacsplusplus_state(int * Nacs, Pdynarray table_acs) {
  Sacs tmpacs,acs;
  int ind,k;

 
  
  if (*Nacs==table_acs->nb_item-1) {

    tmpacs=(Sacs)realloc(table_acs->acs,table_acs->nb_item*2*sizeof(sacs));
    table_acs->nb_item=table_acs->nb_item*2;
    fprintf(stdout, "Realloc general %ld\n",table_acs->nb_item*sizeof(sacs));

    if (!tmpacs)
                {
                        fprintf(stderr, "ERROR: Couldn't realloc memory!\n");
                        exit(-1);
                }
    else { table_acs->acs=tmpacs; acs=tmpacs; }
    
    for (ind=*Nacs+1;ind<table_acs->nb_item;ind++) {
      acs[ind].s=INDEFINI;
      acs[ind].lg=0;
      acs[ind].pos=INDEFINI;
      for (k=0; k<MAXALPHA; k++) {
	acs[ind].arc[k]=INDEFINI;
	acs[ind].larc[k]=INDEFINI;
      }
    }
  }

  (*Nacs)++;
}




/*************************************************************************/
/***************** Construction de l'ACS *********************************/
/*************************************************************************/


/* label comparison */
int arc_cmp (char *S, int i, int j, int lg)
{int k=0;
  while (k<lg && S[i++]==S[j++]) k++;
  return k;
}


/* initialisation of the acs */
int init_acs (char *S, int lgrS, Pdynarray table_acs, int Nacs)
{int i, j,k;

 Sacs acs = table_acs->acs;
 k= table_acs->nb_item;

 for (i=0; i<k; i++) {
    acs[i].s=INDEFINI;
    acs[i].lg=0;
    acs[i].pos=INDEFINI;
    for (j=0; j<MAXALPHA; j++) {
       acs[i].arc[j]=INDEFINI;
       acs[i].larc[j]=INDEFINI;
       }
    }
 acs[INITIAL].s=INDEFINI;
 acs[INITIAL].lg=0;
 acs[INITIAL].pos=0;
 acs[INITIAL].arc[lettoind(S[0])]=FINAL;
 acs[INITIAL].larc[lettoind(S[0])]=lgrS;
 acs[FINAL].s=INITIAL;
 acs[FINAL].lg=lgrS;

 /* modif ici acs[FINAL].pos=lgrS-1 */
 acs[FINAL].pos=lgrS-1;
 return FINAL;
}



/* clone state (E,j) in position Nacs ***/
void dup_etat (char *S, Pdynarray table_acs, int *Nacs, int E, int j)
{int k,  coup, Ej;
  Sacs acs = table_acs->acs;
  
  coup=acs[E].larc[j],
  Ej=acs[E].arc[j];

  nacsplusplus_state(Nacs,table_acs);
  acs= table_acs->acs;
  
 
  for (k=0; k<MAXALPHA; k++) {
    acs[*Nacs].arc[k]=acs[Ej].arc[k];
    acs[*Nacs].larc[k]=acs[Ej].larc[k];
    }
 acs[*Nacs].lg=acs[E].lg+acs[E].larc[j];
 acs[*Nacs].pos=acs[Ej].pos;
 acs[*Nacs].s=acs[Ej].s;
 acs[Ej].s=*Nacs;
 acs[FINAL].s=*Nacs;
 acs[E].arc[j]=*Nacs;

 /* Calcul du suffixe de Nacs */
 while (E!=INITIAL || coup>1) {
    if (E!=INITIAL)
       E=acs[E].s;
    else {
       coup--;
       j=lettoind(S[acs[*Nacs].pos-coup+1]);
       }
    while (acs[E].larc[j] < coup) {         /* Fastfind */
       coup=coup-acs[E].larc[j];
       E=acs[E].arc[j];
       j=lettoind(S[acs[*Nacs].pos-coup+1]);
       }
    if (acs[E].lg+acs[E].larc[j] < acs[acs[E].arc[j]].lg) 
       acs[E].arc[j]=*Nacs;
    else 
       break;
    }

}





/********************************************/
/*** Creation de l'Etat Nacs a la coupure ***/
/* Et s(F)=Nacs  */
void cree_etat(char *S, int lgrS, int i, Pdynarray table_acs, int *Nacs,
               int E, int j, int coup)
{int ind;
  Sacs acs = table_acs->acs;

   nacsplusplus_state(Nacs,table_acs);
  acs= table_acs->acs;

  acs[*Nacs].lg=acs[E].lg+coup; 
 acs[*Nacs].pos=acs[acs[E].arc[j]].pos-acs[E].larc[j]+coup;
 if (i<lgrS) {
    ind=lettoind(S[i]);
    acs[*Nacs].arc[ind]=FINAL;   /* Nouvel arc Nacs->FINAL */
    acs[*Nacs].larc[ind]=lgrS-i; /* de i a la fin          */
    }
 ind=lettoind(S[acs[acs[E].arc[j]].pos-acs[E].larc[j]+coup+1]);
 acs[*Nacs].arc[ind]=acs[E].arc[j];        /* Ancien arc Nacs->(E,j) */
 acs[*Nacs].larc[ind]=acs[E].larc[j]-coup;
 acs[E].arc[j]=*Nacs;                      /* Ancien arc E->Nacs */
 acs[E].larc[j]=coup;
}




/*  SlowFind */
int slowfind(char *S, int lgrS, int *i, Pdynarray table_acs, int *Nacs, int *E, int *j)
{int Ej, LEj,      /* Ej = etat pointe par l'arc E-j->, LEj = lgr de l'arc */
     coup=0;

 Sacs acs = table_acs->acs;

 while (*i<lgrS) {
    *j=lettoind(S[*i]);
    Ej=acs[*E].arc[*j];
    if (Ej == INDEFINI) return 0;
    LEj=acs[*E].larc[*j];
    coup=arc_cmp(S, *i, acs[Ej].pos-LEj+1, LEj);
    (*i)+=coup;
    if (coup < LEj) return coup;
    if (acs[*E].lg+LEj < acs[Ej].lg && Ej!=FINAL) {
       dup_etat(S, table_acs, Nacs, *E, *j); /* Dup(E,j)->Nacs, on repart a Nacs */
       acs = table_acs->acs;
    }
    *E=acs[*E].arc[*j];
    }
 return 0;
}


/*****************************/
/*** Construction de l'ACS ***/
Sacs construit_acs (char *S, int lgrS,long int Nmax, int *Nacs, Pdynarray table_acs)
{int i, j,  /* i=indice de parcours de S , j=arc de E */
     coup,  /* coup = lgr de la coupure sur l'arc */
     E;     /* E=etat de parcours de l'ACS */
  Sacs acs;


 /* Initialisation avec I, F, et la source sur l'arc I->F */
 *Nacs=init_acs(S, lgrS, table_acs, Nmax);
 
 acs = table_acs->acs;

 E=INITIAL;
 i=1;
 
 while (i<lgrS) {
    /* Reperage de la coupure (coup) : SlowFind      */
    /* coup==0 => (E,j) n'existe pas                 */
    /* sinon : (E,j) existe et != en S[i] et coup+1  */
    coup=slowfind(S, lgrS, &i, table_acs, Nacs, &E, &j);
    acs = table_acs->acs;

    /* nouvelle lettre ou dernier mot lu est un suffixe => s(FINAL)=E */
    if (!coup) {
       acs[FINAL].s=E;
       if (E==INITIAL) { /* Nouvelle lettre : nouvel arc I->F */
          acs[INITIAL].arc[j]=FINAL;
          acs[INITIAL].larc[j]=lgrS-i;
          i++;
          }
       }

    else { /* Creation du nouvel etat Nacs et calcul du suffixe */
       cree_etat(S, lgrS, i, table_acs, Nacs, E, j, coup);
       acs= table_acs->acs;
       acs[FINAL].s=*Nacs;            /* s(FINAL) = nouvel etat */

       /* Calcul du suffixe de Nacs */
       while (coup) {
          if (E!=INITIAL)
             E=acs[E].s;
          else {
             if (coup==1) break; /* s(Nacs)=INITIAL */
             coup--;
             j=lettoind(S[acs[*Nacs].pos-coup+1]);
             }
          while (coup && acs[E].larc[j]<=coup) { /* Fastfind */
             coup=coup-acs[E].larc[j];
             E=acs[E].arc[j];
             j=lettoind(S[acs[*Nacs].pos-coup+1]);
             }
          if (!coup) break;  /* Suffixe de Nacs trouve */
	  /* Pas de redirection sur un arc non modifiable => cree etat */
          if (acs[E].lg+acs[E].larc[j]==acs[acs[E].arc[j]].lg) {
             cree_etat(S, lgrS, i, table_acs, Nacs, E, j, coup);
	     acs= table_acs->acs;
             acs[(*Nacs)-1].s=*Nacs;
             }
          else { /* sous- suffixe => redirection de l'arc vers Nacs */
             acs[E].arc[j]=*Nacs;
             acs[E].larc[j]=coup;
             }
          }
       acs[*Nacs].s=E;
       }
       
    /* Ajout des arcs inexistants pour le nouveau suffixe */
    if (i<lgrS) {
       j=lettoind(S[i]);
       while (E!=INITIAL && acs[E].arc[j]==INDEFINI) {
          acs[E].arc[j]=FINAL;
          acs[E].larc[j]=lgrS-i;
          E=acs[E].s;
          }
       }
 }

 //nacsplusplus_state(Nacs,table_acs);
 (*Nacs)++;
 return acs;
}

/***************************************/
/*** Calcul du nombre d'arc de l'ACS ***/
int nb_arcs (Pdynarray table_acs, int Nacs)
{int i, j, nb=0;
 Sacs acs = table_acs->acs;

 for (i=0; i<Nacs; i++)
    for (j=0; j<MAXALPHA; j++)
       if (acs[i].arc[j]!=INDEFINI) nb++;
 return nb;
}



/***********************************************************************/
/*************** Construction du fichier .dot **************************/
/***********************************************************************/


/**********************/
/* Affichage d'un arc */
void aff_arc (char *S, Pdynarray table_acs, int E, int arc, FILE *f)
{int i,k;
  Sacs acs = table_acs->acs;

     k=acs[acs[E].arc[arc]].pos-acs[E].larc[arc]+1;
 for (i=0; i<acs[E].larc[arc]; i++) fputc(S[i+k],f);
}


/********************************/
/* Construction du fichier .dot */
void cree_dotfile(char *S, Pdynarray table_acs, int Nacs, char *ref)
{
 int i,j;
 FILE *f;
 Sacs acs = table_acs->acs;

 f=Fopen(ref, (char *)"w"); 

 fputs("digraph G {\nsize = \"7.5,7.5\";\nratio=compress\n",f);
 fputs("node [shape=circle];\nrankdir=LR;\n",f);
 /* Definition des etats */
 for (i=0; i<Nacs; i++) fprintf(f,"{ rank= same; \"%d\"; }\n",i);

 for (i=0; i<Nacs; i++) fprintf(f,"%d [label=\"%d-p%d\"]\n",i,i,acs[i].pos);
 /* Notation des etats terminaux */
 i=FINAL;
 do {
   fprintf(f,"%d [shape=doublecircle];\n",i);
   i=acs[i].s;
 }
 while (i!=INDEFINI);

 /* lien suffixes */
 for (i=1; i<Nacs; i++)
   fprintf(f,"%d -> %d [style=dotted,weight=0];\n",i,acs[i].s);

 /* arcs */
 for (i=0; i<Nacs; i++)
   for (j=0; j<MAXALPHA; j++)
     if (acs[i].arc[j]!=INDEFINI) {
       fprintf(f,"%d -> %d [weight=10000,label=", i, acs[i].arc[j]);
       aff_arc (S, table_acs, i, j, f);
       if(acs[i].lg+acs[i].larc[j]==acs[acs[i].arc[j]].lg)
	 fputs(",style=bold];\n",f);
       else
	 fputs("];\n",f);
     }
 fputs("}\n",f);
 fclose(f);
}


/* Stats on the automaton du fichier .dot */
void count_dotfile(char *S, Pdynarray table_acs, int Nacs, char *ref,int taille)
{
 int i,j;
 long int nbarc,nbreal;

 fprintf(stdout,"Number of states: %d\n",Nacs);

 nbarc=0;
 nbreal=0;
 for (i=0; i<Nacs; i++)
   for (j=0; j<MAXALPHA; j++)
     if ( table_acs->acs[i].arc[j]!=INDEFINI) {
       nbreal++;
       if ( table_acs->acs[table_acs->acs[i].arc[j]].lg!=taille)
	 nbarc++;

     }

 fprintf(stdout,"Number of transitions: %ld total, %ld non vers puit\n",nbreal,nbarc);
}


// code a node and its transition in a buffer, considering that all the nodes target of a  transition 
// have already been coded and that their position in the main char coding is in asc[node].pos, remplacing
// the former one (this saves one int per node)
int code_node(uchar * Buff, Pdynarray table_acs, int Nacs, int node,int taille,int succ=0) {
  int size=1,j,k,toadd;
  int trans[MAXALPHA];
  int codage_trans=0;
  Sacs acs = table_acs->acs;
  /* codage des transitions sur un octet: 0 -> Pas de transitions, 1 -> transitions dans l'arbre (longue max), 2 -> transition autre */ 

  for (j=0; j<MAXALPHA; j++) {
    if (acs[node].arc[j]!=INDEFINI) {
      k = acs[node].arc[j];
      if (acs[node].lg+acs[node].larc[j]==acs[k].lg)  trans[j]=1;  /* strong transition */
      else trans[j]=2;                  /* normal transition */
    } else {
      trans[j]=0;    /* no transition */
    }
    codage_trans+=trans[j]*tobasethree[j];    /* coding in base 3 */
  }

#ifdef DEBUG
  fprintf(stdout,"Coding node %d :",node);
  
  for (j=0; j<MAXALPHA; j++) {
    fprintf(stdout,"[%c - %d]",indtolet(j),trans[j]);
  }
  fprintf(stdout," res: %d\n",codage_trans);
#endif

  if (codage_trans>255) {
    fprintf(stdout,"\nCoding error, trans > 255: %d\n",codage_trans);
    codage_trans=255;
  } 

  if (acs[node].s==INDEFINI) /* terminal node */
    {
      Buff[0]=(uchar) 255;
      Buff[1]=codage_trans;
      size++;
    } else {

  /* first octet: coding the transitions */
  Buff[0]=(uchar) codage_trans;
  }


  /* first octet: coding the target of the transitions if not to the sink */
  for (j=0; j<MAXALPHA; j++) {
    switch (trans[j]) {
    case 1:
    case 2:
      if (acs[acs[node].arc[j]].lg==taille) toadd=0;  /* to the sink */
      else toadd=1;                                   /* not to the think */

      #ifdef DEBUG
      if (toadd==1)
      fprintf(stdout,"trans %c -> toadd %d, pos target %d",indtolet(j),toadd,acs[acs[node].arc[j]].pos);
      else fprintf(stdout,"trans %c -> toadd %d\n",indtolet(j),toadd);
      #endif

      if (succ!=0) {
       if (trans[j]==2)
	size+=all_encode_int((acs[node].larc[j]<<1)+toadd,Buff+size);
	else 
        size+=all_encode_int(acs[node].larc[j],Buff+size);
      } 
      else size+=all_encode_int((acs[node].larc[j]<<1)+toadd,Buff+size);

 #ifdef DEBUG
      fprintf(stdout," length: %d, encoded %d \n",acs[node].larc[j], (acs[node].larc[j]<<1)+toadd);
 #endif
      if (toadd==1 && !(succ==1 && trans[j]==1)) {  /* if once has to code the number of the target */
	size+=all_encode_int(acs[acs[node].arc[j]].pos,Buff+size);  /* not to the think -> encoding the number of the node */

#ifdef DEBUG
      fprintf(stdout," pos target encoded:  %d \n",acs[acs[node].arc[j]].pos);
 #endif

      }
      break;
    }
  }

  return size;  /* returns the number of bytes of the transition */  
}


/***************************/
/*** FONCTION PRINCIPALE ***/
int main( int c, char* v[])
{char *S;
 long int lgrS=0;
 Pdynarray table_acs;
 int Nacs=0,k,tmppos,i;
 int * tabnodes,tmpsize;
 uchar * buf_coding, * petit_buf;
 int size_coding=INITCODING;
 FILE *binFile;

 if (c<2) { 
   fprintf(stderr,"program %s need 1 arguments: <input file> \n", v[0]); 
   exit(-1);
 }

 S=read_source(v[1], &lgrS);

 table_acs = (dynArray *) malloc(sizeof(dynArray));
 table_acs->acs=(Sacs ) malloc (NBITEM*sizeof(sacs));
 table_acs->nb_item=NBITEM;

 

 //table_acs = new dynamic_array<sacs>(NBITEM);

#ifdef DEBUGSEARCH
 fprintf(stdout,"length: %d\n",lgrS);
#endif


 // table_acs->acs=construit_acs(S, lgrS, lgrS+1, &Nacs, table_acs);
 construit_acs(S, lgrS, lgrS+1, &Nacs, table_acs);
 



 if (lgrS<=30)  {
   std::string dotname(v[1]);
   dotname+=".dot";
   cree_dotfile(S, table_acs, Nacs, (char *)dotname.c_str());
 }

  count_dotfile(S, table_acs, Nacs, v[1],lgrS);


  //acs = table_acs->pointcont();

/********************************************************************************************************/
/*************************************  CLASSIC CODING **************************************************/
/********************************************************************************************************/

   
 /* changing suffixes for the coding -1 -> final, to encode */
 k=FINAL;
 do {
   tmppos=table_acs->acs[k].s;
   table_acs->acs[k].s=INDEFINI;
   k=tmppos;
 }
 while (k!=INDEFINI);

 fprintf(stdout,"Sorting the nodes accordingly to their greatest position:\n");
 tabnodes=(int *)malloc(sizeof(int)*Nacs);
 for (k=0;k<Nacs;k++) tabnodes[k]=k;
 for (k=1;k<Nacs;k++) table_acs->acs[k].pos=table_acs->acs[k].pos+1;
 quickSort(tabnodes,0,Nacs-1,table_acs->acs);
 fprintf(stdout,"...done\n");
 
 

#ifdef DEBUG
 fprintf(stdout,"tri: \n");
 for (k=0;k<Nacs;k++) fprintf(stdout,"[%d (%d)]",tabnodes[k],table_acs->acs[tabnodes[k]].pos);
 fprintf(stdout,"\n");
#endif
 
 /* from here, nodes are sorted according to their depth */

 
  
 petit_buf= (uchar* )malloc(SIZEBUFFUCHAR*sizeof(uchar));  /* buffer to code only one node with transitions */
 buf_coding=(uchar* )malloc(size_coding*sizeof(uchar));  /* dynamic buffer for coding the automaton */

 tmppos=POSFINALNODE;     /* 6 bytes first to store the position of the initial node of the automaton */
 tmppos+=all_encode_int(lgrS,buf_coding+tmppos); /* coding the total length of the sequence  as the first node !!*/

 for (k=1;k<Nacs;k++) {
   
   #ifdef DEBUG
   fprintf(stdout,"Coding node %d :",tabnodes[k]);
   #endif


   tmpsize=code_node(petit_buf,table_acs,Nacs,tabnodes[k],lgrS,0);
   chartable_space(&size_coding,&buf_coding,tmppos+tmpsize);   /* increase the final buffer if necessary */
   memcpy(buf_coding+tmppos,petit_buf,tmpsize);   /* copy the result in the buffer */
   table_acs->acs[tabnodes[k]].pos=tmppos;  /* WARNING:  replace in the automata the position, not sure it is a good idea !!  */
   tmppos+=tmpsize;

   #ifdef DEBUG
   fprintf(stdout,"%d bytes\n",tmpsize);
   #endif

 }
 

 /* initial node */
 tmpsize=all_encode_int(table_acs->acs[tabnodes[Nacs-1]].pos,buf_coding);
 if (tmpsize>6) {
   fprintf(stderr,"encoding failed beacuse of the size of the automaton, abort\n");
   exit(-1);
 }



 /************************************************************************************************/
 /*************************************  SAVING COMP FILE ****************************************/
 /************************************************************************************************/



 std::string compname(v[1]);
 compname+=".cdawg-comp";

 std::cout << "Writing the file: " <<  compname << std::endl;

 binFile = fopen((char*)compname.c_str(),"wb");
 fwrite(buf_coding,sizeof(unsigned char),tmppos,binFile);
 fclose(binFile);


 fprintf(stdout,"..done\n");
 

 #ifdef DEBUG
 
   fprintf(stdout,"text compressed : ");
 
 for (int j=0;j<tmppos;j++) {
   fprintf(stdout,"[%d->%d]",j,(int)buf_coding[j]);
 }

fprintf(stdout,"\n------------------------------\n");
#endif


/*******************************************************************************************************/
/*************************************  SUCCINT CODING **************************************************/
/********************************************************************************************************/

 uchar * buf_coding2;
 int * tabnodes2;
 int * children,*coding_order, last_bit;
 int nb_state_stack,nb_int,target,nodetop,nb_in_coding;
 //sdsl::bit_vector thebitvector;
 int length_init=0;
/* renumbering of the nodes*/

/* WARNING ! no more pos in table, will be used to code the new numbers of nodes for their coding*/
  for (k=0;k<Nacs;k++) table_acs->acs[k].pos=0;

  sdsl::bit_vector thebitvector(2*Nacs,0);  /* bit vector initialisation */
  last_bit=0;

  coding_order=(int *)malloc(sizeof(int)*Nacs); /* the order to code the node at the end of this balanced parenthesis process*/
  coding_order[0]=0;
  nb_in_coding=1;

  tabnodes2=(int *)malloc(sizeof(int)*Nacs);
  children=(int *)malloc(sizeof(int)*Nacs);
  tabnodes2[0]=0; /* initial node */
  children[0]=-1;
  nb_state_stack=1;
  nb_int=1;
  table_acs->acs[0].pos=1;
  table_acs->acs[0].pos=0;


   while (nb_state_stack>0) {  // while there is a remaining state in the stack
     /* */
    

      if (children[nb_state_stack-1]<MAXALPHA-1) {

	if (children[nb_state_stack-1]==-1) {
#ifdef DEBUGPARENTHESIS
      fprintf(stdout, "([%d]",tabnodes2[ nb_state_stack-1]);
#endif
      thebitvector[last_bit++]=1;
	  
	}

	children[nb_state_stack-1]+=1;
	nodetop=tabnodes2[nb_state_stack-1];
	target= table_acs->acs[nodetop].arc[children[nb_state_stack-1]];
	if (target!=INDEFINI) { /* there is a transition */
	  if(table_acs->acs[nodetop].lg+table_acs->acs[nodetop].larc[children[nb_state_stack-1]]
	     ==table_acs->acs[target].lg) { /* strong transition */
	    table_acs->acs[target].pos=nb_int++;
	    tabnodes2[nb_state_stack]=target;
	    table_acs->acs[target].pos=nb_in_coding;
	    coding_order[nb_in_coding++]=target;
	    children[nb_state_stack]=-1; 
	    nb_state_stack++; 
	  } 
	} /* there was a transition */

      } else { /* end of there is node */
      
        /* moving one node back from */
#ifdef DEBUGPARENTHESIS
	fprintf(stdout, "[%d]$pos %d)",tabnodes2[ nb_state_stack-1],table_acs->acs[tabnodes2[ nb_state_stack-1]].pos);
#endif
	last_bit++;
       nb_state_stack--;
      }


   } /* end of while stack*/

#ifdef DEBUGPARENTHESIS


      fprintf(stdout, "\n--------------------------------\n");

      fprintf(stdout, "bitvector = ");
	
	for(size_type ii=0; ii<2*Nacs; ++ii)
	  std::cout << thebitvector[ii];
      
      fprintf(stdout, "\n--------------------------------\n");
      for (k=0;k<nb_in_coding;k++) {
	fprintf(stdout, "-%d",coding_order[k]);
      }
      
      fprintf(stdout, "\n--------------------------------\n");



#endif

      // saving the parenthesis in a file
      std::string test(v[1]);
      
      test+=".cdawg-scc-bp";
      std::cout << "Writing the file: " <<  test << std::endl;

      sdsl::store_to_file(thebitvector,test);
    


      /* coding the node in a Elias-Fano structure, at the same time coding the node in a char tables */
      /* reusing children for the position of the nodes (already allocated) */

      petit_buf= (uchar* )malloc(SIZEBUFFUCHAR*sizeof(uchar));  /* buffer to code only one node with transitions */
      buf_coding2=(uchar* )malloc(size_coding*sizeof(uchar));  /* dynamic buffer for coding the automaton */
  

    /* beginning of the buffer is the length of the sequence ! */
    length_init=all_encode_int(lgrS,buf_coding2);
    /* number of the last node !! */
    length_init+=all_encode_int(table_acs->acs[FINAL].pos,buf_coding2+length_init);
    

    tmppos=0;

    for (k=0;k<nb_in_coding;k++) {
	
#ifdef DEBUGSUCCINT
	fprintf(stdout,"Succint node coding %d :",coding_order[k]);
#endif

       tmpsize=code_node(petit_buf,table_acs,Nacs,coding_order[k],lgrS,1); /* succint coding of the node */

       chartable_space(&size_coding,&buf_coding2,tmppos+tmpsize+length_init);   /* increase the final buffer if necessary */
       memcpy(buf_coding2+tmppos+length_init,petit_buf,tmpsize);   /* copy the result in the buffer */
       
       children[k]=tmppos;
       tmppos+=tmpsize;
       
       #ifdef DEBUGSUCCINT
   fprintf(stdout,"%d bytes\n",tmpsize);
   #endif

      }

#ifdef DEBUGSUCCINT
        fprintf(stdout,"positions  :",coding_order[k]);
	
	for (k=0;k<nb_in_coding;k++) {
	  fprintf(stdout,"pos %d :",children[k]);
	}
#endif

	int bit_tot=1;  /* first 0 */
	for (k=1;k<nb_in_coding;k++) {
	  bit_tot+=children[k]-children[k-1];
	}
    
#ifdef DEBUGSUCCINT
        fprintf(stdout,"\n bit-total  : %d\n",bit_tot);
#endif


 sdsl::bit_vector bidiffvctor(bit_tot,0);  /* bit vector initialisation */
  last_bit=1;
  bidiffvctor[0]=1;
  int nb_bits;

  for (k=1;k<nb_in_coding;k++) {
	  nb_bits=children[k]-children[k-1]-1;
	  for(i=0;i<nb_bits;i++) bidiffvctor[last_bit+i]=0;
	  last_bit+=nb_bits;
	  bidiffvctor[last_bit++]=1;
	}


#ifdef DEBUGPARENTHESIS
      fprintf(stdout, "\n--------------------------------\n");

      fprintf(stdout, "bitvector elias = ");
	
	for(size_type ii=0; ii<bit_tot; ++ii)
	  std::cout << bidiffvctor[ii];
      
      fprintf(stdout, "\n--------------------------------\n");
#endif

      sdsl::sd_vector<>  sdb(bidiffvctor);

      std::string test2(v[1]);
      
      test2+=".cdawg-scc-ps";

      std::cout << "Writing the file: " <<  test2 << std::endl;
      
      sdsl::store_to_file(sdb,test2);
      
#ifdef DEBUGSUCCINT
 
   fprintf(stdout,"text compressed succint : ");
 
   for (int j=length_init;j<tmppos+length_init;j++) {
     fprintf(stdout,"[%d->%d]",j,(int)buf_coding2[j]);
   }
#endif

 /***  SAVING SUCC ***/

 std::string succname(v[1]);
 succname+=".cdawg-succ";

 std::cout << "Writing the file: " <<  succname << std::endl;

 binFile = fopen((char *)succname.c_str(),"wb");
 fwrite(buf_coding2,sizeof(unsigned char),tmppos,binFile);
 fclose(binFile);

 fprintf(stdout,"..done\n");
 

 free(table_acs->acs);
 free(table_acs);
 free(tabnodes);
 free(tabnodes2);
 free(children);

}
