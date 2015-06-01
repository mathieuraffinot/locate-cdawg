#ifndef NEW_CDAWGUTILS_H_
#define NEW_CDAWGUTILS_H_


extern "C" {
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

}

#include "integer_coding.h"
#include "dynamic.h"
#include <sys/types.h>
#include <unistd.h>
#include <iostream>


#define INDEFINI -1
#define INITIAL 0
#define FINAL 1
#define MAXALPHA 5 /* alphabet total */
#define TRUE 1
#define FALSE 0

#define NBITEM 10000
#define NBOCC 1000     /* initial number of occurrences for the dynamical list */

#define INITCODING    10000    /* first value of the number of char of the coding buffer (dynamic buffer) */
#define SIZEBUFFUCHAR  200    

#define POSFINALNODE 6        /* beginning of the coding of the final node of the automation in BYTES */



// #define DEBUG 1 
// #define DEBUGSEARCH 1
// #define DEBUGOCC 1
// #define DEBUGPARENTHESIS 1
// #define DEBUGSUCCINT 1
// #define DEBUGSEARCHSUCC 1
// #define DEBUGOCCSUCC 1
// #define DEBUGSEARCHTEST 1

using ulint = uint64_t;
typedef unsigned char uchar;


#define fatal_error(msg) \
   { printf (msg); \
     exit(1); \
   }


/***************************/
/**** Structures du acs ***/
typedef struct {
    int arc[MAXALPHA],
        larc[MAXALPHA],
        s,
        lg,
        pos;
    } sacs;
typedef sacs *Sacs;


/** Dynamic array of states */
typedef struct {
  int nb_item;
  Sacs acs;
} dynArray, * Pdynarray;


/** table of nodes to sort them by pos **/
typedef struct {
  int  nb_nodes;
  int * nodes;
} nodesTable, * PodesTable;


/** Dynamic array of pattern occurrences */
typedef struct {
  int nb_item;
  int nb_occ;
  int * tabocc;
} occArray, * Poccarray;

/** Dynamic array for the stack of states in the search */
typedef struct {
  int nb_in;
  int nb_tot;
  Sacs acs;
  int * children;
  int * lengths;
} dynStack, * Pdynstack;

typedef struct {
  sacs acs;
  int children;
  int lengths;
} locStack, * Plocstack;




/* quicksort on nodes regarding to their pos*/
void echanger(int tableau[], int a, int b);

void quickSort(int tableau[], int debut, int fin, Sacs acs);

char indtolet (int i);

int lettoind (char c);

int lettoindpos (int i, char c);

int fullofn(char * pattern,int size);

char *Malloc(long int taille);

int get_number_of_patterns_bis(std::string header);

int get_patterns_length_bis(std::string header);

void replaceExt(std::string& s, const std::string& newExt);

// Fonction d'ouverture de fichier (Ecr ou Lec mais pas les 2) 
int Open(char *ref, int mode);
FILE *Fopen(char * ref,char * mode);


// dynamical increment of the number of char in the resulting file
  void chartable_space(int * size, uchar **Buff,int nb_bytes);

// Chargement d'une source S qui est renvoyee 
char *read_source (char *v,long int *lgrS);


void header_error_bis();


#endif /* NEW_CDAWGUTILS_H_ */
