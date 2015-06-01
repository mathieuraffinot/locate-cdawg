
#include "cdawg_utils.h"


/* for base 3 */
int tobasethree[7]={1,3,9,27,81,243,729};


/* quicksort on nodes regarding to their pos*/
void echanger(int tableau[], int a, int b)

{
    int temp = tableau[a];
    tableau[a] = tableau[b];
    tableau[b] = temp;
}

void quickSort(int tableau[], int debut, int fin, Sacs acs)

{
    int gauche = debut-1;
    int droite = fin+1;
    const int pivot = tableau[debut];

    if(debut >= fin)
        return;

    while(1)
    {
        do droite--; while(acs[tableau[droite]].pos < acs[pivot].pos);
        do gauche++; while(acs[tableau[gauche]].pos > acs[pivot].pos);
        if(gauche < droite)
            echanger(tableau, gauche, droite);
        else break;
    }

    quickSort(tableau, debut, droite,acs);
    quickSort(tableau, droite+1, fin,acs);
}

char indtolet (int i)
{
  switch (i) {
  case 0: return 'a';
  case 1: return 'c';
  case 2: return 'g';
  case 3: return 't';
  }
  return 'n';
}

int lettoind (char c)
{ 
  switch (c) {
  case 'a':
  case 'A': return 0;
  case 'b':
  case 'c':
  case 'C': return 1;
  case 'd':
  case 'g':
  case 'G': return 2;
  case 'e':
  case 't':
  case 'T': return 3;
  }
 
  return 4;
}

/* not used any more */
int lettoindpos (int i, char c)
{
  // if (c>='a' &&  c<='z') return(c-'a');

 switch (c) {
  case 'a':
  case 'A': return 0;
  case 'c':
  case 'C': return 1;
  case 'g':
  case 'G': return 2;
  case 'e':
  case 't':
  case 'T': return 3;
  }
 
 return 4;

 printf("\nErreur lettre_to_indice %d, %c, à la position %d\n",c,c,i); 
 exit(1);
}

/* test if a pattern is full of n */
int fullofn(char * pattern,int size) {
  int k=0;
  while (k<size && lettoind(pattern[k])==4) k++;
  if (k==size) return 1;
  return 0;
}

/***************************************/
/***  Fonction d'allocation memoire ***/
char *Malloc(long int taille)
{char *tmp;
  //  fprintf(stdout,"memoire de %ld\n",taille);
  if ((tmp = (char *) malloc ((int) taille)) == NULL)
    { fprintf(stderr,"ERREUR Malloc d'allocation memoire sur %ld\n",taille); exit(-1); }
  return(tmp);
}

//parse pizza&chilli patterns header:
void header_error_bis(){
  std::cout << "Error: malformed header in patterns file" << std::endl;
  std::cout << "Take a look here for more info on the file format: http://pizzachili.dcc.uchile.cl/experiments.html" << std::endl;
	exit(0);
}



int get_number_of_patterns_bis(std::string header){

	size_t start_pos = header.find("number=");
	if (start_pos == std::string::npos or start_pos+7>=header.size())
		header_error_bis();

	start_pos += 7;

	size_t end_pos = header.substr(start_pos).find(" ");
	if (end_pos == std::string::npos)
		header_error_bis();

	int n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

	return n;

}

int get_patterns_length_bis(std::string header){

	size_t start_pos = header.find("length=");
	if (start_pos == std::string::npos or start_pos+7>=header.size())
		header_error_bis();

	start_pos += 7;

	size_t end_pos = header.substr(start_pos).find(" ");
	if (end_pos == std::string::npos)
		header_error_bis();

	int n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

	return n;

}

void replaceExt(std::string& s, const std::string& newExt) {

  std::string::size_type i = s.rfind('.', s.length());

  if (i != std::string::npos) {
      s.replace(i+1, newExt.length(), newExt);
   }
}

/********************************************************************/
/*** Fonction d'ouverture de fichier (Ecr ou Lec mais pas les 2)  ***/
int Open(char *ref, int mode) 
{int d;
 if ((d=open(ref, (mode==1 ? mode|O_CREAT : mode), 0600)) ==-1)
    { fprintf(stderr,"erreur open fichier %s \n",ref); exit(1); }
 return d;
}

FILE *Fopen(char * ref,char * mode) 
{FILE *f;
 if ((f=fopen(ref, mode)) == NULL)
    { fprintf(stderr,"erreur fopen fichier %s \n",ref); exit(1); }
 return f;
}



// dynamical increment of the number of char in the resulting file
void chartable_space(int * size, uchar **Buff,int nb_bytes) {
  uchar * tmpbuff;

  if (nb_bytes>=*size-1) {
    *size=(int)(1.5*(*size));
    tmpbuff=(uchar *)realloc(*Buff,(*size)*sizeof(uchar));
    //   fprintf(stdout, "Realloc char result %d\n",*size);

     if (!tmpbuff)
                {
                        fprintf(stderr, "ERROR char buffer: Couldn't realloc memory!\n");
                        exit(-1);
                }
    else {*Buff=tmpbuff;}
  }

  return;
}



/**************************************************/
/*** Chargement d'une source S qui est renvoyee ***/
/* lgrS!=0 ==> on ne prend que lgrS bases de la source */
char *read_source (char *v,long int *lgrS)
{char *S;
  int desc,i;

 desc=Open(v,0);
 if (*lgrS==0) *lgrS=lseek(desc, 0L, SEEK_END);
 S = (char *) Malloc (((*lgrS)+1)*sizeof(char));
 lseek(desc, 0L, SEEK_SET);
 read(desc, S, *lgrS);
 while (*lgrS>0 && ((S[*lgrS-1]=='\0') || (S[*lgrS-1]=='\n') ))  (*lgrS)--;
 S[*lgrS]='\0';
 close(desc);

 for (i=0;i<*lgrS-1;i++) lettoindpos(i,S[i]);

 return S;
}
