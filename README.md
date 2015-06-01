

Locate-cdawg: replacing sampling by CDAWG localisation in BWT indexing approaches.  

===============
Main author: Mathieu Raffinot (m_raffinot@yahoo.com)

Notes: * original CDAWG building code has been written by Renaud
Vérin. The building algorithm is from the paper:

REF1 [Maxime Crochemore and Renaud Vérin. 1997. Direct Construction of
Compact Directed Acyclic Word Graphs. In Proceedings of the 8th Annual
Symposium on Combinatorial Pattern Matching (CPM '97), Alberto
Apostolico and Jotun Hein (Eds.). Springer-Verlag, London, UK, UK,
116-129.]

* Makefile and directory structure has been written by Nicola Prezza
Also this code can be complementary of the LZ-RLCSA code of Nicola Prezza
see: https://github.com/nicolaprezza/lz-rlcsa

This archive is an implementation of the cdawg succint algorithms of the paper: 

REF2: [Djamal Belazzougui, Fabio Cunial, Travis Gagie, Nicola Prezza, Mathieu Raffinot:
Composite repetition-aware data structures. CoRR abs/1502.05937 (2015)]

#####################################################################################
### Brief description


The CDAWG  data structure (see REF1) is a compact automaton built on an
input text and which permits fast pattern localisation. However, coded
directly, the multicharacter transitions requires a huge memory space
which makes it unpractical in practice. One idea of paper REF2 is to
only keep of the CDAWG the structure that allows to find the positions
of a given pattern ** KNOWING ** that the pattern is in the text. This is
garantied by using a BWT first. From the CDAWG, it only remains to code
each first character of each transition, which now permits to encode
the structure in a relatively small amount of space. The localisation
is done very efficiently in worst time O(M+NOCC), where M is the size
of the pattern and NOCC the number of occurrences in the text.
 
Two CDAWG [lite] encodings are proposed: 
* a "simple" BYTE coding
[.cdawg-comp file], which is bigger than the next one, but allows
faster localisations; 
* a more complex succint coding [.cdawg-succ,
.cdawg-scc, .cdawg-ps files] which total sizes are smaller than the
first one, but which localisation remains very fast.

The search code requires that RCLCSA files (for the BWT) have been
built on the original file and are in the same directory.


The RLCSA is based on the following libraries: RLCSA
(https://github.com/adamnovak/rlcsa), sdsl
(https://github.com/simongog/sdsl-lite), and BWTIL
(https://github.com/nicolaprezza/BWTIL)

#####################################################################################
### Download

To clone the repository, use the --recursive option:

> git clone --recursive http://github.com/mathieuraffinot/locate-cdawg

### Compile

The library has been tested under linux using gcc 4.9.2. You need the SDSL library installed on your system (https://github.com/simongog/sdsl-lite).

To compile, just execute

> make

The above command compiles the rlcsa library (this is done only the first time make is executed) and then all executables.

#####################################################################################
### Run

After compiling, run 

>  cdawg_build.x <input.txt>


This command will create the files <input.text>.cdawg-comp, <input.text>.cdawg-succ, <input.text>.cdawg-scc and  <input.text>.cdawg-ps.

Note: when used on a file with less than 30 letters, the program
outputs a .dot file of tha CDAWG automation, which can be visualize
using the [dot] software.


Run

> classic_search_dawg.x <index_basename> <patterns_file>   [BYTE CODE CDAWG]
or
> succ_search_dawg.x <index_basename> <patterns_file>	   [SUCCINT CODDE CDAWG]

to align a set of patterns in pizza&chilli format (http://pizzachili.dcc.uchile.cl/experiments.html) using the CDAWG index.


#####################################################################################
### Examples


# Example on CERE file:

ORIGINAL FILE

440M 	 cere_0.fasta

CDAWG files

138M 	 cere_0.fasta.cdawg-comp

2,4M 	 cere_0.fasta.cdawg-scc-bp
7,6M 	 cere_0.fasta.cdawg-scc-ps
98M 	 cere_0.fasta.cdawg-succ


RLCA Files (RATE 4)

34M 	 cere_0.fasta.rlcsa.array
451M 	 cere_0.fasta.rlcsa.sa_samples

RLCA Files (RATE 8)

34M 	 cere_0.fasta.rlcsa.array
234M 	 cere_0.fasta.rlcsa.sa_samples


RLCA Files (RATE 32)

34M 	 RATE-32/cere_0.fasta.rlcsa.array
62M 	 RATE-32/cere_0.fasta.rlcsa.sa_samples

RLCA Files (RATE 128)

34M 	 RATE-128/cere_0.fasta.rlcsa.array
17M 	 RATE-128/cere_0.fasta.rlcsa.sa_samples



PATTERNS FILES [using genpattern on cere_0.fasta file]

351 	 cere-patterns-3-100     3: length, 100; number of patterns
551 	 cere-patterns-5-100     5: length, 100; number of patterns
751 	 cere-patterns-7-100     7: length, 100; number of patterns
951 	 cere-patterns-9-100     9: length, 100; number of patterns

-------------------------------------------------------------------------------------------
EXAMPLE TIME RESULTS (ON PATTERNS OF LENGTHS 5)

On an:
Intel(R) Pentium(R) CPU G2030 @ 3.00GHz
RAM 16 GB, DDR3 
Running linux kernel 3.13.0-24 (generic)


** classic_search_dawg.x cere_0.fasta cere-patterns-5-100

[Classic] Loading in time (with RCCSA) 0.582655
Ave. nb occur: 595213.437500; total time: 8.822793; ave. time: 0.088228

** succ_search_dawg.x cere_0.fasta cere-patterns-5-100

[Succint] Loading in time 0.628873
Ave. nb occur: 595213.437500; total time: 12.506963; ave. time: 0.125070


RCLCA:


** ../../lz-rlcsa/rlcsa-search cere_0.fasta
RATE-128/PATTERNS/cere-patterns-5-100 - WARNING: SAMPLE_RATE = 128

Load time : 0 seconds.
Search time : 434 seconds

** rlcsa-search RATE-32/cere_0.fasta cere-patterns-5-100 - WARNING:
SAMPLE_RATE = 32 (from https://github.com/nicolaprezza/lz-rlcsa)

Search time : 91 seconds. (1m 31s)

** ../../lz-rlcsa/rlcsa-search cere_0.fasta
RATE-8/PATTERNS/cere-patterns-5-100 - WARNING: SAMPLE_RATE = 8

Search time : 26 seconds.

** ../../lz-rlcsa/rlcsa-search cere_0.fasta
RATE-4/PATTERNS/cere-patterns-5-100 - WARNING: SAMPLE_RATE = 4

Load time : 4 seconds.
Search time : 18 seconds