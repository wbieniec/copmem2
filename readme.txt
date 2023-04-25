copMEM2
============
copMEM2 is an improved version of copMEM (see https://github.com/wbieniec/copmem).
Compared to the original version (=copMEM), the following improvements have been made to improve the performance of the program, especially in the case of a large output file, i.a.:
* mulithreaded processing,
* improved sorting algorithm for MEM's,
* much faster processing when both genomes are very similar (e.g., two individuals of the same species),
* increasing the sampling steps in some cases, which reduces the hash table and decreases the number of iteration steps,
* quicker search for the beginnings of sequences,
* improved structures for frugal storing of MEM's in memory.


Compilation
===========
Linux systems: g++ is required. The compilation command is:
make all

Windows systems: Visual Studio 2017 (2022 is recommended).
Open copmem.sln, switch to Release mode (if needed), and build.


Usage:
======
Linux:
./copmem2 -o <MEMs_file> [options] <R> <Q>

Windows:
copmem2.exe -o <MEMs_file> [options] <R> <Q>

where the obligatory parameters are:
  -o <MEMs_file> output file with the found MEM information,
  <R>            reference file with one or more FASTA sequences,
  <Q>            query file with one or more FASTA sequences,

and the optional ones are:
  -l n           minimal MEM length (the default value is 100). This value must be >= 50,
  -v | -q        verbose mode (prints more details) | quiet mode,
  -r             calculates reverse complement matches. Does not go with -f switch,
  -b             calculates forward and reverse complement matches. Does not go with -f switch,
  -K 36|44|56    determines length of k-mer for hash calculation. Default is 44,
  -e             forces k2 == 1,
  -H 1|2|3|4|5   selects hash function. Default is 1.
     1: maRushPrime1HashSimplified (based on http://www.amsoftware.narod.ru/algo2.html, but simplified);
     2: xxHash32 by Yann Collet (https://github.com/Cyan4973/xxHash);
     3: xxHash64 by Yann Collet (https://github.com/Cyan4973/xxHash);
     4: MetroHash64 by J. Andrew Rogers (https://github.com/jandrewrogers/MetroHash);
     5: CityHash64 from Google (https://github.com/google/cityhash).
  -t N - runs software on N threads including the main thread (default is 1, maximum is 64),
  -mf - switches to the memory-frugal mode. It reduces the number of bits holding the hash from 29 to 28, the textual buffer size is reduced from 2^24 to 2^22 bytes and may reduce the seed size (K) compared to copmem (see the function processCmd in copMEM2.cpp for details).

Note that the switches -t and -mf were not available in copmem!
  -ilb} -- ignore lowercase bases. Not used by default (i.e., lowercase bases are NOT ignored). When set, all lowercase symbols will be treated as N
  
The list of all parameters is available after ./copmem2 -h. Some of those, not listed above, may be used as diagnostic and for tuning the application in terms of speed for some datasets, yet are not normally needed.


Credits:
========
* fmt library (https://github.com/fmtlib/fmt, maintained by Victor Zverovich (vitaut) and Jonathan Müller (foonathan)),
* xxHash32 and xxHash64 (Yann Collet),
* MetroHash64 (J. Andrew Rogers),
* CityHash64 (Google),
* kxsort (https://github.com/voutcn/kxsort, Dinghua Li (voutcn)),
* Daniel Lemire (our LSD-radix based hybrid sort is partly based on https://github.com/lemire/Code-used-on-Daniel-Lemire-s-blog/tree/master/2021/04/09).


MemoryFill:
-----------
MemoryFill is an auxiliary program that enables "cold start", especially for several runs with the same input files.
Usage:

./memoryfill <amount_of_memory>

Use with <amount_of_memory> close to the amount of the installed physical memory. One can use values succeeded by M or G, like 14G.
M = 2^20, G = 2^30.


Demo:
-----
Download and ungzip the genome archives: 
Homo sapiens https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
Mus musculus https://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz

Ungzip files into separate directories and concatenate single sequence files info one file, e.g.:
In Windows:
type mus\*.fa > mus.all.fa
type hum\*.fa > hum.all.fa
In Linux:
cat mus/*.fa > mus.all.fa
cat hum/*.fa > hum.all.fa

Run demo.cmd or demo.sh depending on the operating system.
