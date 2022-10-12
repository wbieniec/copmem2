/* ================================================================ *
*  CopMEM2.cpp : Main file                                          *
*                                                                   *
*  copMEM is a program for efficient computation of MEMs            *
*  (Maximal Exact Matches) in a pair of genomes.                    *
*  Its main algorithmic idea requires that two internal parameters  *
*  (k1 and k2) are coprime, hence the name.                         *
*                                                                   *
*                                                                   *
*  Copyright (c) 2018-2022 Szymon Grabowski and Wojciech Bieniecki  *
*  All rights reserved                                              *
*                                                                   *
*  This program is free software: you can redistribute it and/or    *
*  modify it under the terms of the GNU General Public License as   *
*  published by the Free Software Foundation, either version 3 of   *
*  the License, or (at your option) any later version.              *
*                                                                   *
*  This program is distributed in the hope that it will be useful,  *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of   *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
*  GNU General Public License for more details.                     *
*                                                                   *
*  You should have received a copy of the GNU General Public        *
*  License along with this program.                                 *
*                                                                   *
*  This file is subject to the terms and conditions defined in the  *
*  file 'license', which is part of this source code package.       *
* ================================================================= */

using namespace std;



#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>
#include <cassert>
#include <cmath>
#include <thread>
#include <mutex>

#include "../common/StringUtils.h"
#include "../common/kxsort.h"
#include "../fmt/format.h"


#if (defined(linux) || defined(__linux) || defined(__linux__))
#define _prefetch(x,y) __builtin_prefetch(x,1,(4-y))
#else
#include <xmmintrin.h>
#define _prefetch(x,y) _mm_prefetch(x,y)
#endif

#include "../common/StopWatch.h"
#include "../common/Hashes.h"
#include "../common/StringUtils.h"

// those define's are only for development and debugging
#define radixsort 1
#define stdfmt 0
#define predsv 1


/////////////////// OWN TYPES ///////////////////////////
typedef pair<string, size_t> SequenceItem;
typedef vector<SequenceItem> SequenceVector;
typedef vector<uint32_t> SequenceVectorR;
typedef tuple<string, char*, size_t> SequenceItem2;
typedef vector<SequenceItem2> SequenceVector2;

typedef tuple <size_t, size_t, vector<size_t>> BlockItem;
typedef vector<BlockItem> BlockVector;

template <class MyUINT> using MatchTuple = tuple <uint32_t, MyUINT, uint32_t>; //{gen2_pos, gen1_pos, length}

typedef tuple<size_t, char*> GenomeData; //<size of buffer, memory pointer>

template<class MyUINT1, class MyUINT2> using HashBuffer = pair<MyUINT1*, MyUINT2* >;

enum class verbosity { v0, v1, v2 };
enum class genometype { reference, query };
enum class RC { no, yes, both };
enum class RefSize { small, big, huge };


struct RadixTraits8 {
	//Tuple <4B, 4B, 4B>
	static const int nBytes = 8;

	int kth_byte(const MatchTuple<uint32_t>& x, int k) {
		switch (k) {
		case 0: return get<1>(x) & 0xFF;
		case 1: return get<1>(x) >> 8 & 0xFF;
		case 2: return get<1>(x) >> 16 & 0xFF;
		case 3: return get<1>(x) >> 24 & 0xFF;
		case 4: return get<0>(x) & 0xFF;
		case 5: return get<0>(x) >> 8 & 0xFF;
		case 6: return get<0>(x) >> 16 & 0xFF;
		case 7: return get<0>(x) >> 24 & 0xFF;
		default: return -1;
		}
		return -1;
	}

	bool compare(const MatchTuple<uint32_t>& x, const MatchTuple<uint32_t>& y) {
		return x < y;
	}
};


struct RadixTraits12 {
	//Tuple <4B, 8B, 4B>
	static const int nBytes = 12;

	int kth_byte(const MatchTuple<uint64_t>& x, int k) {
		switch (k) {
		case 0:  return get<1>(x) & 0xFFULL;
		case 1:  return get<1>(x) >> 8 & 0xFFULL;
		case 2:  return get<1>(x) >> 16 & 0xFFULL;
		case 3:  return get<1>(x) >> 24 & 0xFFULL;
		case 4:  return get<1>(x) >> 32 & 0xFFULL;
		case 5:  return get<1>(x) >> 40 & 0xFFULL;
		case 6:  return get<1>(x) >> 48 & 0xFFULL;
		case 7:  return get<1>(x) >> 56 & 0xFFULL;
		case 8:  return get<0>(x) & 0xFFULL;
		case 9:  return get<0>(x) >> 8 & 0xFFULL;
		case 10: return get<0>(x) >> 16 & 0xFFULL;
		case 11: return get<0>(x) >> 24 & 0xFFULL;
		default: return -1;
		}
		return -1;
	}

	bool compare(const MatchTuple<uint64_t>& x, const MatchTuple<uint64_t>& y) {
		return x < y;
	}
};


class NullBuffer : public streambuf {
public:
	int overflow(int c) { return c; }
};


/////////////////// FUNCTIONS ///////////////////////////
void initHashFuncMatrix();
void displayHelp(const char* progName);
void displayHelp2();
void displayParams();
void initGlobals();
void initDefaults();
void processCmd(int argc, char* argv[]);
GenomeData readMultiFasta(string fn, const char paddingChar, bool removeNs, genometype seqType);
void createBlockBuffer(vector<size_t>& seqStarts);
void deleteBlockBuffer();
SequenceVector2 readBlock(int number, ifstream& f, BlockItem& bi, const char paddingChar, bool removeNs);
void deleteReading(GenomeData& r);
template <class MyUINT1, class MyUINT2>
void deleteHashBuffer(HashBuffer<MyUINT1, MyUINT2>& buf);
template <class MyUINT>
void dumpMEM(string& matchesBuffer, SequenceItem& item1, MatchTuple<MyUINT>& match);
void displayMatchInfo(string& name, size_t count);
#if predsv == 1
SequenceItem findSeqDesc(size_t index, SequenceVectorR& seq);
#else
SequenceItem findSeqDesc(size_t index, SequenceVector& seq);
bool lowerBoundComp(SequenceItem lhs, SequenceItem rhs);
#endif


template <class MyUINT1, class MyUINT2>
void cumThread(MyUINT1* sampledPositions, MyUINT2* cum, char* gen, size_t beg, size_t end, bool last);

//////////////////// GLOBAL CONSTS //////////////////////////
constexpr int MAX_THREADS = 64;
constexpr uint32_t MAX_MULTI = 512;
constexpr char REF_PADDING_CHAR = 123;
constexpr char Q_PADDING_CHAR = 125;
constexpr size_t DEF_MATCH_BLOCK_SIZE = 1 << 21;
constexpr size_t MAX_BLOCK_SIZE = 1ULL << 31;
constexpr size_t MATCHES_SORT_T1 = 1ULL << 10; //up to MATCHES_SORT_T1 items are sorted with std::sort


//////////////////// GLOBAL SEMI CONSTS ////////////////////////////
uint32_t MULTI1; //used for prefetching in creation a Hash table
uint32_t MULTI2; //used for prefetching in creation a Hash table
uint32_t MULTI;  //used for prefetching in processQuery
size_t MATCHES_BUFFER_SIZE; //buffer for sorting matches in query phase
uint32_t LONG_MEM;

//////////////////// GLOBAL VARS ////////////////////////////

uint32_t K;  //length of hash
uint32_t H;  //hash type
uint32_t HS; //hash size exponent (0 -> 28, 1 -> 29, 2 -> 30)
uint32_t HASH_SIZE; //hash length
uint32_t L;  //minimal length of MEM
uint32_t k1; //Reference genome step 
uint32_t k2; //Query genome step 
size_t nThreads; //number of threads
mutex mtx;   //lock for shared arrays
RefSize refSize; //computed Reference size: small, big, huge
RefSize forceBR; //forced Reference size

NullBuffer null_buffer;
ostream null_stream(&null_buffer);
ostream* v1logger, * v2logger;

string matchesFN;
string R_FN;
string Q_FN;
ofstream f_matches[MAX_THREADS]; //files 
size_t MATCH_BLOCK_SIZE[MAX_THREADS];
bool blockSortBreakFlag[MAX_THREADS];
vector< MatchTuple<uint32_t>> matchesCopy[MAX_THREADS];

verbosity isVerbose; //quiet, normal or verbose mode
RC isRC;
bool isMemFru;
vector<char*> blockBuffer; //buffer used in buffered reading of some Query sequences (vector of these for supporting threads)
BlockVector blockVector[MAX_THREADS];
SequenceVector r_SV;
SequenceVector q_SV;
#if predsv == 1
SequenceVectorR r_SV_pred;
size_t min_diff;
size_t min_diff_bitshift;
#endif


//////////////////// IMPLEMENTATION ////////////////////////////

void initDefaults() {
	MULTI1 = 128; //used for creation a Hash table
	MULTI2 = 128;
	MULTI = 256; //used in processQuery
	MATCHES_BUFFER_SIZE = 1ULL << 24; //16MB
	LONG_MEM = 4096;
	K = 44;
	H = 1;
	HS = 1;
	L = 100;
	k1 = 8;
	k2 = 7;
	nThreads = 1;
	isRC = RC::no;
	isMemFru = false;
	forceBR = RefSize::small;
	refSize = RefSize::small;
	isVerbose = verbosity::v1;
	v1logger = &cout;
	v2logger = &null_stream;
	//ios_base::sync_with_stdio(false); //may increase I/O operations speed

	//copy_.resize(DEF_MATCH_BLOCK_SIZE);  // 17.08.2022
}


void initGlobals() {
	HASH_SIZE = 1U << (28 + HS);
	hashFunc = hashFuncMatrix[K][H][HS];
	if (isMemFru)
		MATCHES_BUFFER_SIZE = 1ULL << 22;
}


void displayHelp(const char* progName) {
	cout << "copMEM 2.0, by Szymon Grabowski and Wojciech Bieniecki, October 2022." << endl;
	cout << "Usage: " << progName << " [-l n] [-H l] [-K n] [-e] [-t T] [-v|-q] [-mf]|[-b]|[-r] <-o MEMs_file> <Ref_genome> <Query_genome>\n";
	cout << "Attention: -o is a required parameter. l is optional (default: 100).\n";
	cout << "-o MEMs_file - the output file with matches.\n";
	cout << "-v - verbose mode. Display more details.\n";
	cout << "-q - quiet mode. No screen output.\n";
	cout << "-mf - memory frugal mode.\n";
	cout << "-b - compute forward and reverse-complement matches. Not available with -f.\n";
	cout << "-r - compute only reverse-complement matches. Not available with -f.\n";
	cout << "-l n - minimal length of matches. Default value is 100.\n";
	cout << "-H n - hash function. 1: maRushPrime1HashSimplified (default), 2: xxhash32, 3: xxhash64, 4: metroHash64, 5: cityHash64.\n";
	cout << "-K 32|36|40|44|48|52|56|60|64|68|72|76|80|84|88|92|96 Default K is 44.\n";
	cout << "-e - simulates E-MEM (forces k_2 = 1). If you don't understand it, don't use!\n";
	cout << "-t T - maximum number of threads. Default is 1.\n";
	displayHelp2();
}

void displayHelp2() {
	cout << "-k1 n - manually set k_1. Use with -k2 option.\n";
	cout << "-k2 n - manually set k_2. Use with -k1 option.\n";
	cout << "-multi1 n - manually set a buffer size for HT creation. Range: 64 - " << MAX_MULTI << ". Should be a power of two.\n";
	cout << "-multi2 n - manually set a buffer size for HT creation. Range: 64 - " << MAX_MULTI << ". Should be a power of two.\n";
	cout << "-multi n - manually set a buffer size for processing query. Range: 64 - " << MAX_MULTI << ". Should be a power of two.\n";
	cout << "-hash_bits n - manually set a length of the hash. Valid values: 28, 29, 30.\n";
	cout << "-fbr 1|2 - manually force big Ref - long datatypes for processing. 1 - big, 2 - huge.\n";
	cout << "-lm n - long MEM threshold. Default is " << LONG_MEM << ".\n";
}


void displayParams() {
	cout << "CopMEM 2.0 \nPARAMETERS" << endl;
	cout << "Reference filename: " << R_FN << endl;
	cout << "Query filename: " << Q_FN << endl;
	cout << "l = " << L << endl;
	cout << "K = " << K << endl;
	cout << "HASH_SIZE = " << HASH_SIZE << endl;
	cout << "k1 = " << k1 << endl;
	cout << "k2 = " << k2 << endl;
	cout << "Memory frugal mode = " << (isMemFru ? "Yes" : "No") << endl;
	cout << "Reverse Complement = " << ((isRC == RC::both) ? "Both" : (isRC == RC::yes) ? "Yes" : "No") << endl;
	cout << "Threads = " << nThreads << endl;
	cout << "long mem detection = " << LONG_MEM << endl;

	cout << "Hash function: ";
	for (int k = 32; k <= 96; k++) {
		if (hashFunc == hashFuncMatrix[k][1][HS]) { cout << "maRushPrime1HashSimplified\n"; break; }
		if (hashFunc == hashFuncMatrix[k][2][HS]) { cout << "xxhash32\n"; break; }
		if (hashFunc == hashFuncMatrix[k][3][HS]) { cout << "xxhash64\n"; break; }
		if (hashFunc == hashFuncMatrix[k][4][HS]) { cout << "metroHash64\n"; break; }
		if (hashFunc == hashFuncMatrix[k][5][HS]) { cout << "cityHash64\n"; break; }
	}
}

int _gcd(int n, int m) {
	if (n % m == 0)
		return m;
	else
		return _gcd(m, n % m);
}

int gcd(int n, int m) {
	if (n < 1 || m < 1) {
		cerr << "gcd error: the arguments should be positive.\n";
		exit(1);
	}
	if (n < m) {
		int tmp = n;
		n = m;
		m = tmp;
	}
	return _gcd(n, m);
}

void processCmd(int argc, char* argv[]) {
	bool isOset = false;
	bool isEmemLike = false;
	bool isHashBitsSet = false;
	bool k1set = false;
	bool k2set = false;
	bool kset = false;

	const char* incompleteCmd = " option requires one integer argument.\n";
	for (int i = 1; i < argc - 2; ++i) {
		string arg = argv[i];
		if (arg == "-o") {
			if (i + 1 < argc) {
				matchesFN = argv[++i];
				isOset = true;
			}
			else {
				cerr << "-o requires file name.\n";
				exit(1);
			}
		}
		if (arg == "-v") {
			if (isVerbose == verbosity::v0) {
				cerr << "-v and -q parameters are mutually exclusive.\n";
				exit(1);
			}
			isVerbose = verbosity::v2;
			v2logger = &cout;
		}
		if (arg == "-e") {
			isEmemLike = true;
		}
		if (arg == "-q") {
			if (isVerbose == verbosity::v2) {
				cerr << "-v and -q parameters are mutually exclusive.\n";
				exit(1);
			}
			isVerbose = verbosity::v0;
			v1logger = &null_stream;
		}

		if (arg == "-mf") {
			isMemFru = true;
		}

		if (arg == "-b") {
			if (isRC == RC::yes) {
				cerr << "-b and -r parameters are mutually exclusive.\n";
				exit(1);
			}
			isRC = RC::both;
		}
		if (arg == "-r") {
			if (isRC == RC::both) {
				cerr << "-b and -r parameters are mutually exclusive.\n";
				exit(1);
			}
			isRC = RC::yes;
		}
		if (arg == "-h") {
			displayHelp("");
			exit(0);
		}
		if (arg == "-l") {
			if (i + 1 < argc) {
				L = atoi(argv[++i]);
				if (L < 50) {
					cerr << "Incorrect L value (must be >= 50).\n";
					exit(1);
				}
			}
			else {
				cerr << "L" << incompleteCmd;
				exit(1);
			}
		}
		if (arg == "-K") {
			if (i + 1 < argc) {
				K = atoi(argv[++i]);
				if ((K % 4) || (K < 32) || (K > 96)) {
					cerr << "Incorrect K value (must be 32 to 96, and multiple of 4).\n";
					exit(1);
				}
				kset = true;
			}
			else {
				cerr << "K" << incompleteCmd;
				exit(1);
			}
		}
		if (arg == "-H") {
			if (i + 1 < argc) {
				H = atoi(argv[++i]);
				if (H < 1 || H > 5) {
					cerr << "Incorrect H value (must be from 1 to 5).\n";
					exit(1);
				}
			}
			else {
				cerr << "H" << incompleteCmd;
				exit(1);
			}
		}
		if (arg == "-t") {
			if (i + 1 < argc) {
				nThreads = atoi(argv[++i]);
				if (nThreads < 1 || nThreads > MAX_THREADS) {
					cerr << "Incorrect t value (must be from 1 to " << MAX_THREADS << ").\n";
					exit(1);
				}
			}
			else {
				cerr << "t" << incompleteCmd;
				exit(1);
			}
		}
		if (arg == "-k1") {
			if (i + 1 < argc) {
				k1 = atoi(argv[++i]);
				if (k1 == 0) {
					cerr << "k1" << incompleteCmd;
					exit(1);
				}
				k1set = true;
			}
			else {
				cerr << "k1" << incompleteCmd;
				exit(1);
			}
		}
		if (arg == "-k2") {
			if (i + 1 < argc) {
				k2 = atoi(argv[++i]);
				if (k2 == 0) {
					cerr << "k2" << incompleteCmd;
					exit(1);
				}
				k2set = true;
			}
			else {
				cerr << "k2" << incompleteCmd;
				exit(1);
			}
		}
		if (arg == "-multi1") {
			if (i + 1 < argc) {
				MULTI1 = atoi(argv[++i]);
				if (MULTI1 < 64 || MULTI1 > MAX_MULTI) {
					cerr << "MULTI1" << incompleteCmd;
					exit(1);
				}
			}
			else {
				cerr << "MULTI1" << incompleteCmd;
				exit(1);
			}
		}
		if (arg == "-multi2") {
			if (i + 1 < argc) {
				MULTI2 = atoi(argv[++i]);
				if (MULTI2 < 64 || MULTI2 > MAX_MULTI) {
					cerr << "MULTI2" << incompleteCmd;
					exit(1);
				}
			}
			else {
				cerr << "MULTI2" << incompleteCmd;
				exit(1);
			}
		}
		if (arg == "-multi") {
			if (i + 1 < argc) {
				MULTI = atoi(argv[++i]);
				if (MULTI < 64 || MULTI > MAX_MULTI) {
					cerr << "MULTI" << incompleteCmd;
					exit(1);
				}
			}
			else {
				cerr << "MULTI" << incompleteCmd;
				exit(1);
			}
		}
		if (arg == "-lm") {
			if (i + 1 < argc) {
				LONG_MEM = atoi(argv[++i]);
				if (LONG_MEM < 400 || LONG_MEM > 10000) {
					cerr << "lm" << incompleteCmd;
					exit(1);
				}
			}
			else {
				cerr << "MULTI" << incompleteCmd;
				exit(1);
			}
		}


		if (arg == "-hash_bits") {
			if (i + 1 < argc) {
				int hb = atoi(argv[++i]);
				if (hb < 28 || hb > 30) {
					cerr << "-hash_bits" << incompleteCmd;
					exit(1);
				}
				isHashBitsSet = true;
				HS = hb - 28;
			}
			else {
				cerr << "-hash_bits" << incompleteCmd;
				exit(1);
			}

		}
		if (arg == "-fbr") {
			if (i + 1 < argc) {
				int fbr = atoi(argv[++i]);
				switch (fbr) {
				case 1:forceBR = RefSize::big; break;
				case 2:forceBR = RefSize::huge; break;
				default:cerr << "-fbr bad value" << incompleteCmd; exit(1);
				}
			}
			else {
				cerr << "fbr" << incompleteCmd;
				exit(1);
			}
		}

	}

	if (isOset == false) {
		cerr << "-o not given or specified correctly.\n";
		exit(1);
	}


	//touching files:
	R_FN = argv[argc - 2];
	Q_FN = argv[argc - 1];
	ifstream f;

	f.open(R_FN);
	if (f.fail()) {
		cerr << "\nReference file '" << R_FN << "' does not exist. Quit.";
		exit(1);
	}
	f.close();
	f.open(Q_FN);
	if (f.fail()) {
		cerr << "\nQuery file '" << Q_FN << "' does not exist. Quit.";
		exit(1);
	}
	f.close();
	f_matches[0].open(matchesFN, ios::binary);
	if (f_matches[0].fail()) {
		cerr << "\nFAILED. The -o parameter specifies a file that cannot be created.\n";
		exit(1);
	}
	f_matches[0].close();

	/* optimizing K */
	if (kset) {
		if (L - K + 1 <= 0) {
			cerr << "\nL and K mismatch.\n";
			exit(1);
		}
	}
	else {
		if (!isMemFru) {
			if (L >= 200)
				K = 56;
			else
				if (L >= 80)
					K = 44;
				else
					K = 36;
		}
		else {
			if (L >= 200)
				K = 44;
			else
				K = 36;
		}
	}
	if (!isHashBitsSet)
		if (isMemFru)
			HS = 0;
		else
			HS = 1;



	int tempVar = L - K + 1;

	/* setting k1 and k2 */
	if ((k1set && !k2set) || (!k1set && k2set)) {
		cerr << "\nk1 and k2 parameters must come together.\n";
		exit(1);
	}
	if (k1set && k2set) {
		if (isEmemLike) {
			cerr << "\nk1, k2 and e parameters cannot come together.\n";
			exit(1);
		}
		if (k1 * k2 > (uint32_t)tempVar) {
			cerr << "\nk1 and k2 parameters mismatch. k1*k2 > L-K+1\n";
			exit(1);
		}
		if (gcd(k1, k2) != 1) {
			cerr << "\nk1 and k2 parameters mismatch. k1 and k2 are NOT relatively prime!\n";
			exit(1);
		}
	}
	else {
		k1 = (int)(pow(tempVar, 0.5)) + 1;
		k2 = k1 - 1;
		if (k1 * k2 > (uint32_t)tempVar) {
			--k1;
			--k2;
		}
		uint32_t temp_k1 = k1;
		uint32_t best_k1 = k1;
		while (true) {
			++temp_k1;
			bool condition1 = (temp_k1 * k2 <= (uint32_t)tempVar);
			bool condition2 = (gcd(temp_k1, k2) == 1);
			if (condition1 && condition2)
				best_k1 = temp_k1;
			if (!condition1)
				break;
		}
		k1 = best_k1;
	}

	if (isEmemLike == true) {
		k1 = L - K + 1;
		k2 = 1;
	}
}

template <class MyUINT2>
void gencumTh(MyUINT2* cum, char* gen, size_t beg, size_t end, bool last) {
	const size_t k1MULTI1 = k1 * MULTI1;

	uint32_t hashPositions[MAX_MULTI];
	size_t i;

	for (i = beg; i + K + k1MULTI1 < end; i += k1MULTI1) {
		char* tempPointer = gen + i;
		for (size_t temp = 0; temp < MULTI1; ++temp) {
			hashPositions[temp] = hashFunc(tempPointer) + 2;
			tempPointer += k1;
			_prefetch((char*)(cum + hashPositions[temp]), 1);
		}
		mtx.lock();
		for (size_t temp = 0; temp < MULTI1; ++temp) {
			++cum[hashPositions[temp]];
		}
		mtx.unlock();
	}
	//////////////////// processing the end part of R  //////////////////////
	if (last)
		end = end - K + 1;
	mtx.lock();
	for (; i < end; i += k1) {
		uint32_t h = hashFunc(gen + i) + 2;
		++cum[h];
	}
	mtx.unlock();
}


template <class MyUINT2>
void gencum(GenomeData& genome, MyUINT2* cum) {

	size_t N = get<0>(genome);
	char* gen = get<1>(genome);

	fill(cum, cum + HASH_SIZE + 2, 0);

	size_t N1 = N / (nThreads * k1) * k1;

	thread* th = new thread[nThreads - 1];

	size_t beg, end;
	size_t t;
	if (refSize == RefSize::huge) {
		for (t = 0; t < nThreads - 1; t++) {
			beg = t * N1; end = (t + 1) * N1;
			th[t] = thread(gencumTh<uint64_t>, (uint64_t*)cum, gen, beg, end, false);
		}
		//last thread == current thread
		beg = (nThreads - 1) * N1; end = N;
		gencumTh<MyUINT2>(cum, gen, beg, end, true);
	}
	else {
		for (t = 0; t < nThreads - 1; t++) {
			beg = t * N1; end = (t + 1) * N1;
			th[t] = thread(gencumTh<uint32_t>, (uint32_t*)cum, gen, beg, end, false);
		}
		beg = (nThreads - 1) * N1; end = N;
		gencumTh<MyUINT2>(cum, gen, beg, end, true);
	}

	for (t = 0; t < nThreads - 1; t++)
		th[t].join();

	//////////////////// processing the end part of R //////////////////////
	partial_sum(cum, cum + HASH_SIZE + 2, cum);
}


template <class MyUINT1, class MyUINT2>
void cumThread(MyUINT1* sampledPositions, MyUINT2* cum, char* gen, size_t beg, size_t end, bool last) {
	const unsigned int k1MULTI2 = k1 * MULTI2;

	uint32_t hashPositions[MAX_MULTI];
	size_t i1;
	for (i1 = beg; i1 + K + k1MULTI2 < end; i1 += k1MULTI2) {
		char* tempPointer = gen + i1;
		for (unsigned int temp = 0; temp < MULTI2; ++temp) {
			hashPositions[temp] = hashFunc(tempPointer) + 1;
			tempPointer += k1;
			_prefetch((char*)(cum + hashPositions[temp]), 1);
		}

		for (unsigned int temp = 0; temp < MULTI2; ++temp) {
			_prefetch((char*)(sampledPositions + *(cum + hashPositions[temp])), 1);
		}

		size_t i2 = i1;
		mtx.lock();

		if (refSize == RefSize::big) {
			MyUINT1 i2_over_k1 = (MyUINT1)(i2 / k1);
			for (unsigned int temp = 0; temp < MULTI2; ++temp) {
				sampledPositions[cum[hashPositions[temp]]] = i2_over_k1;
				++cum[hashPositions[temp]];
				//i2 += k1;
				++i2_over_k1;
			}
		}
		else {
			for (unsigned int temp = 0; temp < MULTI2; ++temp) {
				sampledPositions[cum[hashPositions[temp]]] = (MyUINT1)i2;
				++cum[hashPositions[temp]];
				i2 += k1;
			}
		}
		mtx.unlock();
	}
	//////////////////// processing the end part of R //////////////////////
	if (last)
		end = end - K + 1;
	mtx.lock();

	if (refSize == RefSize::big) {
		for (; i1 < end; i1 += k1) {
			uint32_t h = hashFunc(gen + i1) + 1;
			sampledPositions[cum[h]] = (MyUINT1)(i1 / k1);
			++cum[h];
		}
	}
	else {
		for (; i1 < end; i1 += k1) {
			uint32_t h = hashFunc(gen + i1) + 1;
			sampledPositions[cum[h]] = (MyUINT1)i1;
			++cum[h];
		}
	}
	mtx.unlock();
}


template<class MyUINT1, class MyUINT2>
HashBuffer<MyUINT1, MyUINT2> processRef(GenomeData& rGenome) {
	CStopWatch stopwatch;
	stopwatch.start();
	size_t N = get<0>(rGenome);
	char* gen = get<1>(rGenome);

	const size_t hashCount = (N - K + 1) / k1;
	*v1logger << "Hash count = " << hashCount << endl;

	MyUINT1* sampledPositions = new MyUINT1[hashCount + 2];
	MyUINT2* cum = new MyUINT2[HASH_SIZE + 2];


	*v2logger << "\tprocessRef: init = " << stopwatch.stop() << endl;
	stopwatch.resume();
	gencum(rGenome, cum);

	*v2logger << "\tprocessRef: cum(1) = " << stopwatch.stop() << endl;
	stopwatch.resume();


	size_t N1 = N / (nThreads * k1) * k1;

	thread* th = new thread[nThreads - 1];
	int bigref = 2;
	if (sizeof(MyUINT2) == sizeof(uint32_t))
		bigref = 1;
	if (sizeof(MyUINT1) == sizeof(uint32_t))
		bigref = 0;
	size_t beg, end;
	int t;
	if (bigref == 2) {
		for (t = 0; t < nThreads - 1; t++) {
			beg = t * N1;  end = (t + 1) * N1;
			th[t] = thread(cumThread<uint64_t, uint64_t>, (uint64_t*)sampledPositions, (uint64_t*)cum, gen, beg, end, false);
		}
		//last thread == current thread
		beg = (nThreads - 1) * N1; end = N;
		cumThread(sampledPositions, cum, gen, beg, end, true);
	}
	if (bigref == 1) {
		for (t = 0; t < nThreads - 1; t++) {
			beg = t * N1; end = (t + 1) * N1;
			th[t] = thread(cumThread<uint32_t, uint32_t>, (uint32_t*)sampledPositions, (uint32_t*)cum, gen, beg, end, false);
		}
		beg = (nThreads - 1) * N1;
		end = N;
		cumThread(sampledPositions, cum, gen, beg, end, true);
	}
	if (bigref == 0) {
		for (t = 0; t < nThreads - 1; t++) {
			beg = t * N1; end = (t + 1) * N1;
			th[t] = thread(cumThread<uint32_t, uint32_t>, (uint32_t*)sampledPositions, (uint32_t*)cum, gen, beg, end, false);
		}
		beg = (nThreads - 1) * N1;
		end = N;
		cumThread(sampledPositions, cum, gen, beg, end, true);
	}

	for (int t = 0; t < nThreads - 1; t++) th[t].join();
	delete[]th;
	*v2logger << "\tprocessRef: cum(2) = " << stopwatch.stop() << endl;
	return { sampledPositions, cum };
}


#if stdfmt == 1

template <class MyUINT>
void dumpMEM(string& matchesBuffer, SequenceItem& item1, MatchTuple<MyUINT>& match) {
	MyUINT baseindex1 = get<1>(match) - (MyUINT)item1.second;
	MyUINT baseindex2 = get<0>(match);
	matchesBuffer.append(" ");
	matchesBuffer.append(item1.first);
	matchesBuffer.append("\t");
	matchesBuffer.append(to_string(baseindex1));
	matchesBuffer.append("\t");
	matchesBuffer.append(to_string(baseindex2));
	matchesBuffer.append("\t");
	matchesBuffer.append(to_string(get<2>(match)));
	matchesBuffer.append("\n");
}

#else
template <class MyUINT>
void dumpMEM(string& matchesBuffer, SequenceItem& item1, MatchTuple<MyUINT>& match) {
	MyUINT baseindex1 = get<1>(match) - (MyUINT)item1.second;
	MyUINT baseindex2 = get<0>(match);
	MyUINT len = get<2>(match);
	matchesBuffer.append(" ");
	matchesBuffer.append(item1.first);
	matchesBuffer.append("\t");
	matchesBuffer.append(fmt::format_int(baseindex1).c_str());
	matchesBuffer.append("\t");
	matchesBuffer.append(fmt::format_int(baseindex2).c_str());
	matchesBuffer.append("\t");
	matchesBuffer.append(fmt::format_int(len).c_str());
	matchesBuffer.append("\n");
}
#endif

#if predsv == 1
SequenceItem findSeqDesc(size_t index, SequenceVectorR& seq) {
	size_t tmp_index = index >> min_diff_bitshift;

	if ((tmp_index < seq.size() - 1) && (r_SV[seq[tmp_index + 1]].second <= index))
		return r_SV[seq[tmp_index + 1]];
	else
		return r_SV[seq[tmp_index]];
}
#else
bool lowerBoundComp(SequenceItem lhs, SequenceItem rhs) {
	return lhs.second < rhs.second;
}
SequenceItem findSeqDesc(size_t index, SequenceVector& seq) {
	SequenceItem dummySequenceItem = { "", index };
	SequenceItem item = seq[0];
	auto lower = lower_bound(seq.begin(), seq.end(), dummySequenceItem, lowerBoundComp);
	size_t diff = lower - seq.begin();
	return seq[diff - 1];
}
#endif

void displayMatchInfo(string& name, size_t count) {
	switch (count) {
	case 0:  *v2logger << name << ": no matches.\n"; break;
	case 1:  *v2logger << name << ": 1 match.\n"; break;
	default: *v2logger << name << ": " << count << " matches.\n"; break;
	}
}

void writeSeqHeader(string& matchesBuffer, string& seqName) {
	matchesBuffer.append("> ");
	matchesBuffer.append(seqName);
	matchesBuffer.append("\n");
}

void flushMatchesBuffer(ofstream& fmatches, string& matchesBuffer) {
	fmatches << matchesBuffer;
	matchesBuffer.clear();
}

template <class MyUINT>
void LSDRadixSort_StdSort_tuple(size_t thID, vector<MatchTuple<MyUINT>>& data, size_t offset) {

	vector< MatchTuple<uint32_t>>& copy_ = matchesCopy[thID];  // 3.08.2022

	size_t length = data.size();
	if (copy_.size() < length)
		copy_.resize(length);

	uint32_t level2[257];
	uint32_t level3[257];
	uint32_t level4[257];
	uint32_t level5[257];
	uint32_t level6[257];
	uint32_t level7[257];

	fill(level2, level2 + 257, 0);
	fill(level3, level3 + 257, 0);
	fill(level4, level4 + 257, 0);
	fill(level5, level5 + 257, 0);
	fill(level6, level6 + 257, 0);
	fill(level7, level7 + 257, 0);

	for (size_t i = offset; i < length; ++i) {
		uint32_t value0 = get<0>(data[i]);
		uint32_t value1 = get<1>(data[i]);

		level2[((value1 >> 16) & 0xFF) + 1]++;
		level3[((value1 >> 24) & 0xFF) + 1]++;

		level4[((value0) & 0xFF) + 1]++;
		level5[((value0 >> 8) & 0xFF) + 1]++;
		level6[((value0 >> 16) & 0xFF) + 1]++;
		level7[((value0 >> 24) & 0xFF) + 1]++;
	}

	for (size_t i = 1; i < 257; ++i) {
		level2[i] += level2[i - 1];
		level3[i] += level3[i - 1];
		level4[i] += level4[i - 1];
		level5[i] += level5[i - 1];
		level6[i] += level6[i - 1];
		level7[i] += level7[i - 1];
	}
	for (size_t i = offset; i < length; ++i)
		copy_[offset + level2[(get<1>(data[i]) >> 16) & 0xFF]++] = data[i];
	for (size_t i = offset; i < length; ++i)
		data[offset + level3[(get<1>(copy_[i]) >> 24) & 0xFF]++] = copy_[i];

	for (size_t i = offset; i < length; ++i)
		copy_[offset + level4[(get<0>(data[i])) & 0xFF]++] = data[i];
	for (size_t i = offset; i < length; ++i)
		data[offset + level5[(get<0>(copy_[i]) >> 8) & 0xFF]++] = copy_[i];
	for (size_t i = offset; i < length; ++i)
		copy_[offset + level6[(get<0>(data[i]) >> 16) & 0xFF]++] = data[i];
	for (size_t i = offset; i < length; ++i)
		data[offset + level7[(get<0>(copy_[i]) >> 24) & 0xFF]++] = copy_[i];

	size_t beg = offset;
	size_t cur = beg + 1;
	while (cur < length) {
		if (((get<1>(data[cur]) & 0xFFFF0000) != (get<1>(data[beg]) & 0xFFFF0000)) || (get<0>(data[cur]) != get<0>(data[beg]))) {
			if (cur - beg == 1) {
				beg = cur;
				continue;
			}
			std::sort(data.begin() + beg, data.begin() + cur);
			beg = cur;
			++cur;
		}
		else
			++cur;
	}
	std::sort(data.begin() + beg, data.begin() + length);
}


template <class MyUINT>
void postProcess(int thID, vector<MatchTuple<MyUINT>>& matches, string& matchesBuffer, string seqName) {
	//Warning - first match is a guard match. Do not touch!

	if (matches.size() == 1) [[unlikely]] {
		displayMatchInfo(seqName, 0);
		return;
	}

		if (matches.size() > 2) {
			auto itFirstMatch = matches.begin() + 1;
			auto itLastMatch = matches.end();

#if radixsort == 1
			if (sizeof(MyUINT) == sizeof(uint32_t)) {
				if (matches.size() > MATCHES_SORT_T1)
					if (nThreads <= 2)
						LSDRadixSort_StdSort_tuple(thID, matches, 1);
					else
						kx::radix_sort(itFirstMatch, itLastMatch, RadixTraits8());
				else
					std::sort(itFirstMatch, itLastMatch);
			}
			else
			{
				if (matches.size() > MATCHES_SORT_T1)
					kx::radix_sort(itFirstMatch, itLastMatch, RadixTraits12());
				else
					std::sort(itFirstMatch, itLastMatch);
			}
#else
			std::sort(itFirstMatch, itLastMatch);
#endif
		}
	//now matches are sorted. If within OVERLAP there is any less than the guard match, they cannot be properly sorted. Procedure must terminate
//Enough to examine the first one.
	if (matches[0] >= matches[1]) {
		blockSortBreakFlag[thID] = true;
		return;
	}

	SequenceItem seq1item;
	MatchTuple<MyUINT> prev_match = matches[1];
#if predsv == 1
	seq1item = findSeqDesc(get<1>(prev_match), r_SV_pred);
#else
	seq1item = findSeqDesc(get<1>(prev_match), r_SV);
#endif
	dumpMEM(matchesBuffer, seq1item, prev_match);
	uint64_t matchcount = 1ULL;
	for (auto itMatch = matches.begin() + 1; itMatch != matches.end(); ++itMatch) {
		if (prev_match != *itMatch) {
			if (get<1>(prev_match) != get<1>(*itMatch)) {
#if predsv == 1
				seq1item = findSeqDesc(get<1>(*itMatch), r_SV_pred);
#else
				seq1item = findSeqDesc(get<1>(*itMatch), r_SV);
#endif
			}
			dumpMEM(matchesBuffer, seq1item, *itMatch);
			++matchcount;
			if (matchesBuffer.size() > (size_t)(MATCHES_BUFFER_SIZE * 0.95))
				flushMatchesBuffer(f_matches[thID], matchesBuffer);
		}
		else {
		}
		prev_match = *itMatch;
	}
	displayMatchInfo(seqName, matchcount);
	matches.erase(matches.begin() + 1, matches.end());
	matches[0] = { 0,0,0 };
}


template <class MyUINT>
void postProcessPartial(int thID, vector<MatchTuple<MyUINT>>& matches, MatchTuple<MyUINT> match, string& matchesBuffer, string seqName) {
	if (matches.back() == match) {
		return;
	}
	matches.push_back(match);

	if (matches.size() < MATCH_BLOCK_SIZE[thID]) [[likely]]
		return;
	if (blockSortBreakFlag[thID] == true)
		return;

	size_t OVERLAP = 3ULL * MATCH_BLOCK_SIZE[thID] / 4ULL;
	auto itFirstMatch = matches.begin() + 1;
	auto itLastMatch = matches.end();
#if radixsort == 1
	if (sizeof(MyUINT) == sizeof(uint32_t)) {
		if (matches.size() > MATCHES_SORT_T1)
			if (nThreads <= 2)
				LSDRadixSort_StdSort_tuple(thID, matches, 1);
			else
				kx::radix_sort(itFirstMatch, itLastMatch, RadixTraits8());
		else
			std::sort(itFirstMatch, itLastMatch);
	}
	else
	{
		if (matches.size() > MATCHES_SORT_T1)
			kx::radix_sort(itFirstMatch, itLastMatch, RadixTraits12());
		else
			std::sort(itFirstMatch, itLastMatch);
	}
#else
	std::sort(itFirstMatch, itLastMatch);
#endif

	//now matches are sorted. If within OVERLAP there is any less than the guard match, they cannot be properly sorted. Procedure must terminate
	//Enough to examine the first one.
	if (matches[0] >= matches[1]) {
		blockSortBreakFlag[thID] = true;
		return;
	}
#if predsv == 1
	SequenceItem seq1item = findSeqDesc(get<1>(matches[1]), r_SV_pred);
#else
	SequenceItem seq1item = findSeqDesc(get<1>(matches[1]), r_SV);
#endif
	dumpMEM(matchesBuffer, seq1item, matches[1]);
	size_t i = 2;
	for (; i < OVERLAP; ++i) {
		if (matches[i - 1] != matches[i]) {
			if (get<1>(matches[i - 1ULL]) != get<1>(matches[i])) {
#if predsv == 1
				seq1item = findSeqDesc(get<1>(matches[i]), r_SV_pred);
#else
				seq1item = findSeqDesc(get<1>(matches[i]), r_SV);
#endif
			}
			dumpMEM(matchesBuffer, seq1item, matches[i]);
			if (matchesBuffer.size() > (size_t)(MATCHES_BUFFER_SIZE * 0.95))
				flushMatchesBuffer(f_matches[thID], matchesBuffer);
		}
		else {
		}
	}

	for (; (matches[i - 1] == matches[i]) && (i < MATCH_BLOCK_SIZE[thID]); ++i); //look for possible duplicates at the end

	matches[0] = matches[i - 1]; //update guard match. On next part, the last erased should not be greater as the guard.
	matches.erase(matches.begin() + 1, matches.begin() + i);
}


template<class MyUINT1, class MyUINT2>
//MyUINT1 - size of sampledPositions element: small:32b, big:32b, huge:64b
//MyUINT2 - size of cum element: small:32b, big:32b, huge:64b
void processQueryTh(int thID, char* gen1, char* gen2, MyUINT2* cum, MyUINT1* sampledPositions, size_t beg, size_t end, string& matchesBuffer, vector<MatchTuple<uint32_t>>& matches12, vector<MatchTuple<uint64_t>>& matches16, string seqName, bool last) {
	if (end - beg < L)
		return;

	size_t i1;
	const unsigned int k2MULTI = k2 * MULTI;

	char* curr1, * curr2;
	uint32_t hArray[MAX_MULTI];
	MyUINT2 posArray[MAX_MULTI * 2];

	const int L_PLUS_ONE = L + 1;
	const int LK2 = (L - K) / 2;
	const int LK2_MINUS_4 = LK2 - 4;
	const int K_PLUS_LK24 = K + LK2_MINUS_4;
	uint32_t l1, l2, r1, r2;

	bool isLongRightMatch = false;
	ptrdiff_t difference;
	char* rightQlm = nullptr;

	size_t _k1 = (refSize == RefSize::big ? k1 : 1); //multiplier for big entries in sampledPositions

	for (i1 = beg; i1 + K + k2MULTI < end + 1; i1 += k2MULTI) {
		if (blockSortBreakFlag[thID])
			return;

		size_t i1temp = i1;
		curr2 = gen2 + i1;
		for (size_t i2 = 0; i2 < MULTI; ++i2) {
			hArray[i2] = hashFunc(curr2);
			curr2 += k2;
		}
		for (size_t i2 = 0; i2 < MULTI; ++i2) {
			memcpy(posArray + i2 * 2, cum + hArray[i2], sizeof(MyUINT2) * 2);
		}

		curr2 = gen2 + i1;
		for (size_t i2 = 0; i2 < MULTI; ++i2) {
			if (posArray[i2 * 2] == posArray[i2 * 2 + 1]) {
				curr2 += k2;
				continue;
			}

			memcpy(&l2, curr2 - LK2, sizeof(uint32_t));
			memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(uint32_t));

			for (MyUINT2 j = posArray[i2 * 2]; j < posArray[i2 * 2 + 1]; ++j) {
				curr1 = gen1 + sampledPositions[j] * _k1;

				if (isLongRightMatch == false) {
					memcpy(&l1, curr1 - LK2, sizeof(uint32_t));
					memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(uint32_t));
					if (r1 == r2 || l1 == l2) {
						char* p1 = curr1 + K;
						char* p2 = curr2 + K;
						while (*p1 == *p2) {
							++p1;
							++p2;
						}
						char* rightR = p1;
						char* rightQ = p2;
						p1 = curr1 - 1;
						p2 = curr2 - 1;
						while (*p1 == *p2) {
							--p1;
							--p2;
						}
						if (rightR - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
							if (refSize == RefSize::small) {
								MatchTuple<uint32_t> tempMatch{ (uint32_t)(p2 + 1 - gen2),(uint32_t)(p1 + 1 - gen1),(uint32_t)(rightR - p1 - 1) };
								postProcessPartial<uint32_t>(thID, matches12, tempMatch, matchesBuffer, seqName);
							}
							else {
								MatchTuple<uint64_t> tempMatch{ (uint32_t)(p2 + 1 - gen2),(uint64_t)(p1 + 1 - gen1),(uint32_t)(rightR - p1 - 1) };
								postProcessPartial<uint64_t>(thID, matches16, tempMatch, matchesBuffer, seqName);
							}

							if (rightQ - curr2 - 2 >= LONG_MEM) {
								isLongRightMatch = true;
								difference = p1 - p2;
								rightQlm = rightQ;
							}
							else {
								isLongRightMatch = false;
								rightQlm = nullptr;
							}
						}
						else {
						}
					}
				}
				else {
					if (curr2 + K <= rightQlm) {
						if (curr1 - curr2 == difference) {
							continue;
						}
					}
					else {
						isLongRightMatch = false;
						rightQlm = nullptr;
					}
					memcpy(&l1, curr1 - LK2, sizeof(uint32_t));
					memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(uint32_t));
					if (r1 == r2 || l1 == l2) {
						char* p1 = curr1 + K;
						char* p2 = curr2 + K;
						while (*p1 == *p2) {
							++p1;
							++p2;
						}
						char* rightR = p1;
						p1 = curr1 - 1;
						p2 = curr2 - 1;
						while (*p1 == *p2) {
							--p1;
							--p2;
						}
						if (rightR - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
							if (refSize == RefSize::small) {
								MatchTuple<uint32_t> tempMatch{ (uint32_t)(p2 + 1 - gen2),(uint32_t)(p1 + 1 - gen1),(uint32_t)(rightR - p1 - 1) };
								postProcessPartial<uint32_t>(thID, matches12, tempMatch, matchesBuffer, seqName);
							}
							else {
								MatchTuple<uint64_t> tempMatch{ (uint32_t)(p2 + 1 - gen2),(uint64_t)(p1 + 1 - gen1),(uint32_t)(rightR - p1 - 1) };
								postProcessPartial<uint64_t>(thID, matches16, tempMatch, matchesBuffer, seqName);
							}
						}
						else {
							//nop
						}
					}
				}
			}
			curr2 += k2;
		}
	}
	//////////////////// processing the end part of Q  //////////////////////
	if (last)
		end = end - K + 1;
	for (; i1 < end; i1 += k2) {
		curr2 = gen2 + i1;
		memcpy(posArray, cum + hashFunc(curr2), sizeof(MyUINT2) * 2);

		if (posArray[0] == posArray[1]) {
			continue;
		}

		memcpy(&l2, curr2 - LK2, sizeof(uint32_t));
		memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(uint32_t));

		for (MyUINT2 j = posArray[0]; j < posArray[1]; ++j) {
			curr1 = gen1 + sampledPositions[j] * _k1;
			memcpy(&l1, curr1 - LK2, sizeof(uint32_t));
			memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(uint32_t));

			if (r1 == r2 || l1 == l2) {
				char* p1 = curr1 + K;
				char* p2 = curr2 + K;

				while (*p1 == *p2) {
					++p1;
					++p2;
				}
				char* right = p1;

				p1 = curr1 - 1;
				p2 = curr2 - 1;
				while (*p1 == *p2) {
					--p1;
					--p2;
				}

				if (right - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
					if (refSize == RefSize::small) {
						MatchTuple<uint32_t> tempMatch{ (uint32_t)(p2 + 1 - gen2),(uint32_t)(p1 + 1 - gen1),(uint32_t)(right - p1 - 1) };
						postProcessPartial<uint32_t>(thID, matches12, tempMatch, matchesBuffer, seqName);
					}
					else {
						MatchTuple<uint64_t> tempMatch{ (uint32_t)(p2 + 1 - gen2),(uint64_t)(p1 + 1 - gen1),(uint32_t)(right - p1 - 1) };
						postProcessPartial<uint64_t>(thID, matches16, tempMatch, matchesBuffer, seqName);
					}
				}
				else {
					//nop
				}
			}
		}
		curr2 += k2;
	}
}


template<class MyUINT1, class MyUINT2>
//MyUINT1 - size of sampledPositions element: small:32b, big:32b, huge:64b
//MyUINT2 - size of cum element: small:32b, big:32b, huge:64b
void processQueryWorker(int thID, string ofn, size_t N1, char* start1, MyUINT2* cum, MyUINT1* sampledPositions) {
	string matchesBuffer;
	matchesBuffer.reserve(MATCHES_BUFFER_SIZE);
	vector <MatchTuple<uint32_t>> matches12;
	vector <MatchTuple<uint64_t>> matches16;
	matches12.push_back({ 0,0,0 }); // first element is a guard value protect from sorting errors
	matches16.push_back({ 0,0,0 });

	ifstream qf(Q_FN, ios::binary);  //each threads opens query file
	if (qf.fail()) {
		cerr << "Thread " << thID << ": Error opening input file " << Q_FN << "\n";
		exit(-1);
	}

	f_matches[thID].open(ofn, ios::binary); //matches file - output.
	if (f_matches[thID].fail()) {
		cerr << "Thread " << thID << ": Error opening output file " << ofn << "\n";
		exit(-1);
	}

	size_t beg, end;
	GenomeData rGenome = { N1,start1 };
	BlockVector& bv = blockVector[thID];

	for (auto bi : bv) {
		if (blockSortBreakFlag[thID])
			break; //leave outer loop

		mtx.lock();
		SequenceVector2 sv2 = readBlock(thID, qf, bi, Q_PADDING_CHAR, true);
		mtx.unlock();
		for (auto si3 : sv2) {
			size_t N2 = get<2>(si3);
			char* start2 = get<1>(si3);
			if (isRC != RC::yes) {
				beg = 0;
				end = N2;

				if (matchesBuffer.size() > (size_t)(MATCHES_BUFFER_SIZE * 0.95))
					flushMatchesBuffer(f_matches[thID], matchesBuffer);

				writeSeqHeader(matchesBuffer, get<0>(si3));
				processQueryTh<MyUINT1, MyUINT2>(thID, start1, start2, cum, sampledPositions, beg, end, matchesBuffer, matches12, matches16, get<0>(si3), true);
				if (refSize == RefSize::small)
					postProcess<uint32_t>(thID, matches12, matchesBuffer, get<0>(si3));
				else
					postProcess<uint64_t>(thID, matches16, matchesBuffer, get<0>(si3));
			}
			if (isRC != RC::no) {
				char complement[256];
				prepareComplement(complement, Q_PADDING_CHAR);
				reverseComplement(start2, complement, N2);
				beg = 0;
				end = N2;

				if (matchesBuffer.size() > (size_t)(MATCHES_BUFFER_SIZE * 0.95))
					flushMatchesBuffer(f_matches[thID], matchesBuffer);

				string seqNameRev = get<0>(si3) + " Reverse";
				writeSeqHeader(matchesBuffer, seqNameRev);
				processQueryTh<MyUINT1, MyUINT2>(thID, start1, start2, cum, sampledPositions, beg, end, matchesBuffer, matches12, matches16, get<0>(si3), true);
				if (refSize == RefSize::small)
					postProcess<uint32_t>(thID, matches12, matchesBuffer, get<0>(si3) + " Reverse");
				else
					postProcess<uint64_t>(thID, matches16, matchesBuffer, get<0>(si3) + " Reverse");
			}
			if (blockSortBreakFlag[thID])
				break; //leave inner loop
		}
	}
	flushMatchesBuffer(f_matches[thID], matchesBuffer);
	matches12.clear();
	matches16.clear();
	matchesCopy[thID].clear();
	f_matches[thID].close();
	qf.close();

	if (blockSortBreakFlag[thID]) { //sorting issue occured. Do once again with full sorting.
		cout << "Thread " << thID << ": block sort fail, proceed with full sort.\n";
		blockSortBreakFlag[thID] = false;
		MATCH_BLOCK_SIZE[thID] = MAX_BLOCK_SIZE;
		processQueryWorker(thID, ofn, N1, start1, cum, sampledPositions);
	}
}



template<class MyUINT1, class MyUINT2>
//MyUINT1 - size of sampledPositions element: small:32b, big:32b, huge:64b
//MyUINT2 - size of cum element: small:32b, big:32b, huge:64b
void processQueryMT(HashBuffer<MyUINT1, MyUINT2> buffer, GenomeData& rGenome) {
	//0. unpack parameters
	const streamsize READBUFSIZE = 1LL << 21;
	size_t N1 = get<0>(rGenome);
	char* start1 = get<1>(rGenome);
	MyUINT1* sampledPositions = buffer.first;
	MyUINT2* cum = buffer.second;

	//1. create thread array
	thread* th = new thread[nThreads + 1]; //thread number starts of 1

	for (int i = 0; i < nThreads; i++) {
		MATCH_BLOCK_SIZE[i] = DEF_MATCH_BLOCK_SIZE;
		blockSortBreakFlag[i] = false;
	}

	//2. thread loop
	for (int i = 1; i < nThreads; i++) {
		//temporary output files outputFN + "_th" + nTh
		string fn = matchesFN + "_th_" + to_string(i);
		th[i] = thread(processQueryWorker<MyUINT1, MyUINT2>, i, fn, N1, start1, cum, sampledPositions);
	}

	processQueryWorker<MyUINT1, MyUINT2>(0, matchesFN, N1, start1, cum, sampledPositions); //first part is done in main thread
	CStopWatch sw;
	sw.start();
	sw.stop();
	f_matches[0].open(matchesFN, ios::app | ios::binary);
	char* inbuf = new char[READBUFSIZE];
	for (int i = 1; i < nThreads; i++) {
		//3. join threads
		th[i].join();
		string fn = matchesFN + "_th_" + to_string(i);
		ifstream f;
		f.open(fn, ios::binary);
		if (f.fail()) {
			cerr << "Error opening temporary file " << fn << "\n";
			exit(-1);
		}
		//4. merge files and delete
		sw.resume();
		streamsize gcount = 0;
		while (true) {
			f.read(inbuf, READBUFSIZE);
			gcount = f.gcount();
			if (0 == gcount)
				break;
			f_matches[0].write(inbuf, gcount);
		}
		f.close();
		remove(fn.c_str());
		sw.stop();
	}
	sw.resume();
	f_matches[0].close();
	delete[] inbuf;
	delete[] th;
	sw.stop();
	cout << "Time of output merging " << sw.totalTime() << "\n";
}

#if predsv == 1
/**
Function used for fast search in sequence map
*/
size_t findMinDiffBitShift(SequenceVector& seq, size_t N) {
	size_t d = SIZE_MAX;
	size_t last_value = seq.back().second;
	for (size_t i = 1; i < seq.size(); ++i) {
		size_t currd = seq[i].second - seq[i - 1].second;
		if (currd < d)
			d = currd;
	}
	if (N - last_value < d)
		d = N - last_value;

	size_t bitshift = 0;
	size_t cur = 1;
	while (cur <= d) {
		cur *= 2;
		++bitshift;
	}
	min_diff = 1 << (bitshift - 1);
	return bitshift - 1;
}

/**
Function used for fast search in sequence map
*/
SequenceVectorR createPredSV(SequenceVector& seq, size_t N) {
	if (seq.size() >= UINT32_MAX) {
		cerr << "Too many sequences! Exit.\n";
		exit(1);
	}
	uint32_t S = (uint32_t)seq.size();

	SequenceVectorR predSV;
	uint32_t prev = 0;
	size_t cur = 0;
	for (uint32_t i = 0; i < seq.size(); ++i) {
		size_t num = seq[i].second;
		while (num > cur) {
			predSV.push_back(prev);
			cur += min_diff;
		}
		prev = i;
	}
	for (; cur < N; cur += min_diff)
		predSV.push_back(S - 1);

	return predSV;
}
#endif

GenomeData readMultiFasta(string fn, const char paddingChar, bool removeNs, genometype seqType) {
	CStopWatch stopWatch;
	stopWatch.start();
	const char beginChar = '>';
	const char terminatorChar = 0;
	const char eolChar1 = 10;
	const char eolChar2 = 13;
	const char spaceChar = ' ';
	const char filterChar = (char)0xDF;
	const int paddingSize = L + 1;

	if (seqType == genometype::query)
		*v1logger << "Reading query genome ...\n";
	else
		*v1logger << "Reading reference genome ...\n";

	//create a buffer for the whole file + padding at left and right
	ifstream f(fn, ios::ate | ios::binary);
	if (f.fail()) {
		cerr << "\nFile '" << fn << "' does not exist. Quit.";
		exit(1);
	}
	size_t N = f.tellg();
	f.seekg(0, ios::beg);
	char* buf1 = new char[N + 2 * paddingSize];
	if (buf1 == nullptr) {
		cerr << "\nFile '" << fn << "' is too large. Quit.";
		exit(1);
	}
	memset(buf1, paddingChar, paddingSize);
	f.read(buf1 + paddingSize, N);
	*v2logger << "\treadMultiFasta: Reading file from disk " << stopWatch.stop() << "\n";
	stopWatch.resume();

	buf1[paddingSize + N] = terminatorChar; // null-terminate the string
	memset(buf1 + paddingSize + N + 1, paddingChar, paddingSize - 1);
	memset(buf1 + paddingSize + N, terminatorChar, 10);  // >= sizeof(uint_64) is enough
	f.close();

	char* gen = buf1 + paddingSize;
	SequenceVector& seq = r_SV;
	if (seqType == genometype::query)
		seq = q_SV;

	char* dst = gen;
	char* src = gen;

	char tempLine[512];  // FASTA lines shouldn't be that long (even 128 should be ok)

	while (true) {
		if (*src == beginChar) {
			size_t idx = 0;
			while (*src != eolChar1 && *src != eolChar2 && *src != ' ' && *src != terminatorChar) {
				tempLine[idx] = *src++;
				++idx;
			}
			tempLine[idx] = 0;
			seq.push_back({ tempLine + 1, (dst - gen) });  // + 1, as we omit the starting '>'
														   //search for EOL
			while (*src != eolChar1 && *src != eolChar2)
				src++;

			*dst++ = paddingChar;
		}
		else {
			while (*src == eolChar1 || *src == eolChar2) {
				++src;
			}
			if (*src == beginChar)
				continue;
			if (*src == terminatorChar)
				break;

			uint64_t temp2;
			memcpy(&temp2, src, 8);
			while ((temp2 & 0x4040404040404040) == 0x4040404040404040) {
				temp2 = temp2 & 0xDFDFDFDFDFDFDFDF;
				memcpy(dst, &temp2, 8);
				dst += 8;
				src += 8;
				memcpy(&temp2, src, 8);
			}
			while (((*src) & (char)0x40) == (char)0x40) {
				*dst++ = (*src++) & (char)0xDF;
			}
		}
	}

	memset(dst, paddingChar, N + paddingSize - (dst - gen));

	*v2logger << "\treadMultiFasta: Before replaceBadSymbol " << stopWatch.stop() << "\n";
	stopWatch.resume();
	char* beg = gen, * end = gen;

	size_t nThreads2 = (nThreads == 1 ? 1 : 2);  // protects from memory congestion (?)
	const size_t part_size = (dst - gen) / nThreads2;
	thread* th = new thread[nThreads2];
	for (size_t t = 0; t < nThreads2 - 1; t++) {
		beg = gen + t * part_size;
		end = gen + (t + 1ULL) * part_size;
		th[t] = thread(replaceBadSymbol, beg, end, 'N', paddingChar);
	}
	replaceBadSymbol(end, dst, 'N', paddingChar);
	for (size_t t = 0; t < nThreads2 - 1; t++)
		th[t].join();

	for (size_t t = 0; t < nThreads2 - 1; t++) {
		beg = gen + t * part_size;
		end = gen + (t + 1ULL) * part_size;
		th[t] = thread(replaceBadSymbol, beg, end, 'n', paddingChar);
	}
	replaceBadSymbol(end, dst, 'n', paddingChar);
	for (size_t t = 0; t < nThreads2 - 1; t++)
		th[t].join();
	delete[] th;

	for (auto s : seq) {
		*v2logger << s.first << "\t" << s.second << endl;
	}

	*v2logger << "\treadMultiFasta: After replaceBadSymbol " << stopWatch.stop() << "\n";
	stopWatch.resume();

	*v2logger << "\treadMultiFasta: Analysing file " << stopWatch.stop() << "\n";

	*v2logger << "\treadMultiFasta: Genome data size " << (dst - gen) << endl;
	*v2logger << "\treadMultiFasta: Sequences " << seq.size() << endl;

	return { (dst - gen), gen };
}


SequenceVector2 readBlock(int thID, ifstream& f, BlockItem& bi, const char paddingChar, bool removeNs) {
	const char eolChar1 = 10;
	const char eolChar2 = 13;
	const char beginChar = '>';
	const int paddingSize = L + 1;
	size_t N = get<1>(bi);
	vector<size_t>& sv2 = get<2>(bi);

	char* gen = nullptr;
	char* genTemp = blockBuffer[thID] + paddingSize;
	//reading buffer to memory
	f.seekg(get<0>(bi), ios::beg);
	f.read(genTemp, get<1>(bi));
	memset(genTemp + get<1>(bi), 0, paddingSize);

	gen = genTemp;
	char* dst = gen;
	char* src = gen;
	SequenceVector2 sv3;

	char tempLine[512];  // FASTA lines shouldn't be that long (even 128 should be ok)

	//processing sequences: removing header, new lines
	for (size_t i = 0; i < sv2.size(); i++) {
		size_t startHeader = sv2[i];
		size_t nextSeqence = (sv2.size() > (i + 1) ? sv2[i + 1] : get<1>(bi));

		string seqName;
		size_t seqStart;
		size_t seqSize;
		//check, if the sequence starts from '>'
		src = gen + startHeader;
		char* lastChar = gen + (nextSeqence);
		if (*src != beginChar) {
			cerr << "Invalid fasta file.\n";
			exit(0);
		}
		size_t idx = 0;
		while (*src != eolChar1 && *src != eolChar2 && *src != ' ' && src != lastChar) { //reading header line
			tempLine[idx] = *src++;
			++idx;
		}
		tempLine[idx] = 0;
		seqName = tempLine + 1;
		while (*src != eolChar1 && *src != eolChar2) //search for EOL after header
			src++;
		seqStart = dst - gen;
		char* startPtr = dst;
		*dst++ = paddingChar;  //padding char at the begin to make indexing matches from 1

		while (true) { //scanning the sequence
			while (*src == eolChar1 || *src == eolChar2) {
				if (src == lastChar)
					break;
				++src;
			}
			if (src == lastChar)
				break;

			uint64_t temp2; //scan 8 bytes at once
			memcpy(&temp2, src, 8);
			while ((temp2 & 0x4040404040404040) == 0x4040404040404040) {
				temp2 = temp2 & 0xDFDFDFDFDFDFDFDF;
				memcpy(dst, &temp2, 8);
				dst += 8;
				src += 8;
				memcpy(&temp2, src, 8);
			}
			while (((*src) & (char)0x40) == (char)0x40) {
				*dst++ = (*src++) & (char)0xDF;
			}
		}
		seqSize = dst - startPtr;
		sv3.push_back({ seqName, startPtr, seqSize });
	}
	*dst = paddingChar;
	if (removeNs == true) {
		replaceBadSymbol(gen, dst, 'N', paddingChar);
		replaceBadSymbol(gen, dst, 'n', paddingChar);
	}

	return sv3;
}


void createBlockBuffer(vector<size_t>& seqStarts) {
	const size_t MIN_BLOCK_LEN = 1ULL << 25; //block length not less than 32 MB
	//0. Check for the number of threads
	if (nThreads > seqStarts.size() - 1) {
		nThreads = (int)seqStarts.size() - 1;
		*v1logger << "Warning: the number of threads (nThreads) is limited to the number of sequences in Q, i.e. " << nThreads << "\n";
	}

	// 1. Given vector of sequences, obtain best ranges of sequences
	size_t genSize = seqStarts.back() - seqStarts.front();

	size_t* bestDiv = new size_t[nThreads + 1];
	bestDiv[0] = 0;
	for (int t = 0; t < nThreads; t++) {
		bestDiv[t] = t * genSize / nThreads;
	}
	bestDiv[nThreads] = genSize;

	//2. Assign slices of the sequence vector to thread sequence vectors
	vector<size_t>* thSeqStarts = new vector<size_t>[nThreads];
	int th = 0;
	thSeqStarts[th].push_back(seqStarts[0]);
	for (int i = 1; i < seqStarts.size(); ++i) {
		if (seqStarts[i] > bestDiv[th + 1]) {
			thSeqStarts[th].push_back(seqStarts[i]);
			++th;
		}
		thSeqStarts[th].push_back(seqStarts[i]);
	}
	for (th = 0; th < nThreads; th++) {
		// 3. For each vector, obtain maximum sequence size
		size_t maxSeqLen = MIN_BLOCK_LEN;
		for (size_t i = 0; i < thSeqStarts[th].size() - 1; ++i) {
			size_t le = thSeqStarts[th][i + 1] - thSeqStarts[th][i];
			if (maxSeqLen < le)
				maxSeqLen = le;
		}
		// obtain blocks
		// each block consists of:
		// starting ptr in file
		// length of the block
		// list of sequences
		// pointers in sequences must be shifted to beginning of the block
		BlockItem bi;
		get<0>(bi) = thSeqStarts[th][0];
		for (size_t i = 0; i < thSeqStarts[th].size() - 1; ++i) {
			size_t le = thSeqStarts[th][i + 1] - thSeqStarts[th][i];
			if (get<1>(bi) + le > maxSeqLen) {  //If not fit -- push block to the list and start a new block
				blockVector[th].push_back(bi);
				get<0>(bi) = thSeqStarts[th][i];
				get<1>(bi) = le;
				get<2>(bi).clear();
			}
			else
				get<1>(bi) += le;
			get<2>(bi).push_back(thSeqStarts[th][i] - get<0>(bi));
		}
		blockVector[th].push_back(bi);
		*v2logger << "Thread " << th << endl;
		for (auto b : blockVector[th]) {
			*v2logger << get<0>(b) << "\tBlock " << get<1>(b) << "\n";
			for (auto li : get<2>(b)) {
				*v2logger << "\t\tSequence " << li << "\n";
			}
		}

		thSeqStarts[th].clear();
		//Assign memory for block buffers
		size_t blockBufferSize = maxSeqLen + 2 * (L + 1U);
		char* tmpBlock = new char[blockBufferSize];
		memset(tmpBlock, 0, blockBufferSize);
		blockBuffer.push_back(tmpBlock);
	}
	seqStarts.clear();
	delete[]thSeqStarts;
	delete[]bestDiv;
}


void deleteBlockBuffer() {
	blockVector[0].clear();
	for (auto bb : blockBuffer)
		delete[] bb;
}


void scanMultiFasta(string fn) {
	*v1logger << "Scanning Query genome ...\n";
	CStopWatch stopWatch;
	stopWatch.start();
	const char beginChar = '>';
	const size_t scanBufSize = 1ULL << 22;//26
	char* buf1 = new char[scanBufSize];
	if (buf1 == nullptr) {
		cerr << "\nFile '" << fn << "' is too large. Quit.";
		exit(1);
	}

	ifstream f(fn, ios::binary);
	if (f.fail()) {
		cerr << "\nFile '" << fn << "' does not exist. Quit.";
		exit(1);
	}
	size_t bytesRead = 0;
	size_t currPos = 0;
	vector<size_t> seqStarts;
	while (f) {
		f.read(buf1, scanBufSize);
		if (f)
			bytesRead = scanBufSize;
		else
			bytesRead = f.gcount();
		char* foundPos = buf1;
		while (true) {
			foundPos = find(foundPos, buf1 + bytesRead, beginChar);
			if (foundPos != (buf1 + bytesRead)) {
				seqStarts.push_back(currPos + (foundPos - buf1));
				foundPos++;
			}
			else
				break;
		}
		currPos += bytesRead;
	}
	seqStarts.push_back(currPos); //pointer to the end of file
	f.close();
	delete[] buf1;
	createBlockBuffer(seqStarts);
	seqStarts.clear();

	*v2logger << "\tscanMultiFasta: Analysing file " << stopWatch.stop() << "\n";
	*v2logger << "\tscanMultiFasta: Sequence blocks " << blockVector[0].size() << endl;
}


template <class MyUINT1, class MyUINT2>
void deleteHashBuffer(HashBuffer<MyUINT1, MyUINT2>& buf) {
	delete[] buf.first;
	delete[] buf.second;
}


void deleteReading(GenomeData& r) {
	delete[](get<1>(r) - (L + 1));
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		displayHelp(argv[0]);
		return 1;
	}
	initHashFuncMatrix();
	initDefaults();
	processCmd(argc, argv);
	initGlobals();

	if (isVerbose > verbosity::v0) {
		displayParams();
		cout.setf(ios_base::unitbuf);
	}
	CStopWatch stopwatch;
	stopwatch.start();
	size_t maxSeqSize = 0;
	size_t N1 = 0;
	size_t N2 = 0;
	GenomeData  rGenome, qGenome;
	scanMultiFasta(Q_FN);
	rGenome = readMultiFasta(R_FN, REF_PADDING_CHAR, true, genometype::reference);
	*v1logger << "Time of I/O = " << stopwatch.stop() << endl;
	stopwatch.resume();
#if predsv == 1
	min_diff_bitshift = findMinDiffBitShift(r_SV, get<0>(rGenome));
	r_SV_pred = createPredSV(r_SV, get<0>(rGenome));
#endif
	refSize = RefSize::small;  // small Reference
	if ((get<0>(rGenome)) >= (1ULL << 32)) {
		refSize = RefSize::big;  // large Reference
		if ((get<0>(rGenome)) / k1 >= (1ULL << 32))
			refSize = RefSize::huge;  // huge Reference
	}
	if (forceBR > refSize)
		refSize = forceBR;

	if (refSize == RefSize::huge) {
		*v1logger << "WARNING - LARGE reference file (SIZE / k1 > 4GB), all 64-bit arrays.\n";
		pair<uint64_t*, uint64_t*> buffer = processRef<uint64_t, uint64_t>(rGenome);
		*v1logger << "Time of processRef = " << stopwatch.stop() << endl;
		processQueryMT<uint64_t, uint64_t>(buffer, rGenome);
		*v1logger << "Time of processQuery = " << stopwatch.stop() << endl;
		stopwatch.resume();
		deleteBlockBuffer();
		deleteHashBuffer(buffer);
	}
	else
		if (refSize == RefSize::big) {
			*v1logger << "WARNING - BIG reference file (>4GB), some 64-bit arrays.\n";
			pair<uint32_t*, uint32_t*> buffer = processRef<uint32_t, uint32_t>(rGenome);
			*v1logger << "Time of processRef = " << stopwatch.stop() << endl;
			stopwatch.resume();
			processQueryMT<uint32_t, uint32_t>(buffer, rGenome);
			*v1logger << "Time of processQuery = " << stopwatch.stop() << endl;
			stopwatch.resume();
			deleteBlockBuffer();
			deleteHashBuffer(buffer);
		}
		else { //refSize == small
			pair<uint32_t*, uint32_t*> buffer = processRef<uint32_t, uint32_t>(rGenome);
			*v1logger << "Time of processRef = " << stopwatch.stop() << endl;
			stopwatch.resume();
			processQueryMT<uint32_t, uint32_t>(buffer, rGenome);
			*v1logger << "Time of processQuery = " << stopwatch.stop() << endl;
			stopwatch.resume();
			deleteBlockBuffer();
			deleteHashBuffer(buffer);
		}
	deleteReading(rGenome);
	*v1logger << "Time of deleting = " << stopwatch.stop() << "\nTotal time = " << stopwatch.totalTime() << "\nFINISHED\n\n";
	return 0;
}
