/* ================================================================ *
*  Hashes.h : header file for hash definition                       *
*                                                                   *
*  copMEM2 is a program for efficient computation of MEMs           *
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

#pragma once
#ifndef _HASHES_H
#define _HASHES_H
#include <cstdint>
#define XXH_INLINE_ALL
#define XXH_PRIVATE_API
#include "xxhash.h"
#include "metrohash64.h"
#include "city.h"

template <size_t K, std::uint32_t HASH_SIZE_MINUS_ONE>
inline std::uint32_t xxhash32(const char* key) {
	XXH32_hash_t result = XXH32((const void*)key, K, 1286608618UL);
	return (std::uint32_t)(result & HASH_SIZE_MINUS_ONE);
}

template <size_t K, std::uint32_t HASH_SIZE_MINUS_ONE>
inline std::uint32_t xxhash64(const char* key) {
	XXH64_hash_t result = XXH64((const void*)key, K, 9876543210ULL);
	return (std::uint32_t)(result & HASH_SIZE_MINUS_ONE);
}

unsigned long long result = 0ULL;
template <size_t K, std::uint32_t HASH_SIZE_MINUS_ONE>
inline std::uint32_t metroHash64(const char* key) {
	MetroHash64::Hash((const uint8_t*)key, K, (uint8_t*)(&result), 9876543210ULL);
	return (std::uint32_t)(result & HASH_SIZE_MINUS_ONE);
}

template <size_t K, std::uint32_t HASH_SIZE_MINUS_ONE>
inline std::uint32_t cityHash64(const char* key) {
	return (std::uint32_t)(CityHash64(key, K) & HASH_SIZE_MINUS_ONE);
}

// based on http://www.amsoftware.narod.ru/algo2.html
/*
template <size_t K>
inline std::uint32_t maRushPrime1HashSimplified(const char *str) {
	std::uint64_t hash = K;
	for (std::uint32_t j = 0; j < K/4; ) {
		std::uint32_t k;
		memcpy(&k, str, 4);
		k += j++;
		hash ^= k;
		hash *= 171717;
		str += 4;
	}
	return (std::uint32_t)(hash & HASH_SIZE_MINUS_ONE);
}
*/
template <size_t K, std::uint32_t HASH_SIZE_MINUS_ONE>
inline std::uint32_t maRushPrime1HashSimplified(const char* str) {
	std::uint32_t hash = (uint32_t)K;
	for (std::uint32_t j = 0; j < K / 4; j++) {
		std::uint32_t k;
		memcpy(&k, str, 4);
		hash += k;
		hash *= 11;
		//hash *= 171717;
		str += 4;
	}
	return hash & HASH_SIZE_MINUS_ONE;
}





inline std::uint32_t maRushPrime1HashSimplified8(const char *str) {
	std::uint64_t hash = 8;   // 0x3B5E;  
	for (std::uint32_t j = 0; j < 2; ) {
		std::uint32_t k;
		memcpy(&k, str, 4);
		k += j++;
		hash ^= k;
		hash *= 171717;
		str += 4;
	}
	//return (std::uint32_t)(hash & HASH_SIZE_MINUS_ONE);
	return (std::uint32_t)(hash);
}
inline std::uint64_t maRushPrime1HashSimplified12(const std::uint64_t *tab) {
	std::uint64_t hash = 12;   // 0x3B5E;  
	for (std::uint32_t j = 0; j < 2;++j ) {
		hash ^= tab[j];
		hash *= 171717;
	}
	return hash;
}

std::uint32_t(*hashFunc)(const char*);
std::uint32_t(*hashFuncMatrix[97][6][3])(const char*);

void initHashFuncMatrix() {
	hashFuncMatrix[32][1][0] = maRushPrime1HashSimplified<32, (1 << 28) - 1>;
	hashFuncMatrix[32][2][0] = xxhash32<32, (1 << 28) - 1>;
	hashFuncMatrix[32][3][0] = xxhash64<32, (1 << 28) - 1>;
	hashFuncMatrix[32][4][0] = metroHash64<32, (1 << 28) - 1>;
	hashFuncMatrix[32][5][0] = cityHash64<32, (1 << 28) - 1>;
	hashFuncMatrix[36][1][0] = maRushPrime1HashSimplified<36, (1 << 28) - 1>;
	hashFuncMatrix[36][2][0] = xxhash32<36, (1 << 28) - 1>;
	hashFuncMatrix[36][3][0] = xxhash64<36, (1 << 28) - 1>;
	hashFuncMatrix[36][4][0] = metroHash64<36, (1 << 28) - 1>;
	hashFuncMatrix[36][5][0] = cityHash64<36, (1 << 28) - 1>;
	hashFuncMatrix[40][1][0] = maRushPrime1HashSimplified<40, (1 << 28) - 1>;
	hashFuncMatrix[40][2][0] = xxhash32<40, (1 << 28) - 1>;
	hashFuncMatrix[40][3][0] = xxhash64<40, (1 << 28) - 1>;
	hashFuncMatrix[40][4][0] = metroHash64<40, (1 << 28) - 1>;
	hashFuncMatrix[40][5][0] = cityHash64<40, (1 << 28) - 1>;
	hashFuncMatrix[44][1][0] = maRushPrime1HashSimplified<44, (1 << 28) - 1>;
	hashFuncMatrix[44][2][0] = xxhash32<44, (1 << 28) - 1>;
	hashFuncMatrix[44][3][0] = xxhash64<44, (1 << 28) - 1>;
	hashFuncMatrix[44][4][0] = metroHash64<44, (1 << 28) - 1>;
	hashFuncMatrix[44][5][0] = cityHash64<44, (1 << 28) - 1>;
	hashFuncMatrix[48][1][0] = maRushPrime1HashSimplified<48, (1 << 28) - 1>;
	hashFuncMatrix[48][2][0] = xxhash32<48, (1 << 28) - 1>;
	hashFuncMatrix[48][3][0] = xxhash64<48, (1 << 28) - 1>;
	hashFuncMatrix[48][4][0] = metroHash64<48, (1 << 28) - 1>;
	hashFuncMatrix[48][5][0] = cityHash64<48, (1 << 28) - 1>;
	hashFuncMatrix[52][1][0] = maRushPrime1HashSimplified<52, (1 << 28) - 1>;
	hashFuncMatrix[52][2][0] = xxhash32<52, (1 << 28) - 1>;
	hashFuncMatrix[52][3][0] = xxhash64<52, (1 << 28) - 1>;
	hashFuncMatrix[52][4][0] = metroHash64<52, (1 << 28) - 1>;
	hashFuncMatrix[52][5][0] = cityHash64<52, (1 << 28) - 1>;
	hashFuncMatrix[56][1][0] = maRushPrime1HashSimplified<56, (1 << 28) - 1>;
	hashFuncMatrix[56][2][0] = xxhash32<56, (1 << 28) - 1>;
	hashFuncMatrix[56][3][0] = xxhash64<56, (1 << 28) - 1>;
	hashFuncMatrix[56][4][0] = metroHash64<56, (1 << 28) - 1>;
	hashFuncMatrix[56][5][0] = cityHash64<56, (1 << 28) - 1>;
	hashFuncMatrix[60][1][0] = maRushPrime1HashSimplified<60, (1 << 28) - 1>;
	hashFuncMatrix[60][2][0] = xxhash32<60, (1 << 28) - 1>;
	hashFuncMatrix[60][3][0] = xxhash64<60, (1 << 28) - 1>;
	hashFuncMatrix[60][4][0] = metroHash64<60, (1 << 28) - 1>;
	hashFuncMatrix[60][5][0] = cityHash64<60, (1 << 28) - 1>;
	hashFuncMatrix[64][1][0] = maRushPrime1HashSimplified<64, (1 << 28) - 1>;
	hashFuncMatrix[64][2][0] = xxhash32<64, (1 << 28) - 1>;
	hashFuncMatrix[64][3][0] = xxhash64<64, (1 << 28) - 1>;
	hashFuncMatrix[64][4][0] = metroHash64<64, (1 << 28) - 1>;
	hashFuncMatrix[64][5][0] = cityHash64<64, (1 << 28) - 1>;
	hashFuncMatrix[68][1][0] = maRushPrime1HashSimplified<68, (1 << 28) - 1>;
	hashFuncMatrix[68][2][0] = xxhash32<68, (1 << 28) - 1>;
	hashFuncMatrix[68][3][0] = xxhash64<68, (1 << 28) - 1>;
	hashFuncMatrix[68][4][0] = metroHash64<68, (1 << 28) - 1>;
	hashFuncMatrix[68][5][0] = cityHash64<68, (1 << 28) - 1>;
	hashFuncMatrix[72][1][0] = maRushPrime1HashSimplified<72, (1 << 28) - 1>;
	hashFuncMatrix[72][2][0] = xxhash32<72, (1 << 28) - 1>;
	hashFuncMatrix[72][3][0] = xxhash64<72, (1 << 28) - 1>;
	hashFuncMatrix[72][4][0] = metroHash64<72, (1 << 28) - 1>;
	hashFuncMatrix[72][5][0] = cityHash64<72, (1 << 28) - 1>;
	hashFuncMatrix[76][1][0] = maRushPrime1HashSimplified<76, (1 << 28) - 1>;
	hashFuncMatrix[76][2][0] = xxhash32<76, (1 << 28) - 1>;
	hashFuncMatrix[76][3][0] = xxhash64<76, (1 << 28) - 1>;
	hashFuncMatrix[76][4][0] = metroHash64<76, (1 << 28) - 1>;
	hashFuncMatrix[76][5][0] = cityHash64<76, (1 << 28) - 1>;
	hashFuncMatrix[80][1][0] = maRushPrime1HashSimplified<80, (1 << 28) - 1>;
	hashFuncMatrix[80][2][0] = xxhash32<80, (1 << 28) - 1>;
	hashFuncMatrix[80][3][0] = xxhash64<80, (1 << 28) - 1>;
	hashFuncMatrix[80][4][0] = metroHash64<80, (1 << 28) - 1>;
	hashFuncMatrix[80][5][0] = cityHash64<80, (1 << 28) - 1>;
	hashFuncMatrix[84][1][0] = maRushPrime1HashSimplified<84, (1 << 28) - 1>;
	hashFuncMatrix[84][2][0] = xxhash32<84, (1 << 28) - 1>;
	hashFuncMatrix[84][3][0] = xxhash64<84, (1 << 28) - 1>;
	hashFuncMatrix[84][4][0] = metroHash64<84, (1 << 28) - 1>;
	hashFuncMatrix[84][5][0] = cityHash64<84, (1 << 28) - 1>;
	hashFuncMatrix[88][1][0] = maRushPrime1HashSimplified<88, (1 << 28) - 1>;
	hashFuncMatrix[88][2][0] = xxhash32<88, (1 << 28) - 1>;
	hashFuncMatrix[88][3][0] = xxhash64<88, (1 << 28) - 1>;
	hashFuncMatrix[88][4][0] = metroHash64<88, (1 << 28) - 1>;
	hashFuncMatrix[88][5][0] = cityHash64<88, (1 << 28) - 1>;
	hashFuncMatrix[92][1][0] = maRushPrime1HashSimplified<92, (1 << 28) - 1>;
	hashFuncMatrix[92][2][0] = xxhash32<92, (1 << 28) - 1>;
	hashFuncMatrix[92][3][0] = xxhash64<92, (1 << 28) - 1>;
	hashFuncMatrix[92][4][0] = metroHash64<92, (1 << 28) - 1>;
	hashFuncMatrix[92][5][0] = cityHash64<92, (1 << 28) - 1>;
	hashFuncMatrix[96][1][0] = maRushPrime1HashSimplified<96, (1 << 28) - 1>;
	hashFuncMatrix[96][2][0] = xxhash32<96, (1 << 28) - 1>;
	hashFuncMatrix[96][3][0] = xxhash64<96, (1 << 28) - 1>;
	hashFuncMatrix[96][4][0] = metroHash64<96, (1 << 28) - 1>;
	hashFuncMatrix[96][5][0] = cityHash64<96, (1 << 28) - 1>;

	hashFuncMatrix[32][1][1] = maRushPrime1HashSimplified<32, (1 << 29) - 1>;
	hashFuncMatrix[32][2][1] = xxhash32<32, (1 << 29) - 1>;
	hashFuncMatrix[32][3][1] = xxhash64<32, (1 << 29) - 1>;
	hashFuncMatrix[32][4][1] = metroHash64<32, (1 << 29) - 1>;
	hashFuncMatrix[32][5][1] = cityHash64<32, (1 << 29) - 1>;
	hashFuncMatrix[36][1][1] = maRushPrime1HashSimplified<36, (1 << 29) - 1>;
	hashFuncMatrix[36][2][1] = xxhash32<36, (1 << 29) - 1>;
	hashFuncMatrix[36][3][1] = xxhash64<36, (1 << 29) - 1>;
	hashFuncMatrix[36][4][1] = metroHash64<36, (1 << 29) - 1>;
	hashFuncMatrix[36][5][1] = cityHash64<36, (1 << 29) - 1>;
	hashFuncMatrix[40][1][1] = maRushPrime1HashSimplified<40, (1 << 29) - 1>;
	hashFuncMatrix[40][2][1] = xxhash32<40, (1 << 29) - 1>;
	hashFuncMatrix[40][3][1] = xxhash64<40, (1 << 29) - 1>;
	hashFuncMatrix[40][4][1] = metroHash64<40, (1 << 29) - 1>;
	hashFuncMatrix[40][5][1] = cityHash64<40, (1 << 29) - 1>;
	hashFuncMatrix[44][1][1] = maRushPrime1HashSimplified<44, (1 << 29) - 1>;
	hashFuncMatrix[44][2][1] = xxhash32<44, (1 << 29) - 1>;
	hashFuncMatrix[44][3][1] = xxhash64<44, (1 << 29) - 1>;
	hashFuncMatrix[44][4][1] = metroHash64<44, (1 << 29) - 1>;
	hashFuncMatrix[44][5][1] = cityHash64<44, (1 << 29) - 1>;
	hashFuncMatrix[48][1][1] = maRushPrime1HashSimplified<48, (1 << 29) - 1>;
	hashFuncMatrix[48][2][1] = xxhash32<48, (1 << 29) - 1>;
	hashFuncMatrix[48][3][1] = xxhash64<48, (1 << 29) - 1>;
	hashFuncMatrix[48][4][1] = metroHash64<48, (1 << 29) - 1>;
	hashFuncMatrix[48][5][1] = cityHash64<48, (1 << 29) - 1>;
	hashFuncMatrix[52][1][1] = maRushPrime1HashSimplified<52, (1 << 29) - 1>;
	hashFuncMatrix[52][2][1] = xxhash32<52, (1 << 29) - 1>;
	hashFuncMatrix[52][3][1] = xxhash64<52, (1 << 29) - 1>;
	hashFuncMatrix[52][4][1] = metroHash64<52, (1 << 29) - 1>;
	hashFuncMatrix[52][5][1] = cityHash64<52, (1 << 29) - 1>;
	hashFuncMatrix[56][1][1] = maRushPrime1HashSimplified<56, (1 << 29) - 1>;
	hashFuncMatrix[56][2][1] = xxhash32<56, (1 << 29) - 1>;
	hashFuncMatrix[56][3][1] = xxhash64<56, (1 << 29) - 1>;
	hashFuncMatrix[56][4][1] = metroHash64<56, (1 << 29) - 1>;
	hashFuncMatrix[56][5][1] = cityHash64<56, (1 << 29) - 1>;
	hashFuncMatrix[60][1][1] = maRushPrime1HashSimplified<60, (1 << 29) - 1>;
	hashFuncMatrix[60][2][1] = xxhash32<60, (1 << 29) - 1>;
	hashFuncMatrix[60][3][1] = xxhash64<60, (1 << 29) - 1>;
	hashFuncMatrix[60][4][1] = metroHash64<60, (1 << 29) - 1>;
	hashFuncMatrix[60][5][1] = cityHash64<60, (1 << 29) - 1>;
	hashFuncMatrix[64][1][1] = maRushPrime1HashSimplified<64, (1 << 29) - 1>;
	hashFuncMatrix[64][2][1] = xxhash32<64, (1 << 29) - 1>;
	hashFuncMatrix[64][3][1] = xxhash64<64, (1 << 29) - 1>;
	hashFuncMatrix[64][4][1] = metroHash64<64, (1 << 29) - 1>;
	hashFuncMatrix[64][5][1] = cityHash64<64, (1 << 29) - 1>;
	hashFuncMatrix[68][1][1] = maRushPrime1HashSimplified<68, (1 << 29) - 1>;
	hashFuncMatrix[68][2][1] = xxhash32<68, (1 << 29) - 1>;
	hashFuncMatrix[68][3][1] = xxhash64<68, (1 << 29) - 1>;
	hashFuncMatrix[68][4][1] = metroHash64<68, (1 << 29) - 1>;
	hashFuncMatrix[68][5][1] = cityHash64<68, (1 << 29) - 1>;
	hashFuncMatrix[72][1][1] = maRushPrime1HashSimplified<72, (1 << 29) - 1>;
	hashFuncMatrix[72][2][1] = xxhash32<72, (1 << 29) - 1>;
	hashFuncMatrix[72][3][1] = xxhash64<72, (1 << 29) - 1>;
	hashFuncMatrix[72][4][1] = metroHash64<72, (1 << 29) - 1>;
	hashFuncMatrix[72][5][1] = cityHash64<72, (1 << 29) - 1>;
	hashFuncMatrix[76][1][1] = maRushPrime1HashSimplified<76, (1 << 29) - 1>;
	hashFuncMatrix[76][2][1] = xxhash32<76, (1 << 29) - 1>;
	hashFuncMatrix[76][3][1] = xxhash64<76, (1 << 29) - 1>;
	hashFuncMatrix[76][4][1] = metroHash64<76, (1 << 29) - 1>;
	hashFuncMatrix[76][5][1] = cityHash64<76, (1 << 29) - 1>;
	hashFuncMatrix[80][1][1] = maRushPrime1HashSimplified<80, (1 << 29) - 1>;
	hashFuncMatrix[80][2][1] = xxhash32<80, (1 << 29) - 1>;
	hashFuncMatrix[80][3][1] = xxhash64<80, (1 << 29) - 1>;
	hashFuncMatrix[80][4][1] = metroHash64<80, (1 << 29) - 1>;
	hashFuncMatrix[80][5][1] = cityHash64<80, (1 << 29) - 1>;
	hashFuncMatrix[84][1][1] = maRushPrime1HashSimplified<84, (1 << 29) - 1>;
	hashFuncMatrix[84][2][1] = xxhash32<84, (1 << 29) - 1>;
	hashFuncMatrix[84][3][1] = xxhash64<84, (1 << 29) - 1>;
	hashFuncMatrix[84][4][1] = metroHash64<84, (1 << 29) - 1>;
	hashFuncMatrix[84][5][1] = cityHash64<84, (1 << 29) - 1>;
	hashFuncMatrix[88][1][1] = maRushPrime1HashSimplified<88, (1 << 29) - 1>;
	hashFuncMatrix[88][2][1] = xxhash32<88, (1 << 29) - 1>;
	hashFuncMatrix[88][3][1] = xxhash64<88, (1 << 29) - 1>;
	hashFuncMatrix[88][4][1] = metroHash64<88, (1 << 29) - 1>;
	hashFuncMatrix[88][5][1] = cityHash64<88, (1 << 29) - 1>;
	hashFuncMatrix[92][1][1] = maRushPrime1HashSimplified<92, (1 << 29) - 1>;
	hashFuncMatrix[92][2][1] = xxhash32<92, (1 << 29) - 1>;
	hashFuncMatrix[92][3][1] = xxhash64<92, (1 << 29) - 1>;
	hashFuncMatrix[92][4][1] = metroHash64<92, (1 << 29) - 1>;
	hashFuncMatrix[92][5][1] = cityHash64<92, (1 << 29) - 1>;
	hashFuncMatrix[96][1][1] = maRushPrime1HashSimplified<96, (1 << 29) - 1>;
	hashFuncMatrix[96][2][1] = xxhash32<96, (1 << 29) - 1>;
	hashFuncMatrix[96][3][1] = xxhash64<96, (1 << 29) - 1>;
	hashFuncMatrix[96][4][1] = metroHash64<96, (1 << 29) - 1>;
	hashFuncMatrix[96][5][1] = cityHash64<96, (1 << 29) - 1>;

	hashFuncMatrix[32][1][2] = maRushPrime1HashSimplified<32, (1 << 30) - 1>;
	hashFuncMatrix[32][2][2] = xxhash32<32, (1 << 30) - 1>;
	hashFuncMatrix[32][3][2] = xxhash64<32, (1 << 30) - 1>;
	hashFuncMatrix[32][4][2] = metroHash64<32, (1 << 30) - 1>;
	hashFuncMatrix[32][5][2] = cityHash64<32, (1 << 30) - 1>;
	hashFuncMatrix[36][1][2] = maRushPrime1HashSimplified<36, (1 << 30) - 1>;
	hashFuncMatrix[36][2][2] = xxhash32<36, (1 << 30) - 1>;
	hashFuncMatrix[36][3][2] = xxhash64<36, (1 << 30) - 1>;
	hashFuncMatrix[36][4][2] = metroHash64<36, (1 << 30) - 1>;
	hashFuncMatrix[36][5][2] = cityHash64<36, (1 << 30) - 1>;
	hashFuncMatrix[40][1][2] = maRushPrime1HashSimplified<40, (1 << 30) - 1>;
	hashFuncMatrix[40][2][2] = xxhash32<40, (1 << 30) - 1>;
	hashFuncMatrix[40][3][2] = xxhash64<40, (1 << 30) - 1>;
	hashFuncMatrix[40][4][2] = metroHash64<40, (1 << 30) - 1>;
	hashFuncMatrix[40][5][2] = cityHash64<40, (1 << 30) - 1>;
	hashFuncMatrix[44][1][2] = maRushPrime1HashSimplified<44, (1 << 30) - 1>;
	hashFuncMatrix[44][2][2] = xxhash32<44, (1 << 30) - 1>;
	hashFuncMatrix[44][3][2] = xxhash64<44, (1 << 30) - 1>;
	hashFuncMatrix[44][4][2] = metroHash64<44, (1 << 30) - 1>;
	hashFuncMatrix[44][5][2] = cityHash64<44, (1 << 30) - 1>;
	hashFuncMatrix[48][1][2] = maRushPrime1HashSimplified<48, (1 << 30) - 1>;
	hashFuncMatrix[48][2][2] = xxhash32<48, (1 << 30) - 1>;
	hashFuncMatrix[48][3][2] = xxhash64<48, (1 << 30) - 1>;
	hashFuncMatrix[48][4][2] = metroHash64<48, (1 << 30) - 1>;
	hashFuncMatrix[48][5][2] = cityHash64<48, (1 << 30) - 1>;
	hashFuncMatrix[52][1][2] = maRushPrime1HashSimplified<52, (1 << 30) - 1>;
	hashFuncMatrix[52][2][2] = xxhash32<52, (1 << 30) - 1>;
	hashFuncMatrix[52][3][2] = xxhash64<52, (1 << 30) - 1>;
	hashFuncMatrix[52][4][2] = metroHash64<52, (1 << 30) - 1>;
	hashFuncMatrix[52][5][2] = cityHash64<52, (1 << 30) - 1>;
	hashFuncMatrix[56][1][2] = maRushPrime1HashSimplified<56, (1 << 30) - 1>;
	hashFuncMatrix[56][2][2] = xxhash32<56, (1 << 30) - 1>;
	hashFuncMatrix[56][3][2] = xxhash64<56, (1 << 30) - 1>;
	hashFuncMatrix[56][4][2] = metroHash64<56, (1 << 30) - 1>;
	hashFuncMatrix[56][5][2] = cityHash64<56, (1 << 30) - 1>;
	hashFuncMatrix[60][1][2] = maRushPrime1HashSimplified<60, (1 << 30) - 1>;
	hashFuncMatrix[60][2][2] = xxhash32<60, (1 << 30) - 1>;
	hashFuncMatrix[60][3][2] = xxhash64<60, (1 << 30) - 1>;
	hashFuncMatrix[60][4][2] = metroHash64<60, (1 << 30) - 1>;
	hashFuncMatrix[60][5][2] = cityHash64<60, (1 << 30) - 1>;
	hashFuncMatrix[64][1][2] = maRushPrime1HashSimplified<64, (1 << 30) - 1>;
	hashFuncMatrix[64][2][2] = xxhash32<64, (1 << 30) - 1>;
	hashFuncMatrix[64][3][2] = xxhash64<64, (1 << 30) - 1>;
	hashFuncMatrix[64][4][2] = metroHash64<64, (1 << 30) - 1>;
	hashFuncMatrix[64][5][2] = cityHash64<64, (1 << 30) - 1>;
	hashFuncMatrix[68][1][2] = maRushPrime1HashSimplified<68, (1 << 30) - 1>;
	hashFuncMatrix[68][2][2] = xxhash32<68, (1 << 30) - 1>;
	hashFuncMatrix[68][3][2] = xxhash64<68, (1 << 30) - 1>;
	hashFuncMatrix[68][4][2] = metroHash64<68, (1 << 30) - 1>;
	hashFuncMatrix[68][5][2] = cityHash64<68, (1 << 30) - 1>;
	hashFuncMatrix[72][1][2] = maRushPrime1HashSimplified<72, (1 << 30) - 1>;
	hashFuncMatrix[72][2][2] = xxhash32<72, (1 << 30) - 1>;
	hashFuncMatrix[72][3][2] = xxhash64<72, (1 << 30) - 1>;
	hashFuncMatrix[72][4][2] = metroHash64<72, (1 << 30) - 1>;
	hashFuncMatrix[72][5][2] = cityHash64<72, (1 << 30) - 1>;
	hashFuncMatrix[76][1][2] = maRushPrime1HashSimplified<76, (1 << 30) - 1>;
	hashFuncMatrix[76][2][2] = xxhash32<76, (1 << 30) - 1>;
	hashFuncMatrix[76][3][2] = xxhash64<76, (1 << 30) - 1>;
	hashFuncMatrix[76][4][2] = metroHash64<76, (1 << 30) - 1>;
	hashFuncMatrix[76][5][2] = cityHash64<76, (1 << 30) - 1>;
	hashFuncMatrix[80][1][2] = maRushPrime1HashSimplified<80, (1 << 30) - 1>;
	hashFuncMatrix[80][2][2] = xxhash32<80, (1 << 30) - 1>;
	hashFuncMatrix[80][3][2] = xxhash64<80, (1 << 30) - 1>;
	hashFuncMatrix[80][4][2] = metroHash64<80, (1 << 30) - 1>;
	hashFuncMatrix[80][5][2] = cityHash64<80, (1 << 30) - 1>;
	hashFuncMatrix[84][1][2] = maRushPrime1HashSimplified<84, (1 << 30) - 1>;
	hashFuncMatrix[84][2][2] = xxhash32<84, (1 << 30) - 1>;
	hashFuncMatrix[84][3][2] = xxhash64<84, (1 << 30) - 1>;
	hashFuncMatrix[84][4][2] = metroHash64<84, (1 << 30) - 1>;
	hashFuncMatrix[84][5][2] = cityHash64<84, (1 << 30) - 1>;
	hashFuncMatrix[88][1][2] = maRushPrime1HashSimplified<88, (1 << 30) - 1>;
	hashFuncMatrix[88][2][2] = xxhash32<88, (1 << 30) - 1>;
	hashFuncMatrix[88][3][2] = xxhash64<88, (1 << 30) - 1>;
	hashFuncMatrix[88][4][2] = metroHash64<88, (1 << 30) - 1>;
	hashFuncMatrix[88][5][2] = cityHash64<88, (1 << 30) - 1>;
	hashFuncMatrix[92][1][2] = maRushPrime1HashSimplified<92, (1 << 30) - 1>;
	hashFuncMatrix[92][2][2] = xxhash32<92, (1 << 30) - 1>;
	hashFuncMatrix[92][3][2] = xxhash64<92, (1 << 30) - 1>;
	hashFuncMatrix[92][4][2] = metroHash64<92, (1 << 30) - 1>;
	hashFuncMatrix[92][5][2] = cityHash64<92, (1 << 30) - 1>;
	hashFuncMatrix[96][1][2] = maRushPrime1HashSimplified<96, (1 << 30) - 1>;
	hashFuncMatrix[96][2][2] = xxhash32<96, (1 << 30) - 1>;
	hashFuncMatrix[96][3][2] = xxhash64<96, (1 << 30) - 1>;
	hashFuncMatrix[96][4][2] = metroHash64<96, (1 << 30) - 1>;
	hashFuncMatrix[96][5][2] = cityHash64<96, (1 << 30) - 1>;

}

#endif //!_HASHES_H
