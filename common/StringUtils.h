/* ================================================================ *
*  StringUtils.cpp : helper header file                             *
*                                                                   *
*  copMEM2 is a program for efficient computation of MEMs           *
*  (Maximal Exact Matches) in a pair of genomes.                    *
*  Its main algorithmic idea requires that two internal parameters  *
*  (k1 and k2) are coprime, hence the name.                         *
*                                                                   *
*                                                                   *
*  Copyright (c) 2018-2022, Szymon Grabowski and Wojciech Bieniecki *
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
void reverseComplement(char* start, const char* lut, const size_t N);
void prepareComplement(char* complement, const char paddingChar);
void replaceBadSymbol(char* gen, char* dst, const char symbol, const char paddingChar);
void replaceBadSymbols(char* gen, char* dst, std::vector<char> symbols, const char paddingChar);
uint32_t int2cstr(char* tmpbuf, char* buf, uint32_t n);
char* int2cstr_c(char* tmpbuf, char* buf, uint32_t n);
char* int2cstr_t(char* tmpbuf, char* buf, uint32_t n);
char* int2cstr_t_n(char* tmpbuf, char* buf, uint32_t n);
