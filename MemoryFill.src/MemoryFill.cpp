/* ================================================================= *
 *  MemoryFill.cpp : Main file                                       *
 *                                                                   *
 *  MemoryFill is an auxiliary program that enables "cold start",    *
 *  especially for several runs with the same input files.           *
 *                                                                   *
 *                                                                   *
 *  Copyright (c) 2018, Szymon Grabowski and Wojciech Bieniecki      *
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
 
 
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <chrono>
#include <iostream>
#include <fstream>

void displayHelp(char* name) {
	std::cout << "Memory fill, by Szymon Grabowski and Wojciech Bieniecki, May 2018.\n Use: " << name << " <bytes to fill\n";
	std::cout << name << " 1024\n";
	std::cout << name << " 53k\n";
	std::cout << name << " 53M\n";
	std::cout << name << " 53G\n";
}

int main(int ac, char* av[])
{
	std::cout.setf(std::ios_base::unitbuf);
	std::cout << "Memory fill to disable disk caches (--> 'cold start')\n";
	if (ac == 1) {
		displayHelp(av[0]);
		return 0;
	}

	size_t l = strlen(av[1]);
	size_t mult = 1;
	switch (av[1][l - 1]) {
	case 'k':
		mult = 1 << 10;
		break;
	case 'M':
		mult = 1 << 20;
		break;
	case 'G':
		mult = 1 << 30;
		break;
	}
	size_t N = mult * std::atoll(av[1]);
	if (N == 0) {
		displayHelp(av[0]);
		return 0;
	}
	std::cout << N << " bytes\n";
	const size_t num = 100;
	const size_t fillSize = N / num;
	char** arr = new char*[num];
	for (size_t i = 0; i < num; i++) {
		std::cout << ".";
		arr[i] = new char[fillSize];
		if (arr[i] == 0) {
			std::cerr << "out of memory";
			return -1;
		}
		memset(arr[i], 170, fillSize);
	}
	for (size_t i = 0; i < num; i++)
		delete [] arr[i];
	std::this_thread::sleep_for(std::chrono::milliseconds(5000));
	std::cout << " Done!\n";
	return 0;
}
