/* ================================================================ *
*  StringUtils.cpp : helper file                                    *
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
#include <algorithm>
#include <cstring>
#include <string>


void reverseComplement(char* start, const char* lut, const size_t N) {
	char* left = start + 1; // sequence starts from paddingChar
	char* right = start + N - 1;
	while (right > left) {
		char tmp = lut[*left];
		*left = lut[*right];
		*right = tmp;
		++left;
		--right;
	}
	if (left == right)
		*left = lut[*left];
}


void prepareComplement(char* complement, const char paddingChar) {
	memset(complement, paddingChar, 256);
	complement['A'] = 'T';
	complement['a'] = 'T';
	complement['C'] = 'G';
	complement['c'] = 'G';
	complement['G'] = 'C';
	complement['g'] = 'C';
	complement['T'] = 'A';
	complement['t'] = 'A';
	complement['N'] = 'N';
	complement['n'] = 'N';

}


void replaceBadSymbol(char* gen, char* dst, const char symbol, const char paddingChar) {
	char* movingPtr = gen;
	while (true) {
		char* tempPtr = std::find(movingPtr, dst, symbol);
		while (*tempPtr == symbol) {
			*tempPtr = paddingChar;
			++tempPtr;
		}
		if (tempPtr == dst)
			break;
		movingPtr = tempPtr + 1;
	}
}


const char* len_lut = "000001002003004005006007008009010011012013014015016017018019020021022023024025026027028029030031032033034035036037038039040041042043044045046047048049050051052053054055056057058059060061062063064065066067068069070071072073074075076077078079080081082083084085086087088089090091092093094095096097098099100101102103104105106107108109110111112113114115116117118119120121122123124125126127128129130131132133134135136137138139140141142143144145146147148149150151152153154155156157158159160161162163164165166167168169170171172173174175176177178179180181182183184185186187188189190191192193194195196197198199200201202203204205206207208209210211212213214215216217218219220221222223224225226227228229230231232233234235236237238239240241242243244245246247248249250251252253254255256257258259260261262263264265266267268269270271272273274275276277278279280281282283284285286287288289290291292293294295296297298299300301302303304305306307308309310311312313314315316317318319320321322323324325326327328329330331332333334335336337338339340341342343344345346347348349350351352353354355356357358359360361362363364365366367368369370371372373374375376377378379380381382383384385386387388389390391392393394395396397398399400401402403404405406407408409410411412413414415416417418419420421422423424425426427428429430431432433434435436437438439440441442443444445446447448449450451452453454455456457458459460461462463464465466467468469470471472473474475476477478479480481482483484485486487488489490491492493494495496497498499500501502503504505506507508509510511512513514515516517518519520521522523524525526527528529530531532533534535536537538539540541542543544545546547548549550551552553554555556557558559560561562563564565566567568569570571572573574575576577578579580581582583584585586587588589590591592593594595596597598599600601602603604605606607608609610611612613614615616617618619620621622623624625626627628629630631632633634635636637638639640641642643644645646647648649650651652653654655656657658659660661662663664665666667668669670671672673674675676677678679680681682683684685686687688689690691692693694695696697698699700701702703704705706707708709710711712713714715716717718719720721722723724725726727728729730731732733734735736737738739740741742743744745746747748749750751752753754755756757758759760761762763764765766767768769770771772773774775776777778779780781782783784785786787788789790791792793794795796797798799800801802803804805806807808809810811812813814815816817818819820821822823824825826827828829830831832833834835836837838839840841842843844845846847848849850851852853854855856857858859860861862863864865866867868869870871872873874875876877878879880881882883884885886887888889890891892893894895896897898899900901902903904905906907908909910911912913914915916917918919920921922923924925926927928929930931932933934935936937938939940941942943944945946947948949950951952953954955956957958959960961962963964965966967968969970971972973974975976977978979980981982983984985986987988989990991992993994995996997998999";


uint32_t int2cstr(char* tmpbuf, char* outbuf, uint32_t n) {
	char* reversed_tmp = outbuf;

	size_t tmp_pos = 0;
	do {
		uint32_t cur = n % 1000;
		n /= 1000;
		memcpy(tmpbuf + tmp_pos, len_lut + 3LL * cur, 3);
		tmp_pos += 3;
	} while (n > 0);

	size_t tmp_pos_end = tmp_pos;
	tmp_pos -= 3;
	while (*(tmpbuf + tmp_pos) == '0')
		++tmp_pos;

	size_t reversed_tmp_pos = 0;

	// copying the first (which is last in the tmp buffer!) radix-1000 digit (1..3 decimal digits)
	memcpy(reversed_tmp + reversed_tmp_pos, tmpbuf + tmp_pos, tmp_pos_end - tmp_pos);

	reversed_tmp_pos += (tmp_pos_end - tmp_pos);

	tmp_pos_end -= 3;
	while (tmp_pos_end > 0) {
		tmp_pos_end -= 3;
		memcpy(reversed_tmp + reversed_tmp_pos, tmpbuf + tmp_pos_end, 3);
		reversed_tmp_pos += 3;
	}
	//reversed_tmp[reversed_tmp_pos] = 0;
	return reversed_tmp_pos;
}

char* int2cstr_c(char* tmpbuf, char* outbuf, uint32_t n) {
	char* reversed_tmp = outbuf;

	if (n == 0)  [[unlikely]]
	{
		reversed_tmp[0] = '0';
		reversed_tmp[1] = 0;  // terminator
		return reversed_tmp;
	}

	size_t tmp_pos = 0;
	do {
		uint32_t cur = n % 1000;
		n /= 1000;
		memcpy(tmpbuf + tmp_pos, len_lut + 3LL * cur, 3);
		tmp_pos += 3;
	} while (n > 0);

	size_t tmp_pos_end = tmp_pos;
	tmp_pos -= 3;
	while (*(tmpbuf + tmp_pos) == '0')
		++tmp_pos;

	size_t reversed_tmp_pos = 0;

	// copying the first (which is last in the tmp buffer!) radix-1000 digit (1..3 decimal digits)
	memcpy(reversed_tmp + reversed_tmp_pos, tmpbuf + tmp_pos, tmp_pos_end - tmp_pos);

	reversed_tmp_pos += (tmp_pos_end - tmp_pos);

	tmp_pos_end -= 3;
	while (tmp_pos_end > 0) {
		tmp_pos_end -= 3;
		memcpy(reversed_tmp + reversed_tmp_pos, tmpbuf + tmp_pos_end, 3);
		reversed_tmp_pos += 3;
	}
	reversed_tmp[reversed_tmp_pos] = 0;
	return reversed_tmp;
}


char* int2cstr_t(char* tmpbuf, char* outbuf, uint32_t n) {
	// append \t before the number!
	outbuf[0] = '\t';
	char* reversed_tmp = outbuf + 1;

	int tmp_pos = 0;
	while (true) {
		int cur = n % 1000;
		n /= 1000;
		memcpy(tmpbuf + tmp_pos, len_lut + 3 * cur, 3);
		tmp_pos += 3;
		if (n == 0)
			break;
	}
	int tmp_pos_end = tmp_pos;
	tmp_pos -= 3;
	while (*(tmpbuf + tmp_pos) == '0')
		++tmp_pos;

	int reversed_tmp_pos = 0;

	// copying the first (which is last in the tmp buffer!) radix-1000 digit (1..3 decimal digits)
	memcpy(reversed_tmp + reversed_tmp_pos, tmpbuf + tmp_pos, tmp_pos_end - tmp_pos);

	reversed_tmp_pos += (tmp_pos_end - tmp_pos);

	tmp_pos_end -= 3;
	while (tmp_pos_end > 0) {
		tmp_pos_end -= 3;
		memcpy(reversed_tmp + reversed_tmp_pos, tmpbuf + tmp_pos_end, 3);
		reversed_tmp_pos += 3;
	}
	reversed_tmp[reversed_tmp_pos] = 0;  // terminator
	return outbuf;
}


char* int2cstr_t_n(char* tmpbuf, char* outbuf, uint32_t n){
	outbuf[0] = '\t';
	char* reversed_tmp = outbuf + 1;

	int tmp_pos = 0;
	while (true) {
		int cur = n % 1000;
		n /= 1000;
		memcpy(tmpbuf + tmp_pos, len_lut + 3 * cur, 3);
		tmp_pos += 3;
		if (n == 0)
			break;
	}
	int tmp_pos_end = tmp_pos;
	tmp_pos -= 3;
	while (*(tmpbuf + tmp_pos) == '0')
		++tmp_pos;

	int reversed_tmp_pos = 0;

	// copying the first (which is last in the tmp buffer!) radix-1000 digit (1..3 decimal digits)
	memcpy(reversed_tmp + reversed_tmp_pos, tmpbuf + tmp_pos, tmp_pos_end - tmp_pos);

	reversed_tmp_pos += (tmp_pos_end - tmp_pos);

	tmp_pos_end -= 3;
	while (tmp_pos_end > 0) {
		tmp_pos_end -= 3;
		memcpy(reversed_tmp + reversed_tmp_pos, tmpbuf + tmp_pos_end, 3);
		reversed_tmp_pos += 3;
	}

	reversed_tmp[reversed_tmp_pos] = '\n';
	reversed_tmp[reversed_tmp_pos + 1] = 0;  // terminator

	return outbuf;
}
