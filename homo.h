
/*
 * homo.h for PolyScape
 * Copyright (c) 2013 Beifang Niu && Kai Ye WUGSC All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _HOMO_H_
#define _HOMO_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "param.h"

// homopolymer site
class HomoSite {
public:
    HomoSite();
    ~HomoSite();
    // length of homo
    bit8_t  typeLen;
    // homo or microsate 
    // A/C/G/T/AC/AGC
    bit16_t homoType;
    bit16_t length;
    bit16_t frontKmer;
    bit16_t endKmer;
    
    std::string chr;
    // readable
    std::string transfer;

    std::string bases;
    std::string fbases;
    std::string ebases;

    // location
    int location;

    // distribution
    unsigned short **normalDis;
    unsigned short **tumorDis;

    inline void InitType(){  };

    void TransferString();
    void InitialDis();
    void OutputDis();
    void ReleaseMemory();
    void PouroutDis(std::ofstream &fout);

    protected:
        // xxx
};

#endif //_HOMO_H_

