/*
 * homo.cpp for PolyScape 
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

#include <iostream>
#include <sstream>
#include <bitset>
#include <omp.h>
#include "homo.h"
#include "polyscan.h"

extern Param param;
extern Param paramd;
extern PolyScan polyscan;
extern char homo_code[];

HomoSite::HomoSite()
    : homoType(0)
    , length(0)
    , frontKmer(0)
    , endKmer(0)
    , typeLen(0)
    , chr("")
    , transfer("")
    , bases("")
    , fbases("")
    , ebases("")
    , location(0)
    , normalDis(NULL)
    , tumorDis(NULL)
{
    InitType();
};

HomoSite::~HomoSite() {
    // xxxxx
};

// transfer binary to string
void HomoSite::TransferString() {
    // assign memory
    //chr.resize(MAX_TRANSFER_LINE_LENGTH);
    std::stringstream ss;
    bit8_t tch = 0;
    char tchbuff[MAX_FLANK_REGION];
    for (int i=typeLen-1; i>=0; i--) {
        tch = 0;
        tch = homoType&3;
        tchbuff[i] = homo_code[tch];
        homoType = homoType>>2;
    }
    // load type
    for (int i=0; i<typeLen; i++) {
        ss<<tchbuff[i];
    }
    ss << '\t';
    for (int i=param.ContextLength-1; i>=0; i--) {
        tch = 0;
        tch = frontKmer&3;
        tchbuff[i] = homo_code[tch];
        frontKmer = frontKmer>>2;
    }
    // load front flank region
    for (int i=0; i<param.ContextLength; i++) {
        ss<<tchbuff[i];
    }
    ss << '\t';
    for (int i=param.ContextLength-1; i>=0; i--) {
        tch = 0;
        tch = endKmer&3;
        tchbuff[i] = homo_code[tch];
        endKmer = endKmer>>2;
    }
    // load end flank region
    for (int i=0; i<param.ContextLength; i++) {
        ss<<tchbuff[i];
    }
    transfer = ss.str();
    
    ss.clear();
    ss.str("");
}

// initial distribution
void HomoSite::InitialDis() {
    normalDis = new unsigned short *[polyscan.totalBamPairsNum];
    tumorDis  = new unsigned short *[polyscan.totalBamPairsNum];
    for (unsigned int j=0; j<polyscan.totalBamPairsNum; j++) {
        normalDis[j] = new unsigned short [paramd.s_dispots];
        tumorDis[j]  = new unsigned short [paramd.s_dispots];
        for (unsigned int k=0; k<paramd.s_dispots; k++) {
            normalDis[j][k] = 0;
            tumorDis[j][k]  = 0;
        }
    }
}

// Out distribution
void HomoSite::OutputDis() {
    for (unsigned int j=0; j<polyscan.totalBamPairsNum; j++) {
        for (unsigned int k=0; k<paramd.s_dispots; k++) {
            std::cout<<normalDis[j][k];
        }
        for (unsigned int k=0; k<paramd.s_dispots; k++) {
            std::cout<<tumorDis[j][k];
        }
    }
}

// Pourout distribution
void HomoSite::PouroutDis(std::ofstream &fout) {
    fout << chr << " "
         << location << " "
         << fbases << " "
         << length << "["
         << bases  << "] "
         << ebases <<"\t";
    for (unsigned int j=0; j<polyscan.totalBamPairsNum; j++) {
        fout << polyscan.totalBamPairs[j].sName << "|";
        for (unsigned int k=0; k<paramd.s_dispots; k++) {
            fout <<normalDis[j][k] << " ";
        }
        fout << "|";
        for (unsigned int k=0; k<paramd.s_dispots; k++) {
            fout<<tumorDis[j][k] << " ";
        }
        fout << "\t";
    }
    fout << "\n";

}

// release memory
void HomoSite::ReleaseMemory() {
    for (unsigned int k=0; k<polyscan.totalBamPairsNum; k++) {
        delete [] normalDis[k];
        delete []  tumorDis[k]; 
    }
    delete [] normalDis;
    delete [] tumorDis;
}

