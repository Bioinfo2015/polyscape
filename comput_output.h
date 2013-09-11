
/*
 * comput_output.h for PolyScape
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


#ifndef _COMPUT_OUTPUT_H_
#define _COMPUT_OUTPUT_H_

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include "param.h"
#include "chi.h"

class DistributionPerSample {
public:
    DistributionPerSample () {
	SampleName = "";
	NormalCoverage = 0;
	TumorCoverage = 0;
	bool WithSufficientCoverage = false;
	bool NormalWithSufficientCoverage = false;
	double Difference = 0.0;
	PValue = 0.0;
	somatic = false;
	Genotype[0] = -2;
	Genotype[1] = -2;
	WithGenotype = false;
    }
    std::string SampleName;
    unsigned int NormalReadCount[100];
    unsigned int TumourReadCount[100];    
    unsigned NormalCoverage;
    unsigned TumorCoverage;
    bool WithSufficientCoverage;
    bool NormalWithSufficientCoverage;
    double Difference;
    double PValue;
    bool somatic;
    bool WithGenotype;
    int Genotype[2];
    //void CalculateDifference(); // difference between Normal and Tumor distribution if WithSufficientCoverage is true
    void UpdateAll();
    void ComputeGenotype();
    bool WhetherHomSiteWithLength(unsigned int Length);
    int GetLength() {
	if (WithGenotype == false) return -1;
	if (Genotype[0] == Genotype[1]) return 	Genotype[0];
	else return -1;
    }
};



// homopolymer site
class SiteInOneSample {
public:
    SiteInOneSample() {
	ChrName = "";
	pos = 0;
	front_kmer = "";
	end_kmer = "";
	repeat = "";
	NumberOfRepeats = 0;	
    }

    ~SiteInOneSample() {
	AllSamples.clear();
    }


    // common features
    std::string ChrName;
    unsigned int pos;
    std::string front_kmer;
    std::string end_kmer;
    std::string repeat;
    unsigned int NumberOfRepeats;

    std::vector< DistributionPerSample > AllSamples;

};


double DistanceBetweenTwo(unsigned int * FirstOriginal, unsigned int * SecondOriginal)
{

double SmallDouble = 0.0000000001;

/*
    std::cout << "in DistanceBetweenTwo" << std::endl;
    for (unsigned int i; i < dispots; i++) std::cout << FirstOriginal[i] << " ";
    std::cout << std::endl;
    for (unsigned int i; i < dispots; i++) std::cout << SecondOriginal[i] << " ";
    std::cout << std::endl;
*/     
    // delare
    double *Min;
    double *Max;

    double *FirstNormalized;
    double *SecondNormalized;

    unsigned int dispots = 100;

    FirstNormalized   = new double [dispots + 2];
    SecondNormalized  = new double [dispots + 2];

    unsigned int sumFirst = 0;
    unsigned int sumSecond = 0;



    // sum 
    for (unsigned i = 0; i< dispots; i++)
    {
        sumFirst  += FirstOriginal[i];
        sumSecond += SecondOriginal[i];
    }

    // normalization
    for (unsigned i = 0; i < dispots; i++)
    {	if (FirstOriginal[i] < 1) FirstNormalized[i]  = 0.0;
        else FirstNormalized[i]  = (FirstOriginal[i] / (double)sumFirst);

	if (SecondOriginal[i] < 1) SecondNormalized[i] = 0.0;
        else SecondNormalized[i] = SecondOriginal[i] / (double)sumSecond;
	//std::cout << "First  " << i << "\t" << FirstNormalized[i] << "\t" << FirstOriginal[i] << "\t" << sumFirst << std::endl;;
	//std::cout << "Second " << i << "\t" << SecondNormalized[i] << "\t" << SecondOriginal[i] << "\t" << sumSecond << std::endl;;
    }

    FirstNormalized[dispots]  = 0.0;
    SecondNormalized[dispots] = 0.0;

    Min  = new double [dispots + 2];
    Max  = new double [dispots + 2];

    // initialization 
    for (unsigned i=0; i<=dispots; i++)
    {
        Min[i] = 0;
        Max[i] = 0;
    }

    // get min and max per position
    for (unsigned i=0; i<=dispots; i++)
    {
        //cout << i << " " << First.Normalized[i] << " " << Second.Normalized[i] << endl;
        if (FirstNormalized[i] <= SecondNormalized[i])
        {
            Min[i] = FirstNormalized[i];
            Max[i] = SecondNormalized[i];
        }
        else
        {
            Min[i] = SecondNormalized[i];
            Max[i] = FirstNormalized[i];
        }
	if (Min[i] < SmallDouble) Min[i] = 0.0;
	if (Max[i] < SmallDouble) Max[i] = 0.0;
	
    }

    // caculate area
    double AreaMin = 0.0;
    double AreaMax = 0.0;

    for (unsigned i=0; i< dispots; i++)
    {
        AreaMin += (Min[i] + Min[i + 1]);
        AreaMax += (Max[i] + Max[i + 1]);
    }
    AreaMin = AreaMin;
    AreaMax = AreaMax;

    //std::cout << AreaMax << "\t" << AreaMin << "\t" << (AreaMax - AreaMin)/AreaMax << std::endl;
    if (AreaMax < AreaMin) {
	std::cout << "something is wrong with distance calculation " << AreaMax << " " << AreaMin << std::endl;
	for (unsigned index = 0; index < 100; index++) {
		std::cout << index << "\t" << FirstOriginal[index] << "\t" << SecondOriginal[index] << "\t" << sumFirst << "\t" << sumSecond << "\t" << FirstNormalized[index] << "\t" << SecondNormalized[index] << "\t" << Min[index] << "\t" << Max[index] << std::endl;
	}
    }


    delete [] Min;
    delete [] Max;

    delete [] FirstNormalized;
    delete [] SecondNormalized;

    return (AreaMax - AreaMin)/AreaMax;

}

bool DistributionPerSample::WhetherHomSiteWithLength(unsigned int Length) {
	if (Genotype[0] == -2 || Genotype[1] == -2) 
		ComputeGenotype();
	if (Genotype[0] == Length && Genotype[1] == Length)
		return true;
	else return false;
}

void DistributionPerSample::ComputeGenotype() {
//void GetGenotype(const std::vector< unsigned int > & dis, Genotype & genotype) {
    unsigned int Offset = 1;
    unsigned int CoverageCutoff = 20;
    unsigned int first  = 0;
    unsigned int second = 0;
    unsigned int Sum    = 0;

    //std::cout << "one site " << std::endl;

    // find the largest number 
    for (unsigned int pos_index = 0; pos_index < 100; pos_index++) { // NormalReadCount
	//if (pos_index < 30)
	//	std::cout << NormalReadCount[pos_index] << " ";
       Sum += NormalReadCount[pos_index]; 
       if (NormalReadCount[pos_index] > NormalReadCount[first]) {
          first = pos_index;
       }
    }
    //std::cout << std::endl << "Sum " << Sum << " first " << first << std::endl;
    if (Sum < CoverageCutoff) {
	Genotype[0] = -1;
	Genotype[1] = -1;
	WithGenotype = false;
	//std::cout << WithGenotype << "\t" << Genotype[0] << " " << Genotype[1] << std::endl;
	return;
    }
    if (first == 0) second = 1;
    for (unsigned int pos_index = 0; pos_index < 100; pos_index++) {
       if (pos_index == first) continue; 
       if (NormalReadCount[pos_index] > NormalReadCount[second]) {
          second = pos_index;
       }
    }
    float first_ratio  = NormalReadCount[first]  / (float) Sum;
    float second_ratio = NormalReadCount[second] / (float) Sum;
    //std::cout << first_ratio << " " <<  second_ratio << std::endl;
    if (first_ratio > 0.7) {
        Genotype[0] = first + Offset; 
        Genotype[1]  = first + Offset;
        WithGenotype = true; 
    }
    else if (first_ratio > 0.3 && second_ratio < first_ratio * 0.66) {
        Genotype[0] = first + Offset; 
        Genotype[1]  = first + Offset;
        WithGenotype = true; 
    }
    else if (first_ratio > 0.3 && second_ratio >= first_ratio * 0.66) {
        Genotype[0] = first + Offset; 
        Genotype[1]  = second + Offset; 
        WithGenotype = true; 
    }
    //std::cout << WithGenotype << "\t" << Genotype[0] << " " << Genotype[1] << std::endl;
    //std::cout << "end of one site \n" << std::endl;
    return;

} 

void DistributionPerSample::UpdateAll() {
	NormalCoverage = 0;
	TumorCoverage = 0;
	unsigned int cutoff = 20;
	for (unsigned int index = 0; index < 100; index++) {
		NormalCoverage += NormalReadCount[index];
		TumorCoverage += TumourReadCount[index];
	}
	Difference = 0.0;
	if (NormalCoverage >= 20) 
		NormalWithSufficientCoverage = true;
	else NormalWithSufficientCoverage = false;
	if (NormalCoverage >= 20 && TumorCoverage >= 20) {
		WithSufficientCoverage = true;
		Difference = DistanceBetweenTwo(NormalReadCount, TumourReadCount);
		//PValue = get_chisqr_p(NormalReadCount, TumourReadCount);
		PValue = X2BetweenTwo(NormalReadCount, TumourReadCount, 100);
		//PValue = 0;
		if (PValue < 0.001) somatic = true;
		//if (somatic) {
			/*
			std::cout << std::endl;
			for (unsigned index = 0; index < 100; index++) {
				std::cout << NormalReadCount[index] << " ";
			}
			std::cout << std::endl;
			for (unsigned index = 0; index < 100; index++) {
				std::cout << TumourReadCount[index] << " ";
			}
			std::cout << std::endl;

			std::cout << PValue << std::endl;
			*/
		//}
	}
	else {
		WithSufficientCoverage = false;
		Difference = -1.0;
		PValue = 1;
	}
	ComputeGenotype();
}



#endif //_COMPUT_OUTPUT_H_

