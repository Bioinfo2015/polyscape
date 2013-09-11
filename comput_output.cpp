/*
 * comput_output.cpp for PolyScape
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

// System header files
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>
#include <unistd.h>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>

// Static function declaration
#include "param.h"
#include "polyscan.h"
#include "distribution.h"
#include "utilities.h"
#include "comput_output.h"
#include "perrepeatsize.h"
#include "chi.h"

// branch 
#include "cmds.h"

//Param paramd;
//PolyScan polyscan;

std::string DisFile;
std::string OutputFile;
std::string PairFile;
//std::string disFile;

//std::ifstream finH;
//std::ifstream finM;
//std::ifstream finB;
//std::ofstream foutD;

//std::string one_region;

void ComUsage(void) {
    std::cerr<<"\nUsage:  plolyscape com [options] \n\n"

        <<"       -i   <string>   homopolymer and microsates distribution file\n"
        <<"       -p   <string>   pairing file: two bam files per line\n"
        <<"       -o   <string>   output file\n\n"
        <<"       \n"
 
       <<"       -h   help\n\n"
        << std::endl;
    exit(1);
}


int cGetOptions(int rgc, char *rgv[]) {
    int i;
    for (i=1; i<rgc; i++) {
        if (rgv[i][0] != '-') return i;
        switch(rgv[i][1]) {
            case 'i': DisFile = rgv[++i]; break;
            case 'p': PairFile = rgv[++i]; break;
            case 'o': OutputFile  = rgv[++i]; break;
            break;
            case 'h':ComUsage();
            case '?':ComUsage();    
        }
    }
    return i;
}

void printdistrubition(unsigned int * input) {
	for (unsigned i = 0; i < 30; i++) {
		std::cout << input[i] << " ";
	}
	std::cout << std::endl;
}

void UpdateClusterAssignment(std::vector <double> & AllDistances, std::vector <double> & Centers, std::vector <unsigned> & ClusterAssignment) {
	unsigned int TempClusterID;
	double TempDistance, TempDistanceMin;
	std::cout << "center: ";
	for (unsigned index = 0; index < Centers.size(); index++) std::cout << Centers[index] << " ";
	std::cout << std::endl;
	for (unsigned i = 0; i < AllDistances.size(); i++) {
		TempClusterID = 0;
		TempDistanceMin = AllDistances[i] - Centers[0];
		if (TempDistanceMin < 0) TempDistanceMin = TempDistanceMin * (-1);  
		//std::cout << 0 << " " << TempDistanceMin << std::endl;
		for (unsigned int CenterIndex = 0; CenterIndex < Centers.size(); CenterIndex++) {
			
			TempDistance = AllDistances[i] - Centers[CenterIndex];
			if (TempDistance < 0) TempDistance = TempDistance * (-1);
			//std::cout << CenterIndex << "\t" << TempDistance << std::endl;
			if (TempDistance < TempDistanceMin) {
				TempClusterID = CenterIndex;
				TempDistanceMin = TempDistance;
				//if (TempDistance > 0.5) std::cout << TempDistance << "\t" << TempClusterID << "\t" << TempDistanceMin << std::endl; 
			}
		}
		ClusterAssignment[i] = TempClusterID;
		//std::cout << AllDistances[i] << ":" << TempClusterID << " ";
	}
	//std::cout << std::endl;
}

void UpdateClusterCenter(std::vector <double> & AllDistances, std::vector <double> & Centers, std::vector <unsigned> & ClusterAssignment) {
	
	std::vector <double> SumOfDistances;
	std::vector <unsigned int> NumberOfInstances;
	for (unsigned CenterIndex = 0; CenterIndex < Centers.size(); CenterIndex++) {
		SumOfDistances.push_back(0.0);
		NumberOfInstances.push_back(0);
	}
	for (unsigned int index = 0; index < AllDistances.size(); index++) {
		SumOfDistances[ClusterAssignment[index]] +=  AllDistances[index];
		NumberOfInstances[ClusterAssignment[index]]++;
	}
	for (unsigned CenterIndex = 0; CenterIndex < Centers.size(); CenterIndex++) {
		Centers[CenterIndex] = SumOfDistances[CenterIndex] / NumberOfInstances[CenterIndex];
		std::cout << CenterIndex << "\t" << SumOfDistances[CenterIndex] << "\t" << NumberOfInstances[CenterIndex] << "\t" << Centers[CenterIndex] << std::endl;
	}
}

void GetClusterCenter(std::vector <double> & AllDistances, std::vector <double> & Centers) {
	std::sort(AllDistances.begin(), AllDistances.end());
	//initialize center as the beginning, middle and end, and so on.
	unsigned int NumberOfClusters = Centers.size();
	//for (unsigned int index = 0; index < AllDistances.size(); index++) std::cout << AllDistances[index] << " ";
	//std::cout << std::endl;
	unsigned int StepSize = AllDistances.size() / NumberOfClusters;
	double step = 1.0 / (NumberOfClusters + 1);
	for (unsigned i = 0; i < NumberOfClusters; i++) {
		Centers[i] =  i * step;
		//std::cout << "GetClusterCenter" << "\t" << i << "\t" << Centers[i] << std::endl; 
	}
	std::vector <unsigned> ClusterAssignment (AllDistances.size(), 0);
	std::vector <unsigned> PreviousClusterAssignment(ClusterAssignment);
	for (unsigned int IternationCounter = 0; IternationCounter < 20; IternationCounter++) {
		UpdateClusterAssignment(AllDistances, Centers, ClusterAssignment);
		UpdateClusterCenter(AllDistances, Centers, ClusterAssignment);
		//for (unsigned index = 0; index < Centers.size(); index++) {
			//std::cout << "center " << IternationCounter << "\t" << index << "\t" << Centers[index] << "; ";
		//}
		//std::cout << std::endl;
		if (ClusterAssignment == PreviousClusterAssignment) break;
		else PreviousClusterAssignment = ClusterAssignment;
	}
}

void GetCutOffs(std::vector <double> & AllDistances, std::vector <double> & CutOffs) {
	unsigned NumOfCutoffs = CutOffs.size();
	unsigned NumOfClusters = NumOfCutoffs + 1;
	std::vector <double> Centers;
	for (unsigned i = 0; i < NumOfClusters; i++) 
		Centers.push_back(0.0);

	GetClusterCenter(AllDistances, Centers);
	for (unsigned i = 0; i < NumOfClusters - 1; i++) {
		CutOffs[i] = (Centers[i] + Centers[i + 1]) / 2;
		std::cout << Centers[i] << "\t" << CutOffs[i] << std::endl;
	}
	//std::cout << CutOffs[0] << "\t" << CutOffs[1] << std::endl;
}

void Display(unsigned int * ReferenceDistribution, unsigned int * NormalReadCount, unsigned int NumDataPoint, double PValue) {
	std::cout << PValue << std::endl;
	for (unsigned int i = 0; i < 30; i++) {
		std::cout << ReferenceDistribution[i] << "\t";
		//std::cout << NormalReadCount[i] << std::endl;
	}
	std::cout << std::endl;
	for (unsigned int i = 0; i < 30; i++) {
		//std::cout << ReferenceDistribution[i] << "\t";
		std::cout << NormalReadCount[i] << "\t";
	}
	std::cout << std::endl;
}


int HomoAndMicrosateCom(int argc, char *argv[]) {

    if (argc == 1) ComUsage();
    for (int i=0; i<argc; i++) {
        std::cout <<argv[i]<<' ';
    }
    Initial_Time();
    std::cout << "\nStart at:  " << Curr_Time() << std::endl;
    int noptions = cGetOptions(argc, argv);

    std::ifstream input(DisFile.c_str());
    std::ifstream input2(DisFile.c_str());
    std::ifstream pair(PairFile.c_str());
    std::ofstream output(OutputFile.c_str());
    std::ofstream output_p_somatic((OutputFile + "_p_somatic").c_str());
    //std::ofstream output_p_germline((OutputFile + "_p_germline").c_str());
    std::ofstream output_somatic((OutputFile + "_somatic").c_str());
    std::ofstream output_germline((OutputFile + "_germline").c_str());
    std::ofstream output_distribution((OutputFile + "_dis").c_str());
    unsigned NumberOfPairs = 0;
    std::string TempStr, TempSampleName;
    std::vector <std::string> SampleNames;
    while (pair >> TempSampleName) {
	getline(pair, TempStr);
	SampleNames.push_back(TempSampleName);
	NumberOfPairs++;
    }
    std::cout << "There are " << NumberOfPairs << " samples" << std::endl;

    std::vector <double> SumDifference;
    std::vector <unsigned int> NumberOfMSIDataPoints;
    std::vector <unsigned int> NumberOfDataPoints;
    for (unsigned int index = 0; index < NumberOfPairs; index++) {
    	SumDifference.push_back(0.0);
	NumberOfDataPoints.push_back(0);
	NumberOfMSIDataPoints.push_back(0);
    }

    unsigned int tempint;

    std::vector <double> AllDistances;

    unsigned PrecisionNumber = 5;
    std::cout.precision(PrecisionNumber);
    output.precision(2);
    output_p_somatic.precision(PrecisionNumber);
    //output_p_germline.precision(PrecisionNumber);
    output_somatic.precision(PrecisionNumber);
    output_germline.precision(PrecisionNumber);
    //output_distribution.precision(PrecisionNumber);
    // store the consensus
    std::cout << "initialize std::vector <PerRepeatSize> AllRepeatSizes" << std::endl;
    std::vector <PerRepeatSize> AllRepeatSizes;
    for (unsigned int index = 0; index < 100; index++) { // this 100 should be replaced by the actual maximum length of repeat unit. 
	PerRepeatSize temp_one;
	AllRepeatSizes.push_back(temp_one);
    }

    SiteInOneSample tempone;
    std::vector <SiteInOneSample> AllSites;
    std::cout << "scan for germline and somatic variants" << std::endl;
    unsigned int Count = 0;
    while (input >> tempone.ChrName >> tempone.pos >> tempone.front_kmer >> tempone.NumberOfRepeats >> tempone.repeat >> tempone.end_kmer) {
	if (Count % 1000 == 0 && Count) std::cout << Count << std::endl;
	Count++;
	//std::cout << tempone.ChrName << "\t" << tempone.pos << std::endl;
	tempone.AllSamples.clear();
	unsigned RepeatSize = tempone.repeat.size();


	bool ReportSomatic = false;
	bool ReportGermline = false;
	for (unsigned PairIndex = 0; PairIndex < NumberOfPairs; PairIndex++) {
		
		DistributionPerSample OneSample;
		input >> OneSample.SampleName >> TempStr;
		for (unsigned index = 0; index < 100; index++) {
			input >> tempint;
			OneSample.NormalReadCount[index] = tempint;
			if (index < 6) OneSample.NormalReadCount[index] = 0; // remove this once beifang changed his code
		}
		input >> TempStr;
		for (unsigned index = 0; index < 100; index++) {
                	input >> tempint;
 	                OneSample.TumourReadCount[index] = tempint;
			if (index < 6) OneSample.TumourReadCount[index] = 0; // remove this once beifang changed his code
       		}
		OneSample.UpdateAll();

		if (OneSample.WithSufficientCoverage) {
			NumberOfDataPoints[PairIndex]++;
			AllDistances.push_back(OneSample.Difference);
			output_p_somatic << log10(OneSample.PValue) << std::endl;
			 if (OneSample.somatic) {
				//std::cout << "
				NumberOfMSIDataPoints[PairIndex]++;
				SumDifference[PairIndex] += OneSample.Difference;
			 }
		} 
		// sum up read counts for each repeat size and each number of repeats. 
		//if (OneSample.NormalWithSufficientCoverage && OneSample.GetLength() == tempone.NumberOfRepeats) {
		//	AllRepeatSizes[RepeatSize].Add(tempone.NumberOfRepeats, OneSample);
		//}
		tempone.AllSamples.push_back(OneSample);
		if (OneSample.WithGenotype) ReportGermline = true;
		if (OneSample.somatic) ReportSomatic = true; 
	}
	if (ReportSomatic) {
		output_somatic << tempone.ChrName << "\t" << tempone.pos << "\t" << tempone.front_kmer 
			<< "\t" << tempone.NumberOfRepeats << "\t" << tempone.repeat << "\t" << tempone.end_kmer;
		for (unsigned PairIndex = 0; PairIndex < NumberOfPairs; PairIndex++) {
			output_somatic << "\t" << std::fixed << tempone.AllSamples[PairIndex].Difference;
		}
		output_somatic << std::endl;
	}
	if (ReportGermline) {
		output_germline << tempone.ChrName << "\t" << tempone.pos << "\t" << tempone.front_kmer 
			<< "\t" << tempone.NumberOfRepeats << "\t" << tempone.repeat << "\t" << tempone.end_kmer;
		for (unsigned PairIndex = 0; PairIndex < NumberOfPairs; PairIndex++) {
			output_germline << "\t" << tempone.AllSamples[PairIndex].Genotype[0] << "|" << tempone.AllSamples[PairIndex].Genotype[1];
		}
		output_germline << std::endl;
	}
    }

    std::cout << "output distribution" << std::endl;
    //std::ofstream output_d((OutputFile + "_d").c_str()); 
    for (unsigned int RepeatSize = 0; RepeatSize < 6; RepeatSize++) {
	for (unsigned NumberOfRepeat = 0; NumberOfRepeat < 20; NumberOfRepeat++) {
		output_distribution << RepeatSize << ":" << NumberOfRepeat << "\t";
		AllRepeatSizes[RepeatSize].ReportDistribution(NumberOfRepeat, output_distribution);
	}
    }

/*
    std::cout << "get germline varinats" << std::endl;
    //std::ofstream output_p_germline((OutputFile + "_P_Germline").c_str()); 
    double PValue;
    //ReportNormalizedDistribution
    // scan the file for the second time to call germline
    while (input2 >> tempone.ChrName >> tempone.pos >> tempone.front_kmer >> tempone.NumberOfRepeats >> tempone.repeat >> tempone.end_kmer) {
	output_germline << tempone.ChrName << "\t" << tempone.pos << "\t" << tempone.front_kmer << "\t" 
		<< tempone.NumberOfRepeats << "\t" << tempone.repeat << "\t" << tempone.end_kmer;
	//std::endl;
	//tempone.AllSamples.clear();
	unsigned RepeatSize = tempone.repeat.size();
	unsigned int * ReferenceDistribution = AllRepeatSizes[RepeatSize].GetDistribution(tempone.NumberOfRepeats);
	for (unsigned PairIndex = 0; PairIndex < NumberOfPairs; PairIndex++) {
		
		DistributionPerSample OneSample;
		input2 >> OneSample.SampleName >> TempStr;
		for (unsigned index = 0; index < 100; index++) {
			input2 >> tempint;
			OneSample.NormalReadCount[index] = tempint;
		}
		input2 >> TempStr;
		for (unsigned index = 0; index < 100; index++) {
                	input2 >> tempint;
 	               OneSample.TumourReadCount[index] = tempint;
       		}
		OneSample.UpdateAll();

		// sum up read counts for each repeat size and each number of repeats. 
		if (OneSample.NormalWithSufficientCoverage) {
			//AllRepeatSizes[RepeatSize].Add(tempone.NumberOfRepeats, OneSample);
			PValue = X2BetweenTwo(ReferenceDistribution, OneSample.NormalReadCount, 100);
			if (PValue < 0.001) output_germline << "\t1";
			else output_germline << "\t0";
			output_p_germline << log10(PValue) << std::endl;
			//if (PValue < 0.001) 
			//	Display(ReferenceDistribution, OneSample.NormalReadCount, 100, PValue);
		}
		else output_germline << "\t-1";
	//	tempone.AllSamples.push_back(OneSample);
	}
	output_germline << std::endl;
    }
*/
    std::string ChrName;
    unsigned int pos;
    std::string front_kmer;
    std::string end_kmer;
    std::string repeat;
    unsigned int NumberOfRepeats;

	std::cout << "output difference" << std::endl;  
    output << "SampleID\tTotal_Number_of_Sites\tNumber_of_Somatic_Sites\t%" << std::endl; 
    for (unsigned PairIndex = 0; PairIndex < NumberOfPairs; PairIndex++) {
	output << SampleNames[PairIndex] << "\t" << NumberOfDataPoints[PairIndex] << "\t" << NumberOfMSIDataPoints[PairIndex] << "\t" << std::fixed << (NumberOfMSIDataPoints[PairIndex] / (double)NumberOfDataPoints[PairIndex]) * 100.0 << std::endl;//"\t";// << SumDifference[PairIndex] << "\t"; 
    	//if (NumberOfMSIDataPoints[PairIndex]) output << SumDifference[PairIndex] / NumberOfMSIDataPoints[PairIndex] << std::endl;
	//else output << 0.0 << std::endl;
    }
    
    return 0;
}


