/*
 * PerRepeatSize.h for PolyScape

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


#ifndef _PERREPEATSIZE_H_
#define _PERREPEATSIZE_H_

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include "param.h"
#include "comput_output.h"

const unsigned int BigNumber = 1000000;

class PerSize {
public:
	PerSize () {
		//unsigned int tempint;
		for (unsigned int index = 0;  index < 100; index++) {
			//distribution.push_back(0);
			distribution[index] = 0;
			//TooManyDataPoints = false;
		}
		TooManyDataPoints = false;
	}
	~PerSize () {
		//distribution.clear();
	}
	unsigned int distribution[100];
	double NormalizedDistribution[100];
	unsigned int NormalizedIntDistribution[100];
	bool TooManyDataPoints;

	void Add(const DistributionPerSample & one, unsigned int NumberOfRepeats) {
		unsigned Sum = 0;
		for (unsigned int index = 5; index < 100; index++) Sum += one.NormalReadCount[index];
		
		// if this one is not a pure site, ignore
		if ((one.NormalReadCount[NumberOfRepeats] / (double)Sum) < 0.8) return;

		if (TooManyDataPoints) return;
		for (unsigned int index = 5; index < 100; index++) {
			distribution[index] += one.NormalReadCount[index];
			if (distribution[index] > BigNumber)
				TooManyDataPoints = true;
		}
	}
	unsigned int * GetDistribution() {
		return distribution;
	}
	
	void ReportDistribution(std::ofstream & output) {
		output << TooManyDataPoints << " ";
		for (unsigned int index = 0; index < 100; index++)
			output << "\t" << distribution[index];
		output << std::endl;
	}

	void ReportNormalizedDistribution(std::ofstream & output) {
		output.precision(5);
		output << TooManyDataPoints << " ";
		unsigned int Sum = 0;
		for (unsigned int index = 0; index < 100; index++) {
			Sum += distribution[index];
		}
		double DoubleSum = Sum;
		for (unsigned int index = 0; index < 100; index++) {
			NormalizedDistribution[index] = distribution[index] / DoubleSum;
		}
		for (unsigned int index = 0; index < 100; index++)
			output << "\t" << std::fixed << NormalizedDistribution[index];
		output << std::endl;
	}
	unsigned int * GetNormalizedDistribution() {
		unsigned int MaxValue = 0;
		for (unsigned int index = 0; index < 100; index++) {
			if (distribution[index] > MaxValue) MaxValue = distribution[index]; 
		}
		unsigned int Wei = (unsigned int)log10((double)MaxValue);
		if (Wei <= 4) {
			
			return distribution;
		}
		else {
			unsigned int Divid = pow(10, Wei - 4);
			for (unsigned int index = 0; index < 100; index++) {
				NormalizedIntDistribution[index] = (unsigned int)(distribution[index]/Divid);
			}
			return NormalizedIntDistribution;
		}
	}

protected:
};



class PerRepeatSize { // AT or ATT, 
public:	
	PerRepeatSize() {
		PerSize TempOne;
		for (unsigned int index = 0; index < 100; index++) { // repeated 100 times?
			AllSize.push_back(TempOne);
		}
	}
	~PerRepeatSize() {
		AllSize.clear();
	}

	std::vector <PerSize> AllSize;

	void Add(unsigned int NumberOfRepeats, const DistributionPerSample & one) {
		AllSize[NumberOfRepeats].Add(one, NumberOfRepeats);
	}
	unsigned int * GetDistribution(unsigned int NumberOfRepeats) {
		return AllSize[NumberOfRepeats].GetDistribution();
	}
	unsigned int * GetNormalizedDistribution(unsigned int NumberOfRepeats) {
		return AllSize[NumberOfRepeats].GetNormalizedDistribution();
	}
	void ReportDistribution(unsigned int NumberOfRepeats, std::ofstream & output) {
		AllSize[NumberOfRepeats].ReportDistribution(output);
	}
	void ReportNormalizedDistribution(unsigned int NumberOfRepeats, std::ofstream & output) {
		AllSize[NumberOfRepeats].ReportNormalizedDistribution(output);
	}

protected:
};

#endif //_PERREPEATSIZE_H_
