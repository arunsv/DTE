
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "TimingToSequence.hpp"

#define MIN_TIME_SLICE 1e-6 //Minimum time slice to accomodate boundary cases

TimingToSequence::TimingToSequence()
{
	timingFileName = "";
    uCalcMode = 1;
    dDelay = 0.0;
}

TimingToSequence::~TimingToSequence()
{
    
}

void TimingToSequence::setFileName(std::string fName)
{
	timingFileName = fName;
}

void TimingToSequence::setBinSize(double dBSize)
{
    dBinSize = dBSize;
}

void TimingToSequence::setDelay(double dValue)
{
	dDelay = dValue;
}

void TimingToSequence::setCalcMode(unsigned short int uMode)
{
    uCalcMode = uMode;
}

bool TimingToSequence::getTiming()
{
	bool bFlag = readActivityTimingFromFile();
    
	if(!bFlag) std::cout << "Unable to fetch timing data." << std::endl;
    
    return bFlag;
}

bool TimingToSequence::readActivityTimingFromFile()
{
	bool bFlag=false;
	std::string line;
    
	std::ifstream ipFile(timingFileName.c_str());
    
	if (ipFile.is_open())
	{
		while (std::getline(ipFile,line)) // Get each line
		{
            std::stringstream ss;
            ss.str(line);
            std::vector<double> vdTemp;
            std::string sTemp;
            
            while(std::getline(ss,sTemp,','))
            {
                vdTemp.push_back(atof(sTemp.c_str()));
            }
            
            vTiming.push_back(vdTemp);
            
        }
        
        bFlag = true;
	}
	else
	{
		bFlag=false;
		std::cerr << "Unable to process file" << std::endl;
	}
    
	return bFlag;
    
}

unsigned int TimingToSequence::getTotalTrains()
{
    return (unsigned int) vTiming.size();
}

void TimingToSequence::displayTimingTrainStats()
{
    std::cout << "Number of trains detected : " << vTiming.size() << std::endl;
    
    for (size_t i=0; i < vTiming.size(); ++i)
    {
        std::cout << "Train number : " << i << " Size : " << vTiming[i].size() << std::endl;
        DisplayTrain(i);
    }
    std::cout << std::endl;
}

//Get the necessary sequences

void TimingToSequence::getSequences(unsigned int idxX, unsigned int idxY, std::vector<unsigned int>& seqX,std::vector<unsigned int>& seqY)
{
    
    std::vector<double> vTrainX, vTrainY;
    vTrainX = vTiming[idxX];
    vTrainY = vTiming[idxY];
    
    //Add delay to vTrainX; We then compute the TE from the delayed
    //version of train X to train Y. This is the second method
    if (uCalcMode == 2)
    {
        for (size_t idx=0; idx < vTrainX.size(); ++idx)
        {
            vTrainX[idx] = vTrainX[idx]+dDelay;
        }
    
    }
    
    ExtractSequences(vTrainX,vTrainY,seqX,seqY);
    //DisplaySequences(idxX,idxY,seqX,seqY);
}
void TimingToSequence::ExtractSequences(std::vector<double>& vTX, std::vector<double>& vTY,std::vector<unsigned int>& seqX,std::vector<unsigned int>& seqY)
{
    //First we need to find the total time interval and hence the number of bins
    //Thus we find the minimum starting time and the maximum ending time to cover the whole space
    //This will make the sequence lengths equal
    
    
    double minStartTime = (std::min<double>(vTX.front(),vTY.front()));
    double maxEndTime = (std::max<double>(vTX.back(),vTY.back()));
    unsigned long iSeqLength = (unsigned long) ((maxEndTime-minStartTime)/dBinSize) + 1;
    
    //Resize the vectors and fill with zeros
    seqX.resize(iSeqLength, 0);
    seqY.resize(iSeqLength, 0);
    
    //Form the sequences with zeros and ones
    
    for (size_t j=0; j < vTX.size();++j)
    {
        size_t binNo = (int) ((vTX[j]- minStartTime)/dBinSize);
        if (seqX[binNo] == 0) seqX[binNo] = 1;
    }
    
    for (size_t j=0; j < vTY.size();++j)
    {
        size_t binNo = (int) ((vTY[j]- minStartTime)/dBinSize);
        if (seqY[binNo] == 0) seqY[binNo] = 1;
    }
    
    //Make the delayed version of the X signal here
    // Following method 1
    if (uCalcMode == 1)
    {
        unsigned int dDiscreteDelay = ceil((dDelay / dBinSize));
        std::vector<unsigned int> vTemp(dDiscreteDelay,0);
        std::vector<unsigned int> vTemp1 = vTemp;
        //Construct a vector of zeros of size dDiscreteDelay
        // Append this to seqX from the beginning and to seqY at the end
        vTemp.insert(vTemp.end(),seqX.begin(),seqX.end());
        seqX = vTemp;
        seqY.insert(seqY.end(),vTemp1.begin(),vTemp1.end());
        
    }
    
}

void TimingToSequence::DisplayTrain(size_t idx)
{
    std::cout << "Individual spike times  : " << idx << std::endl;
    std::vector<double> vTemp = vTiming[idx];
    
    for (size_t idx=0; idx < vTemp.size(); ++idx)
    {
        std::cout << vTemp[idx] << "  ";
    }
    std::cout << std::endl;

}

void TimingToSequence::DisplaySequences(unsigned int idX,unsigned int idY,std::vector<unsigned int>& sX,std::vector<unsigned int>& sY)
{
    std::cout << std::endl << "Printing the sequence corresponding to train number : " << idX << std::endl;
    std::cout << "Sequence length : " << sX.size() << std::endl;
    for (size_t i=0; i < sX.size(); ++i)
    {
        std::cout << sX[i] << "  ";
    }
    std::cout << std::endl;
    
    std::cout << std::endl << "Printing the sequence corresponding to train number : " << idY << std::endl;
    std::cout << "Sequence length : " << sY.size() << std::endl;
    for (size_t i=0; i < sY.size(); ++i)
    {
        std::cout << sY[i] << "  ";
    }
    std::cout << std::endl;

}


