#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include "TEDriver.hpp"

TEDriver::TEDriver(TimingToSequence* TS)
{
    theTS = TS;
}

TEDriver::~TEDriver()
{

}

void TEDriver::setOutputFile(std::string sF1)
{
    sTEFile = sF1;
    
}
void TEDriver::computeTransEntropyMain()
{
    
    unsigned int iTotTrains = theTS->getTotalTrains();
    
    for (unsigned int idxX=0; idxX < iTotTrains; ++idxX)
        for (unsigned int idxY=0; idxY < iTotTrains; ++idxY)
        {
            //Skip iteration if source and reciever are the same process
            if (idxX == idxY) continue;
            
            //Get the sequences corresponding to train numbers idxX and idxY
            theTS->getSequences(idxX, idxY,seqX,seqY);
            
            //Calculate transfer entropy and store
            
            TEntropy* myTE = new TEntropy;
            myTE->LoadSequences(seqX,seqY);
            
            try
            {
                //Compute the transfer entropy
                
                myTE->computeTransferEntropy();
                
                //Get transfer entropy and store
                double TExy=0.0;
                TExy = myTE->getTransferEntropy();
                std::pair<unsigned int,unsigned int> vPair = std::make_pair(idxX,idxY);
                
                mTEMap.insert(std::make_pair(vPair,TExy));
                
            }
            catch (int i)
            {
                if ( i == TEntropy::ZeroProbability)
                {
                    std::cerr << std::endl << " Encountered a zero probability sub-sequence. Terminating program. " << std::endl ;
                }
            }
            
            delete myTE;
        }
        
}

void TEDriver::writeInfo()
{
	std::stringstream ss;
    
    std::ofstream ofs;
    ofs.open (sTEFile.c_str(), std::ofstream::out);
    if (ofs.is_open())
    {
        ss << "# Delayed Transfer Entropy values in bits" << std::endl;
        std::map<std::pair<unsigned int,unsigned int>,double>::const_iterator itrTE = mTEMap.begin();
        
        for (; itrTE != mTEMap.end(); ++itrTE)
        {
            ss << itrTE->first.first << " " << itrTE->first.second << " " << itrTE->second << std::endl;
        }
        ss << std::endl;
        ofs << ss.rdbuf();
        
    }
    else std::cout << "Unable to write the transfer entropy output file." << std::endl;
    ofs.close();
    
}
