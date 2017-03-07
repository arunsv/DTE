
#ifndef TEDRIVER_HPP_
#define TEDRIVER_HPP_

#include "TimingToSequence.hpp"
#include "TEntropy.hpp"

class TEDriver{
    
public:
    TEDriver(TimingToSequence*);
    ~TEDriver();
    
    void setOutputFile(std::string);
    void computeTransEntropyMain();
    void writeInfo();
    
private:
    
    std::string sTEFile;
    TimingToSequence* theTS;
    std::vector<unsigned int> seqX, seqY;
    std::map<std::pair<unsigned int,unsigned int>,double> mTEMap;
    
};


#endif /* TEDRIVER_HPP_ */
