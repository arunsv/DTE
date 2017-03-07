
#ifndef TIMINGTOSEQUENCE_HPP_
#define TIMINGTOSEQUENCE_HPP_

#include <string>
#include <vector>

class TimingToSequence{
    
public:
    
	TimingToSequence();
	~TimingToSequence();
    
	void setFileName(std::string );
    void setBinSize(double);
    void setDelay(double);
    void setCalcMode(unsigned short int);
    unsigned int getTotalTrains();
	bool getTiming();
    void getSequences(unsigned int, unsigned int, std::vector<unsigned int>&,std::vector<unsigned int>&);
    
private:
    
    //Private members
    double dBinSize;
    unsigned short int uCalcMode;
    double dDelay;
    
    std::string timingFileName;
    std::vector<std::vector<double> > vTiming;
    
    //Private methods
    bool readActivityTimingFromFile();
    void displayTimingTrainStats();
    void ExtractSequences(std::vector<double>&,std::vector<double>&,std::vector<unsigned int>&,std::vector<unsigned int>&);
    
    //Mainly for debug purposes
    void DisplaySequences(unsigned int,unsigned int,std::vector<unsigned int>&,std::vector<unsigned int>&);
    void DisplayTrain(size_t);

};


#endif /* TIMINGTOSEQUENCE_HPP_ */
