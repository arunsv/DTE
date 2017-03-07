
#ifndef TEENTROPY_HPP_
#define TEENTROPY_HPP_

#include <vector>
#include <map>

class TEntropy {

public:
	/* Constructor and Destructor*/
	TEntropy();
	~TEntropy();

	/*Sequence Related */
	void LoadSequences(std::vector <unsigned int>& ,std::vector <unsigned int>& );
	void ClearSequences();

	/* Get/Set methods */
	void SetWindowSize(unsigned int winSize);
	unsigned int GetWindowSize();
	double getTransferEntropy();

    /*Computation related */
	void computeTransferEntropy();
	void caculateMutualInformation( );
    
    //Static member variable for tracking zeroprobabilities
	static const int ZeroProbability = 10;


private:

	// Private member variables

	std::vector<unsigned int> seqX;
	std::vector<unsigned int> seqY;

	//Sequence vectors
	//Transfer Entropy Calculations
	std::vector<unsigned int> vXh,vYh,vXhXn1,vYhYn1,vXhYh,vYhXh,vXhYhYn1,vYhXhXn1;
	//Entropy and Mutual information calculations
	std::vector<unsigned int>vX,vY,vXY,vYX;

	//Probability maps
	typedef std::map<std::vector<unsigned int>,double > pMap;
	typedef std::map<std::vector<unsigned int>,double >::iterator pMapIt;
	//Transfer Entropy Calculations
	pMap pXh,pXhXn1,pYh,pYhYn1,pXhYh,pYhXh,pXhYhYn1,pYhXhXn1;
	//Entropy and Mutual information calculations
	pMap pX,pY,pXY,pYX;

	unsigned int wS; //Window size
	double tranEntropy;

	// Private member functions
	void setupSamples(std::vector <unsigned int>& ,std::vector <unsigned int>& );

	void updateAllProbMaps(bool bPartial,bool bBlock);
	void updateProbMap(std::vector<unsigned int>&, pMap&);
	void clearAllSubSeqVectors();
	void calculateAllJointProbabilities();
	void calculateJointProbability(pMap&);
	double caculateTransferEntropyForSeqPair(pMap&, pMap&, pMap&, pMap&);
	double computeOneShannonEntropy(pMap&);
	double getProbSubSeq(pMap&,std::vector<unsigned int>&);
    

    //Debug related
	void printSeqProbability(pMap&,std::string);
	void printVector(std::vector<unsigned int>& , std::string );
};


#endif /* TEENTROPY_HPP_ */
