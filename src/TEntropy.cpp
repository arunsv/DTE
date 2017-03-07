

#include <iostream>
#include <cmath>
#include <sstream>
#include "TEntropy.hpp"

TEntropy::TEntropy()
{
    wS=1;
    tranEntropy = 0.0;
}

TEntropy::~TEntropy()
{

}

void TEntropy::LoadSequences(std::vector <unsigned int>& series1,std::vector <unsigned int>& series2)
{
	//Error condition

	seqX = series1;
	seqY = series2;
}

void TEntropy::ClearSequences()
{
	seqX.clear();
	seqY.clear();
}


void TEntropy::SetWindowSize(unsigned int winSize)
{
	wS = winSize;
}

unsigned int TEntropy::GetWindowSize()
{
	return wS;
}


void TEntropy::computeTransferEntropy()
{
	setupSamples(seqX,seqY);
	calculateAllJointProbabilities();

	//Actual calculation of  Transfer Entropy
	try {		tranEntropy = caculateTransferEntropyForSeqPair(pXhYhYn1, pXhYh, pYhYn1, pYh);		}
	catch (int i)  {		if (i ==ZeroProbability )		throw ZeroProbability	;		}

}

void TEntropy::setupSamples(std::vector <unsigned int>& seqX,std::vector <unsigned int>& seqY)
{

	//First extract the various vectors

	unsigned long sL=seqX.size();

	for (size_t i=0;i < sL; ++i)
	{

		if (i < (sL-wS)) // Transfer entropy stuff
		{
			//Initialize the vXh and VYh vectors.
			vXh.assign(wS,0);
			vYh.assign(wS,0);

			for (unsigned int j=0; j < wS; ++j)
			{
				vXh[j]= seqX[i+j];
				vYh[j]= seqY[i+j];

			} // Index j

			//Do all the other insertions

			vXhYh=vXh; vXhYh.insert(vXhYh.end(),vYh.begin(),vYh.end());
			vYhXh=vYh;	vYhXh.insert(vYhXh.end(),vXh.begin(),vXh.end());
			vXhXn1=vXh;	vXhXn1.push_back(seqX[i+wS]);
			vYhYn1=vYh;	vYhYn1.push_back(seqY[i+wS]);
			vXhYhYn1=vXhYh; vXhYhYn1.push_back(seqY[i+wS]);
			vYhXhXn1=vYhXh; vYhXhXn1.push_back(seqX[i+wS]);

			//Now update the frequencies
			updateAllProbMaps(false,true); //Update all maps

		}

		if (i < sL) // Entropy and mutual information stuff
		{
			//Initialize the vXh and VYh vectors.
			vX.assign(1,0);
			vY.assign(1,0);

			//vX and vY can hold just 1 symbol as opposed to the others
			vX[0]=seqX[i];
			vY[0]=seqY[i];
			vXY=vX;
			vXY.insert(vXY.end(),vY.begin(),vY.end());
			vYX=vY;
			vYX.insert(vYX.end(),vX.begin(),vX.end());

			//Now update the frequencies
			updateAllProbMaps(false,false); //Update all symbol maps
		}

		clearAllSubSeqVectors();

	} // Index i

	//Finally update only the Xh,Yh,XhYh,YhXh vectors with the last sub-sequence
	//Remember that previously we neglected the last sub-sequence by looping only till i < (sL-wS)
	//This was because we needed the last one for the sequences with (n+1)th bit and we stopped there for all sub-sequence vectors
	//Now update the frequencies

	size_t i= (sL-wS);
	//Initialize the vX and VY vectors.
	vXh.assign(wS,0);vYh.assign(wS,0);

	for (unsigned int j=0; j < wS; ++j)
	{
		vXh[j]= seqX[i+j];
		vYh[j]= seqY[i+j];

	} // Index j
	vXhYh=vXh; vXhYh.insert(vXhYh.end(),vYh.begin(),vYh.end());
	vYhXh=vYh;	vYhXh.insert(vYhXh.end(),vXh.begin(),vXh.end());
	updateAllProbMaps(true,true); // Update only the partial maps
	clearAllSubSeqVectors();

} // End of function

void TEntropy::updateAllProbMaps(bool bPartial, bool bBlock)
{
	if (bBlock) // Update block probabilities
	{
		updateProbMap(vXh,pXh);
		updateProbMap(vYh,pYh);
		updateProbMap(vXhYh,pXhYh);
		updateProbMap(vYhXh,pYhXh);

		if (!bPartial) //If bPartial is false then update all block probabilities
		{
			updateProbMap(vXhXn1,pXhXn1);
			updateProbMap(vYhYn1,pYhYn1);
			updateProbMap(vXhYhYn1,pXhYhYn1);
			updateProbMap(vYhXhXn1,pYhXhXn1);
		}
	}
	else //Update symbol probabilities
	{
		updateProbMap(vX,pX);
		updateProbMap(vY,pY);
		updateProbMap(vXY,pXY);
		updateProbMap(vYX,pYX);
	}

}

void TEntropy::updateProbMap(std::vector<unsigned int>& v, pMap& p)
{

	if (p.find(v) == p.end())
	{
		p[v]=1;
	}
	else
	{
		p[v]=p[v]+1;
	}
}

void TEntropy::clearAllSubSeqVectors()
{
	vXh.clear();
	vYh.clear();
	vXhXn1.clear();
	vYhYn1.clear();
	vXhYh.clear();
	vYhXh.clear();
	vXhYhYn1.clear();
	vYhXhXn1.clear();

	vX.clear();
	vY.clear();
	vXY.clear();
	vYX.clear();

}

void TEntropy::calculateAllJointProbabilities()
{
	calculateJointProbability(pX);
	calculateJointProbability(pY);
	calculateJointProbability(pXY);
	calculateJointProbability(pYX);

	calculateJointProbability(pXh);
	calculateJointProbability(pYh);
	calculateJointProbability(pXhYh);
	calculateJointProbability(pYhXh);

	calculateJointProbability(pXhXn1);
	calculateJointProbability(pYhYn1);
	calculateJointProbability(pXhYhYn1);
	calculateJointProbability(pYhXhXn1);

}

void TEntropy::calculateJointProbability(pMap& p)
{

	pMapIt it1;
	double count=0.0;

	for (it1 = p.begin(); it1 != p.end(); ++it1)
	{
		count = count+it1->second;
	}

	for (it1 = p.begin(); it1 != p.end(); ++it1)
	{
		it1->second = it1->second/count;
	}
}

//Compute the transfer entropy from the probability  maps
double TEntropy::caculateTransferEntropyForSeqPair(pMap& Pxhyhyn1, pMap& Pxhyh, pMap& Pyhyn1, pMap&  Pyh)
{
	double Probxhyhyn1,Probxhyh,Probyh,Probyhyn1,Txy=0.0;
	pMapIt itrM,itrM1;

	//Loop over the full symbol
	for (itrM=Pxhyhyn1.begin(); itrM!=Pxhyhyn1.end(); ++itrM )
	{
		std::vector<unsigned int> Vxhyhyn1=itrM->first;
		Probxhyhyn1=itrM->second;

		//Extract the sub-sequence symbols of interest and then read-off their probabilities

		std::vector<unsigned int>::const_iterator itrV = Vxhyhyn1.begin();
		std::vector<unsigned int>Vxhyh = std::vector<unsigned int>(itrV,itrV+2*wS);
		std::vector<unsigned int>Vyh = std::vector<unsigned int>(itrV+wS,itrV+2*wS);
		std::vector<unsigned int>Vyhyn1 = std::vector<unsigned int>(itrV+wS,itrV+2*wS+1);

		try {

			Probxhyh = getProbSubSeq(Pxhyh,Vxhyh);
			Probyh = getProbSubSeq(Pyh,Vyh);
			Probyhyn1 = getProbSubSeq(Pyhyn1,Vyhyn1);

			double txy=0.0; // Local contribution to transfer entropy

			txy =   log(Probxhyhyn1)/log(2.00);
			txy -=  log(Probxhyh)/log(2.00);
			txy -=  log(Probyhyn1)/log(2.00);
			txy += log(Probyh)/log(2.00);

			txy *= Probxhyhyn1;
			Txy += txy; // Adding up the local contributions of the transfer entropy
		}
		catch (int i)
		{
			if (i ==ZeroProbability )	throw ZeroProbability	;
		}

	}
	return Txy;
}


double TEntropy::getProbSubSeq(pMap& p,std::vector<unsigned int>& v)
{
	pMapIt itr = p.find(v);
	double prob=0.0;

	if (itr != p.end())
	{
			prob = itr->second;
	}
	else
	{
		prob = 0.0;
		throw ZeroProbability;
	}

	return prob;

}


double TEntropy::computeOneShannonEntropy(pMap& p)
{
	pMapIt itrM=p.begin();
	double ent=0.0,prob = 0.0;
	//Loop over the symbol set
	for (; itrM!=p.end(); ++itrM )
	{
		prob=itrM->second;
		//Here there is no danger of encountering a zero probability since we don't do symbol splitting
		// with the possibility of a sub-symbol not existing
		ent += prob*((log(prob)/log(2.00)));
	}

	return -ent;
}

double TEntropy::getTransferEntropy()
{
	return tranEntropy;
}

//Debug related

void TEntropy::printSeqProbability(pMap& p,std::string mName)
{

	pMapIt it1;

	std::cout << "Printing sequence and probability. for map : " << mName << std::endl;
	std::stringstream ss;

	for (it1 = p.begin(); it1 != p.end(); ++it1)
	{
		std::vector<unsigned int> v=it1->first;

		for (unsigned int i=0; i < v.size(); ++i)
		{
			ss << v[i] << " ";
		}
		ss << "Probability : " << it1->second << "\n";
	}

	std::cout << ss.str() << std::endl;

}

void TEntropy::printVector(std::vector<unsigned int>& v, std::string vName)
{
	std::cout << "Printing vector : " << vName << std::endl;
	for (unsigned int i=0; i< v.size(); ++i)
	{
			std::cout << v[i] << " ";
	}

	std::cout << std::endl;

}
