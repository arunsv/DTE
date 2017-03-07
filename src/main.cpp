//
//  main.cpp
//  TimingToTransferEntropy
//
//  Created by Arun Sathanur on 5/10/14.
//  Copyright (c) 2014 ACE Lab, UWEE. All rights reserved.
//

/* Standard includes*/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

/*Specific includes*/
#include "TimingToSequence.hpp"
#include "TEDriver.hpp"


int main(int argc, char** argv)
{
    
	if (argc == 6) // Exactly five arguments allowed
	{
		//Read in the input and output files
		std::string sWorkingDir(argv[1]);
        std::string sTimingFile(argv[2]);
        double dBinSize=atof(argv[3]);
        double dDelay = atof(argv[4]);
        unsigned short int uCalcMode = atoi(argv[5]);
        
        std::string sTEFile = sWorkingDir + std::string("/TE.txt");
        
        //Read in the spike trains from the file
    
		TimingToSequence* myTS = new TimingToSequence;
        myTS->setFileName(sTimingFile);
        myTS->setBinSize(dBinSize);
        myTS->setDelay(dDelay);
        myTS->setCalcMode(uCalcMode);
        
		bool bFlag = myTS->getTiming();

        if (bFlag)
        {
            TEDriver* myTEDriver = new TEDriver(myTS);
            
            myTEDriver->setOutputFile(sTEFile);
            myTEDriver->computeTransEntropyMain();
            myTEDriver->writeInfo();
            
            delete myTEDriver;
        }
        
        delete myTS;
	}
	else
	{
		std::cerr << "The program supports exactly five arguments. " << std::endl;
		std::cerr << "Usage : " << std::endl;
   		std::cerr << "./DTE workingDir timingFile binResolution delay TE_Calculation_Mode" << std::endl;
	}
    

	return 0;
    
}

