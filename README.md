# DTE
Compute the delayed transfer entropy (DTE); recovery of delay profile between two signals via computation of DTE

Description of scripts, prerequisites and  how to run the scripts


The DTE package comprises a core C++ program and a Python driver script. The C++ code is a general transfer entropy computation framework that reads in  time series data from a timing file and computes and outputs  transfer entropy values between all possible pairs. Even though the code can handle an arbitrary number of timing trains, our focus will be between two causally related signals with various delay models as driven by the Python script. Note that the delay profiles can be arbitrary; we use well known models for verification. 

The Python driver script helps to generate two causally relate time series according to known forward delay models. For illustration purposes, three different delay models are considered namely a constant delay, a uniform delay distribution and a Gaussian delay distribution. Once the two time series are generated, the delayed transfer entropy is computed by sweeping the delay parameter in equation 3. The resulting profile is seen to provide evidence of the forward delay model which is unseen by the transfer entropy engine. The technical details are explained in the accompanying paper [1].

Description of the C++ engine:


To build the program just run make within the “src” sub-directory. The program supports exactly five arguments and the usage is as follows:
./DTE workingDir timingFile binResolution delay TE_Calculation_Mode

DTE is the executable
workingDir is the working directory where the timing file can be found and the outputs are written.
timingFile is a simple text file with the timing information. Each line corresponds to one node and the timing of each nodal activity is given as times of tweet/re-tweet separated by a comma.
binResolution is the size of the bin used for discretization in the same units as the timing in the timing file.
delay is the swept-delay in the DTE profile computation 
TE_Calculation_Mode refers to the manner in which DTE computation is carried out and can take on the values 1 or 2. Current default is 1 while 2 refers to an “experimental” mode that is still under development.


Description of the Python driver script

The driver python script needs numpy and matplotlib.pyplot libraries to work. Two classes are defined within the script to make data movement easy. The first one is dParams and it controls the delay model. The second is the TEParams that stores the parameters related to process generation and DTE computation. 
	
  The script first generates a Poisson process X according to a given rate. Then using a pre-defined delay model, it generates a process Y causally dependent on X according to the defined delay model, Under the “Constant” model, each spike in X undergoes the same delay to become process Y. Under the “Uniform” model, each spike in X undergoes a delay drawn from a pre-defined uniform distribution to become process Y. Finally, under the “Gaussian” model, each spike in X undergoes a delay drawn from a pre-defined Gaussian distribution to become process Y. Then using only these two timing trains, the swept DTE profile is calculated by sweeping the deay as a parameter.  The shape of this profile will point to the actual delay models. Further details are available in the accompanying paper [1]. 

Three sample plots obtained from the script are included in the “samplePlots” sub-directory for illustration purposes. In all three cases process Y was causally generated from process X with different delay models (Constant, Uniform and Gauss). As seen from the plots, the swept DTE profile clearly points to the underlying forward model. Note that the reverse DTE (from Y to X) is always very small as it should be. Note that the plots are based on a single realization of the processes involved. A smoother behavior can be expected by running a Monte Carlo analysis with say 30-50 runs and averaging. 

Located in the /Python sub-directory, the calling syntax of the script is 
python DTEDriver.py dirName progMode delayMode

Examples : 
python DTEDriver.py testUnif save Uniform
python DTEDriver.py testUnif load

dirName is the working directory name and a new directory will be created under ../workingDir. If it already exists, the user will be prompted to choose a different name.
progMode is the program running mode. The choices are save and load. Under save, the program does the computation, saves the swept DTE data in the working directory in the numpy compressed format and plots the results. Under the load option, the program does not do any computation, instead it loads a previously saved data in the numpy compressed format and does the plotting.
delayMode is the delay distribution to use. Presently the choices are constant, uniform and Gaussian. This argument is not needed when the progMode argument is load.
