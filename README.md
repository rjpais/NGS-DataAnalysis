# Pipeline Developer versions for MiniON data analysis using MEDAKA METHOD 

## In this branch repository contains multiple versions of pipelines developed for bioinformatic analysis of NGS data from oxoford nanopore thecnologies. Pipelines are alpha versions that were developed in Python 3 and tested on a ubunto os. Parameterization is only possible withing variable changing inside each script. Documentation is provided also provided in each script. Here, we include pipelines that where fully personalized versions of medaka and and systematic adaptation of ARTIC.      


## Pre-instalation software requisites for most scripts:
* anaconda python 3.7 distribution or later 
* mathplotlib
* biopython 
* anaconda python distribution 
* Medaka 1.2.1 (anaconda forge installation)
* Nanofilt 2.6.0  
* NanoStats 1.4.0 
* jupyter notebook  

## Running instructions:

1) On the command line using anaconda environment, activate medaka environment by typing  


2) Activate the medaka/artic environment on anaconda depending on the script your using.
 
     (base) $ conda activate medaka
    
          or 

     (base) $ conda activate artic
     

3) On the command line, go to the path where the script is located and run the following command: 

      (medaka) $ python SCRIPT_NAME.py
      
       or
      
      (artic) $ python SCRIPT_NAME.py  


Some input variables may be asked on the command line console and file choose GUI will pop up. 
If all works, you will see a folders for each sample and outputs files within will be generated systematically. 

