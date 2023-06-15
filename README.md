# Experimental VAMOS++ AGATA SPIDER Implementation On Nuclear Spectroscopy

User guide of EVASIONS code

Author: C. Fougères (2023)
Contact: cfougeres@anl.gov

References
 - *Understanding the cosmic abundance of 22Na: lifetime measurements in 23Mg* (C. Fougères, Thesis, Normandie University, https://hal.science/tel-03768524v1)
	--> chapter 4
 - *Search for 22Na in novae supported by a novel method for measuring femtosecond nuclear lifetimes* (C. Fougères et al., arXiv preprint, https://doi.org/10.48550/arXiv.2212.06302)

/#################/

Requirements
1. linux based operating system or MacOS
2.  CERN/ROOT version 6  (see https://root.cern.ch/).
3. `g++` compiler 
4. `git` 
5. `make`

(Optional PDF reader)

/#################/


Installation
1. Get code at:
> git clone https://github.com/CFougeres/EVASIONScode
2. In the main directory 'EVASIONcode/', compile:
> make

/#################/

Run 
> ./EVASIONSsim +'#5 inputs''

	input_1: if codeA for for gamma emission used ==1 (else ==0)
	input_2: if codeB for for particle emission used ==1 (else ==0)
	input_3: if width determination ==1 (else ==0)
	input_4: if saving TTrees ==1 (else ==0)
	input_5: if live plots ==1 (else ==0)			--> once opened, click on TCanvas to launch the plot app.
							  	    to stop Enter Ctrl+C in terminal

2 code versions available (see doc/MCalgo.pdf)

- 'CodeA gamma' == {analysis of gamma-rays -- recoil events, velocities profiles, lifetime and total width measurements ...}  
- 'CodeB_gparticle' == {analysis of unbound_particle -- recoil events, excitation functions, partial width measurements ...}  

/#################/

USER INPUTS
- Experiment and physics case general inputs

	Location: 'inputs.dat'

	Rem: parameters names are explicit

- Decay data to be analyzed (based on AGATA-VAMOS-SPIDER experiment)

        Location "data/"

        Formalism TChain with "TreeMaster", 
                
                Branches:
                        
                        Brho, Theta, Phi, ThetaLdeg (from VAMOS)
                        
                        CodeA  = trackE, trackX1, trackY1, trackZ1  (from AGATA)
                        
                        CodeB = "particle Si_dE, Si_Eres, Si_theta"  (from SPIDER)	    

        Example of extraction root script: "src/functions/Backup/extractionData.C"
		
- Stopping power tables from SRIM [J. F. Ziegler et al., SRIM Co., United States of America 6th (2013)] with S.P. in keV.micron

        Location "src/functions/SP_SRIM/"
	
        5 files in SRIM .txt formalism (see "examples/")
                
                (1) 'beam_in_matter.dat'
	        
                (2) 'recoil_in_matter.dat'
		
                (3) 'ejectil_in_matter.dat'
		
                (4) 'particle_in_matter.dat'
		
                (5) 'particle_in_silicon.dat'


/#################/

Saved data

Location "tree_results/"


/#################/

Source codes

- Main	

	Location "src/"

- Auxialliaries	

	Location "src/functions/"
       
...


/#################/

Compiled version

	Location "lib/"
        
