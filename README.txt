QuickPIC scripts introduction
Erik Adli, Dec 13 2011

OVERVIEW :
**********

The main steps for running a QPIC job are :
  I) Prepare the input file  (rpinput)
 II) Run QPIC with the prepared input file on UCLA cluster
III) Analyze the QPIC output

A set of Matlab scripts and functions have been written to aid in steps I) and III).  To get started the relevant Matlab scripts to run are :
	quickpicsim_param2qp.m   	for rpinput generation
	quickpicsim_ana.m   		for analysis of output




DETAILS :
**********


  I) Prepare the input file
---------------------------

All the input for the plasma PIC simulation performed by QPIC is provided by a single file named "rpinput".  The rpinput contains detailed information about :
* beam parameters
* plasma parameters
* output options 
* simulation technical details

The parameters are often interlinked so care must be taken when modifying an rpinput.   A function (my_gen_rpinput.m) has been written to create a consistent rpinput from physical beam and plasma parameters. The function reads a template rpinput file, changes the relevant parameters and writes an updated rpinput.  The script :
	quickpicsim_param2qp.m
contains three working examples on how to use this function, including 1) generate a single gaussian beam, 2) generate two gaussian beams and 3) generate a single beam with an arbitrary z-distribution

Notes:
1) The rpinput (see for example the "rpinput_template") contains quite thorough documentation of all the input options.  This file is the best documentation of the input option available (at least for the author of this readme).




 II) Run QPIC with the prepared input file on UCLA cluster
----------------------------------------------------------


Example commands to login and use and to the UCLA cluster :

login :
	ssh -Y eadli@hoffman2.idre.ucla.edu
go to your scratch folder:
        cd $SCRATCH
make a simulation folder :
	mkdir testrun
get the last version of qpic, for example from :
	cp ~eadli/QPIC_EXE/qpic.e /u/scratch/e/eadli/testrun/qpic.e
From your the computer where you generated your rpinput, copy your rpinput to the cluster :
	scp rpinput eadli@hoffman2.idre.ucla.edu:/u/scratch/e/eadli/testrun/


Example of how to create a script for QPIC job submission :

Run the executable :
	mpi.q
Enter the following commands (described as you go along) :
	B
	qpic.e	
	1024
	8
	y
	msa
	8
	<return>
	n
	Q

Submit a job :
	qsub qpic.e.cmd	

Check job status :
	qstat -u eadli

Copy suggested tar script :
	cp ~eadli/tar_qp.sh .

Tar output results :
	sh ~/tar_qp.sh testrun

On your local machine :
	scp eadli@hoffman2.idre.ucla.edu:scratch/testrun.tgz . )
Untar the output folder :
	tar xvzg testrun.tgz




III) Analyze the QPIC output
----------------------------


The output of QPIC is in HDF-format, and is split into a number of fields for the fields, beam/plasma densities and 6D beam phase space.   A script has been written to automatically parse the QPIC output, convert the data into SI-units and store it in a defined data structure.  The script name is :
	quickpicsim_ana.m

The first section contains user parameters.  In particular.  The quickpic output folder name must be defined in the script :
	datadir = '~/mydatafolder/testrun/'

The script runs through each saved timestep and stores the data in a structure called "qp".  Examples of how to retrieve data from the data structure :
a) Plot the stored field slice FEZ after the first saved time step :
	pcolor(qp(1).FEZ')
b) Get the stored sigma_x for timestep 2, beam 1
	qp(2).PP(1).sigma_x
c) Get the 6D phase space for timestep 1, bean :
	qp(1).PP(1).BEAM
The 6D phase space is stored in the following format :
Col #	
1	x  [um]
2	y  [um]
3	z  [um]
4	x' [urad]
5	y' [urad]
6	E  [GeV]
	

Notes:
1) The supported stored physical quantities for the moment are QEB, QEP, FEZ and beam phase space PHA-BEAM.  Other quantities can be easily added if needed.
2) The script contains several plot routines which illustrates how to use the data.  The plotting is separated from the data processing, and is activated by the flag "do_plot = 1"









