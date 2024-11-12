# Brief description

py-ncsm is is a simple NCSM code for solving the A=3 nucleon system written
in python. The code is not optimized or very fast, but istead written in a 
simple and understandable manner for ease of use. 

In the directory 'deuteron/' there is also a script to solve for the deuteron 
binding energy in a harmonic oscillator basis.

License: See LICENSE.txt file.

Written by: 
Oliver Thim,
Chalmers University of Technology
v1.0 Autumn 2024

## Important things to note about the code:

- This code is primarely written to consider the triton (^3H) in the limit of
  isospin symmetry. This is reminded through the flag:

  isospin_sym = True

  If isospin_sym is this is set to false, the code will terminate since that is 
  not implemented.

## External dependencies:

- FORTRAN code for computing Moshinsky brackets. (TODO)

- Code to compute Wigner symbols: wigxjpf (installable via pip)

# Installing and running the code:

- Create the conda environment using the environment.yml file and activate 
  the environment:

  $ conda env-create -f environment.yml
  $ conda activate py-ncsm-env

- Install the package wigxjpf, you can install it via pip.

- Install the package gmosh using f2py by going to the directory src/external
  and follow the instruction in that directorys README file.

- The code needs pre-calculated 3N basis states and the transformation matrix
  between the partially antisymmetrized and fully antisymmetrized 3N basis.
  To pre-calculate and store states for a given Nmax run the command: 

  $ src/states/compute_save_basis.py --Nmax=<Nmax_list>

  With e.g. <Nmax_list> =0,2,4.

  This generates a directory in src/states/Nmax=<Nmax>_data/ that contains the 
  pre-calculated data. Note that the storage space needed might not negligable
  for Nmax>~40. NOTE that you need to set (J,T,pi) values in the file
  compute_save_basis.py. The default is the triton channel (1,1,1).

- To run the code you just need to specify an input file and run the code.
  
  If an input file is NOT provided the variables are set in the default section 
  of the code and the output will by default be printed to screen. This mode
  can be useful e.g. for testing. E.g. 

  $ python run-py-ncsm.py

  An input file is provided as by running:

  $ python run-py-nscm.py example_input_file_1.txt

  ----------------------------- input_file.txt --------------------------------
  [settings]
  nmax_arr         = [0,2,4]
  hbar_omega       = 24    # hbar_omega in MeV.
  isospin_sym      = True
  interaction_file = interactions/idaho_n3lo_nmax_40_hw_24.txt
  output_file      = example-inout-files/out_small.txt
  ------------------------- END input_file.txt --------------------------------


  This will generate the corresponding output file <outdir>/<out_file_name>

- The output file has the following general structure:

  ----------------------------- output_file.txt -------------------------------
  <version_of_this_code>
  infile=<input_file_name>
  <copy_of_input_parameters>
  <run_information_from_the_code>
  
  -------------------------------------------------------------------------------
                                  *** OUTPUT ***
  -------------------------------------------------------------------------------
  Nmax     E 	 	  dim (Gamma) dim (alpha)
  0        12.16523   1        	  2       
  2        4.88032    5        	  14      
  4        1.82073    15       	  44      

  ------------------------- END output_file.txt -------------------------------

- To benchmark the code, you can verify that the output of the code is the same
  and the corresponing example input/output files in example-inout-files/

- Have fun!
