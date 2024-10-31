# Brief description

py-ncsm is is a simple NCSM code for solving the A=2,3 nucleon system written
in python. The code is not optimized or very fast, but istead written in a 
simple and understandable manner for ease of use.

<License>

Written by: 
Oliver Thim,
Chalmers University of Technology
v1.0 Autum 2024

## Important things to note about the code:

- This code is primarely written to consider the triton (^3H) in the limit of
  isospin symmetry.

  ISOSPIN_SYM = True

  If this is this is set to false, the code will curretly crash since that is 
  not implemented.

## External dependencies:

- FORTRAN code for computing Moshinsky brackets. (TODO)

- Code to compute Wigner symbols: wigxjpf (installable via pip)

# Installing and running the code:

- Create the conda environment using the environment.yml file and activate 
  the environment.

- Install the package wigxjpf. You can install it via pip.

- Install the package gmosh using f2py by going to the directory src/external
  and follow the instruction in that directorys README file.

- The code needs pre-calculated 3N basis states and the transformation matrix
  between the partially antisymmetrized and fully antisymmetrized 3N basis.
  To pre-calculate and store states for a given Nmax run the command: 

  $ src/states/compute_save_basis.py -Nmax <Nmax_list>

  With e.g. <Nmax_list> = [0,2,4].

  This generates a directory in src/states/Nmax=<Nmax>_data/ that contains the 
  pre-calculated data. Note that the storage space needed might not negligable
  for Nmax>~40.

- To run the code you just need to specify an input file and run the code.
  
  If an input file is NOT provided the variables are set in the default section 
  of the code and the output will by default be printed to screen. This mode
  can be useful e.g. for testing. E.g. 

  $ python run-py-ncsm.py

  An input file is provided as by running:

  $ python run-py-nscm.py < example_input_file_1.txt

  ------------------------ example_input_file_1.txt ---------------------------


  -----------------------------------------------------------------------------


  This will generate the corresponding output file <outdir>/<out_file_name>

- The output file has the following general structure:

  --------------------------- output_file.txt ---------------------------------
  <version_of_this_code>
  infile=<input_file_name>
  <copy_of_input_file>

  <run_information_from_the_code>
  
  *** OUTPUT ***
  -----------------------------------------------------------------------------
  <output_data>
  -----------------------------------------------------------------------------

- To benchmark the code, you can verify that the output of the code is the same
  and the corresponing example input/output files in example-inout-files/

- Have fun!
