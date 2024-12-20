# Brief description

py-ncsm is a simple NCSM code for solving the A=3 nucleon system written
in python. The code is not optimized, or very flexible --- but written in an 
understandable manner. This code is the result of the author reading 
a PhD course in Computational Nuclear Physics at Chalmers 2024.

The code is currently set up only to solve the ^3H system using an effective 
isoscalar approximation. Note that the code is used at your own risk and there
can be bugs in the code. If you find a bug, or 
have suggestions for improvement please contact the author, or make a github 
issue.

In the directory 'deuteron/' there is also a script to solve for the deuteron 
binding energy in a harmonic oscillator basis.

License: See LICENSE.txt file.

Written by: 
Oliver Thim,
oliver.thim@chalmers.se
Chalmers University of Technology
v1.0 Autumn 2024

## Important things to note about the code:

- This code is primarily written to consider the triton (^3H). If 

  isospin_sym = False

  the function get_two_body_HO_potential_el(...) in setup_Hamiltonian.py
  will average the np and nn potentials as it should be treated in the
  triton. If

  isospin_sym = True

  the function get_two_body_HO_potential_el(...) will just take the np 
  potential element.


## External dependencies:

- FORTRAN code for computing Moshinsky brackets.

- Code to compute Wigner symbols: pywigxjpf (installable via pip)

# Installing and running the code:

- Create the conda environment using the environment.yml file and activate 
  the environment:

  $ conda env create -f environment.yml
  $ conda activate py-ncsm-env

- Install the package wigxjpf, you can install it via pip.

- Install the package gmosh using f2py by going to the directory src/external
  and follow the instructions in the README file in that directory.

- The code needs pre-calculated 3N basis states and the transformation matrix
  between the partially antisymmetrized and fully antisymmetrized 3N basis.
  To pre-calculate and store states for a given Nmax run the command: 

  $ python src/states/compute_save_basis.py --NMAX=<Nmax_list>

  With e.g. <Nmax_list> =0,2,4.

  This generates a directory in src/states/Nmax=<Nmax>_data/ that contains the 
  pre-calculated data. Note that the storage space needed might not be 
  negligable for Nmax>~40. NOTE that you need to set (2*J,2*T,pi) values in the file
  compute_save_basis.py. The default is the triton channel (1,1,1).

- To run the code you just need to specify an input file.
  
  If an input file is NOT provided the variables are set in the default section 
  of the code and the output will by default be printed to stdout. This mode
  can be useful e.g. for testing:

  $ python run-py-ncsm.py

  An input file is provided by running:

  $ python run-py-nscm.py example_inout_files/in_small.txt

  ------------------------------- in_small.txt --------------------------------
  [settings]
  nmax_arr         = [0,2,4,10,20]
  hbar_omega       = 24    # hbar_omega in MeV.
  isospin_sym      = False
  fast_comp        = False # Not implemented
  dim_lanczos      = 2000
  interaction_file = interactions/idaho_n3lo_nmax_40_hw_24_Np_80_finite.txt
  output_file      = example-inout-files/out_small.txt
  --------------------------- END in_small.txt --------------------------------


  This will generate the corresponding output file <outdir>/<out_file_name>.

- The output file has the following general structure:

  ----------------------------- out_small.txt -------------------------------
  <version_of_this_code>
  infile=<input_file_name>
  <copy_of_input_parameters>
  <run_information_from_the_code>
  
  -----------------------------------------------------------------------------
                                  *** OUTPUT ***
  -----------------------------------------------------------------------------
  Nmax  E                    dim (Gamma)  dim (alpha)  Comp. H (ms)  Diag. H (ms)  method
  0     12.1652310190506654  1            2            0             0             exact   
  2     4.8803190292006242   5            14           1             0             exact   
  4     1.8207328044041804   15           44           3             0             exact   
  10    -5.0530076036149802  108          322          135           1             exact   
  20    -7.7727667430020322  632          1892         4385          26            exact   

  --------------------------- END out_small.txt -------------------------------
- The format needed in the interaction files is described in the README file 
  in the directory 'interactions/'.

- To benchmark the code, you can verify that the output from running the 
  example input files are the same as in 'example-inout-files/out_benchmark/'.

- Have fun!
