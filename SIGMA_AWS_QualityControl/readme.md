SIGMA-AWS-QC program

The SIGMA-AWS-QC program is a quality control (QC) program for meteorological data observed by two automatic weather stations (SIGMA-A and SIGMA-B) installed in Northwest Greenland since 2012.
The dataset before QC processing is stored in ./Lv1.1_to_1.2/input.
The data after processing will be stored in ./Lv1.1_to_1.2/output.

To run this program, please execute the commands in the order of;
./Lv1.1_to_Lv.1.2
$make
$cd run
$make 
$sh run_make_L1.2.sh

However, it is necessary to edit "site_flag" in run_make_L1.2.sh before execution.
Define the name of the site for data processing as "site_flag" before execution.

The main program is make_Level1.2data.f90.
The module related to data loading is data_IO.f90
QC_L1.1.f90 is a module for subroutines related to QualityControl
sol_info.f90 is a module for calculating solar zenith and azimuth angles
The make file in the source directory is a makefile for compilation.

In creating Level 1.3 dataset from Level 1.2 dataset, the procedure is basically the same as above.
However, the input data for this process is set to refer to Lv1.1_to_1.2/output.
Subroutines_Lv1.2.f90 used in this process for subroutines to calculate other physical quantities and processes required by QC, which is a module that compiles subroutines required for QC and subroutines to calculate other physical quantities.

Those programs have been tested on CentOS7, and gfortran and Intel fortran compiler.
By editing the FC and FCFLAGS variables in the makefile, the compiler can be selected according to the runtime environment of the executor.

Author information
Name: Motoshi Nishimura
Affiriate: National Institute of Polar Research

License: CC-BY 4.0