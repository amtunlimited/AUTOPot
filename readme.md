#AUTOPot

####A wrapping library for the POTLIB2001 standard to *Mathematica*

####Written By Aaron Tagliaboschi and Dr. Jeremy Maddox

####Department of Chemistry, Western Kentucky University

***

##Abstract

Chemical reactions may be described theoretically using Born-Oppenheimer
potential energy surfaces. For many chemical systems, these have been
incorporated within an on-line repository of FORTRAN codes (POTLIB). The goal
of this project is to be able to install and utilize POTLIB subroutines in
higher level languages, such as *Mathematica*. To do this we have utilized the
GNU C compiler and FORTRAN compiler to link the object files along with
MathLink drivers written in C. We have found this to be an effective method and
can now successfully generate MathLink executables and call the POTLIB
functions from within *Mathematica* without a serious cost to runtime. This is
helpful for analyzing reactions in *Mathematica* and using some of its higher-
level numerical functions. Thus far, we have been able to make graphs of the
potential energy surface in internal coordinates and Jacobi coordinates, as
well as apply numerical methods for finding the minimum energy reaction
coordinate and tracing Newtonian trajectories.

##Installation

No installation is necessary beyond simply saving the folder to the computer.
To use the included scripts and to build the libraries, you will need GCC,
GFORTRAN, and perl installed on your system. If you are using Linux, these are
normally installed by default. Compiling on/building for Windows is possible,
but you will to install Cygwin. [See below](#cygwin) for more details.

##Usage Instruction

To use this library, you have to build the POTLIB source library against pot.c
and utility.f along with generating code for MathLink using the mprep utility
that comes with the MathLink libraries. To make things simpler, two scripts are
included with this package to automatically generate everything needed to
create a build package, including all generated code, the library file itself,
and a makefile.

The first script `configure` is a shell script for single packages and is
the backbone of the batch script as well. The call format is:

    ./configure POTLIBfile platform
    
The first argument `POTLIBfile` is the actual POTLIB FORTRAN file to be used
and `platform` is the platform you are building for, either 'cygwin', 'linux',
or 'mac'. The script copies the library to the AUTOPot folder and seds through
the templates for the makefile and mprep MathLink template files. To package,
just make sure to copy the src, tmp, and makefile in the same folder.

**Note:** The makefile assumes the mprep utility, ML32i3 library, and
mathlink.h header file are installed in the execution path, where your computer
and GCC defaults to. The files are found in:

    [Mathematica install folder]/SystemFiles/Links/MathLink/DeveloperKit/[Operating System]/CompilerAdditions

The second script `batch-config` compiles batches of MathLink excecutables It's
format is:

    ./batch-config POTfolder binfolder platform
    
The first argument `POTfolder` is the folder containing the POTLIB FORTRAN files
ready for compiling. The second argument `binfolder` is the folder the resulting
executables will go. The last argument `platform` is for the `configure` script.

This script builds each .f file in `POTfolder` using the `configure` script.

##Cygwin

To compile MathLink executables on Windows, this package uses *Cygwin*, a Linux
emulator, and *MinGW*, a minimalist GCC compiler for Windows. In Cygwin, you
need to have the following packages installed:

* mingw-gcc
* mingw-gfortran
* make
* perl

You will also have to install the mprep executable, ML32i3 library, and
mathlink.h in the correct locations. Other than that, it works exactly as in
Linux. The final product is a MathLink executable.

**Note:** Because of licensing of `cygwin1.dll`, if you distribute Cygwin-built
MathLink executables you much release the source code under the GPLv3 or a
compatible license.

##License

This software is licensed under the GPLv3. For more information, please see
`COPYING`