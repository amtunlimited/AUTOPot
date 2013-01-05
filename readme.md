#AUTOPot

##A wrapping library for the POTLIB2001 standard to Mathematica

###Written By Aaron Tagliaboschi and Dr. Jeremy Maddox

###Department of Chemistry, Western Kentucky University

***

##Abstract

Chemical reactions may be described theoretically using Born-Oppenheimer potential energy
surfaces. For many chemical systems, these have been incorporated within an on-line
repository of FORTRAN codes (POTLIB). The goal of this project is to be able to install and
utilize POTLIB subroutines in higher level languages, such as Mathematica. To do this we
have utilized the GNU C compiler and FORTRAN compiler to link the object files along with
MathLink drivers written in C. We have found this to be an effective method and can now
successfully generate MathLink executables and call the POTLIB functions from within
Mathematica without a serious cost to runtime. This is helpful for analyzing reactions in
Mathematica and using some of its higher-level numerical functions. Thus far, we have been
able to make graphs of the potential energy surface in internal coordinates and Jacobi
coordinates, as well as apply numerical methods for finding the minimum energy reaction
coordinate and tracing Newtonian trajectories.

##Installation

No installation is necessary beyond simply saving the folder to the computer.

##Use


