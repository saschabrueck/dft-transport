About
-----
A massively parallel quantum mechanical device simulator based on DFT (Density Functional Theory) and NEGF (Non-Equilibrium Green's Function) formalisms capable of simulating realistically large devices consisting of several thousands of atoms. For more information see the references below.

Installation
------------
An installer script is distributed with the source code that installs CP2K, OMEN, and all the solvers/libraries that may be employed for calculations. In case manual installation is preferred, the script is well-documented with the steps that are needed to be taken.

Compile and Run
---------------
In order to compile the code, follow the following steps:

* Compile CP2K  
    * Compile CP2K with target `libcp2k`, for example:  
`> make -j N ARCH=Linux-x86-64-gfortran VERSION=popt libcp2k`
* Compile OMEN
    * cd to `makefiles/` and using a sample `.mk` file (e.g. `arch1.mk`) write a `.mk` file (say, `myarch.mk`) according to your local installations.
    * cd to the source directory, `src/`, and run configure:  
`> ./configure --with-arch=myarch`
    * Various solvers can be enabled by passing one or multiple of the following options to the configure script:  
`--with-pardiso --with-mumps --with-superlu --with-pexsi --with-splitsolve --with-fempoisson`
    * Now, running make in the current directory (`src/`) will create an executable named `transport`:  
`> make -j N`  
 
To run:  
`> mpirun -np N ./transport myinput.inp`  
 
The input file, `myinput.inp`, should contain a `&TRANSPORT` section. For more details, see the corresponding page on CP2K's Reference Manual: [TRANSPORT](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/TRANSPORT.html).

Example
-------
An example input file, `gnr.inp` , for a graphene nanoribbon system composed of 96 atoms can be found in the `tests` directory.

References
----------
1. S. Brück, M. Calderara, M. H. Bani-Hashemian, J. VandeVondele, and M. Luisier. *Towards Ab-Initio Simulations of Nanowire Field-Effect Transistors*. Proceedings of the International Workshop on Computational Electronics (IWCE), June 2014, Paris, France. doi: [10.1109/IWCE.2014.6865831](http://doi.org/10.1109/IWCE.2014.6865831).
2. S. Brück, M. Calderara, M. H. Bani-Hashemian, J. VandeVondele, and M. Luisier. *Efficient Algorithms for Large-Scale Quantum Transport Calculations*. The Journal of Chemical Physics 147(7): 074116, Aug. 2017. doi:[10.1063/1.4998421](http://doi.org/10.1063/1.4998421).

Theses
------
1. S. Brück, *Ab-initio Quantum Transport Simulations for Nanoelectronic Devices*. PhD thesis, ETH Zurich, 2017. doi:[10.3929/ethz-b-000226622](http://doi.org/10.3929/ethz-b-000226622).
2. M. Calderara, *SplitSolve, an Algorithm for Ab-Initio Quantum Transport Simulations*. PhD thesis, ETH Zurich, 2016. doi:[10.3929/ethz-a-010781848](https://doi.org/10.3929/ethz-a-010781848).
3. M. H. Bani-Hashemian, *Large-Scale Nanoelectronic Device Simulation from First Principles*. PhD thesis, ETH Zurich, 2016. doi:[10.3929/ethz-a-010811338](http://doi.org/10.3929/ethz-a-010811338).

Other Resources
---------------
1. [Nano-TCAD group](https://nano-tcad.ee.ethz.ch), ETH Zurich.  
2. CP2K user manual: [TRANSPORT](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/TRANSPORT.html).  
3. [CP2K wiki-page](https://www.cp2k.org/howto:cp2k_omen).
