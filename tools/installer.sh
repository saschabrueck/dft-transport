#!/bin/bash -e

# /*
# Copyright (c) 2017 ETH Zurich
# Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# */

# =========================================================================
# \brief    script for installing CP2K-OMEN 
# This script installs CP2K-OMEN and compiles the code in 4 steps:
# STEP 1: installing CP2K;
# STEP 2: installing solvers/libraries that may be used by OMEN;
# STEP 3: generating a .mk file (machine specific Makefile) for OMEN; 
# STEP 4: compiling CP2K-OMEN.
#
# Although all the steps can be done manually, for convenience and to have 
# a well-matched installation, it is advised to use the script for at least 
# the first two steps. 
#
# \history  created: 22-07-2017 
# \author   Mohammad Hossein Bani-Hashemian 
# =========================================================================

# * Set the following paths *******
TOPDIR=~/Projects/cp2k-omen
cp2kDIR=${TOPDIR}/cp2k
omenDIR=${TOPDIR}/omen

installDIR=${omenDIR}/install
buildDIR=${installDIR}/build
libsDIR=${installDIR}/libs

# * Preparation for STEP 1 ********
# Install either the current development version (dev) of CP2K or a release version (release) (>=4.1)
cp2kDEVorRELEASE="dev"
cp2kRELEASEVER="7.1"

# * Preparation for STEP 2 ********
# Specify the versions of the solvers that you would like to install:
# An unset version would mean the solver will not be installed.  
parmetisVER="4.0.3"
hypreVER="2.11.2"
SuiteSparseVER="4.5.5"
qhullVER="2015-src-7.2.0"
mumpsVER="5.1.1"

# If you would like to use SplitSolve, set magmaVER. 
# Note that SplitSolve uses MAGMA and MAGMA ver>=2.0 requires CUDA ver>=5.0
cudaDIR="/usr/local/cuda-7.5"
magmaVER="2.2.0"

# If you have downloaded the shared library file of PARDISO, specify the name 
# of the file (pardisoSOFILE), the path to the .so file (pardisoLIBDIR), and
# the path to the license file (pardisoLICPATH). 
# If you do not wish to use PARDISO, leave pardisoLIBDIR unset. 
#pardisoSOFILE="libpardiso500-MPI-GNU472-X86-64.so"
#pardisoLIBDIR=${libsDIR}/pardiso/lib
#pardisoLICPATH=${libsDIR}/pardiso/lib
pardisoSOFILE=
pardisoLIBDIR=
pardisoLICPATH=

# * Preparation for STEP 3 ********
# This step is optional
generate_makefile="yes"
# Specify the name of your machine:
machine="mont-fort1"

# * Preparation for STEP 4 ********
# This step is optional
compile_cp2komen="yes"
# Specify the number of processes for compilation:
Nproc=
# Specify the cp2k VERSION for compilation:
cp2k_target="popt"
# Set the OMEN's configure flags, i.e. solvers to be enabled:
omenCONFIGURE_ARGS="--with-pexsi --with-mumps"

# * You should not need to change anything below this line, unless you would like to modify the steps.

# STEP 1: ******************************************************************************************

# install CP2K: ===========================================================
# Here the script uses the CP2K toolchain script to install CP2K together
# with PEXSI and SuperLU_DIST that all are used by OMEN as well.
# -------------------------------------------------------------------------
echo "Installing CP2K ==========================================="
cd ${TOPDIR} 

if [ "${cp2kDEVorRELEASE}" = "dev" ] ; then
   git clone --recursive https://github.com/cp2k/cp2k.git cp2k
elif [ "${cp2kDEVorRELEASE}" = "release" ] ; then
   if [ -z ${cp2kRELEASEVER} ]; then
      echo "CP2K release version (cp2kRELEASEVER) needs to be set."
      exit 1
   fi
   git clone -b support/v${cp2kRELEASEVER} https://github.com/cp2k/cp2k.git cp2k
fi

cp2k_toolchainDIR=${cp2kDIR}/tools/toolchain

cd ${cp2k_toolchainDIR}
./install_cp2k_toolchain.sh \
                 --with-gcc \
                 --with-mpich \
                 --with-pexsi \
                 --no-check-certificate 

cp ${cp2k_toolchainDIR}/install/arch/* ${cp2kDIR}/arch/  

# IMPORTANT: Set up environment variables before moving on:
source ${cp2k_toolchainDIR}/install/setup

# define lib path for .mk file (STEP 3)
lib_cp2k="-L${cp2kDIR}/lib/local/${cp2k_target} -lcp2k -lxsmmf -lxsmm -lxcf03 -lxc -lint2 -lsymspg -lelpa -lfftw3"
lib_dbcsr="-L${cp2kDIR}/lib/local/${cp2k_target}/exts/dbcsr -ldbcsr"

# STEP 2: ******************************************************************************************

# install OMEN solvers: ===================================================
# the following are installed by CP2K toolchain
# get versions of libraries that have been just installed by CP2K toolchain
scalapackVER=$(awk -F'"' '/^scalapack_ver=/ {print $2}' ${cp2k_toolchainDIR}/scripts/install_scalapack.sh)
mpichVER=$(awk -F'"' '/^mpich_ver=/ {print $2}' ${cp2k_toolchainDIR}/scripts/install_mpich.sh)
superluVER=$(awk -F'"' '/^superlu_ver=/ {print $2}' ${cp2k_toolchainDIR}/scripts/install_superlu.sh)
pexsiVER=$(awk -F'"' '/^pexsi_ver=/ {print $2}' ${cp2k_toolchainDIR}/scripts/install_pexsi.sh)
openblasVER=$(awk -F'"' '/^openblas_ver=/ {print $2}' ${cp2k_toolchainDIR}/scripts/install_openblas.sh)

# -------------------------------------------------------------------------
mkdir -p ${installDIR}
mkdir -p ${libsDIR}
mkdir -p ${buildDIR}

# script for setting environment variables
cp2komenENVSETUP="${installDIR}/cp2komen_envsetup"
cat <<EOF > ${cp2komenENVSETUP}
#!/bin/bash
EOF

echo "Installing OMEN solvers ==================================="
echo " "

# PEXSI -------------------------------------------------------------------
# copy PEXSI headers that are not copied via CP2K toolchain script 
cp -r ${cp2k_toolchainDIR}/build/pexsi_v${pexsiVER}/include/pexsi/ ${cp2k_toolchainDIR}/install/pexsi-${pexsiVER}/include/

# define include and lib paths for .mk file (STEP 3)
inc_sludist="-I\$(TOOLCHAIN)/install/superlu_dist-${superluVER}/include/"
lib_sludist="-L\$(TOOLCHAIN)/install/superlu_dist-${superluVER}/lib -lsuperlu_dist"
inc_pexsi="-I\$(TOOLCHAIN)/install/pexsi-${pexsiVER}/include/"
lib_pexsi="-L\$(TOOLCHAIN)/install/pexsi-${pexsiVER}/lib -lpexsi"

# ParMETIS ----------------------------------------------------------------
if [ ! -z ${parmetisVER} ]; then
   echo "Installing ParMETIS ======================================="
   mkdir -p ${libsDIR}/parmetis-${parmetisVER}
   mkdir -p ${libsDIR}/parmetis-${parmetisVER}/lib
   cd ${buildDIR}

   wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-${parmetisVER}.tar.gz
   tar -xf parmetis-${parmetisVER}.tar.gz

   cd ${buildDIR}/parmetis-${parmetisVER}

   make config \
        cc=${MPICC:-mpicc} \
        cxx=${MPICXX:-mpic++} \
        prefix=${buildDIR}/parmetis-${parmetisVER} > configure.log 2>&1
   make > make.log 2>&1
   make install > install.log 2>&1

   cd metis
   make config \
        cc=${MPICC:-mpicc} \
        cxx=${MPICXX:-mpic++} \
        prefix=${buildDIR}/parmetis-${parmetisVER} > configure.log 2>&1
   make > make.log 2>&1
   make install > install.log 2>&1

   find . -name "*.a" -type f -exec cp {} ${libsDIR}/parmetis-${parmetisVER}/lib \;

   # define include and lib paths for .mk file (STEP 3)
   lib_parmetis="-L\$(LIB_TOP)/parmetis-${parmetisVER}/lib -lparmetis -lmetis"
fi

# HYPRE -------------------------------------------------------------------
if [ ! -z ${hypreVER} ]; then
   echo "Installing hypre =========================================="
   mkdir -p ${libsDIR}/hypre
   mkdir -p ${libsDIR}/hypre/bin
   mkdir -p ${libsDIR}/hypre/libexec
   mkdir -p ${libsDIR}/hypre/lib
   mkdir -p ${libsDIR}/hypre/include
   cd ${buildDIR}

   wget https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-${hypreVER}.tar.gz
   tar -xf hypre-"${hypreVER}".tar.gz

   cd ${buildDIR}/hypre-${hypreVER}/src

   ./configure -q --bindir=${libsDIR}/hypre/bin \
                  --libexecdir=${libsDIR}/hypre/libexec \
                  --libdir=${libsDIR}/hypre/lib \
                  --includedir=${libsDIR}/hypre/include \
                  --without-superlu

   make install > install.log 2>&1

   # define include and lib paths for .mk file (STEP 3)
   inc_hypre="-I\$(LIB_TOP)/hypre/include/"
   lib_hypre="-L\$(LIB_TOP)/hypre/lib -lHYPRE"
fi

# SuiteSparse -------------------------------------------------------------
if [ ! -z ${SuiteSparseVER} ]; then
   echo "Installing SuiteSparse ===================================="
   mkdir -p ${libsDIR}/SuiteSparse
   mkdir -p ${libsDIR}/SuiteSparse/static
   cd ${buildDIR}

   wget http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-${SuiteSparseVER}.tar.gz
   tar -xf SuiteSparse-${SuiteSparseVER}.tar.gz

   cd ${buildDIR}/SuiteSparse

   make install INSTALL_LIB=${libsDIR}/SuiteSparse/lib \
                INSTALL_INCLUDE=${libsDIR}/SuiteSparse/include \
                BLAS=${cp2k_toolchainDIR}/install/openblas-${openblasVER}/lib/libopenblas.a  \
                LAPACK=${cp2k_toolchainDIR}/install/scalapack-${scalapackVER}/lib/libscalapack.a \
                AUTOCC=no \
                CUDA=no \
                MY_METIS_LIB=${buildDIR}/parmetis-${parmetisVER}/lib/libmetis.a \
                MY_METIS_INC=${buildDIR}/parmetis-${parmetisVER}/include \
                > install.log 2>&1

   find . -name "*.a" -type f -exec cp {} ${libsDIR}/SuiteSparse/static \;

   # define include and lib paths for .mk file (STEP 3)
   inc_ssparse="-I\$(LIB_TOP)/SuiteSparse/include/"
   lib_ssparse="-L\$(LIB_TOP)/SuiteSparse/static -lumfpack -lamd -lcholmod -lcolamd -lccolamd -lcamd -lsuitesparseconfig"
fi

# Qhull -------------------------------------------------------------------
if [ ! -z ${qhullVER} ]; then
   echo "Installing qhull =========================================="
   mkdir -p ${libsDIR}/qhull
   cd ${buildDIR}

   wget http://www.qhull.org/download/qhull-${qhullVER}.tgz
   mkdir -p ${buildDIR}/qhull
   tar -xf qhull-${qhullVER}.tgz -C qhull/ --strip-components 1

   cd ${buildDIR}/qhull

   mv ${buildDIR}/qhull/Makefile ${buildDIR}/qhull/Makefile.orig
   sed -e "s|\(DESTDIR *=\).*|\1 ${libsDIR}/qhull|g" \
       -e "s|\(CC *=\).*|\1 gcc -w|g" ${buildDIR}/qhull/Makefile.orig > ${buildDIR}/qhull/Makefile
   make all > make.log 2>&1 
   make install > install.log 2>&1

   # define include and lib paths for .mk file (STEP 3)
   inc_qhull="-I\$(LIB_TOP)/qhull/include/libqhull/"
   lib_qhull="-L\$(LIB_TOP)/qhull/lib -lqhullstatic"
fi

# MUMPS -------------------------------------------------------------------
if [ ! -z ${mumpsVER} ]; then
   echo "Installing MUMPS =========================================="
   mkdir -p ${libsDIR}/MUMPS
   cd ${buildDIR}

   wget http://mumps.enseeiht.fr/MUMPS_${mumpsVER}.tar.gz
   mkdir -p ${buildDIR}/MUMPS
   tar -xf MUMPS_${mumpsVER}.tar.gz -C MUMPS/ --strip-components 1

   cd ${buildDIR}/MUMPS

   sed -e "/^#.*LMETISDIR *=/s/^#//" \
       -e "/^#.*IMETIS *=/s/^#//" \
       -e "s|\(LMETISDIR *=\).*|\1 ${buildDIR}/parmetis-${parmetisVER}/lib|g" \
       -e "s|\(IMETIS *=\).*|\1 -I${buildDIR}/parmetis-${parmetisVER}/include|g" \
       -e "/^#.*-lparmetis/s/^#//" \
       -e "s|\(CC *=\).*|\1 gcc|g" \
       -e "s|\(FC *=\).*|\1 mpif90|g" \
       -e "s|\(FL *=\).*|\1 mpif90|g" \
       -e "s|\(LAPACK *=\).*|\1 -L/${cp2k_toolchainDIR}/install/scalapack-${scalapackVER}/lib -lscalapack|g" \
       -e "s|\(SCALAP *=\).*|\1 -L/${cp2k_toolchainDIR}/install/scalapack-${scalapackVER}/lib -lscalapack|g" \
       -e "s|\(INCPAR *=\).*|\1 -I/${cp2k_toolchainDIR}/install/mpich-${mpichVER}/include|g" \
       -e "s|\(LIBPAR *=\).*|\1 \$(SCALAP) \$(LAPACK) -L/${cp2k_toolchainDIR}/install/mpich-${mpichVER}/lib -lmpi -lgfortran|g" \
       -e "s|\(LIBBLAS *=\).*|\1 -L/${cp2k_toolchainDIR}/install/openblas-${openblasVER}/lib -lopenblas|g" \
       -e "s|\(OPTF *=\).*|\1 -O -w|g" \
       -e "s|\(OPTC *=\).*|\1 -O3 -w|g" \
          Make.inc/Makefile.inc.generic > Makefile.inc

   make alllib > install.log 2>&1
   cp -r ${buildDIR}/MUMPS/lib/ ${buildDIR}/MUMPS/include/ ${libsDIR}/MUMPS

   # define include and lib paths for .mk file (STEP 3)
   inc_mumps="-I\$(LIB_TOP)/MUMPS/include/"
   lib_mumps="-L\$(LIB_TOP)/MUMPS/lib -lzmumps -ldmumps -lmumps_common -lpord -lscalapack"
fi

# MAGMA -------------------------------------------------------------------
if [ ! -z ${magmaVER} ]; then
   echo "Installing MAGMA =========================================="
   mkdir -p ${libsDIR}/magma
   cd ${buildDIR}

   wget http://icl.cs.utk.edu/projectsfiles/magma/downloads/magma-${magmaVER}.tar.gz
   mkdir -p ${buildDIR}/magma
   tar -xf magma-${magmaVER}.tar.gz -C magma/ --strip-components 1

   cd ${buildDIR}/magma

   cp ${buildDIR}/magma/make.inc-examples/make.inc.openblas ${buildDIR}/magma/

   if [ -z ${cudaDIR} ]; then
      echo "Path to the CUDA Toolkit (cudaDIR) needs to be set."
      exit 1
   fi

   sed -e "/^#.*OPENBLASDIR *?=/s/^#//" \
       -e "/^#.*CUDADIR *?=/s/^#//" \
       -e "s|\(OPENBLASDIR *?=\).*|\1 ${cp2k_toolchainDIR}/install/openblas-${openblasVER}|g" \
       -e "s|\(CUDADIR *?=\).*|\1 ${cudaDIR}|g" \
       -e "s|\(NVCC *=\).*|\1 ${cudaDIR}/bin/nvcc -Xcompiler=--std=gnu++98 -D__GNUC__=4 -D__GNUC_MINOR__=9|g" \
       -e "/^NVCCFLAGS/ s/$/ -w/" \
       -e "/^-.*include/s/^-/#-/" \
          make.inc.openblas > make.inc

   make install prefix=${libsDIR}/magma > install.log 2>&1

   # define include and lib paths for .mk file (STEP 3)
   nvcc_path="${cudaDIR}/bin/nvcc"
   nvcc_flags="-Xcompiler=--std=gnu++98 -D__GNUC__=4 -D__GNUC_MINOR__=9 -w"
   inc_cuda="-I${cudaDIR}/include/"
   lib_cuda="-L${cudaDIR}/lib64/ -lcudart -lcublas -lcusparse -lblas"
   inc_magma="-I\$(LIB_TOP)/magma/include/"
   lib_magma="-L\$(LIB_TOP)/magma/lib -lmagma"

# set CUDA related environment variables 
cat <<EOF >> ${cp2komenENVSETUP}
LD_LIBRARY_PATH=${cudaDIR}/lib64:\$LD_LIBRARY_PATH
PATH=${cudaDIR}/bin:\$PATH
EOF

# add path to the MAGMA shared library files to LD_LIBRARY_PATH
cat <<EOF >> ${cp2komenENVSETUP}
LD_LIBRARY_PATH=${libsDIR}/magma/lib:\$LD_LIBRARY_PATH
EOF
fi

# PARDISO -----------------------------------------------------------------
# Downloading PARDISO needs registration on their website. Therefore, this 
# script does not automatically install PARDISO.  

# define include and lib paths for .mk file (STEP 3)
if [ ! -z ${pardisoLIBDIR} ]; then
   pardiso_file=${pardisoSOFILE/"lib"/ -l}
   lib_pardiso="-L${pardisoLIBDIR}/${pardiso_file/".so"/}"

# add path to the PARDISO shared library file to LD_LIBRARY_PATH
# add path to pardiso.lic to PATH
cat <<EOF >> ${cp2komenENVSETUP}
LD_LIBRARY_PATH=${pardisoLIBDIR}:\$LD_LIBRARY_PATH
PARDISO_LIC_PATH=${pardisoLICPATH}
EOF
fi

echo "Done! ====================================================="
echo " "

# STEP 3: ******************************************************************************************

# generate a .mk file:
# -------------------------------------------------------------------------
if [ "${generate_makefile}" = "yes" ] ; then
   echo "Generate a .mk file ======================================="
   cd ${omenDIR}

   topdir=${installDIR}
   libtop="\$(TOP_DIR)/libs"

   sed -e "s|\(source \).*|\1 ${cp2k_toolchainDIR}/install/setup ; source ${cp2komenENVSETUP}|g" \
       -e "s|\(TOP_DIR *=\).*|\1 ${topdir}|g" \
       -e "s|\(LIB_TOP *=\).*|\1 ${libtop}|g" \
       -e "s|\(TOOLCHAIN *=\).*|\1 ${cp2k_toolchainDIR}|g" \
       -e "s|\(INCSSPARSE *=\).*|\1 ${inc_ssparse}|g" \
       -e "s|\(INCMUMPS *=\).*|\1 ${inc_mumps}|g" \
       -e "s|\(INCHYPRE *=\).*|\1 ${inc_hypre}|g" \
       -e "s|\(INCQHULL *=\).*|\1 ${inc_qhull}|g" \
       -e "s|\(INCSLUDIST *=\).*|\1 ${inc_sludist}|g" \
       -e "s|\(INCPEXSI *=\).*|\1 ${inc_pexsi}|g" \
       -e "s|\(LIBSSPARSE *=\).*|\1 ${lib_ssparse}|g" \
       -e "s|\(LIBMUMPS *=\).*|\1 ${lib_mumps}|g" \
       -e "s|\(LIBHYPRE *=\).*|\1 ${lib_hypre}|g" \
       -e "s|\(LIBQHULL *=\).*|\1 ${lib_qhull}|g" \
       -e "s|\(LIBPARMETIS *=\).*|\1 ${lib_parmetis}|g" \
       -e "s|\(LIBSLUDIST *=\).*|\1 ${lib_sludist}|g" \
       -e "s|\(LIBPEXSI *=\).*|\1 ${lib_pexsi}|g" \
       -e "s|\(LIBDBCSR *=\).*|\1 ${lib_dbcsr}|g" \
       -e "s|\(LIBCP2K *=\).*|\1 ${lib_cp2k}|g" \
       -e "s|\(LIBPARDISO *=\).*|\1 ${lib_pardiso}|g" \
       -e "s|\(NVCC *=\).*|\1 ${nvcc_path}|g" \
       -e "s|\(NVCCFLAGS *=\).*|\1 ${nvcc_flags}|g" \
       -e "s|\(INCCUDA *=\).*|\1 ${inc_cuda}|g" \
       -e "s|\(LIBCUDA *=\).*|\1 ${lib_cuda}|g" \
       -e "s|\(INCMAGMA *=\).*|\1 ${inc_magma}|g" \
       -e "s|\(LIBMAGMA *=\).*|\1 ${lib_magma}|g" \
           ${omenDIR}/makefiles/arch.tmpl > ${omenDIR}/makefiles/${machine}.mk

   echo "Done! ====================================================="
   echo " "
fi

# STEP 4: ******************************************************************************************

# compile CP2K-OMEN:
# -------------------------------------------------------------------------
if [ "${compile_cp2komen}" = "yes" ] ; then
   source ${cp2k_toolchainDIR}/install/setup
   # compile CP2K
   echo "Compiling cp2k with target ${cp2k_target} libcp2k ========="
   cd ${cp2kDIR}
#   make -j ${Nproc} ARCH=local VERSION=${cp2k_target} 
   make -j ${Nproc} ARCH=local VERSION=${cp2k_target} libcp2k

   source ${cp2komenENVSETUP}
   # compile OMEN
   echo "Compiling OMEN ============================================"
   cd ${omenDIR}/src
   ./configure --with-arch=${machine}.mk ${omenCONFIGURE_ARGS}
   make -j ${Nproc} 
fi

cd ${omenDIR}


