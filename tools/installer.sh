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
TOPDIR=/scratch/seyedb
cp2kDIR=${TOPDIR}/cp2k
omenDIR=${TOPDIR}/omen

installDIR=${omenDIR}/install
buildDIR=${installDIR}/build
libsDIR=${installDIR}/libs

# * Preparation for STEP 1 ********
# Install either the current truck version of CP2K or a release version (>=4.1)
cp2kTRUNK="yes"
cp2kRELEASE="no"
cp2kRELEASEVER="4_1"

# * Preparation for STEP 2 ********
# Specify the versions of the solvers that you would like to install:
# An unset version would mean the solver will not be installed.  
hypreVER="2.11.2"
SuiteSparseVER="4.5.5"
qhullVER="2015-src-7.2.0"
mumpsVER="5.1.1"

# If you would like to use SplitSolve, set magmaVER. 
# Note that SplitSolve uses MAGMA and MAGMA ver>=2.0 requires CUDA ver>=5.0
cudaDIR="/usr/local/cuda-7.5"
magmaVER="2.2.0"

# If you have downloaded the shared library file of PARDISO, specify the name 
# of the file and the path to it here, otherwise, leave pardisoDIR unset:
pardisoSOFILE="libpardiso500-MPI-GNU472-X86-64.so"
#pardisoDIR=${libsDIR}/pardiso/lib
pardisoDIR=

# * Preparation for STEP 3 ********
# This step is optional
generate_makefile="yes"
# Specify the name of your machine:
machine="mont-fort1"

# * Preparation for STEP 4 ********
# This step is optional
compile_cp2komen="yes"
# Specify the cp2k VERSION for compilation:
cp2k_target="popt"
# Set the OMEN's configure flags, i.e. solvers to be enabled:
omenCONFIGURE_ARGS="--with-pexsi --with-mumps"

# * You should not need to change anything below this line, unless you would like to modify the steps.

# STEP 1: ******************************************************************************************

# install CP2K: ===========================================================
# Here the script uses the CP2K toolchain script to install CP2K together
# with PEXSI, SuperLU_DIST, and ParMETIS that all are used by OMEN as well.
# -------------------------------------------------------------------------
echo "installing CP2K ==========================================="
cd ${TOPDIR} 

if [ "${cp2kTRUNK}" = "yes" ] ; then

elif [ "${cp2kRELEASE}" = "yes" ] ; then 
   svn checkout http://svn.code.sf.net/p/cp2k/code/branches/cp2k-${cp2kRELEASEVER}-branch cp2k
fi

cp2k_toolchainDIR=${cp2kDIR}/cp2k/tools/toolchain

cd ${cp2k_toolchainDIR}
./install_cp2k_toolchain.sh \
                 --with-gcc \
                 --with-binutils \
                 --with-mpich \
                 --with-pexsi \
                 --no-check-certificate 

cp ${cp2k_toolchainDIR}/install/arch/* ${cp2kDIR}/cp2k/arch/  

# IMPORTANT: Set up environment variables before moving on:
source ${cp2k_toolchainDIR}/install/setup

# get versions of libraries that have been just installed by CP2K toolchain
source ${cp2k_toolchainDIR}/scripts/package_versions.sh

# define lib path for .mk file (STEP 3)
lib_cp2k="-L${cp2kDIR}/cp2k/lib/local/${cp2k_target} -lcp2k -lxsmmf -lxsmm -lderiv -lint -lxcf90 -lxc -lfftw3"

# STEP 2: ******************************************************************************************

# install OMEN solvers: ===================================================
# the following are installed by CP2K toolchain
parmetisVER=${parmetis_ver}
reflapackVER=${reflapack_ver}
scalapackVER=${scalapack_ver}
mpichVER=${mpich_ver}
superluVER=${superlu_ver}
pexsiVER=${pexsi_ver}
openblasVER=${openblas_ver}
# -------------------------------------------------------------------------
mkdir -p ${installDIR}
mkdir -p ${libsDIR}
mkdir -p ${buildDIR}

echo "installing OMEN solvers ==================================="
echo " "

# PEXSI -------------------------------------------------------------------
# copy PEXSI headers that are not copied via CP2K toolchain script 
cp -r ${cp2k_toolchainDIR}/build/pexsi_v${pexsiVER}/include/pexsi/ ${cp2k_toolchainDIR}/install/pexsi-${pexsiVER}/include/

# define include and lib paths for .mk file (STEP 3)
lib_parmetis="-L\$(TOOLCHAIN)/install/parmetis-${parmetisVER}/lib -lparmetis -lmetis"
inc_sludist="-I\$(TOOLCHAIN)/install/superlu_dist-${superluVER}/include/"
lib_sludist="-L\$(TOOLCHAIN)/install/superlu_dist-${superluVER}/lib -lsuperlu_dist"
inc_pexsi="-I\$(TOOLCHAIN)/install/pexsi-${pexsiVER}/include/"
lib_pexsi="-L\$(TOOLCHAIN)/install/pexsi-${pexsiVER}/lib -lpexsi"

# HYPRE -------------------------------------------------------------------
if [ ! -z ${hypreVER} ]; then
   echo "installing hypre =========================================="
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
   
   make install > install.log

   # define include and lib paths for .mk file (STEP 3)
   inc_hypre="-I\$(LIB_TOP)/hypre/include/"
   lib_hypre="-L\$(LIB_TOP)/hypre/lib -lHYPRE"
fi

# SuiteSparse -------------------------------------------------------------
if [ ! -z ${SuiteSparseVER} ]; then
   echo "installing SuiteSparse ===================================="
   mkdir -p ${libsDIR}/SuiteSparse
   mkdir -p ${libsDIR}/SuiteSparse/static
   cd ${buildDIR}
   
   wget http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-${SuiteSparseVER}.tar.gz
   tar -xf SuiteSparse-${SuiteSparseVER}.tar.gz
   
   cd ${buildDIR}/SuiteSparse
   
   make install INSTALL_LIB=${libsDIR}/SuiteSparse/lib \
                INSTALL_INCLUDE=${libsDIR}/SuiteSparse/include \
                AUTOCC=no \
                CUDA=no \
                MY_METIS_LIB=${cp2k_toolchainDIR}/install/parmetis-${parmetisVER}/lib/libmetis.a \
                MY_METIS_INC=${cp2k_toolchainDIR}/install/parmetis-${parmetisVER}/include \
                > install.log
   
   find . -name "*.a" -type f -exec cp {} ${libsDIR}/SuiteSparse/static \;

   # define include and lib paths for .mk file (STEP 3)
   inc_ssparse="-I\$(LIB_TOP)/SuiteSparse/include/"
   lib_ssparse="-L\$(LIB_TOP)/SuiteSparse/static -lumfpack -lamd -lcholmod -lcolamd -lccolamd -lcamd -lsuitesparseconfig"
fi

# Qhull -------------------------------------------------------------------
if [ ! -z ${qhullVER} ]; then
   echo "installing qhull =========================================="
   mkdir -p ${libsDIR}/qhull
   cd ${buildDIR}
   
   wget http://www.qhull.org/download/qhull-${qhullVER}.tgz
   mkdir -p ${buildDIR}/qhull
   tar -xf qhull-${qhullVER}.tgz -C qhull/ --strip-components 1
   
   cd ${buildDIR}/qhull
   
   mv ${buildDIR}/qhull/Makefile ${buildDIR}/qhull/Makefile.orig
   sed -e "s|\(DESTDIR *=\).*|\1 ${libsDIR}/qhull|g" \
       -e "s|\(CC *=\).*|\1 gcc -w|g" ${buildDIR}/qhull/Makefile.orig > ${buildDIR}/qhull/Makefile
   make all > install.log 
   make install

   # define include and lib paths for .mk file (STEP 3)
   inc_qhull="-I\$(LIB_TOP)/qhull/include/libqhull/"
   lib_qhull="-L\$(LIB_TOP)/qhull/lib -lqhullstatic"
fi

# MUMPS -------------------------------------------------------------------
if [ ! -z ${mumpsVER} ]; then
   echo "installing MUMPS =========================================="
   mkdir -p ${libsDIR}/MUMPS
   cd ${buildDIR}
   
   wget http://mumps.enseeiht.fr/MUMPS_${mumpsVER}.tar.gz
   mkdir -p ${buildDIR}/MUMPS
   tar -xf MUMPS_${mumpsVER}.tar.gz -C MUMPS/ --strip-components 1
   
   cd ${buildDIR}/MUMPS
   
   sed -e "/^#.*LMETISDIR *=/s/^#//" \
       -e "/^#.*IMETIS *=/s/^#//" \
       -e "s|\(LMETISDIR *=\).*|\1 ${cp2k_toolchainDIR}/install/parmetis-${parmetisVER}/lib|g" \
       -e "s|\(IMETIS *=\).*|\1 -I${cp2k_toolchainDIR}/install/parmetis-${parmetisVER}/include|g" \
       -e "/^#.*-lparmetis/s/^#//" \
       -e "s|\(CC *=\).*|\1 gcc|g" \
       -e "s|\(FC *=\).*|\1 mpif90|g" \
       -e "s|\(FL *=\).*|\1 mpif90|g" \
       -e "s|\(LAPACK *=\).*|\1 -L/${cp2k_toolchainDIR}/install/lapack-${reflapackVER}/lib -llapack|g" \
       -e "s|\(SCALAP *=\).*|\1 -L/${cp2k_toolchainDIR}/install/lapack-${scalapackVER}/lib -lscalapack|g" \
       -e "s|\(INCPAR *=\).*|\1 -I/${cp2k_toolchainDIR}/install/mpich-${mpichVER}/include|g" \
       -e "s|\(LIBPAR *=\).*|\1 \$(SCALAP) \$(LAPACK) -L/${cp2k_toolchainDIR}/install/mpich-${mpichVER}/lib -lmpi -lgfortran|g" \
       -e "s|\(LIBBLAS *=\).*|\1 -L/${cp2k_toolchainDIR}/install/lapack-${scalapackVER}/lib -lblas|g" \
       -e "s|\(OPTF *=\).*|\1 -O -w|g" \
       -e "s|\(OPTC *=\).*|\1 -O3 -w|g" \
          Make.inc/Makefile.inc.generic > Makefile.inc
   
   make alllib > install.log
   cp -r ${buildDIR}/MUMPS/lib/ ${buildDIR}/MUMPS/include/ ${libsDIR}/MUMPS

   # define include and lib paths for .mk file (STEP 3)
   inc_mumps="-I\$(LIB_TOP)/MUMPS/include/"
   lib_mumps="-L\$(LIB_TOP)/MUMPS/lib -lzmumps -ldmumps -lmumps_common -lpord -lscalapack"
fi

# MAGMA -------------------------------------------------------------------
if [ ! -z ${magmaVER} ]; then
   echo "installing MAGMA =========================================="
   mkdir -p ${libsDIR}/magma
   cd ${buildDIR}
   
   wget http://icl.cs.utk.edu/projectsfiles/magma/downloads/magma-${magmaVER}.tar.gz
   mkdir -p ${buildDIR}/magma
   tar -xf magma-${magmaVER}.tar.gz -C magma/ --strip-components 1
   
   cd ${buildDIR}/magma
   
   cp ${buildDIR}/magma/make.inc-examples/make.inc.openblas ${buildDIR}/magma/
   
   sed -e "/^#.*OPENBLASDIR *?=/s/^#//" \
       -e "/^#.*CUDADIR *?=/s/^#//" \
       -e "s|\(OPENBLASDIR *?=\).*|\1 ${cp2k_toolchainDIR}/install/openblas-${openblasVER}|g" \
       -e "s|\(CUDADIR *?=\).*|\1 ${cudaDIR}|g" \
       -e "s|\(NVCC *=\).*|\1 ${cudaDIR}/bin/nvcc -Xcompiler=--std=gnu++98 -D__GNUC__=4 -D__GNUC_MINOR__=9|g" \
       -e "/^NVCCFLAGS/ s/$/ -w/" \
       -e "/^-.*include/s/^-/#-/" \
          make.inc.openblas > make.inc
   
   make install prefix=${libsDIR}/magma > install.log

   # define include and lib paths for .mk file (STEP 3)
   nvcc_path="${cudaDIR}/bin/nvcc"
   nvcc_flags="-Xcompiler=--std=gnu++98 -D__GNUC__=4 -D__GNUC_MINOR__=9 -w"
   inc_cuda="-I${cudaDIR}/include/"
   lib_cuda="-L${cudaDIR}/lib64/ -lcudart -lcublas -lcusparse -lblas"
   inc_magma="-I\$(LIB_TOP)/magma/include/"
   lib_magma="-L\$(LIB_TOP)/magma/lib -lmagma"
fi

# PARDISO -----------------------------------------------------------------
# Downloading PARDISO needs registration on their website. Therefore, this 
# script does not automatically install PARDISO.  

# define include and lib paths for .mk file (STEP 3)
if [ ! -z ${pardisoDIR} ]; then
   pardiso_file=${pardisoSOFILE/"lib"/ -l}
   lib_pardiso="-L${pardisoDIR}/${pardiso_file/".so"/}"
fi

echo "Done! ====================================================="
echo " "

# STEP 3: ******************************************************************************************

# generate a .mk file:
# -------------------------------------------------------------------------
if [ "${generate_makefile}" = "yes" ] ; then
   echo "generate a .mk file ======================================="
   cd ${omenDIR}
   
   topdir=${installDIR}
   libtop="\$(TOP_DIR)/libs"
   
   sed -e "s|\(source \).*|\1 ${cp2k_toolchainDIR}/install/setup|g" \
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
   # compile CP2K
   echo "compiling cp2k with target ${cp2k_target} libcp2k ========="
   cd ${cp2kDIR}/cp2k/src
   ln -sf ../makefiles/Makefile .
   make -j ARCH=local VERSION=${cp2k_target} 
   make -j ARCH=local VERSION=${cp2k_target} libcp2k
   
   # compile OMEN
   echo "compiling OMEN ============================================"
   cd ${omenDIR}/src
   ./configure --with-arch=${machine}.mk ${omenCONFIGURE_ARGS}
   make -j 
fi

cd ${omenDIR}


