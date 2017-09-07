#!/bin/bash -e

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

# * set the following paths *******
TOPDIR=/scratch/seyedb
cp2kDIR=${TOPDIR}/cp2k
omenDIR=${TOPDIR}/omen

installDIR=${omenDIR}/install
buildDIR=${installDIR}/build
libsDIR=${installDIR}/libs

# * preparation for STEP 1 ********
# install either the current truck version of CP2K or a release version (>=4.1)
cp2kTRUNK="yes"
cp2kRELEASE="no"
cp2kRELEASEVER="4_1"

# * preparation for STEP 2 ********
# with or without CUDA (needed for SplitSolve and MAGMA):
# Note that MAGMA v >2.0 requires CUDA v >5.0
withCUDA="no"
cudaDIR="/usr/local/cuda-7.5"

# Specify the versions of the solvers that you would like to install:
hypreVER="2.11.2"
SuiteSparseVER="4.5.5"
qhullVER="2015-src-7.2.0"
mumpsVER="5.1.1"
magmaVER="2.2.0"

# If you have downloaded the shared library file of PARDISO, specify the name 
# of the file and the path to it here, otherwise, leave pardisoDIR unset:
pardisoSOFILE="libpardiso500-MPI-GNU472-X86-64.so"
#pardisoDIR=${libsDIR}/pardiso/lib
pardisoDIR=

# * preparation for STEP 3 ********
# specify the name of your machine:
machine="mont-fort1"

# * preparation for STEP 4 ********
# specify the cp2k VERSION for compilation:
cp2k_target="popt"
# set the OMEN's configure flags, i.e. solvers to be enabled:
omenCONFIGURE_ARGS="--with-pexsi --with-mumps"

# * you should not need to change anything below this line, unless you would like to 
#   edit some steps or skip them by commenting them out.

# STEP 1: ******************************************************************************************

# install CP2K: ===========================================================
# Here the script uses the CP2K toolchain script to install CP2K together
# with PEXSI, SuperLU_DIST, and ParMETIS that all are used by OMEN as well.
# -------------------------------------------------------------------------
echo "installing CP2K ==========================================="
cd ${TOPDIR} 

if [ "${cp2kTRUNK}" = "yes" ] ; then
   svn checkout http://svn.code.sf.net/p/cp2k/code/trunk cp2k
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

source ${cp2k_toolchainDIR}/scripts/common_vars.sh

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

# HYPRE -------------------------------------------------------------------
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

# SuiteSparse -------------------------------------------------------------
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

# Qhull -------------------------------------------------------------------
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

# MUMPS -------------------------------------------------------------------
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

# MAGMA -------------------------------------------------------------------
if [ "${withCUDA}" = "yes" ] ; then
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
else
   echo "MAGMA will not be installed. In order to install MAGMA, please set withCUDA="yes"."
fi

# PARDISO -----------------------------------------------------------------
# Downloading PARDISO needs registration on their website. Therefore, this 
# script does not automatically install PARDISO.  

echo "Done! ====================================================="
echo " "

# STEP 3: ******************************************************************************************

# generate a .mk file:
# -------------------------------------------------------------------------
echo "generate a .mk file ======================================="
cd ${omenDIR}

sed -e "s|\(source \).*|\1 ${cp2k_toolchainDIR}/install/setup|g" \
    -e "s|\(TOP_DIR *=\).*|\1 ${installDIR}|g" \
    -e "s|\(LIB_TOP *=\).*|\1 \$(TOP_DIR)/libs|g" \
    -e "s|\(TOOLCHAIN *=\).*|\1 ${cp2k_toolchainDIR}|g" \
        ${omenDIR}/makefiles/arch.tmpl > ${omenDIR}/makefiles/${machine}.mk

sed -i \
    -e "s|\(INCSSPARSE *=\).*|\1 -I\$(LIB_TOP)/SuiteSparse/include/|g" \
    -e "s|\(INCMUMPS *=\).*|\1 -I\$(LIB_TOP)/MUMPS/include/|g" \
    -e "s|\(INCHYPRE *=\).*|\1 -I\$(LIB_TOP)/hypre/include/|g" \
    -e "s|\(INCQHULL *=\).*|\1 -I\$(LIB_TOP)/qhull/include/libqhull/|g" \
    -e "s|\(INCSLUDIST *=\).*|\1 -I\$(TOOLCHAIN)/install/superlu_dist-${superluVER}/include/|g" \
    -e "s|\(INCPEXSI *=\).*|\1 -I\$(TOOLCHAIN)/install/pexsi-${pexsiVER}/include/|g" \
        ${omenDIR}/makefiles/${machine}.mk

sed -i \
    -e "s|\(LIBSSPARSE *=\).*|\1 -L\$(LIB_TOP)/SuiteSparse/static -lumfpack -lamd -lcholmod -lcolamd -lccolamd -lcamd -lsuitesparseconfig |g" \
    -e "s|\(LIBMUMPS *=\).*|\1 -L\$(LIB_TOP)/MUMPS/lib -lzmumps -ldmumps -lmumps_common -lpord -lscalapack |g" \
    -e "s|\(LIBHYPRE *=\).*|\1 -L\$(LIB_TOP)/hypre/lib -lHYPRE |g" \
    -e "s|\(LIBQHULL *=\).*|\1 -L\$(LIB_TOP)/qhull/lib -lqhullstatic |g" \
    -e "s|\(LIBPARMETIS *=\).*|\1 -L\$(TOOLCHAIN)/install/parmetis-${parmetisVER}/lib -lparmetis -lmetis |g" \
    -e "s|\(LIBSLUDIST *=\).*|\1 -L\$(TOOLCHAIN)/install/superlu_dist-${superluVER}/lib -lsuperlu_dist |g" \
    -e "s|\(LIBPEXSI *=\).*|\1 -L\$(TOOLCHAIN)/install/pexsi-${pexsiVER}/lib -lpexsi |g" \
    -e "s|\(LIBCP2K *=\).*|\1 -L${cp2kDIR}/cp2k/lib/local/${cp2k_target} -lcp2k -lxsmmf -lxsmm -lderiv -lint -lxcf90 -lxc -lfftw3 |g" \
        ${omenDIR}/makefiles/${machine}.mk

if [ ! -z ${pardisoDIR} ]; then
   pardiso_lib=${pardisoSOFILE/"lib"/ -l}
   sed -i \
       -e "s|\(LIBPARDISO *=\).*|\1 -L${pardisoDIR}/${pardiso_lib/".so"/ } |g" \
           ${omenDIR}/makefiles/${machine}.mk
fi

if [ "${withCUDA}" = "yes" ] ; then
   sed -i \
       -e "s|\(NVCC *=\).*|\1 ${cudaDIR}/bin/nvcc|g" \
       -e "s|\(NVCCFLAGS *=\).*|\1 -Xcompiler=--std=gnu++98 -D__GNUC__=4 -D__GNUC_MINOR__=9 -w |g" \
       -e "s|\(INCCUDA *=\).*|\1 -I${cudaDIR}/include/ |g" \
       -e "s|\(LIBCUDA *=\).*|\1 -L${cudaDIR}/lib64/ -lcudart -lcublas -lcusparse -lblas |g" \
       -e "s|\(INCMAGMA *=\).*|\1 -I\$(LIB_TOP)/magma/include/|g" \
       -e "s|\(LIBMAGMA *=\).*|\1 -L\$(LIB_TOP)/magma/lib -lmagma |g" \
           ${omenDIR}/makefiles/${machine}.mk
fi

echo "Done! ====================================================="
echo " "

# STEP 4: ******************************************************************************************
# Feel free to comment out the following lines and do this step manually.
# 
# compile CP2K-OMEN:
# -------------------------------------------------------------------------
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

cd ${omenDIR}


