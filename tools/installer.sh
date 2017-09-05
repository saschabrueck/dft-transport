#!/bin/bash -e

# ========================================================================
# \brief    script for installing CP2K-OMEN 
# \history  created: 22-07-2017 
# \author   Mohammad Hossein Bani-Hashemian 
# ========================================================================

TOPDIR=/scratch/seyedb
cp2kDIR=${TOPDIR}/cp2k
omenDIR=${TOPDIR}/omen

installDIR=${omenDIR}/install
buildDIR=${installDIR}/build
libsDIR=${installDIR}/libs

machine="mont-fort1"

# with or without CUDA (needed for SplitSolve and MAGMA)  
# Note: MAGMA v >2.0 requires CUDA v >5.0
withCUDA="no"
cudaDIR="/usr/local/cuda-7.5"

# specify the cp2k VERSION for the last step (compilation):
cp2k_target="popt"

# STEP 1: ******************************************************************************************

# install CP2K: ===========================================================
# Note that, to have a well-matched installation, this script uses the CP2K
# toolchain script to install CP2K with PEXSI, SuperLU_DIST, and ParMETIS. 
# -------------------------------------------------------------------------
echo "installing CP2K ==========================================="
cd ${TOPDIR} 

# install either the currnet trunk or the latest release  
svn checkout http://svn.code.sf.net/p/cp2k/code/trunk cp2k
#svn checkout http://svn.code.sf.net/p/cp2k/code/branches/cp2k-4_1-branch cp2k

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

## STEP 2: ******************************************************************************************
#
## install OMEN solvers: ===================================================
## specify the versions here 
hypreVER="2.11.2"
SuiteSparseVER="4.5.5"
qhullVER="2015-src-7.2.0"
mumpsVER="5.1.1"
magmaVER="2.2.0"

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

# copy PEXSI headers that are not copied via CP2K toolchain script 
# -------------------------------------------------------------------------
cp -r ${cp2k_toolchainDIR}/build/pexsi_v${pexsiVER}/include/pexsi/ ${cp2k_toolchainDIR}/install/pexsi-${pexsiVER}/include/

# install hypre:
# -------------------------------------------------------------------------
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

# install SuiteSparse:
# -------------------------------------------------------------------------
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

# install qhull:
# -------------------------------------------------------------------------
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

# install MUMPS:
# -------------------------------------------------------------------------
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

# install MAGMA:
# -------------------------------------------------------------------------
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

# install PARDISO:
# -------------------------------------------------------------------------
# Downloading PARDISO needs registration on their website. Therefore, this 
# script does not automatically install PARDISO. If you have the .so file, 
# specify the name of the shared object file and the path to the files here, 
# otherwise leave the following two variables unset:

#pardisoSOFILE="libpardiso500-MPI-GNU472-X86-64.so"
#pardisoDIR=${libsDIR}/pardiso/lib
pardisoSOFILE=
pardisoDIR=

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
    -e "s|\(INCSSPARSE *=\).*|\1 \$(LIB_TOP)/SuiteSparse/include/|g" \
    -e "s|\(INCMUMPS *=\).*|\1 \$(LIB_TOP)/MUMPS/include/|g" \
    -e "s|\(INCHYPRE *=\).*|\1 \$(LIB_TOP)/hypre/include/|g" \
    -e "s|\(INCQHULL *=\).*|\1 \$(LIB_TOP)/qhull/include/libqhull/|g" \
    -e "s|\(INCSLUDIST *=\).*|\1 \$(TOOLCHAIN)/install/superlu_dist-${superluVER}/include/|g" \
    -e "s|\(INCPEXSI *=\).*|\1 \$(TOOLCHAIN)/install/pexsi-${pexsiVER}/include/|g" \
        ${omenDIR}/makefiles/${machine}.mk

sed -i \
    -e "s|\(LIBSSPARSE *=\).*|\1 \$(LIB_TOP)/SuiteSparse/static -lumfpack -lamd -lcholmod -lcolamd -lccolamd -lcamd -lsuitesparseconfig |g" \
    -e "s|\(LIBMUMPS *=\).*|\1 \$(LIB_TOP)/MUMPS/lib -lzmumps -ldmumps -lmumps_common -lpord -lscalapack |g" \
    -e "s|\(LIBHYPRE *=\).*|\1 \$(LIB_TOP)/hypre/lib -lHYPRE |g" \
    -e "s|\(LIBQHULL *=\).*|\1 \$(LIB_TOP)/qhull/lib -lqhullstatic |g" \
    -e "s|\(LIBPARMETIS *=\).*|\1 \$(TOOLCHAIN)/install/parmetis-${parmetisVER}/lib -lparmetis -lmetis |g" \
    -e "s|\(LIBSLUDIST *=\).*|\1 \$(TOOLCHAIN)/install/superlu_dist-${superluVER}/lib -lsuperlu_dist |g" \
    -e "s|\(LIBPEXSI *=\).*|\1 \$(TOOLCHAIN)/install/pexsi-${pexsiVER}/lib -lpexsi |g" \
    -e "s|\(LIBCP2K *=\).*|\1 ${cp2kDIR}/cp2k/lib/local/${cp2k_target} -lcp2k -lxsmmf -lxsmm -lderiv -lint -lxcf90 -lxc -lfftw3 |g" \
        ${omenDIR}/makefiles/${machine}.mk

if [ ! -z ${pardisoDIR} ]; then
   pardiso_lib=${pardisoSOFILE/"lib"/ -l}
   sed -i \
       -e "s|\(LIBPARDISO *=\).*|\1 ${pardisoDIR}/${pardiso_lib/".so"/ } |g" \
           ${omenDIR}/makefiles/${machine}.mk
fi

if [ "${withCUDA}" = "yes" ] ; then
   sed -i \
       -e "s|\(NVCC *=\).*|\1 ${cudaDIR}/bin/nvcc|g" \
       -e "s|\(NVCCFLAGS *=\).*|\1 -Xcompiler=--std=gnu++98 -D__GNUC__=4 -D__GNUC_MINOR__=9 -w |g" \
       -e "s|\(INCCUDA *=\).*|\1 ${cudaDIR}/include/ |g" \
       -e "s|\(LIBCUDA *=\).*|\1 ${cudaDIR}/lib64/ -lcudart -lcublas -lcusparse -lblas |g" \
       -e "s|\(INCMAGMA *=\).*|\1 \$(LIB_TOP)/magma/include/|g" \
       -e "s|\(LIBMAGMA *=\).*|\1 \$(LIB_TOP)/magma/lib -lmagma |g" \
           ${omenDIR}/makefiles/${machine}.mk
fi

echo "Done! ====================================================="
echo " "

# STEP 4: ******************************************************************************************
# Feel free to comment out the following lines and do this step manually.
# 
# compile CP2K-OMEN:
# -------------------------------------------------------------------------
# example: compile CP2K popt
echo "compiling cp2k with target ${cp2k_target} libcp2k ========="
cd ${cp2kDIR}/cp2k/src
ln -sf ../makefiles/Makefile .
make -j ARCH=local VERSION=${cp2k_target} 
make -j ARCH=local VERSION=${cp2k_target} libcp2k

# example: compile OMEN with PEXSI and MUMPS
echo "compiling OMEN ============================================"
cd ${omenDIR}/src
./configure --with-arch=${machine}.mk --with-pexsi --with-mumps
make -j 

cd ${omenDIR}


