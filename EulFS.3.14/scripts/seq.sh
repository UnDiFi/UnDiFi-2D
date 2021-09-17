cd $1
echo "**********************************************************"
echo Running testcase:
echo $1
echo "**********************************************************"
pwd
cp .petscrc ${HOME}
echo "********************"
echo running ${CURRENT_VERSION} on PETSC_ARCH = ${PETSC_ARCH}
echo "********************"
${CURRENT_VERSION} -itmax ${ITMAX} -linearization ${SNES} -data_dir ${DATA_DIR} > out.${PETSC_ARCH}
cut -c5-9 -c32- convhst.l2 > convhst.${PETSC_ARCH} 
