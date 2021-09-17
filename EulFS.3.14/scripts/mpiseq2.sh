cd $1
echo "**********************************************************"
echo Running testcase:
echo $1
echo "**********************************************************"
pwd
cp .petscrc ${HOME}
echo "********************"
echo running sequential version
echo "********************"
${SEQ_CODE} -itmax ${ITMAX} -linearization ${SNES} -data_dir ${DATA_DIR} > out.seq
cut -c5-9 -c32- convhst.l2 > tmp.seq 
echo "********************"
echo running parallel version
echo "********************"
${MPIRUN} -np ${NUM_PROCS} ${MPI_CODE} -itmax ${ITMAX} -linearization ${SNES} -data_dir ${DATA_DIR}np${NUM_PROCS}/ > out.mpi
echo $1 >> ${HOME}/differences
cut -c5-9 -c32- convhst.l2 > tmp.mpi 
diff tmp.seq tmp.mpi >> ${HOME}/differences
