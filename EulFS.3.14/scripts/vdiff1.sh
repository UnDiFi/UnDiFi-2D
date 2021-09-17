cd ${TCASE}
echo "**********************************************************"
echo Running testcase:
echo ${TCASE}
echo "**********************************************************"
pwd
cp .petscrc ${HOME}
echo "********************"
echo running ${CURRENT_VERSION}
echo "********************"
${CURRENT_VERSION} -itmax ${ITMAX} -linearization ${SNES} -data_dir ${DATA_DIR} > out10.10 
cut -c7- convhst.l2 > tmp.new 
echo "********************"
echo running ${OLD_VERSION}
echo "********************"
${OLD_VERSION} -itmax ${ITMAX} -linearization ${SNES} -data_dir ${DATA_DIR} > out10.9
echo $1 >> ${HOME}/differences
cut -c7- convhst.l2 > tmp.old 
diff tmp.new tmp.old >> ${HOME}/differences
