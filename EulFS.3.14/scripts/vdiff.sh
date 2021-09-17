cd $1
echo "**********************************************************"
echo Running testcase:
echo $1
echo "**********************************************************"
pwd
cp .petscrc ${HOME}
echo "********************"
echo running ${CURRENT_VERSION}
echo "********************"
${CURRENT_VERSION} -itmax ${ITMAX} -linearization ${SNES} -data_dir ${DATA_DIR} > out10.10 
cut -c5-9 -c32- convhst.l2 > tmp.new 
echo "********************"
echo running ${OLD_VERSION}
echo "********************"
${OLD_VERSION} -itmax ${ITMAX} -linearization ${SNES} -data_dir ${DATA_DIR} > out10.9
echo $1 >> ${HOME}/differences
cut -c5-9 -c32- convhst.l2 > tmp.old 
diff tmp.new tmp.old >> ${HOME}/differences
