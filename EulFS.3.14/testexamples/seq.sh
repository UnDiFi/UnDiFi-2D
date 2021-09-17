#include ${PETSC_DIR}/bmake/${PETSC_ARCH}/packages
#
# SEQ_CODE = ${HOME}/bin/${PETSC_ARCH}/eulfs
# if the full path is set for the seq code then the output is sent
# to ${HOME}/bin/${PETSC_ARCH} !!! at least on alpha
#
#SEQ_CODE = eulfs11.8
#MPI_CODE = ${HOME}/bin/${PETSC_ARCH}/peulfs
#
#
#
#
runseq:
	@cp .petscrc ${HOME}
	@echo "running code " ${SEQ_CODE} " with runtime options :" ${FSPL_OPTS}
	@${SEQ_CODE} ${FSPL_OPTS} -data_dir ${DATA_DIR} > out.${PETSC_ARCH}; \
	cut -c1-10,33- convhst.l2 > convhst.seq.${PETSC_ARCH}
diff:
	@if (diff convhst.seq.${ARCH1} convhst.seq.${ARCH2} >> ${LOGFILE} ) \
	then echo "No differences have been found while running" ${TCASE}; else echo "Possible problem with " ${TCASE} ", check for differences"; echo; fi
vdiff:
	@cp .petscrc ${HOME}; rm -f convhst.?.${PETSC_ARCH}
	@echo >> ${LOGFILE}
	@echo
	@echo "running testcase : " ${TCASE} " with: " ${CODE1} >> ${LOGFILE}
	@echo "running testcase with: " ${CODE1} ${FSPL_OPTS}
	@for i in file010.dat file015.dat fold010.dat; do \
	if ([ -e $$i ]) then rm -v $$i; fi; \
	if ([ -e $$i.BAK ]) then cp -vp $$i.BAK $$i; else echo 'There is no' $$i ' to be copied'; fi \
	done
	@if (${CODE1} ${FSPL_OPTS} -data_dir ${DATA_DIR} > out.1.${PETSC_ARCH}) \
	then (cut -c1-10,33- convhst.l2 > convhst.1.${PETSC_ARCH}) \
	else echo; echo ${CODE1} has failed on ${TCASE}; echo; echo "error:" ${CODE1} has failed on ${TCASE} >> ${LOGFILE}; fi
	@echo
	@echo "running testcase with: " ${CODE2} ${FSPL_OPTS}
	@echo "running testcase : " ${TCASE} " with: " ${CODE2} >> ${LOGFILE}
	@for i in file010.dat file015.dat fold010.dat; do \
	if ([ -e $$i ]) then rm -v $$i; fi; \
	if ([ -e $$i.BAK ]) then cp -vp $$i.BAK $$i; else echo 'There is no' $$i 'to be copied'; fi \
	done
	@if (${CODE2} ${FSPL_OPTS} -data_dir ${DATA_DIR} > out.2.${PETSC_ARCH}) \
	then (cut -c1-10,33- convhst.l2 > convhst.2.${PETSC_ARCH}) \
	else echo; echo ${CODE2} has failed on ${TCASE}; echo; echo "error:" ${CODE2} has failed on ${TCASE} >> ${LOGFILE};  fi
	@if (diff convhst.1.${PETSC_ARCH} convhst.2.${PETSC_ARCH} > scratchfile ) \
	then echo "No differences have been found while running" ${TCASE}; echo "No differences have been found while running" ${TCASE} >> ${LOGFILE}; \
	else echo "Possible problem with " ${TCASE} ", check for differences"; echo "Possible problem with " ${TCASE} ", differences are listed below" >> ${LOGFILE}; cat scratchfile >> ${LOGFILE}; fi

pdiff:
	@cp .petscrc ${HOME}
	@echo "running testcase : " ${TCASE} " with: " ${SEQ_CODE} >> ${LOGFILE}
	@echo "running testcase with: " ${SEQ_CODE}
	@${SEQ_CODE} ${FSPL_OPTS} -data_dir ${DATA_DIR} > out1.${PETSC_ARCH}; \
	cut -c1-10,33- convhst.l2 > convhst.1.${PETSC_ARCH}
	@echo "running testcase with: " ${MPI_CODE} "on " ${NPES} "processors"
	@echo "running testcase : " ${TCASE} " with: " ${MPI_CODE} "on " ${NPES} "processors" >> ${LOGFILE}
	@mpirun \-np ${NPES} ${MPI_CODE} ${FSPL_OPTS} -data_dir ${DATA_DIR}/np${NPES} > out.${PETSC_ARCH}; \
	cut -c1-10,33- convhst.l2 > convhst.${NPES}.${PETSC_ARCH}
	@if (diff convhst.1.${PETSC_ARCH} convhst.${NPES}.${PETSC_ARCH} >> ${LOGFILE} ) \
	then true; else echo "Possible problem, check for differences"; fi
clean:	
	@rm -f convhst.* output.??? out* fort.* log* slip-free* no-slip* timing.??? mytime.??? fspl.out
