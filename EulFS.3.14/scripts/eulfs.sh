#
# A sh script to compare different versions of the code
# differences are sent to ${HOME}/differences
#
#
#
export ITMAX=5
export SNES=picard
export CURRENT_VERSION=eulfs11
export OLD_VERSION=eulfs3
#
echo running on `uname -a` `date` > ${HOME}/differences
echo checking differences between ${CURRENT_VERSION} and ${OLD_VERSION} > ${HOME}/differences
#
# this is used to test scalar problems (advection-diffusion)
#
export DATA_DIR=${HOME}/grids/2D/scalar/advdiff/21x21/
export NUM_PROCS=4
#
export TCASE=${HOME}/testcases/2D/scalar/tcase3
sh vdiff1.sh ${TCASE}
#
# this is used to test scalar problems (advection)
#
export DATA_DIR=${HOME}/grids/2D/scalar/source3/41x41/
export NUM_PROCS=4
#
export TCASE=${HOME}/testcases/2D/scalar/tcase6
sh vdiff1.sh ${TCASE}
