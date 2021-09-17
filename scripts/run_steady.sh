#!/bin/bash

# ----------------------------------
# Colors
# ----------------------------------
NOCOLOR='\033[0m'
RED='\033[0;31m'
GREEN='\033[0;32m'
ORANGE='\033[0;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
LIGHTGRAY='\033[0;37m'
DARKGRAY='\033[1;30m'
LIGHTRED='\033[1;31m'
LIGHTGREEN='\033[1;32m'
YELLOW='\033[1;33m'
LIGHTBLUE='\033[1;34m'
LIGHTPURPLE='\033[1;35m'
LIGHTCYAN='\033[1;36m'
WHITE='\033[1;37m'
bold=$(tput bold)
normal=$(tput sgr0)

# ---------------------------------------------------------------------------------------------------------------------------------
helpFunction()
{
   echo ""
   echo -e "${bold} Usage: $0 ${normal}-s ${RED}Solver${NOCOLOR} -m ${BLUE}Mode${NOCOLOR} -f ${GREEN}Flow${NOCOLOR} ${normal}"
   echo -e "${RED}\t-s Solver: neo or eulfs ${NOCOLOR}"
   echo -e "${BLUE}\t-m Mode: capturing or fitting ${NOCOLOR}"
   echo -e "${GREEN}\t-f Flow: steady or unsteady ${NOCOLOR}"
   echo ""
   echo ""
   exit 1 # Exit script after printing help
}
# ---------------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------------
while getopts ":s:m:f:" arg
do
  case $arg in
    s) Solver=$(echo $OPTARG | tr '[:upper:]' '[:lower:]');;
    m) Mode=$(echo $OPTARG   | tr '[:upper:]' '[:lower:]');;
    f) Flow=$(echo $OPTARG   | tr '[:upper:]' '[:lower:]');;
    ?) helpFunction;; # Print helpFunction in case parameter is non-existent
  esac
done
echo -e "\n$Solver $Mode $Flow\n"
# ---------------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------------
# Print helpFunction in case parameters are empty
if [ -z "$Solver" ] || [ -z "$Mode" ] || [ -z "$Flow" ]
then
   echo -e "${YELLOW} Some or all of the parameters are empty! ${NOCOLOR}";
   helpFunction
fi
# ---------------------------------------------------------------------------------------------------------------------------------

testname=`basename $PWD`
echo $testname

# TODO: this can be made input arg as well
iters=501
echo $iters

# ---------------------------------------------------------------------------------------------------------------------------------
# Redirection
if [ "$Solver" = "neo" ] && [ "$Mode" = "capturing" ] && [ "$Flow" = "steady" ]; then
# Generate input files for NEO
  ../../bin/triangle2grd << !
  na00.1   # fname
!
  ../../bin/na2vvvv << !
  na00.1   # fname
!

# modify NEO's input file for SC simulation
  cp NEO_data/textinput/inputfile-exp.txt.SC NEO_data/textinput/inputfile-exp.txt

# Write empty NEO_data/input/vel.dat
  echo "" > NEO_data/input/vel.dat

# Launch NEO
  ../../bin/CRD_euler
# ---------------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------------
elif [ "$Solver" = "neo" ] && [ "$Mode" = "capturing" ] && [ "$Flow" = "unsteady" ]; then
# Generate input files for NEO
  ../../bin/triangle2grd << !
  na00.1   # fname
!
  ../../bin/na2vvvv << !
  na00.1   # fname
!

# modify NEO's input file for SC simulation
  cp NEO_data/textinput/inputfile-exp.txt.SC NEO_data/textinput/inputfile-exp.txt

# Write empty NEO_data/input/vel.dat
  echo "" > NEO_data/input/vel.dat

# Launch NEO
  ../../bin/CRD_euler
# ---------------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------------
elif [ "$Solver" = "eulfs" ] && [ "$Mode" = "capturing" ] && [ "$Flow" = "steady" ]; then
# Generate file00[123].dat
  ../../bin/triangle2dat-NEW-x86_64 << !
  na00.1   # fname
  n        # Are there any periodic surfaces?
!
# Be sure that no any .petsrc is present in the ~/home/
  rm ~/home/.petscrc

# Launch EulFS
  ../../bin/EulFS_x86_64 -itmax $iters
# ---------------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------------
elif [ "$Solver" = "eulfs" ] && [ "$Mode" = "capturing" ] && [ "$Flow" = "unsteady" ]; then
# Generate file00[123].dat
  ../../bin/triangle2dat-NEW-x86_64 << !
  na00.1   # fname
  n        # Are there any periodic surfaces?
!
# Be sure that no any .petsrc is present in the ~/home/
  rm ~/home/.petscrc

# Launch EulFS
  ../../bin/EulFS_x86_64 -itmax $iters
# ---------------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------------
elif [ "$Solver" = "neo" ] && [ "$Mode" = "fitting" ] && [ "$Flow" = "steady" ]; then
# modify NEO's input file for SF simulation
  cp NEO_data/textinput/inputfile-exp.txt.SF NEO_data/textinput/inputfile-exp.txt

# Write empty NEO_data/input/vel.dat
  echo "" > NEO_data/input/vel.dat

#                            nbegin, nsteps, eulfs, steady, testcase,   logfile
  ../../bin/UnDiFi-2D_x86_64 0       $iters  false  true    $testname | tee run.log
# ---------------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------------
elif [ "$Solver" = "neo" ] && [ "$Mode" = "fitting" ] && [ "$Flow" = "unsteady" ]; then
# modify NEO's input file for SF simulation
  cp NEO_data/textinput/inputfile-exp.txt.SF NEO_data/textinput/inputfile-exp.txt

  ../../bin/UnDiFi-2D_x86_64 0       $iters  false  false   $testname | tee run.log
# ---------------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------------
elif [ "$Solver" = "eulfs" ] && [ "$Mode" = "fitting" ] && [ "$Flow" = "steady" ]; then
  ../../bin/UnDiFi-2D_x86_64 0       $iters  true   true    $testname | tee run.log
# ---------------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------------
elif [ "$Solver" = "eulfs" ] && [ "$Mode" = "fitting" ] && [ "$Flow" = "unsteady" ]; then
  ../../bin/UnDiFi-2D_x86_64 0       $iters  true   false   $testname | tee run.log
# ---------------------------------------------------------------------------------------------------------------------------------


fi
# ---------------------------------------------------------------------------------------------------------------------------------
