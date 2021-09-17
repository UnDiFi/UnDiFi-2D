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


helpFunction()
{
   echo ""
   echo -e "${bold} Usage: $0 ${normal}-s ${RED}Solver${NOCOLOR} "
   echo -e "${RED}\t-s Solver: neo or eulfs ${NOCOLOR}"
   echo ""
   echo ""
   exit 1 # Exit script after printing help
}


while getopts "s:" arg
do
  case $arg in
    s) Solver=$(echo $OPTARG | tr '[:upper:]' '[:lower:]');;
    ?) helpFunction;; # Print helpFunction in case parameter is non-existent
  esac
done
echo -e "\n$Solver \n"


# Print helpFunction in case parameters are empty
if [ -z "$Solver" ]
then
   echo -e "${YELLOW} Some or all of the parameters are empty! ${NOCOLOR}";
   helpFunction
fi



for dir in `ls -d step?????`

do
   echo $dir
   cp inp $dir
   cd $dir
   if [ "$Solver" = "eulfs" ]
   then
     ../../../bin/fsplot-$HOSTTYPE
      preplot file012.dat $dir-eulfs.plt
###### FOR PARAVIEW USERS #######
#      ../../../bin/dat2paraview
#        cp paraviewE.dat ../$dir-eulfs.dat  
##############################	  
   elif [ "$Solver" = "neo" ]
   then
      preplot vvvv.dat $dir-neo.plt
###### FOR PARAVIEW USERS #######
#      ../../../bin/dat2paraview
#        cp paraviewN.dat ../$dir-neo.dat  
##############################	
    fi

   mv $dir*.plt ..
   cd ..
done

