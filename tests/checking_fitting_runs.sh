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


# Print helpFunction in case parameters are empty
if [ -z "$Solver" ]
then
   echo -e "${YELLOW} Some or all of the parameters are empty! ${NOCOLOR}";
   helpFunction
fi




echo '******************************'
echo 'Checking  the UnFiDi-2D runs '


if [ "$Solver" = "eulfs" ]
then
  echo '    with eulfs'



#  for dir in `ls -d [C][i]*1`
  for dir in `ls -d [QCNMRS]*`

  do
    cd $dir 
#    if   [-e step00501]   &&   [-e convergenza.dat] ; 
    if   [ -e step00501/file010.dat ] &&  [ -e convergenza.dat ]   
    then
     if [ -e checksum_fitting_eulfs ]
     then
      tail -n 1 convergenza.dat > out.log
      diff out.log checksum_fitting_eulfs > out
      if [   $? -eq 0  ]  
      then 
        echo test $dir ok 
      else
       echo test $dir failed 
      fi
      rm out
      rm out.log
     else
      echo test $dir checksum missing
     fi
   else
        echo test $dir not completed or executed 
   fi
   cd ..
  done

  echo '*********************'

  elif [ "$Solver" = "neo" ]
   then
   echo '    with neo'


   for dir in `ls -d [QCNMRS]*`
   do
    cd $dir
    if  [ -e step00501/vvvv.dat ]  &&  [ -e residual_norm.dat ] 
    then

     if [ -e checksum_fitting_neo ]
     then
      tail -n 1 residual_norm.dat > out.log
      diff out.log checksum_fitting_neo > out
      if [   $? -eq 0  ]
      then
        echo test $dir ok
      else
       echo test $dir failed
      fi
      rm out
      rm out.log
     else
      echo test $dir checksum missing 
     fi
   else
        echo test $dir not completed or executed 
   fi
   cd ..
  done
  echo '*********************'

   fi

