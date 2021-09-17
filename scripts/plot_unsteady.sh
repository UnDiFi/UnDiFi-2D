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
   echo -e "${bold} Usage: $0 ${normal}-s ${RED}Solver${NOCOLOR} -m ${BLUE}Mode${NOCOLOR} "
   echo -e "${RED}\t-s Solver: neo or eulfs ${NOCOLOR}"
   echo -e "${BLUE}\t-m Mode: capturing or fitting ${NOCOLOR}"
   echo ""
   echo ""
   exit 1 # Exit script after printing help
}


while getopts "s:m:" arg
do
  case $arg in
    s) Solver=$(echo $OPTARG | tr '[:upper:]' '[:lower:]');;
    m) Mode=$(echo $OPTARG   | tr '[:upper:]' '[:lower:]');;
    ?) helpFunction;; # Print helpFunction in case parameter is non-existent
  esac
done
echo -e "\n$Solver $Mode \n"


# Print helpFunction in case parameters are empty
if [ -z "$Solver" ] || [ -z "$Mode" ]
then
   echo -e "${YELLOW} Some or all of the parameters are empty! ${NOCOLOR}";
   helpFunction
fi

fsplot=../../bin/fsplot-$HOSTTYPE

if [ "$Mode" = "fitting" ]; then
################################

for dir in `ls -d step?????`

do
   echo $dir
   cp inp $dir
   cd $dir
   if [ "$Solver" = "eulfs" ]
   ##########################
   then
      ../../../bin/fsplot-$HOSTTYPE
      preplot file012.dat $dir-eulfs.plt
   ######## FOR PARAVIEW USERS #########
   # ../../../bin/dat2paraview
   # cp paraviewE.dat ../$dir-eulfs.dat
   #####################################
   elif [ "$Solver" = "neo" ]
   ##########################
   then
     preplot vvvv.dat $dir-neo.plt
   ######## FOR PARAVIEW USERS #########
   # ../../../bin/dat2paraview
   # cp paraviewN.dat ../$dir-neo.dat
   #####################################
   fi
   mv $dir*.plt ..
   cd ..
done

elif [ "$Mode" = "capturing" ]; then
####################################

   if [ "$Solver" = "eulfs" ]
   ##########################
   then
      for i in `ls file010_*.dat | cut -c9-14`
        do
            f2=file010_$i.dat
            f3=eulfs_$i.plt
            echo processing $f2 $f3
            ln -sf $f2 file003.dat
            if ($fsplot >& fsplot.log) then true; else echo "problem with " $fsplot; exit 1; fi
            if (preplot file012.dat $f3 >& preplot.log) then true; else echo "problem with " preplot; exit 1; fi
        done
        rm file003.dat
   elif [ "$Solver" = "neo" ]
   ##########################
   then
        # In this case, the .dat files can be directly opened with Tecplot
        # but it is convenient to convert them into binary format so that
        # can be opened also with Paraview or Visit.
        cd NEO_data/output
        for i in `ls vvvv*.dat | cut -c5-10`
            do
                finp=vvvv$i.dat
                fout=neo_$i.plt
                echo processing $finp $fout
                preplot $finp $fout
                cp $fout ../../.
            done
        cd ../..

   fi
fi
