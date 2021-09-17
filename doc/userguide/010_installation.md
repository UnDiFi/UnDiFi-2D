\hypertarget{installation}{}

# Installation \label{chap:installation}
The following chapter describes the installation process on
a Linux machine. This also includes the installation of all
required prerequisites UnDiFi-2D has been developed on GNU/Linux
architectures. Other OS are not supported (and in general there is no
better alternative to GNU/Linux :-)

The provided codes (UnDiFi-2D, NEO, EulFS and auxiliary
packages) have been successfully compiled with the following compilers:

- GNU gfortran (version 4.7.0 or higher)
- Intel Fortran Compiler ifort (version 12.0 or higher)

Other compilers are not supported.

The codes are constituted by several modules, therefore, there are many
dependences on auxiliary sofwares (see Section \ref{sec:auxiliary}).
The easiest way to compile the code is to start with the provided
script **compile_all.sh**. It is then necessary that the system has
"Make" program (preferably GNU make http://www.gnu.org/software/make/).

## Prerequisites
UnDiFi-2D has been tested for various Linux distributions. This includes
Ubuntu 14.04 LTS, 16.04 LTS and 18.04 LTS, 20.04.2 LTS, OpenSUSE 42.1,
Leap 15.1, Leap 15.2, Tumbleweed.

The suggested packages in this section can be replaced by self compiled
versions. The required packages for the Ubuntu Linux distributions are
listed below. Under Ubuntu, they can be obtained using the apt environment:

    sudo apt-get install git

* git
* cmake
* make
* automake
* gfortran
* gcc/g++
* texlive-latex and texlive-lstaddons
* libX11 (see for xorg package)
* LibX11-dev (Ubuntu) LibX11-devel (Opensuse)
* pandoc (version <= 2.9) and pandoc-citeproc
* libkrb5-dev (Ubuntu) / krb5-devel (Opensuse)

Sometimes on some the sytems (based on OpenSuse) the program krb5-config is 
installed in a directory from which is not possible to run it. In this case is sufficient to add the path /usr/bin/mit/bin to the $PATH variable.

On some systems it may be necessary to increase the size of the stack
(part of the memory used to store information about active subroutines)
in order to execute UnDiFi-2D correctly. This is done using the command:

~~~~~~~
ulimit -s unlimited
~~~~~~~

from the command line. For convenience, you can add this line to your
`.bashrc`.

## Obtaining the source \label{sec:download_source}
The UnDiFi-2D repository is available at GitHub. To obtain the most recent
version you have two possibilities:

* Clone the UnDiFi-2D repository from Github

    git clone https://github.com/UnDiFi-2D/UnDiFi-2D.git

* Download UnDiFi-2D from Github:

    wget https://github.com/UnDiFi-2D/UnDiFi-2D/archive/master.tar.gz
    tar xzf master.tar.gz

## Auxiliary softwares \label{sec:auxiliary}
In addition to the shock-capturing gasdynamic solvers, EulFS and
NEO, there are several auxiliary codes, contained in `source_utils`
directory, which are called during the steps of the UnDiFi-2D
algorithm or during the post-processing phase.
They are listed below:

- dat2triangle: reads files in EulFS format then writes two input files
for Triangle, .node and .poly
- dat2paraview: converts both Neo and EulFS output into Paraview files (*.dat)
- fsplot: performs post-processing for EulFS
- Grid_0: creates the neogrid0.grd file
- Na00x2vvvv: writes the file vvvv_input.dat
- NEO2triangle: reads files in NEO format then writes two input files
for Triangle, .node and .poly
- triangle2dat: converts a Triangle mesh to a .dat file format. (The
input file for Triangle can be created using dat2triangle)
- Triangle2grd: converts triangle format to .grd format (NEO)
- Triangle: mesh generator for UnDiFi-2D
- Tecplot: post-processing output generator for UnDiFi-2D (optional)
- Petsc: necessary to use EulFS (ver. 3.14.6)

## Compiling the code \label{sec:compilingthecode}
* Open a terminal
* Go to the UnDiFi-2D directory
* Compile the code by running the script **./compile_all.sh**

For a list of all compiler options visit Section
\ref{sec:compileroptions}. Finally, all the executables are contained
in the UnDiFi-2D `bin` directory.

## Directory paths \label{sec:installation_directory}
In order to run a testcase, there is no need to manually set any
environment variable since they are automatically set during the
compilation. Nevertheless, it may be useful to make the executables
system-wise visible. There are at least two different ways to enable
it:

1. You can add an alias for the path to your executable. Add a command
of the form:

~~~~~~~
alias undifi='$UNDIFIROOT/bin/UnDiFi-2D-2D_x86_64'
~~~~~~~

to the bottom of the file `~/.bashrc`. Source your `~/.bashrc`
afterwards with:

~~~~~~~
. ~/.bashrc
~~~~~~~

2. You can add the UnDiFi-2D binary directory to your `$PATH` environment
variable by adding:

~~~~~~~
export PATH=$PATH:$UNDIFIROOT/bin
~~~~~~~

to the bottom of the file `~/.bashrc` and sourcing your `~/.bashrc`
afterwards.

Now you are ready for the utilization of UnDiFi-2D and all the
auxiliary executables.
