petsc11=$HOME/src/petsc-3.11.4
petsc12=$HOME/src/petsc-3.12.5
diff $petsc11/lib/petsc/conf/variables $petsc12/lib/petsc/conf/variables > log
diff $petsc11/lib/petsc/conf/rules $petsc12/lib/petsc/conf/rules >> log

