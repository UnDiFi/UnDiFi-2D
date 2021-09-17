# "make" builds all

# Switch to build all document parts - or only ready parts
BUILDALL=

SRC=000_userguide.md \
		010_installation.md \
		020_workflow.md \
		040_solver.md \
		050_mesh.md \
		060_tutorials.md \
		070_Eulfs_tutorial.md
#       030_algorithm.md
SRC += references.md

all: userguide
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' SUCCESS: DOCUMENTATION GENERATED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

userguide:
#	pandoc $(SRC) -o userguide.pdf -s -N --citeproc --from markdown+implicit_figures --listings --template=../template.tex
	pandoc $(SRC) -o userguide.pdf -s -N --filter pandoc-citeproc --from markdown+implicit_figures --listings --template=../template.tex

.PHONY: userguide
