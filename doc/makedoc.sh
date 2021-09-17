#doxygen  ../.doxygenconfig
#firefox html/index.html

#sudo apt-get install pandoc
#sudo zypper install pandoc

cd userguide
make
if [ $? -ne 0 ]
then
        echo "Compiling in doc ... failed"
        exit  1
fi

#okular userguide.pdf
