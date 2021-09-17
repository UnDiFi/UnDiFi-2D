code=eulfs3.1.0-$(PETSC_ARCH)
for i in 1 2 3 4
do
echo hello world!
# run one step of the EulFS code
$code
# run the dat 2 angener translator code
if (echo 1 | dat2an ) \
        then echo "dat2an returned ok"; else echo "dat2an has failed"; exit 1; fi
# run the angener code
if (side ) \
        then echo "side returned ok"; else echo "side has failed"; exit 1; fi
done
