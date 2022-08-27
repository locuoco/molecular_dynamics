cd ./math
echo ======== DFT ========
./dft.sh
./dft
echo ======== HELPER ========
./helper.sh
./helper
cd ../physics
echo ======== EWALD ========
./ewald.sh
./ewald
echo ======== MOLECULE ========
./molecule.sh
./molecule
echo ======== PPPM ========
./pppm.sh
./pppm
echo ======== TENSOR ========
./tensor.sh
./tensor
cd ./integrator
echo ======== INTEGRATOR ========
./integrator.sh
./integrator
cd ../../utils
echo ======== PARALLEL_SORT ========
./parallel_sort.sh
./parallel_sort