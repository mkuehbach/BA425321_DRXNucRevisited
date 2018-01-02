# BA425321_DRXNucRevisited
This repository provides supplementary pieces of information to the paper draft by Markus K\"uhbach, Muhammad Imran, Franz Roters, and Markus Bambach on
"Quantitative Implications of Dislocation Density Gradients for SIBM-Based Discontinuous Dynamic Recrystallization Nucleation Models"

Requirements:
-Working Cxx11 compatible C/C++ compiler toolchain (Ubuntu 16.04 and above)<br>
-MessagePassingInterface, either openmpi or mpich<br>
-cmake and make

Check for all in Linux console via:
 cmake --version
 make --version
 gcc --version
 g++ --version
 mpicc --version
 mpicxx --version

How to compile:
-Download tar and unpack, e.g. tar -xvf archivename.tar.gz
-Go into build directory
-Type cmake -DCMAKE_BUILD_TYPE=Release ..
-make

Execute via:
-Modify XML settings file in build folder
-mpiexec -n <nprocesses> ./drxnuc <UnsignedInteger> <NameOfXMLSettingsFile> for example
 mpiexec -n 1 ./drxnuc 1 DRXNUC.Input.xml
