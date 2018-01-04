# BA425321_DRXNucRevisited
This repository provides supplementary pieces of information to the paper draft by Markus Kuehbach, Muhammad Imran, Franz Roters, and Markus Bambach on
"Quantitative Implications of Dislocation Density Gradients for SIBM-Based Discontinuous Dynamic Recrystallization Nucleation Models"

Requirements:<br>
-Working Cxx11 compatible C/C++ compiler toolchain (Ubuntu 16.04 and above)<br>
-MessagePassingInterface, either openmpi or mpich<br>
-cmake and make<br>

Check for existence and functionality of all through Linux shell/console via:<br>
 cmake --version<br>
 make --version<br>
 gcc --version<br>
 g++ --version<br>
 mpicc --version<br>
 mpicxx --version<br>

How to compile:<br>
-Download tar and unpack, e.g. tar -xvf archivename.tar.gz<br>
-Go into build directory<br>
-Type cmake -DCMAKE_BUILD_TYPE=Release ..<br>
-make<br>

Execute via:<br>
-Modify XML settings file in build folder<br>
-mpiexec -n <nprocesses> ./drxnuc <UnsignedInteger> <NameOfXMLSettingsFile><br> 
 for example:    mpiexec -n 1 ./drxnuc 1 DRXNUC.Input.xml<br>
