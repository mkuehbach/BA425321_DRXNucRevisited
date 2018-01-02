cmake_minimum_required(VERSION 3.0)

################################################################################################################
##USER INTERACTION##############################################################################################
################################################################################################################
#please name your project accordingly
set(MYPROJECTNAME "drxnuc")

#choose compiler
set(EMPLOY_INTELCOMPILER OFF)
set(EMPLOY_GNUCOMPILER ON)
set(VERBOSE_VECTORIZATION OFF)

#optimization level -O0, -O1, -O2, or -O3
set(MYOPTLEVEL "-O3")

#additional verbose information


################################################################################################################
##AUTOMATIC SECTION#############################################################################################
################################################################################################################
#user input sanity checks
if(EMPLOY_INTELCOMPILER AND EMPLOY_GNUCOMPILER)
	message([FATAL_ERROR] "You cannot utilize two compiler at the same time!")
endif()


#automatically assign project name and compiler flags
project(${MYPROJECTNAME})
set(CMAKE_BUILD_DIR "build")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MYOPTLEVEL}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MYOPTLEVEL}")

#setting up compiler-specifics
#intel path
if(EMPLOY_INTELCOMPILER)
	message([STATUS] "Employing the Intel compiler!")
	if (VERBOSE_VECTORIZATION)
		set(MYVERBOSE "-qopt-report-phase=vec")
	endif()
	add_definitions("${MYOPTLEVEL} ${MYVERBOSE}")
	message([WARNING] "Currently not SIMD flags provided for Intel compile options")
	add_definitions("-qopenmp -lpthread")
	add_definitions("-std=c++0x")
	add_definitions("-Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat 
		-Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull 
		-Wparentheses -Wpointer-sign -Wreorder -Wreturn-type -Wsign-compare 
		-Wsequence-point -Wtrigraphs -Wunused-variable")
#gcc path
elseif(EMPLOY_GNUCOMPILER)
	message([STATUS] "Employing the GNU compiler!")
	if(VERBOSE_VECTORIZATION)
		set(MYVERBOSE "-fopt-info")
	endif()
	add_definitions("${MYOPTLEVEL} ${MYVERBOSE} ${MYSIMDFLAGS}")
	set(MYOMP "-fopenmp -lpthread")
	add_definitions(${MYOMP})
	add_definitions("-std=c++11")
	add_definitions("-Wall -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat 
		-Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull 
		-Wparentheses -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point 
		-Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable")
else()
	message([FATAL_ERROR] "You have to utilize a compiler!")
endif()

message([STATUS] "Projectname is ${MYPROJECTNAME}")
message([STATUS] "We utilize optimization level ${MYOPTLEVEL}")

#parallelization - MPI process-level
#query location of MPI library
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_executable(${MYPROJECTNAME}
	../src/DRXNUC_Information.cpp
	../src/DRXNUC_Math.cpp
	../src/DRXNUC_Settings.cpp
	../src/DRXNUC_Kernel.cpp
	../src/DRXNUC_Main.cpp
)

#linking process
target_link_libraries(${MYPROJECTNAME} ${MYOMP} ${MPI_LIBRARIES} )

#MPI compilation settings
if(MPI_COMPILE_FLAGS)
  set_target_properties(${MYPROJECTNAME} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${MYPROJECTNAME} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
