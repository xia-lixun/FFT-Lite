#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/CChirpZ.o \
	${OBJECTDIR}/CChirpZd.o \
	${OBJECTDIR}/CCoef.o \
	${OBJECTDIR}/CCoef4.o \
	${OBJECTDIR}/CCoef4d.o \
	${OBJECTDIR}/CCoefd.o \
	${OBJECTDIR}/CFastFourier.o \
	${OBJECTDIR}/CFourier.o \
	${OBJECTDIR}/CSort.o \
	${OBJECTDIR}/CSort4.o \
	${OBJECTDIR}/CTrans.o \
	${OBJECTDIR}/CTrans4.o \
	${OBJECTDIR}/TestMain.o


# C Compiler Flags
CFLAGS=-m64

# CC Compiler Flags
CCFLAGS=-m64 -msse4.1 -O3 -march=native
CXXFLAGS=-m64 -msse4.1 -O3 -march=native

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=--64

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fft-lite

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fft-lite: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fft-lite ${OBJECTFILES} ${LDLIBSOPTIONS} -lfftw3f -lfftw3 -lm -lrt

${OBJECTDIR}/CChirpZ.o: CChirpZ.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CChirpZ.o CChirpZ.cpp

${OBJECTDIR}/CChirpZd.o: CChirpZd.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CChirpZd.o CChirpZd.cpp

${OBJECTDIR}/CCoef.o: CCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CCoef.o CCoef.cpp

${OBJECTDIR}/CCoef4.o: CCoef4.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CCoef4.o CCoef4.cpp

${OBJECTDIR}/CCoef4d.o: CCoef4d.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CCoef4d.o CCoef4d.cpp

${OBJECTDIR}/CCoefd.o: CCoefd.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CCoefd.o CCoefd.cpp

${OBJECTDIR}/CFastFourier.o: CFastFourier.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CFastFourier.o CFastFourier.cpp

${OBJECTDIR}/CFourier.o: CFourier.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CFourier.o CFourier.cpp

${OBJECTDIR}/CSort.o: CSort.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CSort.o CSort.cpp

${OBJECTDIR}/CSort4.o: CSort4.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CSort4.o CSort4.cpp

${OBJECTDIR}/CTrans.o: CTrans.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CTrans.o CTrans.cpp

${OBJECTDIR}/CTrans4.o: CTrans4.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CTrans4.o CTrans4.cpp

${OBJECTDIR}/TestMain.o: TestMain.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/TestMain.o TestMain.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fft-lite

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
