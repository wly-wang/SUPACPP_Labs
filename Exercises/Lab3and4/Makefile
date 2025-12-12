# Compiler and flags
CC     = g++
FLAGS  = -std=c++20 -w
LIBS   = -I ../../GNUplot/ -lboost_iostreams

# Default: build BOTH executables
all: TestFiniteFunctions CustomFunctions

# --------------------------------------------------------
# Build TestFiniteFunctions executable
# --------------------------------------------------------
TestFiniteFunctions: TestFiniteFunctions.o FiniteFunctions.o
	$(CC) $(FLAGS) TestFiniteFunctions.o FiniteFunctions.o $(LIBS) -o TestFiniteFunctions

TestFiniteFunctions.o: TestFiniteFunctions.cxx FiniteFunctions.h
	$(CC) $(FLAGS) $(LIBS) -c TestFiniteFunctions.cxx

# --------------------------------------------------------
# Build CustomFunctions executable (your new main)
# --------------------------------------------------------
CustomFunctions: CustomFunctions.o FiniteFunctions.o
	$(CC) $(FLAGS) CustomFunctions.o FiniteFunctions.o $(LIBS) -o CustomFunctions

CustomFunctions.o: CustomFunctions.cxx CustomFunctions.h FiniteFunctions.h
	$(CC) $(FLAGS) $(LIBS) -c CustomFunctions.cxx

# --------------------------------------------------------
# Base class object file
# --------------------------------------------------------
FiniteFunctions.o: FiniteFunctions.cxx FiniteFunctions.h
	$(CC) $(FLAGS) $(LIBS) -c FiniteFunctions.cxx

# --------------------------------------------------------
# Cleaning rules
# --------------------------------------------------------
clean:
	rm -f *.o *~

cleantarget:
	rm -f TestFiniteFunctions CustomFunctions
