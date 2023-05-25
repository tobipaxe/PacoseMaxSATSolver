This is hopefully a easy to use SAT solver proxy.

To generate a makefile create a build directory (mkdir build).

mkdir build
cd build
cmake ..
make		// generates a bin directory where a library (.a) and a bin/SATSolverProxy binary is produced


How to generate one library:
ar -x libabc.a
ar -x libxyz.a
ar -qc libaz.a  *.o
ar -qc libaz.a  *.or
