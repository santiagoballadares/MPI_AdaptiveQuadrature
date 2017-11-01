# MPI_AdaptiveQuadrature
Simple implementation of the adaptive quadrature algorithm using the bag of tasks pattern with C and MPI.

Compile with: 

      `gcc mpi_quadrature.c stack.h stack.c -lmsmpi -o mpi_quadrature.exe`

Run with: 

      `mpiexec -n 8 mpi_quadrature.exe` 

where `-n x` sets the number of processes

Setup for windows:
1. Install MinGW
2. Install MS-MPI
3. Adapt MS-MPI for MinGW
  - Create the `libmsmpi64.a` library with the MinGW64 tools `gendef` and `dlltool`
      ```
      > gendef msmpi.dll                                  # generate msmpi.def
      > dlltool -d msmpi.def -D msmpi.dll -l libmsmpi.a   # generate the (static) library file libmsmpi.a
      ```
  - Copy the new library to where g++ looks for them, e.g. `/mingw64/lib`
  - Modify the header file `mpi.h`. Add `#include <stdint.h>` above `typedef __int64 MPI_Aint`
  - Copy the modified header file `mpi.h` to the default include folder e.g. `/mingw64/include`


For more details see:
- http://www.math.ucla.edu/~wotaoyin/windows_coding.html
- https://github.com/coderefinery/autocmake/issues/85 
