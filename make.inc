
ifeq ($(HOST),linux-gfortran)

OBJSUF=o
MODSUF=mod
FC=gfortran -c -std=legacy
FFLAGS=-O3
FLINK=gfortran -o $(PROJECT)

endif

ifeq ($(HOST),linux-gfortran-openmp)

OBJSUF=o
MODSUF=mod
FC=gfortran -c -std=legacy
FFLAGS=-O3 --openmp
FLINK=gfortran -o $(PROJECT) --openmp
### export OMP_NUM_THREADS=4
### export OMP_STACKSIZE=1024M

endif

ifeq ($(HOST),macos-gfortran)

OBJSUF=o
MODSUF=mod
FC=gfortran -c -std=legacy
FFLAGS=-O3
FLINK=gfortran -o $(PROJECT) -Wl,-stack_size,0x20000000

endif

ifeq ($(HOST),macos-gfortran-openmp)

OBJSUF=o
MODSUF=mod
FC=gfortran -c -std=legacy
FFLAGS=-O3 --openmp
FLINK=gfortran -o $(PROJECT) --openmp -Wl,-stack_size,0x20000000
### export OMP_NUM_THREADS=4
### export OMP_STACKSIZE=1024M

endif

ifeq ($(HOST),macos-gfortran-quad)

OBJSUF=o
MODSUF=mod
FC=gfortran -c -std=legacy
FFLAGS=-O3 -freal-8-real-16
FLINK=gfortran -o $(PROJECT) -Wl,-stack_size,0x20000000

endif

ifeq ($(HOST),linux-gfortran-openmp-quad)

OBJSUF=o
MODSUF=mod
FC=gfortran -c -std=legacy
FFLAGS=-O3 --openmp -freal-8-real-16
FLINK=gfortran -o $(PROJECT) --openmp

endif

ifeq ($(HOST),macos-gfortran-openmp-quad)

OBJSUF=o
MODSUF=mod
FC=gfortran -c -std=legacy
FFLAGS=-O3 --openmp -freal-8-real-16
FLINK=gfortran -o $(PROJECT) --openmp -Wl,-stack_size,0x20000000
### export OMP_NUM_THREADS=4
### export OMP_STACKSIZE=1024M

endif

