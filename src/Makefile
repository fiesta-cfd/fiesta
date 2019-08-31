default:
		mpicc -c -g mpi_init.c
		mpicc -c -g input.c -llua53 -lm -ldl
		mpicc -g main.c -llua53 -lm -ldl -lcgns input.o mpi_init.o
