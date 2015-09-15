#CC=	/ifs/scratch/c2b2/ip_lab/sy2515/HPC/shuo-gcc-4.9.1/bin/g++
CC=	g++
OPT=	-std=c++11 -static-libstdc++
LIBS=	-lm -lpthread
SRCS=	main.cpp basic.cpp expression.cpp genotype.cpp optimization.cpp parameter_init.cpp parameter_save.cpp opt_subroutine.cpp parameter_test.cpp batch.cpp opt_multi_thread.cpp
OBJS=	main.o basic.o expression.o genotype.o optimization.o parameter_init.o parameter_save.o opt_subroutine.o parameter_test.o batch.o opt_multi_thread.o

EXECUTABLE=	main



all: main clean mv

$(OBJS): $(SRCS)
	$(CC) -c $*.cpp $(OPT)

main: $(OBJS)
	$(CC) -o $(EXECUTABLE) $(OBJS) $(OPT) $(LIBS)

clean:
	-rm -f *.o

mv:
	@chmod 755 $(EXECUTABLE)
#	@mv $(EXECUTABLE_NAIVE) ../IBD_C_upgrade_test/
#	-@../IBD_C_upgrade_test/$(EXECUTABLE_NAIVE)
