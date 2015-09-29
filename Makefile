##
## when "make", the training program will be compiled
## when "make test", the testing program will be compiled
##


#CC=	/ifs/scratch/c2b2/ip_lab/sy2515/HPC/shuo-gcc-4.9.1/bin/g++
CC=	g++
OPT=	-std=c++11 -static-libstdc++
LIBS=	-lm -lpthread
SRCS=	main.cpp basic.cpp expression.cpp genotype.cpp optimization.cpp parameter_init.cpp parameter_save.cpp opt_subroutine.cpp parameter_test.cpp batch.cpp opt_multi_thread.cpp
OBJS=	main.o basic.o expression.o genotype.o optimization.o parameter_init.o parameter_save.o opt_subroutine.o parameter_test.o batch.o opt_multi_thread.o
SRCS_TEST=	test.cpp basic.cpp expression.cpp genotype.cpp batch.cpp test_para_read.cpp test_predict.cpp test_save.cpp
OBJS_TEST=	test.o basic.o expression.o genotype.o batch.o test_para_read.o test_predict.o test_save.o

EXECUTABLE=	main
EXECUTABLE_TEST=	test



### compiling starts from here
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


test: main_test clean_test mv_test

$(OBJS_TEST): $(SRCS_TEST)
	$(CC) -c $*.cpp $(OPT)

main_test: $(OBJS_TEST)
	$(CC) -o $(EXECUTABLE_TEST) $(OBJS_TEST) $(OPT) $(LIBS)

clean_test:
	-rm -f *.o

mv_test:
	@chmod 755 $(EXECUTABLE_TEST)
#	@mv $(EXECUTABLE_NAIVE) ../IBD_C_upgrade_test/
#	-@../IBD_C_upgrade_test/$(EXECUTABLE_NAIVE)
