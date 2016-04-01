##
## when "make train", the training program will be compiled
## when "make test", the testing program will be compiled
##

CC=	g++
#OPT=	-std=c++11 -static-libstdc++ -Wall -pg
OPT=	-std=c++11 -static-libstdc++
LIBS=	-lm -lpthread
SRCS_TRAIN=	main.cpp basic.cpp expression.cpp genotype.cpp optimization.cpp parameter_init.cpp parameter_save.cpp opt_subroutine.cpp batch.cpp opt_multi_thread.cpp opt_nn_acfunc.cpp opt_para_save.cpp opt_debugger.cpp libfunc_matrix.cpp opt_hierarchy.cpp libfunc_matrix_mt.cpp
OBJS_TRAIN=	main.o basic.o expression.o genotype.o optimization.o parameter_init.o parameter_save.o opt_subroutine.o batch.o opt_multi_thread.o opt_nn_acfunc.o opt_para_save.o opt_debugger.o libfunc_matrix.o opt_hierarchy.o libfunc_matrix_mt.o
SRCS_TEST=	test_main.cpp basic.cpp expression.cpp genotype.cpp batch.cpp test_para_read.cpp test_predict.cpp test_save.cpp opt_nn_acfunc.cpp
OBJS_TEST=	test_main.o basic.o expression.o genotype.o batch.o test_para_read.o test_predict.o test_save.o opt_nn_acfunc.o

EXECUTABLE_TRAIN=	train
EXECUTABLE_TEST=	test



all:
	@echo Please choose from \"make train\" and \"make test\" to compile and link...



### compiling starts from here
train: main_train clean mv_train


$(OBJS_TRAIN): $(SRCS_TRAIN)
	$(CC) -c $*.cpp $(OPT)

main_train: $(OBJS_TRAIN)
	$(CC) -o $(EXECUTABLE_TRAIN) $(OBJS_TRAIN) $(OPT) $(LIBS)



test: main_test clean clean


$(OBJS_TEST): $(SRCS_TEST)
	$(CC) -c $*.cpp $(OPT)

main_test: $(OBJS_TEST)
	$(CC) -o $(EXECUTABLE_TEST) $(OBJS_TEST) $(OPT) $(LIBS)




clean:
	-rm -f *.o



mv_train:
	@chmod 755 $(EXECUTABLE_TRAIN)
#	@mv $(EXECUTABLE_NAIVE) ../IBD_C_upgrade_test/
#	-@../IBD_C_upgrade_test/$(EXECUTABLE_NAIVE)

mv_test:
	@chmod 755 $(EXECUTABLE_TEST)
#	@mv $(EXECUTABLE_NAIVE) ../IBD_C_upgrade_test/
#	-@../IBD_C_upgrade_test/$(EXECUTABLE_NAIVE)

