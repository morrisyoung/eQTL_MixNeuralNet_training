CC=	g++
OPT=	-std=c++11
#LIBS=	-lm -lpthread
LIBS=
SRCS=	main.cpp basic.cpp expression.cpp genotype.cpp
OBJS=	main.o basic.o expression.o genotype.o

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
