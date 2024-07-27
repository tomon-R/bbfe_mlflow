CC = mpic++
FC = mpif90

CFLAGS = -g -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I../../submodule/monolis/include
LIBS = -L/usr/local/lib -lm -lstdc++ -L../../submodule/monolis/lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -llapack

INCLUDES_BB = -I../../include
LIBS_BB = -L../../lib -lBB -lBBFE_std -lBBFE_sys -lBBFE_elemmat -lBBFE_manusol

TARGET1 = mlflow_fs_sups
OBJS1 = mlflow_fs_sups.o mlflow_impfunc.o mlflow_elemmat.o fluid_core.o fluid_elemmat.o mlflow_utils.o 

.SUFFIXES: .c .cpp .o

all: $(TARGET1) 

$(TARGET1): $(OBJS1)
	$(FC) $(CFLAGS) $(INCLUDES) $(INCLUDES_BB) $(LDFLAGS) -o $@ $^ $(LIBS_BB) $(LIBS)

clean:
		-rm -f $(OBJS1) $(TARGET1) $(OBJS2) $(TARGET2)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(INCLUDES_BB) $(LIBS_BB)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(INCLUDES_BB)  $(LIBS_BB)
