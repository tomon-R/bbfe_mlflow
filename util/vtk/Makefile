CC       := g++
CFLAGS   := -Wall -Wextra -pedantic -g

INCLUDES = -I/usr/local/include -I../../include
LIBS = -L/usr/local/lib -lm -lstdc++ -L../../lib -lBB -lBBFE_sys 

TARGET1 = vtkddm_hex
OBJ1 = vtkddm_hex.o
TARGET2 = vtk_wireframe_hex
OBJ2 = vtk_wireframe_hex.o
TARGET3 = txt2vtk_hex_nodalvals
OBJ3 = txt2vtk_hex_nodalvals.o

all: build clean_data

build: $(TARGET1) $(TARGET2) $(TARGET3)

$(TARGET1): $(OBJ1)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)
$(TARGET2): $(OBJ2)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)
$(TARGET3): $(OBJ3)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(LIBS) 

.PHONY clean:;
	$(RM) -f $(OBJ1) $(TARGET1)
	$(RM) -f $(OBJ2) $(TARGET2)
	$(RM) -f $(OBJ3) $(TARGET3)

clean_data:
	$(RM) -f $(OBJ1)
	$(RM) -f $(OBJ2)
	$(RM) -f $(OBJ3)
