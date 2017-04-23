cc=g++
LIBS=$(shell root-config --libs)
INCLUDE=$(shell root-config --cflags)

exe=calc
obj=calc.o DMCalc.o
$(exe):$(obj)
	$(cc) $(INCLUDE) -o $(exe) $(obj) $(LIBS)
calc.o:calc.C DMCalc.h
	$(cc) -c $(INCLUDE) calc.C
DMCalc.o: DMCalc.C DMCalc.h
	$(cc) -c $(INCLUDE) DMCalc.C

.PHONY: clean
clean:
	rm *.o $(exe)

