LIBS=-L/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lm -ldl -m64 -I/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/include -lRooFitCore -lRooFit -lrt -lMinuit -lFoam -lMathMore
CFLAGS=-I/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/include
COMILELAG=-O3 -std=c++17
DVAR=-DUSEROOT
all:
	nvcc $(LIBS) $(CFLAGS) -dc -w src/Event.cu -I./include $(DVAR) $(COMILELAG) -o Obj/Event.o
	nvcc $(LIBS) $(CFLAGS) -dc -w src/DPD.cu -I./include $(DVAR) $(COMILELAG) -o Obj/DPD.o
	nvcc $(LIBS) $(CFLAGS) -dc -w src/Amplitude.cu -I./include $(DVAR) $(COMILELAG) -o Obj/Amplitude.o
	nvcc $(LIBS) $(CFLAGS) -dc -w src/NLL_estimator.cu -I./include $(DVAR) $(COMILELAG) -o Obj/NLL_estimator.o
	nvcc $(LIBS) $(CFLAGS) -dc -w src/Minimize.cu -I./include $(DVAR) $(COMILELAG) -o Obj/Minimize.o
	g++ -c ./src/jsoncpp.cpp -o ./Obj/jsoncpp.o
	nvcc $(LIBS) $(CFLAGS) -dc -w ./src/Fit.cu -I./include $(DVAR) $(COMILELAG) -o Obj/Fit.o
	nvcc -w Obj/Event.o Obj/DPD.o Obj/Amplitude.o Obj/NLL_estimator.o Obj/Minimize.o ./Obj/jsoncpp.o Obj/Fit.o $(LIBS) $(CFLAGS) $(DVAR) $(COMILELAG) -o Fit.exe

clean:
	rm Obj/Event.o Obj/DPD.o Obj/Amplitude.o Obj/NLL_estimator.o Obj/Minimize.o ./Obj/jsoncpp.o Obj/Fit.o
