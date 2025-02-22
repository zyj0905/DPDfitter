LIBS=-L/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lm -ldl -m64 -I/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/include -lRooFitCore -lRooFit -lrt -lMinuit -lFoam -lMathMore
CFLAGS=-I/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/includ
COMILELAG=-O3 -std=c++17
DVAR=-DUSEROOT
all:
	g++ $(LIBS) $(CFLAGS) -c src/Event.cxx -I./include $(DVAR) $(COMILELAG) -o Obj/Event.o
	g++ $(LIBS) $(CFLAGS) -c src/DPD.cxx -I./include $(DVAR) $(COMILELAG) -o Obj/DPD.o
	g++ $(LIBS) $(CFLAGS) -c src/Amplitude.cxx -I./include $(DVAR) $(COMILELAG) -o Obj/Amplitude.o
	g++ $(LIBS) $(CFLAGS) -c src/NLL_estimator.cxx -I./include $(DVAR) $(COMILELAG) -o Obj/NLL_estimator.o
	g++ -c ./src/jsoncpp.cpp -o ./Obj/jsoncpp.o
	g++ $(LIBS) $(CFLAGS) -c ./src/SampleMC.cxx -I./include $(DVAR) $(COMILELAG) -o Obj/SampleMC.o
	g++ -w Obj/Event.o Obj/DPD.o Obj/Amplitude.o Obj/NLL_estimator.o ./Obj/jsoncpp.o Obj/SampleMC.o $(LIBS) $(CFLAGS) $(DVAR) $(COMILELAG) -o SampleMC.exe

clean:
	rm Obj/Event.o Obj/DPD.o Obj/Amplitude.o Obj/NLL_estimator.o ./Obj/jsoncpp.o Obj/SampleMC.o
