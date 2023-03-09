########################################################################
#
#
#########################################################################

CXX = g++
CFLAGS = $(shell root-config --cflags) -Isrc/
LIBS = $(shell root-config --glibs) -lGeom -lEve -lRGL
OBJS = lib/main.o lib/CodeB_particle.o lib/CodeA_gamma.o lib/MC_sim.o lib/Extraction_SRIM_stopping_powers.o lib/kinematics.o lib/Tools.o

all: EVASIONSsim

EVASIONSsim: $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LIBS)

lib/main.o: src/main.cpp lib/CodeA_gamma.o lib/CodeB_particle.o
	$(CXX) $(CFLAGS) -c src/main.cpp -o $@

lib/CodeA_gamma.o: src/functions/CodeA_gamma.C lib/MC_sim.o lib/Extraction_SRIM_stopping_powers.o
	$(CXX) $(CFLAGS) -c src/functions/CodeA_gamma.C -o $@

lib/CodeB_particle.o: src/functions/CodeB_particle.C
	$(CXX) $(CFLAGS) -c src/functions/CodeB_particle.C -o $@
	
lib/MC_sim.o: src/functions/MC_sim.C lib/kinematics.o lib/Tools.o
	$(CXX) $(CFLAGS) -c src/functions/MC_sim.C -o $@

lib/Extraction_SRIM_stopping_powers.o: src/functions/Extraction_SRIM_stopping_powers.C
	$(CXX) $(CFLAGS) -c src/functions/Extraction_SRIM_stopping_powers.C -o $@
	
lib/kinematics.o: src/functions/kinematics.C lib/Tools.o
	$(CXX) $(CFLAGS) -c src/functions/kinematics.C -o $@

lib/Tools.o: src/functions/Tools.C
	$(CXX) $(CFLAGS) -c src/functions/Tools.C -o $@
	
clean:
	-rm EVASIONSsim lib/*.o
