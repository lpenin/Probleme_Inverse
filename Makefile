DEBUG_FLAG = -g3 -DDEBUG -std=c++11
OPTIM_FLAG = -O3 -w -DNDEBUG -std=c++11


run : Chaleur.cc fonction.cpp fonction.h spline.h
#g++ $(DEBUG_FLAG)  Chaleur.cc fonction.cpp -o run
	g++ $(OPTIM_FLAG)  Chaleur.cc fonction.cpp -o run


clean :
	rm -f *.o *~ run
