# Compilateur utilisé
CC=g++
CC_para=mpic++


run : Chaleur.cc fonction.cpp fonction.h
	$(CC_para) -O3 -std=c++11  Chaleur.cc fonction.cpp -o run

#si on a des trucs a tester :
test : test.cc
	$(CC) test.cc  $(CXX_FLAGS) -o run_test




# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ run
