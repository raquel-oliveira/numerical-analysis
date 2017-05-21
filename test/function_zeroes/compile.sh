#g++ testes.cpp -o "../../bin/testes" -I "../../include/" -std=c++11 -Wall
#g++ main.cpp -o "../../bin/main" -I "../../include/" -std=c++11 -Wall -O3
#g++ generateMatrix.cpp -o "../../bin/generateMatrix" -I "../../include/" -std=c++11 -Wall
g++ find_zeroes.cpp -o ../../bin/find_zeroes -I "../../include" -std=c++11 -Wall
g++ complex_roots.cpp -o ../../bin/complex_roots -I "../../include" -std=c++14 -Wall
g++ basin_attraction.cpp -o ../../bin/basin_attraction -I "../../include" -std=c++14 `pkg-config opencv --libs --cflags` -I -g
