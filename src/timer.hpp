#include "Kokkos_Core.hpp"
#include <string>

using namespace std;

class fiestaTimer {

public:
    string name;
    fiestaTimer(string n);
    void accumulate();
    double get();
    void clear();
    void reset();
    void start();
    void stop();
    string getf();

private:
    double time;
    Kokkos::Timer timer;
};
