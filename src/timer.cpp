#include "Kokkos_Core.hpp"
#include "timer.hpp"
#include <string>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace std;

fiestaTimer::fiestaTimer(){
    timer = new Kokkos::Timer();
    timer->reset();
    time = 0.0;
    description = "No Description";
}

fiestaTimer::fiestaTimer(string n_) : description(n_) {
    timer = new Kokkos::Timer();
    timer->reset();
    time = 0.0;
}

void fiestaTimer::accumulate(){
    time = time + timer->seconds();
}
double fiestaTimer::get(){
    return time;
}
void fiestaTimer::clear(){
    time = 0.0;
}
void fiestaTimer::reset(){
    timer->reset();
}
void fiestaTimer::start(){
    timer->reset();
    time = 0.0;
}
void fiestaTimer::stop(){
    time = timer->seconds();
}
string fiestaTimer::getf(){
    int d;
    int h;
    int m;
    double s;
    double temp;

    d = 0;
    h = 0;
    m = 0;
    s = 0.0;

    temp = time;
    d = time/864000;

    temp = temp - d*864000;
    h = temp/3600;

    temp = temp - h*3600;
    m = temp/60;

    temp = temp - m*60;
    s = temp;

    stringstream ss;

    ss.precision(2);
    ss.setf(ios::fixed);
    if (d > 0)
        ss << d << "d";
    if (h > 0)
        ss << h << "h";
    if (m > 0)
        ss << m << "m";
    ss << s << "s";

    return ss.str();
}
string fiestaTimer::describe(){
    return description;
}
