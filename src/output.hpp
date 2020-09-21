#include "input.hpp"
#include <string>

using namespace std;

void printSplash(int cFlag);
void printConfig(struct inputConfig cf);

enum Color { RED, GRE, YEL, BLU, MAG, CYA, NON };

string c(int cFlag, Color color);
