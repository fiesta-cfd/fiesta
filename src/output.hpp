#include "input.hpp"
#include <string>

using namespace std;

void printSplash(void);
void printConfig(struct inputConfig cf);

enum Color {
    RED,
    GRE,
    YEL,
    BLU,
    MAG,
    CYA,
    NON
};

string c(Color color);
