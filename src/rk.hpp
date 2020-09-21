#ifndef RK_H
#define RK_H

#include "input.hpp"
#include "rkfunction.hpp"

void rkAdvance(struct inputConfig&, class rk_func*);

#endif
