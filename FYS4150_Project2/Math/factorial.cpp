#include "factorial.h"
#include <cassert>

int factorial(int n) {
    assert(n>=0);
    int fact = 1;
    for (int i = 2; i <= n; i++) {
        fact *= i;
    }
    return fact;
}

