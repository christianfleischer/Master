#include "getcoeffcients.h"
#include <armadillo>

using namespace std;


GetCoeffcients::GetCoeffcients() {
}

void GetCoeffcients::retrieveFromFile(string fileName) {
    mat loadCoefficients;
    loadCoefficients.load(fileName, raw_ascii);

}

