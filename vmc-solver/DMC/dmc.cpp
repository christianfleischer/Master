#include "dmc.h"

DMC::DMC()
{

}

void DMC::runDMC() {
    cycle = 1;
    while (cycle <= m_numberOfEquilibrationSteps) {
        for (int j = 0; j < numberOfWalkers; j++) {

        }

        cycle++;
    }

    cycle = 1;
    while (cycle <= N_c) {
        cycle++;
    }


}

void DMC::setEquilibrationSteps(int equilibration) {
    m_numberOfEquilibrationSteps = equilibration;
}
