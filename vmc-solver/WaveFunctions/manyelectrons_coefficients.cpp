#include "manyelectrons_coefficients.h"
#include <cmath>
#include <cassert>
#include "../InitialStates/randomuniform.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include "../Hamiltonians/hamiltonian.h"
#include <iostream>
#include "Math/factorial.h"

ManyElectronsCoefficients::ManyElectronsCoefficients(System* system, double alpha, double beta, double omega, double C, bool Jastrow) :
        WaveFunction(system) {
    assert(omega > 0);
    m_omega = omega;
    assert(alpha >= 0);
    assert(system->getNumberOfDimensions() > 0 && system->getNumberOfDimensions() <= 3);
    m_numberOfDimensions = system->getNumberOfDimensions();
    //m_alpha = alpha;
    //m_alphaOmega = m_alpha*m_omega;
    m_Jastrow = Jastrow;
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_C = C;
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_halfNumberOfParticles = m_numberOfParticles/2.;

    if (m_numberOfParticles == 1) {
        setUpSlaterDetOneParticle();
    }
    else {
        setUpSlaterDet();
    }

    setUpDistances();
    setUpJastrowMat();
}

double ManyElectronsCoefficients::evaluate(std::vector<class Particle*> particles) {
    // Evaluates the wave function using brute force.

    mat spinUpSlater;
    mat spinDownSlater;
    if (m_numberOfParticles == 1) {
        spinUpSlater = zeros<mat>(m_numberOfParticles, m_numberOfParticles);
        spinDownSlater = zeros<mat>(m_numberOfParticles, m_numberOfParticles);
        spinUpSlater(0,0) = m_SPWFMat(0,0);
    }
    else {
        spinUpSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);
        spinDownSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);

        for (int i=0; i < m_halfNumberOfParticles; i++) {
            //std::vector<double> rSpinUp = particles[i]->getPosition();//m_system->getParticles()[i]->getPosition();
            //double xSpinUp = rSpinUp[0];
            //double ySpinUp = rSpinUp[1];
            //std::vector<double> rSpinDown = particles[i+m_halfNumberOfParticles]->getPosition();//m_system->getParticles()[i+m_halfNumberOfParticles]->getPosition();
            //double xSpinDown = rSpinDown[0];
            //double ySpinDown = rSpinDown[1];

            for (int j=0; j < m_halfNumberOfParticles; j++) {
                //int nx = m_quantumNumbers(j, 0);
                //int ny = m_quantumNumbers(j, 1);
                spinUpSlater(i,j) = m_SPWFMat(i,j);//evaluateSingleParticleWF(nx, ny, xSpinUp, ySpinUp);
                spinDownSlater(i,j) = m_SPWFMat(i+m_halfNumberOfParticles,j);//evaluateSingleParticleWF(nx, ny, xSpinDown, ySpinDown);
            }
        }
    }

    double beta = m_parameters[1];
    double exponent = 0;
    if (m_Jastrow) {
        for (int i=0; i < m_numberOfParticles; i++) {
            //std::vector<double> r_i = particles[i]->getPosition();

            for (int j=i+1; j < m_numberOfParticles; j++) {
                //std::vector<double> r_j = particles[j]->getPosition();
                //double r_ij = (r_i[0] - r_j[0])*(r_i[0] - r_j[0]) + (r_i[1] - r_j[1])*(r_i[1] - r_j[1]);
                double r_ij = m_distances(i,j);//sqrt(r_ij);
                double denom = 1+beta*r_ij;
                exponent += m_a(i,j)*r_ij/denom;
            }
        }
    }

    double waveFunction;
    if (m_numberOfParticles == 1) {
        waveFunction = det(spinUpSlater)*exp(exponent);
    }
    else {
        waveFunction = det(spinDownSlater)*det(spinUpSlater)*exp(exponent);
    }

    return waveFunction;
}

std::vector<double> ManyElectronsCoefficients::computeDerivative(std::vector<class Particle*> particles) {
    //Calculates ∇ψ/ψ for the wave function.

    int i = m_system->getCurrentParticle();
    int numberOfParticles = m_system->getNumberOfParticles();
    std::vector<double> derivative(numberOfParticles*m_numberOfDimensions);

    for (int d = 0; d < m_numberOfDimensions; d++) {
        derivative[i*m_numberOfDimensions+d] = computeSlaterGradient(i)[d];
    }

//    derivative[i*m_numberOfDimensions] = computeSlaterGradient(i)[0]
//                                        ;//+computeJastrowGradient(particles, i)[0];
//    derivative[i*m_numberOfDimensions+1] = computeSlaterGradient(i)[1]
//                                          ;//+computeJastrowGradient(particles, i)[1];

    if (m_Jastrow) {
        for (int d = 0; d < m_numberOfDimensions; d++) {
            derivative[i*m_numberOfDimensions+d] += m_JastrowGrad(i,d);//computeJastrowGradient(particles, i)[0];
            //derivative[i*m_numberOfDimensions+1] += m_JastrowGrad(i,1);//computeJastrowGradient(particles, i)[1];
        }
    }
    return derivative;
    //return 0;
}

std::vector<double> ManyElectronsCoefficients::computeSlaterGradient(/*std::vector<class Particle*> particles, */int i) {
    // Computes the gradient of the Slater part of the wave function.
    std::vector<double> slaterGradient(m_numberOfDimensions);
    for (int d = 0; d < m_numberOfDimensions; d++) {
        slaterGradient[d] = 0;
    }
//    slaterGradient[0] = 0;
//    slaterGradient[1] = 0;
    //double x = particles[i]->getPosition()[0];
    //double y = particles[i]->getPosition()[1];

    if (i < m_halfNumberOfParticles) {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            //int nx = m_quantumNumbers(j, 0);
            //int ny = m_quantumNumbers(j, 1);
            vec SPWFGradient = m_SPWFDMat(i,j);//computeSPWFDerivative(nx, ny, x, y);
            for (int d = 0; d < m_numberOfDimensions; d++) {
                slaterGradient[d] += SPWFGradient[d]*m_spinUpSlaterInverse(j,i);
            }
//            slaterGradient[0] += SPWFGradient[0]*m_spinUpSlaterInverse(j,i);
//            slaterGradient[1] += SPWFGradient[1]*m_spinUpSlaterInverse(j,i);
        }
    }
    else {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            //int nx = m_quantumNumbers(j, 0);
            //int ny = m_quantumNumbers(j, 1);
            vec SPWFGradient = m_SPWFDMat(i,j);//computeSPWFDerivative(nx, ny, x, y);
            for (int d = 0; d < m_numberOfDimensions; d++) {
                slaterGradient[d] += SPWFGradient[d]*m_spinDownSlaterInverse(j, i-m_halfNumberOfParticles);
            }
//            slaterGradient[0] += SPWFGradient[0]*m_spinDownSlaterInverse(j, i-m_halfNumberOfParticles);
//            slaterGradient[1] += SPWFGradient[1]*m_spinDownSlaterInverse(j, i-m_halfNumberOfParticles);
        }
    }

    return slaterGradient;

}

std::vector<double> ManyElectronsCoefficients::computeJastrowGradient(std::vector<class Particle*> particles, int k) {
    // Computes the gradient of the Jastrow part of the wave function.
    std::vector<double> jastrowGradient(m_numberOfDimensions);
    for (int d = 0; d < m_numberOfDimensions; d++) {
        jastrowGradient[d] = 0;
    }

//    jastrowGradient[0] = jastrowGradient[1] = 0;

    double beta = m_parameters[1];
    std::vector<double> r_k = particles[k]->getPosition();

    for (int j=0; j < k; j++) {
        std::vector<double> r_j = particles[j]->getPosition();
        double r_kj = m_distances(k,j);//(r_k[0]-r_j[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_j[1])*(r_k[1]-r_j[1]);
        //r_kj = sqrt(r_kj);
        double denom = 1 + beta*r_kj;
        for (int d = 0; d < m_numberOfDimensions; d++) {
            jastrowGradient[d] += (r_k[d]-r_j[d])/r_kj * m_a(k,j)/(denom*denom);
        }
//        jastrowGradient[0] += (r_k[0]-r_j[0])/r_kj * m_a(k, j)/(denom*denom);
//        jastrowGradient[1] += (r_k[1]-r_j[1])/r_kj * m_a(k, j)/(denom*denom);
    }

    for (int j=k+1; j < m_numberOfParticles; j++) {
        std::vector<double> r_j = particles[j]->getPosition();
        double r_kj = m_distances(k,j);//(r_k[0]-r_j[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_j[1])*(r_k[1]-r_j[1]);
        //r_kj = sqrt(r_kj);
        double denom = 1 + beta*r_kj;
        for (int d = 0; d < m_numberOfDimensions; d++) {
            jastrowGradient[d] += (r_k[d]-r_j[d])/r_kj * m_a(k,j)/(denom*denom);
        }
//        jastrowGradient[0] += (r_k[0]-r_j[0])/r_kj * m_a(k, j)/(denom*denom);
//        jastrowGradient[1] += (r_k[1]-r_j[1])/r_kj * m_a(k, j)/(denom*denom);
    }

    return jastrowGradient;

}

double ManyElectronsCoefficients::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    //Calculates ∇²ψ/ψ for the wave function.

    double slaterLaplacian = 0;
    double jastrowLaplacian = 0;
    double crossTerm = 0;

    if (m_numberOfParticles == 1) {
        slaterLaplacian += m_SPWFDDMat(0,0)*m_spinUpSlaterInverse(0,0);
    }
    else {
        for (int i=0; i < m_halfNumberOfParticles; i++) {
            //std::vector<double> r_i = particles[i]->getPosition();
            //double x = r_i[0];
            //double y = r_i[1];

            for (int j=0; j < m_halfNumberOfParticles; j++) {
                //int nx = m_quantumNumbers(j, 0);
                //int ny = m_quantumNumbers(j, 1);
                slaterLaplacian += m_SPWFDDMat(i,j)//computeSPWFDoubleDerivative(nx, ny, x, y)
                                  *m_spinUpSlaterInverse(j,i);
            }
        }
        for (int i=m_halfNumberOfParticles; i < m_numberOfParticles; i++) {
            //std::vector<double> r_i = particles[i]->getPosition();
            //double x = r_i[0];
            //double y = r_i[1];

            for (int j=0; j < m_halfNumberOfParticles; j++) {
                //int nx = m_quantumNumbers(j, 0);
                //int ny = m_quantumNumbers(j, 1);
                slaterLaplacian += m_SPWFDDMat(i,j)//computeSPWFDoubleDerivative(nx, ny, x, y)
                                  *m_spinDownSlaterInverse(j,i-m_halfNumberOfParticles);
            }
        }
    }

    if (m_Jastrow) {
        double beta = m_parameters[1];
        double jastrowSum1 = 0;
        double jastrowSum2 = 0;
        int dim = m_numberOfDimensions;

        for (int k=0; k < m_numberOfParticles; k++) {
            std::vector<double> r_k = particles[k]->getPosition();

            for (int j=0; j < k; j++) {
                std::vector<double> r_j = particles[j]->getPosition();
                double r_kj = m_distances(k,j);//(r_k[0]-r_j[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_j[1])*(r_k[1]-r_j[1]);
                //r_kj = sqrt(r_kj);
                double denom_kj = (1+beta*r_kj);
                jastrowSum2 += (dim-1)*m_a(k,j)/(r_kj*denom_kj*denom_kj) - 2*m_a(k,j)*beta/(denom_kj*denom_kj*denom_kj);

                for (int i=0; i < k; i++) {
                    std::vector<double> r_i = particles[i]->getPosition();
                    double r_ki = m_distances(k,i);//(r_k[0]-r_i[0])*(r_k[0]-r_i[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_i[1]);
                    //r_ki = sqrt(r_ki);
                    double denom_ki = (1+beta*r_ki);

                    //double factor1 = (r_k[0]-r_i[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_j[1]);
                    double factor1 = 0;
                    for (int d = 0; d < m_numberOfDimensions; d++) {
                        factor1 += (r_k[d]-r_i[d])*(r_k[d]-r_j[d]);
                    }

                    factor1 /= r_ki*r_kj;
                    double factor2 = m_a(k,i)/(denom_ki*denom_ki) * m_a(k,j)/(denom_kj*denom_kj);
                    jastrowSum1 += factor1*factor2;
                }
                for (int i=k+1; i < m_numberOfParticles; i++) {
                    std::vector<double> r_i = particles[i]->getPosition();
                    double r_ki = m_distances(k,i);//(r_k[0]-r_i[0])*(r_k[0]-r_i[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_i[1]);
                    //r_ki = sqrt(r_ki);
                    double denom_ki = (1+beta*r_ki);

                    //double factor1 = (r_k[0]-r_i[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_j[1]);
                    double factor1 = 0;
                    for (int d = 0; d < m_numberOfDimensions; d++) {
                        factor1 += (r_k[d]-r_i[d])*(r_k[d]-r_j[d]);
                    }

                    factor1 /= r_ki*r_kj;
                    double factor2 = m_a(k,i)/(denom_ki*denom_ki) * m_a(k,j)/(denom_kj*denom_kj);
                    jastrowSum1 += factor1*factor2;
                }
            }
            for (int j=k+1; j < m_numberOfParticles; j++) {
                std::vector<double> r_j = particles[j]->getPosition();
                double r_kj = m_distances(k,j);//(r_k[0]-r_j[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_j[1])*(r_k[1]-r_j[1]);
                //r_kj = sqrt(r_kj);
                double denom_kj = (1+beta*r_kj);
                //jastrowSum2 += (d-1)*m_a(k,j)/(r_kj*denom_kj*denom_kj) - 2*m_a(k,j)*beta/(denom_kj*denom_kj*denom_kj);

                for (int i=0; i < k; i++) {
                    std::vector<double> r_i = particles[i]->getPosition();
                    double r_ki = m_distances(k,i);//(r_k[0]-r_i[0])*(r_k[0]-r_i[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_i[1]);
                    //r_ki = sqrt(r_ki);
                    double denom_ki = (1+beta*r_ki);

                    //double factor1 = (r_k[0]-r_i[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_j[1]);
                    double factor1 = 0;
                    for (int d = 0; d < m_numberOfDimensions; d++) {
                        factor1 += (r_k[d]-r_i[d])*(r_k[d]-r_j[d]);
                    }

                    factor1 /= r_ki*r_kj;
                    double factor2 = m_a(k,i)/(denom_ki*denom_ki) * m_a(k,j)/(denom_kj*denom_kj);
                    jastrowSum1 += factor1*factor2;
                }
                for (int i=k+1; i < m_numberOfParticles; i++) {
                    std::vector<double> r_i = particles[i]->getPosition();
                    double r_ki = m_distances(k,i);//(r_k[0]-r_i[0])*(r_k[0]-r_i[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_i[1]);
                    //r_ki = sqrt(r_ki);
                    double denom_ki = (1+beta*r_ki);

                    //double factor1 = (r_k[0]-r_i[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_j[1]);
                    double factor1 = 0;
                    for (int d = 0; d < m_numberOfDimensions; d++) {
                        factor1 += (r_k[d]-r_i[d])*(r_k[d]-r_j[d]);
                    }

                    factor1 /= r_ki*r_kj;
                    double factor2 = m_a(k,i)/(denom_ki*denom_ki) * m_a(k,j)/(denom_kj*denom_kj);
                    jastrowSum1 += factor1*factor2;
                }
            }
        }

        jastrowLaplacian = jastrowSum1 + 2*jastrowSum2;

        for (int d=0; d < m_numberOfDimensions; d++) {
            for (int i=0; i < m_numberOfParticles; i++) {
                crossTerm += computeSlaterGradient(i)[d]*m_JastrowGrad(i,d);//computeJastrowGradient(particles, i)[d];
            }
        }
    }

    double laplacian = slaterLaplacian + jastrowLaplacian + 2*crossTerm;
    return laplacian;
    //return 0;
}

std::vector<double> ManyElectronsCoefficients::computeDerivativeWrtParameters(std::vector<Particle *> particles){
    // Calculates the derivative w.r.t. alpha and beta for the interacting wave function using the analytical expression.

    std::vector<double> derivative(2);
    double slaterUpAlphaDerivative = 0;
    double slaterDownAlphaDerivative = 0;

    for (int i=0; i < m_halfNumberOfParticles; i++) {
        std::vector<double> rSpinUp = particles[i]->getPosition();
//        double xSpinUp = rSpinUp[0];
//        double ySpinUp = rSpinUp[1];
        std::vector<double> rSpinDown = particles[i+m_halfNumberOfParticles]->getPosition();
//        double xSpinDown = rSpinDown[0];
//        double ySpinDown = rSpinDown[1];

        for (int j=0; j < m_halfNumberOfParticles; j++) {
//            int nx = m_quantumNumbers(j, 0);
//            int ny = m_quantumNumbers(j, 1);

            vec n(m_numberOfDimensions);
            for (int d = 0; d < m_numberOfDimensions; d++) {
                n[d] = m_quantumNumbers(j, d);
            }

            slaterUpAlphaDerivative += m_system->getHamiltonian()->computeSPWFAlphaDerivative(n, rSpinUp, j)
                                       *m_spinUpSlaterInverse(j,i);
            slaterDownAlphaDerivative += m_system->getHamiltonian()->computeSPWFAlphaDerivative(n, rSpinDown, j)
                                         *m_spinDownSlaterInverse(j,i);
        }
    }

    double beta = m_parameters[1];
    double exponent = 0;
    double betaDerivative = 0;
    //double r_ij = 0;

    for (int i=0; i < m_numberOfParticles; i++) {
        //std::vector<double> r_i = particles[i]->getPosition();

        for (int j=i+1; j < m_numberOfParticles; j++) {
            //std::vector<double> r_j = particles[j]->getPosition();

            //for (int k=0; k < m_numberOfDimensions; k++) {
            //    r_ij += (r_i[k]-r_j[k])*(r_i[k]-r_j[k]);
            //}
            //r_ij = sqrt(r_ij);
            double r_ij = m_distances(i,j);

            exponent += m_a(i,j)*r_ij/(1. + beta*r_ij);
            double denom = (1+beta*r_ij);
            betaDerivative -= m_a(i,j)*r_ij*r_ij/(denom*denom);
        }
    }

    derivative[0] = (slaterUpAlphaDerivative + slaterDownAlphaDerivative)*evaluate(particles);//*exp(exponent);
    derivative[1] = betaDerivative*evaluate(particles);

    return derivative;
    //return 0;
}

double ManyElectronsCoefficients::computeMetropolisRatio(std::vector<Particle *> particles,
                                            int currentParticle, std::vector<double> positionChange) {
    // Function for calculating the wave function part of the Metropolis ratio,
    // both the Slater part and the Jastrow part.

    //std::vector<double> positionOld = particles[currentParticle]->getPosition();
    m_distancesOld = m_distances;

    for (int i=0; i<m_numberOfDimensions; i++) {
        particles[currentParticle]->adjustPosition(positionChange[i], i);
    }

    //std::vector<double> positionNew = particles[currentParticle]->getPosition();
    m_system->getWaveFunction()->updateDistances(currentParticle);
    m_system->getWaveFunction()->updateSPWFMat(currentParticle);
    m_system->getWaveFunction()->updateJastrow(currentParticle);

    int i = currentParticle;
    double ratioSlaterDet = 0;

    if (i < m_halfNumberOfParticles) {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            //int nx = m_quantumNumbers(j, 0);
            //int ny = m_quantumNumbers(j, 1);
            ratioSlaterDet += m_spinUpSlaterInverse(j,i)
                             *m_SPWFMat(i,j);//evaluateSingleParticleWF(nx, ny, positionNew[0], positionNew[1]);
        }
    }
    else {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            //int nx = m_quantumNumbers(j, 0);
            //int ny = m_quantumNumbers(j, 1);
            ratioSlaterDet += m_spinDownSlaterInverse(j, i-m_halfNumberOfParticles)
                             *m_SPWFMat(i,j);//evaluateSingleParticleWF(nx, ny, positionNew[0], positionNew[1]);
        }
    }

    double beta = m_parameters[1];
    double exponent = 0;

    if (m_Jastrow) {
        for (int j=0; j < i; j++) {
            //double r_ijNew = 0;
            //double r_ijOld = 0;
            //std::vector<double> r_j = particles[j]->getPosition();

            //for (int d=0; d < m_numberOfDimensions; d++) {
                //r_ijNew += (positionNew[d] - r_j[d])
                //          *(positionNew[d] - r_j[d]);
                //r_ijOld += (positionOld[d] - r_j[d])
                //          *(positionOld[d] - r_j[d]);
            //}
            //r_ijNew = sqrt(r_ijNew);
            //r_ijOld = sqrt(r_ijOld);
            double r_ijNew = m_distances(i,j);
            double r_ijOld = m_distancesOld(i,j);

            exponent += m_a(i,j)*r_ijNew / (1. + beta*r_ijNew);
            exponent -= m_a(i,j)*r_ijOld / (1. + beta*r_ijOld);
        }
        for (int j=i+1; j < m_numberOfParticles; j++) {
            //double r_ijNew = 0;
            //double r_ijOld = 0;
            //std::vector<double> r_j = particles[j]->getPosition();

            //for (int d=0; d < m_numberOfDimensions; d++) {
                //r_ijNew += (positionNew[d] - r_j[d])
                //          *(positionNew[d] - r_j[d]);
                //r_ijOld += (positionOld[d] - r_j[d])
                //          *(positionOld[d] - r_j[d]);
            //}
            //r_ijNew = sqrt(r_ijNew);
            //r_ijOld = sqrt(r_ijOld);
            double r_ijNew = m_distances(i,j);
            double r_ijOld = m_distancesOld(i,j);

            exponent += m_a(i,j)*r_ijNew / (1. + beta*r_ijNew);
            exponent -= m_a(i,j)*r_ijOld / (1. + beta*r_ijOld);
        }
    }
    double ratioJastrowFactor = exp(exponent);
    m_ratioSlaterDet = ratioSlaterDet;
    m_metropolisRatio = ratioSlaterDet*ratioJastrowFactor;

    return m_metropolisRatio;

}

void ManyElectronsCoefficients::setUpSlaterDetOneParticle() {

    m_system->retrieveCoefficientsFromFile("../diagonalization/PlotAndData/Coefficients.dat", m_cCoefficients);

    vec cColVec = m_cCoefficients.slice(0).col(0);
    int nMaxCoeff = cColVec.size();
    int numberOfEigstates;

    rowvec cRowVec = m_cCoefficients.slice(0).row(0);
    m_nPrimeMax = cRowVec.size();

    int nMax = 0;

    if (m_numberOfDimensions == 1) {
        nMax = m_halfNumberOfParticles;
    }
    else if (m_numberOfDimensions == 2) {
        for (int n=1; n<=m_halfNumberOfParticles; n++) {
            numberOfEigstates = 0.5*n*n + 0.5*n;
            if (numberOfEigstates == m_halfNumberOfParticles) {
                nMax = n;
                break;
            }
        }
    }
    else if (m_numberOfDimensions == 3) {
        for (int n=1; n<=m_halfNumberOfParticles; n++) {
            numberOfEigstates = (n*n*n + 3*n*n +2*n)/6.;
            if (numberOfEigstates == m_halfNumberOfParticles) {
                nMax = n;
                break;
            }
        }
    }

    if (nMaxCoeff > nMax) {
        nMax = nMaxCoeff;
    }

    if (m_numberOfDimensions == 2) {
        //numberOfEigstates = int(0.5*(nMax+1)*(nMax+2));
        numberOfEigstates = int(0.5*(nMax)*(nMax+1));
    }

    else if (m_numberOfDimensions == 3) {
        //numberOfEigstates = int((nMax+1)*(nMax+2)*(nMax+3)/6.);
        numberOfEigstates = int((nMax)*(nMax+1)*(nMax+2)/6.);
    }

    else { numberOfEigstates = nMax; }

    m_numberOfEigstates = numberOfEigstates;


    m_quantumNumbers = zeros<mat>(numberOfEigstates, m_numberOfDimensions);

    if (m_numberOfDimensions == 1) {
        for (int p = 0; p < numberOfEigstates; p++) {
            m_quantumNumbers(p, 0) = p;
        }
    }

    else if (m_numberOfDimensions == 2) {
        int n = 0;
        int nx = 0;
        int ny = 0;

        for (int p = 0; p < numberOfEigstates; p++) {
            m_quantumNumbers(p, 0) = nx;    m_quantumNumbers(p, 1) = ny;
            if (ny == n) {
                n++;
                nx = n;
                ny = 0;
            }
            else {
                nx--;
                ny++;
            }
        }
    }
    else {
        int j = 0;
        for (int i=0; i<nMax; i++) {
            for (int nx = 0; nx <= i; nx++) {
                for (int ny = 0; ny <= i; ny++) {
                    for (int nz = 0; nz <= i; nz++) {
                        if ( (nx+ny+nz) == i ) {
                            m_quantumNumbers(j,0) = nx;
                            m_quantumNumbers(j,1) = ny;
                            m_quantumNumbers(j,2) = nz;
                            j++;
                        }
                    }
                }
            }
        }

//        int i = 0;

//        for (int nx = 0; nx < 2; nx++) {
//            for (int ny = 0; ny < 2; ny++) {
//                for (int nz = 0; nz < 2; nz++) {
//                    if (nx+ny+nz < 2) {
//                        m_quantumNumbers(i,0) = nx;
//                        m_quantumNumbers(i,1) = ny;
//                        m_quantumNumbers(i,2) = nz;
//                        i++;
//                    }
//                }
//            }
//        }

//        for (int nx = 0; nx < nMax; nx++) {
//            for (int ny = 0; ny < nMax; ny++) {
//                for (int nz = 0; nz < nMax; nz++) {
//                    if (nx+ny+nz < nMax) {
//                        m_quantumNumbersC(i,0) = nx;
//                        m_quantumNumbersC(i,1) = ny;
//                        m_quantumNumbersC(i,2) = nz;
//                        i++;
//                    }
//                }
//            }
//        }
    }

    //----Square well----
    if (m_system->getSquareWellFlag()) {
        for (int d = 0; d < m_numberOfDimensions; d++) {
            m_quantumNumbers(0, d) += 1;
        }
    }
    //----Square well----

    if (m_cCoefficients.slice(0).is_square()) {
        mat cCoeffProd = m_cCoefficients.slice(0);
        for (int d = 1; d < m_numberOfDimensions; d++) {
            cCoeffProd %= m_cCoefficients.slice(d);
        }
        m_cDeterminant = det(cCoeffProd);
    }
    else { m_cDeterminant = 1.; }
    cout << m_cDeterminant << endl;

    // Below m_numberOfParticle instead of m_halfNumberOfParticles. Can't have half of one particle.
    m_spinUpSlater = zeros<mat>(m_numberOfParticles, m_numberOfParticles);
    m_spinDownSlater = zeros<mat>(m_numberOfParticles, m_numberOfParticles);

    m_SPWFMat = zeros<mat>(m_numberOfParticles, m_numberOfParticles);
    m_SPWFDMat = field<vec>(m_numberOfParticles, m_numberOfParticles);
    m_SPWFDDMat = zeros<mat>(m_numberOfParticles, m_numberOfParticles);

    double alpha = m_parameters[0];
    m_system->getHamiltonian()->setAlpha(alpha);

    m_SPWFDMat(0,0) = zeros<vec>(m_numberOfDimensions);

    std::vector<double> r = m_system->getParticles()[0]->getPosition();
    double r2 = 0;

    vec n(m_numberOfDimensions);

    for (int d = 0; d < m_numberOfDimensions; d++) {
        n[d] = m_quantumNumbers(0, d);
        r2 += r[d]*r[d];
    }

    double expFactor = exp(-alpha*m_omega*(r2)*0.5);

    m_system->getHamiltonian()->setExpFactor(expFactor);

    m_expFactorsDim = zeros(m_numberOfDimensions);
    // Precalculate exp factors for each dim to be used in super position of coeffs + eigstates.
    for (int d = 0; d < m_numberOfDimensions; d++) {
        m_expFactorsDim[d] = exp(-alpha*m_omega*(r[d]*r[d])*0.5);
    }

    vec nTemp(m_numberOfDimensions);
//    //m_spinUpSlater(i,j) = m_cDeterminant*evaluateSingleParticleWF(n, rSpinUp);
//    for (int m = 0; m < m_nPrimeMax; m++) {
//        nTemp = conv_to<vec>::from(m_quantumNumbers.row(0));
//        double C = 1;
//        for (int d = 0; d < m_numberOfDimensions; ++d) {
//            C *= m_cCoefficients(nTemp(d), m, d);
//            //nTemp[d] = m_quantumNumbers.col(m);
//        }

//        m_spinUpSlater(0,0) += C*m_system->getHamiltonian()->evaluateSingleParticleWF(nTemp, r, 0);
//        vec temp = m_system->getHamiltonian()->computeSPWFDerivative(nTemp, r, 0);
//        temp *= C;
//        m_SPWFDMat(0,0) += temp;
//        m_SPWFDDMat(0,0) += C*m_system->getHamiltonian()->computeSPWFDoubleDerivative(nTemp, r, 0);
//    }

    for (int eig = 0; eig < m_numberOfEigstates; eig++) {
        double term = 1;
        vec termD(m_numberOfDimensions);
        double termDD = 0;
        double coefficients = 1;
        vec qNums = conv_to<vec>::from(m_quantumNumbers.row(eig));
        for (int d = 0; d < m_numberOfDimensions; d++) {
            term *= harmonicOscillatorBasis(r[d], qNums[d], d); //SWITCH OUT m for n[d].
            termD[d] = harmonicOscillatorBasisDerivative(r, qNums, d);
            termDD += harmonicOscillatorBasisDoubleDerivative(r, qNums, d);
            coefficients *= m_cCoefficients(qNums[d], 0, d);
        }
        m_spinUpSlater(0,0) += coefficients*term;
        m_SPWFDMat(0,0) += coefficients*termD;
        m_SPWFDDMat(0,0) += coefficients*termDD;
    }

    m_SPWFMat(0,0) = m_spinUpSlater(0,0);
    //cout << m_SPWFMat << "     " << r[0] << endl;
    m_spinDownSlater(0,0) = m_spinUpSlater(0,0);

    m_spinUpSlaterInverse = m_spinUpSlater.i();
    m_spinDownSlaterInverse = m_spinDownSlater.i();
}

void ManyElectronsCoefficients::setUpSlaterDet() {
    // Function for setting up the Slater determinant at the begining of the simulation.

    m_system->retrieveCoefficientsFromFile("../diagonalization/PlotAndData/Coefficients.dat", m_cCoefficients);

    vec cColVec = m_cCoefficients.slice(0).col(0);
    int nMaxCoeff = cColVec.size();
    int numberOfEigstates = 0;

    rowvec cRowVec = m_cCoefficients.slice(0).row(0);
    m_nPrimeMax = cRowVec.size();
    assert(m_nPrimeMax >= m_halfNumberOfParticles); //nPrimeMax needs to be larger than half the number of particles.

    int nMax = 0;

    if (m_numberOfDimensions == 1) {
        nMax = m_halfNumberOfParticles;
    }
    else if (m_numberOfDimensions == 2) {
        for (int n=1; n<=m_halfNumberOfParticles; n++) {
            numberOfEigstates = 0.5*n*n + 0.5*n;
            if (numberOfEigstates == m_halfNumberOfParticles) {
                nMax = n;
                break;
            }
        }
    }
    else if (m_numberOfDimensions == 3) {
        for (int n=1; n<=m_halfNumberOfParticles; n++) {
            numberOfEigstates = (n*n*n + 3*n*n +2*n)/6.;
            if (numberOfEigstates == m_halfNumberOfParticles) {
                nMax = n;
                break;
            }
        }
    }

    if (nMaxCoeff > nMax) {
        nMax = nMaxCoeff;
    }

    if (m_numberOfDimensions == 2) {
        //numberOfEigstates = int(0.5*(nMax+1)*(nMax+2));
        numberOfEigstates = int(0.5*(nMax)*(nMax+1));
    }

    else if (m_numberOfDimensions == 3) {
        //numberOfEigstates = int((nMax+1)*(nMax+2)*(nMax+3)/6.);
        numberOfEigstates = int((nMax)*(nMax+1)*(nMax+2)/6.);
    }

    else { numberOfEigstates = nMax; }

    m_numberOfEigstates = numberOfEigstates;


    m_quantumNumbers = zeros<mat>(numberOfEigstates, m_numberOfDimensions);
    m_quantumNumbersDouble = zeros<mat>(numberOfEigstates, m_numberOfDimensions);

    if (m_numberOfDimensions == 1) {
        for (int p = 0; p < numberOfEigstates; p++) {
            m_quantumNumbers(p, 0) = p;
        }
    }

    else if (m_numberOfDimensions == 2) {
        int n = 0;
        int nx = 0;
        int ny = 0;

        for (int p = 0; p < numberOfEigstates; p++) {
            m_quantumNumbers(p, 0) = nx;    m_quantumNumbers(p, 1) = ny;
            if (ny == n) {
                n++;
                nx = n;
                ny = 0;
            }
            else {
                nx--;
                ny++;
            }
        }
    }
    else {
        int j = 0;
        for (int i=0; i<nMax; i++) {
            for (int nx = 0; nx <= i; nx++) {
                for (int ny = 0; ny <= i; ny++) {
                    for (int nz = 0; nz <= i; nz++) {
                        if ( (nx+ny+nz) == i ) {
                            m_quantumNumbers(j,0) = nx;
                            m_quantumNumbers(j,1) = ny;
                            m_quantumNumbers(j,2) = nz;
                            j++;
                        }
                    }
                }
            }
        }

//        int i = 0;

//        for (int nx = 0; nx < 2; nx++) {
//            for (int ny = 0; ny < 2; ny++) {
//                for (int nz = 0; nz < 2; nz++) {
//                    if (nx+ny+nz < 2) {
//                        m_quantumNumbers(i,0) = nx;
//                        m_quantumNumbers(i,1) = ny;
//                        m_quantumNumbers(i,2) = nz;
//                        i++;
//                    }
//                }
//            }
//        }

//        for (int nx = 0; nx < nMax; nx++) {
//            for (int ny = 0; ny < nMax; ny++) {
//                for (int nz = 0; nz < nMax; nz++) {
//                    if (nx+ny+nz < nMax) {
//                        m_quantumNumbersC(i,0) = nx;
//                        m_quantumNumbersC(i,1) = ny;
//                        m_quantumNumbersC(i,2) = nz;
//                        i++;
//                    }
//                }
//            }
//        }
    }

    bool overlappingWells = true;
    if (m_system->getDoubleWellFlag()) {
        for (int d = 0; d < m_numberOfDimensions; d++) {
            if (m_system->getL()(d) != 0) { overlappingWells = false; }
        }
    }

    if (m_system->getDoubleWellFlag() /*&& m_numberOfParticles > 2*/ && !overlappingWells) {
        mat quantumNumbersDoubleWell = zeros<mat>(m_numberOfEigstates, m_numberOfDimensions);
        m_quantumNumbersDouble = m_quantumNumbers;
        for (int p = 0; p < m_numberOfEigstates; p+=2) {
            for (int d = 0; d < m_numberOfDimensions; d++) {
                quantumNumbersDoubleWell(p, d) = m_quantumNumbers(p/2, d);
                if (p+1 < m_numberOfEigstates) {
                    quantumNumbersDoubleWell(p+1, d) = m_quantumNumbers(p/2, d);
                }
            }
        }

        m_quantumNumbers = quantumNumbersDoubleWell;
    }
    else { m_quantumNumbersDouble = m_quantumNumbers; }
    m_quantumNumbers = m_quantumNumbersDouble;
    //m_quantumNumbersDouble = m_quantumNumbers;
    //cout << m_quantumNumbers << endl;
    //cout << m_quantumNumbersDouble << endl;

//    //----Square well----
//    if (m_system->getSquareWellFlag()) {
//        for (int p = 0; p < m_halfNumberOfParticles; p++) {
//            for (int d = 0; d < m_numberOfDimensions; d++) {
//                m_quantumNumbers(p, d) += 1;
//            }
//        }
//    }
//    //----Square well----

    if (m_cCoefficients.slice(0).is_square()) {
        mat cCoeffProd = m_cCoefficients.slice(0);
        for (int d = 1; d < m_numberOfDimensions; d++) {
            cCoeffProd %= m_cCoefficients.slice(d);
        }
        m_cDeterminant = det(cCoeffProd);
    }
    else { m_cDeterminant = 1.; }
    //cout << m_cDeterminant << endl;

    m_a = zeros<mat>(m_numberOfParticles, m_numberOfParticles);
    int half = m_halfNumberOfParticles;
    for (int i=0; i < m_numberOfParticles; i++) {
        for (int j=0; j < m_numberOfParticles; j++) {
            if ( ((i < half) && (j < half)) || ((i >= half) && (j >= half)) ) {
                if (m_numberOfDimensions == 2) { m_a(i,j) = 1./3; }
                else if (m_numberOfDimensions == 3) { m_a(i,j) = 1./4; }
            }
            else {
                if (m_numberOfDimensions == 2) { m_a(i,j) = 1.; }
                else if (m_numberOfDimensions == 3) { m_a(i,j) = 1./2; }
            }
        }
    }

    m_spinUpSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);
    m_spinDownSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);

    m_SPWFMat = zeros<mat>(m_numberOfParticles, m_halfNumberOfParticles);
    m_SPWFDMat = field<vec>(m_numberOfParticles, m_halfNumberOfParticles);
    m_SPWFDDMat = zeros<mat>(m_numberOfParticles, m_halfNumberOfParticles);

    double alpha = m_parameters[0];
    m_system->getHamiltonian()->setAlpha(alpha);

    for (int i=0; i<m_numberOfParticles; i++) {
        for (int j=0; j<m_halfNumberOfParticles; j++) {
            m_SPWFDMat(i,j) = zeros<vec>(m_numberOfDimensions);
        }
    }

    m_expFactorsDim = zeros(m_numberOfDimensions);
    for (int i=0; i < m_halfNumberOfParticles; i++) {
        std::vector<double> rSpinUp = m_system->getParticles()[i]->getPosition();
//        double xSpinUp = rSpinUp[0];
//        double ySpinUp = rSpinUp[1];
        std::vector<double> rSpinDown = m_system->getParticles()[i+m_halfNumberOfParticles]->getPosition();
//        double xSpinDown = rSpinDown[0];
//        double ySpinDown = rSpinDown[1];

        for (int j=0; j < m_halfNumberOfParticles; j++) {
//            nx = m_quantumNumbers(j, 0);
//            ny = m_quantumNumbers(j, 1);
            vec n(m_numberOfDimensions);
            double r2SpinUp = 0;
            double r2SpinDown = 0;
            for (int d = 0; d < m_numberOfDimensions; d++) {
                if (m_system->getDoubleWellFlag()) {
                    n[d] = m_quantumNumbersDouble(j, d);
                    //cout << n[d] << endl;
                }
                else {
                    n[d] = m_quantumNumbers(j, d);
                }


                r2SpinUp += rSpinUp[d]*rSpinUp[d];
                r2SpinDown += rSpinDown[d]*rSpinDown[d];
            }

            double expFactor = exp(-alpha*m_omega*(r2SpinUp)*0.5);
            m_system->getHamiltonian()->setExpFactor(expFactor);

            // Precalculate exp factors for each dim to be used in super position of coeffs + eigstates.
            for (int d = 0; d < m_numberOfDimensions; d++) {
                m_expFactorsDim[d] = exp(-alpha*m_omega*(rSpinUp[d]*rSpinUp[d])*0.5);
            }

//            vec nTemp(m_numberOfDimensions);
//            //m_spinUpSlater(i,j) = m_cDeterminant*evaluateSingleParticleWF(n, rSpinUp);
//            for (int m = 0; m < m_nPrimeMax/*m_numberOfEigstates*/; m++) {
//                nTemp = conv_to<vec>::from(m_quantumNumbers.row(j));
//                double C = 1;
//                for (int d = 0; d < m_numberOfDimensions; ++d) {
//                    C *= m_cCoefficients(nTemp(d), m, d);
//                    //nTemp[d] = m_quantumNumbers.col(m);
//                }

//                m_spinUpSlater(i,j) += C*m_system->getHamiltonian()->evaluateSingleParticleWF(nTemp, rSpinUp, j);
//                vec temp = m_system->getHamiltonian()->computeSPWFDerivative(nTemp, rSpinUp, j);
//                temp *= C;
//                m_SPWFDMat(i,j) += temp;
//                m_SPWFDDMat(i,j) += C*m_system->getHamiltonian()->computeSPWFDoubleDerivative(nTemp, rSpinUp, j);
//            }


//            m_SPWFMat(i,j) = m_spinUpSlater(i,j);

            //m_SPWFDMat(i,j) = computeSPWFDerivative(n, rSpinUp);
            //m_SPWFDMat(i,j) *= m_cDeterminant;

            //m_SPWFDDMat(i,j) = m_cDeterminant*computeSPWFDoubleDerivative(n, rSpinUp);

            for (int eig = 0; eig < m_numberOfEigstates; eig++) {
                double term = 1;
                vec termD(m_numberOfDimensions);
                double termDD = 0;
                double coefficients = 1;
                vec qNums = conv_to<vec>::from(m_quantumNumbers.row(eig));
                for (int d = 0; d < m_numberOfDimensions; d++) {
                    term *= m_cCoefficients(qNums[d], n[d], d)*harmonicOscillatorBasis(rSpinUp[d], qNums[d], d);
                    termD[d] = harmonicOscillatorBasisDerivative(rSpinUp, qNums, d);
                    termDD += harmonicOscillatorBasisDoubleDerivative(rSpinUp, qNums, d);
                    coefficients *= m_cCoefficients(qNums[d], n[d], d);
                }
                m_spinUpSlater(i,j) += term;
                m_SPWFDMat(i,j) += coefficients*termD;
                m_SPWFDDMat(i,j) += coefficients*termDD;
            }

            m_SPWFMat(i,j) = m_spinUpSlater(i,j);


            expFactor = exp(-alpha*m_omega*(r2SpinUp)*0.5);
            m_system->getHamiltonian()->setExpFactor(expFactor);

            // Precalculate exp factors for each dim to be used in super position of coeffs + eigstates.
            for (int d = 0; d < m_numberOfDimensions; d++) {
                m_expFactorsDim[d] = exp(-alpha*m_omega*(rSpinDown[d]*rSpinDown[d])*0.5);
            }

            //m_spinDownSlater(i,j) = m_cDeterminant*evaluateSingleParticleWF(n, rSpinDown);

//            for (int m = 0; m < m_nPrimeMax/*m_numberOfEigstates*/; m++) {
//                nTemp = conv_to<vec>::from(m_quantumNumbers.row(j));
//                double C = 1;
//                for (int d = 0; d < m_numberOfDimensions; ++d) {
//                    C *= m_cCoefficients(nTemp(d), m, d);
//                    //nTemp[d] = m_quantumNumbers.col(m);
//                }
//                m_spinDownSlater(i,j) += C*m_system->getHamiltonian()->evaluateSingleParticleWF(nTemp, rSpinDown, j);
//                vec temp = m_system->getHamiltonian()->computeSPWFDerivative(nTemp, rSpinDown, j);
//                temp *= C;
//                m_SPWFDMat(i+m_halfNumberOfParticles,j) += temp;
//                m_SPWFDDMat(i+m_halfNumberOfParticles,j) += C*m_system->getHamiltonian()->computeSPWFDoubleDerivative(nTemp, rSpinDown, j);
//            }

//            m_SPWFMat(i+m_halfNumberOfParticles, j) = m_spinDownSlater(i,j);

            for (int eig = 0; eig < m_numberOfEigstates; eig++) {
                double term = 1;
                vec termD(m_numberOfDimensions);
                double termDD = 0;
                double coefficients = 1;
                vec qNums = conv_to<vec>::from(m_quantumNumbers.row(eig));
                for (int d = 0; d < m_numberOfDimensions; d++) {
                    term *= harmonicOscillatorBasis(rSpinDown[d], qNums[d], d);
                    termD[d] = harmonicOscillatorBasisDerivative(rSpinDown, qNums, d);
                    termDD += harmonicOscillatorBasisDoubleDerivative(rSpinDown, qNums, d);
                    coefficients *= m_cCoefficients(qNums[d], n[d], d);
                }
                m_spinDownSlater(i,j) += coefficients*term;
                m_SPWFDMat(i+m_halfNumberOfParticles, j) += coefficients*termD;
                m_SPWFDDMat(i+m_halfNumberOfParticles, j) += coefficients*termDD;
            }

            m_SPWFMat(i+m_halfNumberOfParticles, j) = m_spinDownSlater(i,j);

            //m_SPWFDMat(i+m_halfNumberOfParticles,j) = computeSPWFDerivative(n, rSpinDown);
            //m_SPWFDMat(i+m_halfNumberOfParticles,j) *= m_cDeterminant;

            //m_SPWFDDMat(i+m_halfNumberOfParticles,j) = m_cDeterminant*computeSPWFDoubleDerivative(n, rSpinDown);

        }
    }


    //cout << m_spinUpSlater << endl;
    m_spinUpSlaterInverse = m_spinUpSlater.i();
    m_spinDownSlaterInverse = m_spinDownSlater.i();
}

void ManyElectronsCoefficients::setUpDistances() {

    m_distances = zeros<mat>(m_numberOfParticles, m_numberOfParticles);

    for (int i=0; i<m_numberOfParticles; i++) {
        std::vector<double> r_i = m_system->getParticles()[i]->getPosition();

        for (int j=i+1; j<m_numberOfParticles; j++) {
            std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
            double r_ij = 0;

            for (int d = 0; d < m_numberOfDimensions; d++) {
                r_ij += (r_i[d]-r_j[d])*(r_i[d]-r_j[d]);
            }
            m_distances(i,j) = m_distances(j,i) = sqrt(r_ij);
        }
    }
}

void ManyElectronsCoefficients::setUpJastrowMat() {

    m_JastrowMat = zeros<cube>(m_numberOfParticles, m_numberOfParticles, m_numberOfDimensions);
    m_JastrowGrad = zeros(m_numberOfParticles, m_numberOfDimensions);
    double beta = m_parameters[1];

    for (int i=0; i<m_numberOfParticles; i++) {
        std::vector<double> r_i = m_system->getParticles()[i]->getPosition();
        for (int j=0; j < i; j++) {
            std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
            double r_ij = m_distances(i,j);
            double denom = 1 + beta*r_ij;

            for (int d = 0; d < m_numberOfDimensions; d++) {
                m_JastrowGrad(i,d) += (r_i[d]-r_j[d])/r_ij * m_a(i, j)/(denom*denom);
                //m_JastrowGrad(i,1) += (r_i[1]-r_j[1])/r_ij * m_a(i, j)/(denom*denom);
            }
        }

        for (int j=i+1; j < m_numberOfParticles; j++) {
            std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
            double r_ij = m_distances(i,j);
            double denom = 1 + beta*r_ij;

            for (int d = 0; d < m_numberOfDimensions; d++) {
                m_JastrowMat(i,j,d) = (r_i[d]-r_j[d])/r_ij * m_a(i, j)/(denom*denom);
                //m_JastrowMat(i,j,1) = (r_i[1]-r_j[1])/r_ij * m_a(i, j)/(denom*denom);
                m_JastrowMat(j,i,d) = -m_JastrowMat(i,j,d);
                //m_JastrowMat(j,i,1) = -m_JastrowMat(i,j,1);

                m_JastrowGrad(i,d) += m_JastrowMat(i,j,d);
                //m_JastrowGrad(i,1) += m_JastrowMat(i,j,1);
            }
        }
    }
}

void ManyElectronsCoefficients::updateSlaterDet(int currentParticle) {
    // Function for updating the Slater determinant after every accepted metropolis step.
    int i = currentParticle;
    //std::vector<double> r_i = m_system->getParticles()[i]->getPosition();

    if (i < m_halfNumberOfParticles) {
        mat spinUpSlaterInverseOld = m_spinUpSlaterInverse;

        for (int j=0; j < i; j++) {
            double sum = 0;

            for (int l=0; l <m_halfNumberOfParticles; l++) {
                //int nx = m_quantumNumbers(l, 0);
                //int ny = m_quantumNumbers(l, 1);
                sum += m_SPWFMat(i,l)//evaluateSingleParticleWF(nx, ny, r_i[0], r_i[1])
                      *spinUpSlaterInverseOld(l,j);
            }
            for (int k=0; k < m_halfNumberOfParticles; k++) {
                m_spinUpSlaterInverse(k,j) = spinUpSlaterInverseOld(k,j)
                                            -(sum/m_ratioSlaterDet)*spinUpSlaterInverseOld(k,i);
            }
        }
        for (int j=i+1; j < m_halfNumberOfParticles; j++) {
            double sum = 0;

            for (int l=0; l <m_halfNumberOfParticles; l++) {
                //int nx = m_quantumNumbers(l, 0);
                //int ny = m_quantumNumbers(l, 1);
                sum += m_SPWFMat(i,l)//evaluateSingleParticleWF(nx, ny, r_i[0], r_i[1])
                      *spinUpSlaterInverseOld(l,j);
            }
            for (int k=0; k < m_halfNumberOfParticles; k++) {
                m_spinUpSlaterInverse(k,j) = spinUpSlaterInverseOld(k,j)
                                            -(sum/m_ratioSlaterDet)*spinUpSlaterInverseOld(k,i);
            }
        }
        for (int k=0; k < m_halfNumberOfParticles; k++) {
            m_spinUpSlaterInverse(k,i) = spinUpSlaterInverseOld(k,i)/m_ratioSlaterDet;
        }
    }
    else {
        double iHalf = i-m_halfNumberOfParticles;
        mat spinDownSlaterInverseOld = m_spinDownSlaterInverse;

        for (int j=0; j < iHalf; j++) {
            double sum = 0;

            for (int l=0; l < m_halfNumberOfParticles; l++) {
                //int nx = m_quantumNumbers(l, 0);
                //int ny = m_quantumNumbers(l, 1);
                sum += m_SPWFMat(i,l)//evaluateSingleParticleWF(nx, ny, r_i[0], r_i[1])
                      *spinDownSlaterInverseOld(l,j);
            }
            for (int k=0; k < m_halfNumberOfParticles; k++) {
                m_spinDownSlaterInverse(k,j) = spinDownSlaterInverseOld(k,j)
                                              -(sum/m_ratioSlaterDet)
                                               *spinDownSlaterInverseOld(k, iHalf);
            }
        }
        for (int j=iHalf+1; j < m_halfNumberOfParticles; j++) {
            double sum = 0;

            for (int l=0; l < m_halfNumberOfParticles; l++) {
                //int nx = m_quantumNumbers(l, 0);
                //int ny = m_quantumNumbers(l, 1);
                sum += m_SPWFMat(i,l)//evaluateSingleParticleWF(nx, ny, r_i[0], r_i[1])
                      *spinDownSlaterInverseOld(l,j);
            }
            for (int k=0; k < m_halfNumberOfParticles; k++) {
                m_spinDownSlaterInverse(k,j) = spinDownSlaterInverseOld(k,j)
                                              -(sum/m_ratioSlaterDet)
                                               *spinDownSlaterInverseOld(k, iHalf);
            }
        }
        for (int k=0; k < m_halfNumberOfParticles; k++) {
            m_spinDownSlaterInverse(k, iHalf) = spinDownSlaterInverseOld(k, iHalf)/m_ratioSlaterDet;
        }
    }
}

void ManyElectronsCoefficients::updateDistances(int currentParticle) {
    // Function for updating the distances between electrons.
    int i = currentParticle;
    std::vector<double> r_i = m_system->getParticles()[i]->getPosition();

    for (int j=0; j<i; j++) {
        std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
        double r_ij = 0;

        for (int d = 0; d < m_numberOfDimensions; d++) {
            r_ij += (r_i[d]-r_j[d])*(r_i[d]-r_j[d]);
        }
        m_distances(i,j) = m_distances(j,i) = sqrt(r_ij);
    }

    for (int j=i+1; j<m_numberOfParticles; j++) {
        std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
        double r_ij = 0;

        for (int d = 0; d < m_numberOfDimensions; d++) {
            r_ij += (r_i[d]-r_j[d])*(r_i[d]-r_j[d]);
        }
        m_distances(i,j) = m_distances(j,i) = sqrt(r_ij);
    }

}

void ManyElectronsCoefficients::updateSPWFMat(int currentParticle) {

    int i = currentParticle;
    std::vector<double> r_i = m_system->getParticles()[i]->getPosition();
    double alpha = m_parameters[0];

    double r2 = 0;
    for (int d = 0; d < m_numberOfDimensions; d++) {
        r2 += r_i[d]*r_i[d];
    }

    double expFactor = exp(-alpha*m_omega*(r2)*0.5);
    m_system->getHamiltonian()->setExpFactor(expFactor);

    // Precalculate exp factors for each dim to be used in super position of coeffs + eigstates.
    for (int d = 0; d < m_numberOfDimensions; d++) {
        m_expFactorsDim[d] = exp(-alpha*m_omega*(r_i[d]*r_i[d])*0.5);
    }

    if (m_numberOfParticles == 1) {
//        vec n(m_numberOfDimensions);
//        for (int d = 0; d < m_numberOfDimensions; d++) {
//            n[d] = m_quantumNumbers(0, d);
//        }

        m_SPWFMat(0,0) = 0;

        m_SPWFDMat(0,0) = zeros(m_SPWFDMat(0,0).size());

        m_SPWFDDMat(0,0) = 0;

//        vec nTemp(m_numberOfDimensions);
//        for (int m = 0; m < m_nPrimeMax; m++) {
//            nTemp = conv_to<vec>::from(m_quantumNumbers.row(0));
//            double C = 1;
//            for (int d = 0; d < m_numberOfDimensions; ++d) {
//                C *= m_cCoefficients(nTemp(d), m, d);
//                //nTemp[d] = m_quantumNumbers.col(m);
//            }

//            m_SPWFMat(0,0) += C*m_system->getHamiltonian()->evaluateSingleParticleWF(nTemp, r_i, 0);
//            vec temp = m_system->getHamiltonian()->computeSPWFDerivative(nTemp, r_i, 0);
//            temp *= C;
//            m_SPWFDMat(0,0) += temp;
//            m_SPWFDDMat(0,0) += C*m_system->getHamiltonian()->computeSPWFDoubleDerivative(nTemp, r_i, 0);
//        }
        for (int eig = 0; eig < m_numberOfEigstates; eig++) {
            double term = 1;
            vec termD(m_numberOfDimensions);
            double termDD = 0;
            double coefficients = 1;
            vec qNums = conv_to<vec>::from(m_quantumNumbers.row(eig));
            for (int d = 0; d < m_numberOfDimensions; d++) {
                term *= m_cCoefficients(qNums[d], 0, d)*harmonicOscillatorBasis(r_i[d], qNums[d], d);
                termD[d] = harmonicOscillatorBasisDerivative(r_i, qNums, d);
                termDD += harmonicOscillatorBasisDoubleDerivative(r_i, qNums, d);
                coefficients *= m_cCoefficients(qNums[d], 0, d);
                //cout << m_cCoefficients(qNum, m, d) << endl;
            }
            m_SPWFMat(0,0) += term;
            m_SPWFDMat(0,0) += coefficients*termD;
            m_SPWFDDMat(0,0) += coefficients*termDD;
        }
        //cout << m_SPWFMat << "     " << r_i[0] << endl;
        //m_SPWFMat = sqrt(m_SPWFMat%m_SPWFMat/dot(m_SPWFMat,m_SPWFMat));
    }
    else {
        for (int j=0; j<m_halfNumberOfParticles; j++) {
//          int nx = m_quantumNumbers(j, 0);
//          int ny = m_quantumNumbers(j, 1);
            vec n(m_numberOfDimensions);
            for (int d = 0; d < m_numberOfDimensions; d++) {
                if (m_system->getDoubleWellFlag()) {
                    n[d] = m_quantumNumbersDouble(j, d);
                }
                else {
                    n[d] = m_quantumNumbers(j, d);
                }
            }

            m_SPWFMat(i,j) = 0;

            m_SPWFDMat(i,j) = zeros(m_SPWFDMat(i,j).size());

            m_SPWFDDMat(i,j) = 0;

//            vec nTemp(m_numberOfDimensions);
//            for (int m = 0; m < m_nPrimeMax/*m_numberOfEigstates*/; m++) {
//                nTemp = conv_to<vec>::from(m_quantumNumbers.row(j));
//                double C = 1;
//                for (int d = 0; d < m_numberOfDimensions; ++d) {
//                    C *= m_cCoefficients(nTemp(d), m, d);
//                    //nTemp[d] = m_quantumNumbers.col(m);
//                }

//                m_SPWFMat(i,j) += C*m_system->getHamiltonian()->evaluateSingleParticleWF(nTemp, r_i, j);
//                vec temp = m_system->getHamiltonian()->computeSPWFDerivative(nTemp, r_i, j);
//                temp *= C;
//                m_SPWFDMat(i,j) += temp;
//                m_SPWFDDMat(i,j) += C*m_system->getHamiltonian()->computeSPWFDoubleDerivative(nTemp, r_i, j);
//            }


            for (int eig = 0; eig < m_numberOfEigstates; eig++) {
                double term = 1;
                vec termD(m_numberOfDimensions);
                double termDD = 0;
                double coefficients = 1;
                vec qNums = conv_to<vec>::from(m_quantumNumbers.row(eig));
                for (int d = 0; d < m_numberOfDimensions; d++) {
                    term *= harmonicOscillatorBasis(r_i[d], qNums[d], d);
                    termD[d] = harmonicOscillatorBasisDerivative(r_i, qNums, d);
                    termDD += harmonicOscillatorBasisDoubleDerivative(r_i, qNums, d);
                    coefficients *= m_cCoefficients(qNums[d], n[d], d);
                    //cout << m_cCoefficients(qNum, m, d) << endl;
                }
                m_SPWFMat(i,j) += coefficients*term;
                m_SPWFDMat(i,j) += coefficients*termD;
                m_SPWFDDMat(i,j) += coefficients*termDD;
            }
            //cout << m_SPWFMat << "     " << r_i[0] << endl;

//          m_SPWFMat(i,j) = m_cDeterminant*evaluateSingleParticleWF(n, r_i);

//          m_SPWFDMat(i,j) = computeSPWFDerivative(n, r_i);
//          m_SPWFDMat(i,j) *= m_cDeterminant;

//          m_SPWFDDMat(i,j) = m_cDeterminant*computeSPWFDoubleDerivative(n, r_i);
        }
    }

}

void ManyElectronsCoefficients::updateJastrow(int currentParticle) {

    int p = currentParticle;
    std::vector<double> r_p = m_system->getParticles()[p]->getPosition();
    double beta = m_parameters[1];
    m_JastrowMatOld = m_JastrowMat;

    for (int j=0; j<p; j++) {
        std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
        double r_pj = m_distances(p,j);
        double denom = 1 + beta*r_pj;
        for (int d = 0; d < m_numberOfDimensions; d++) {
            m_JastrowMat(p,j,d) = (r_p[d]-r_j[d])/r_pj * m_a(p, j)/(denom*denom);
            //m_JastrowMat(p,j,1) = (r_p[1]-r_j[1])/r_pj * m_a(p, j)/(denom*denom);
            m_JastrowMat(j,p,d) = -m_JastrowMat(p,j,d);
            //m_JastrowMat(j,p,1) = -m_JastrowMat(p,j,1);
        }
    }
    for (int j=p+1; j<m_numberOfParticles; j++) {
        std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
        double r_pj = m_distances(p,j);
        double denom = 1 + beta*r_pj;
        for (int d = 0; d < m_numberOfDimensions; d++) {
            m_JastrowMat(p,j,d) = (r_p[d]-r_j[d])/r_pj * m_a(p, j)/(denom*denom);
            //m_JastrowMat(p,j,1) = (r_p[1]-r_j[1])/r_pj * m_a(p, j)/(denom*denom);
            m_JastrowMat(j,p,d) = -m_JastrowMat(p,j,d);
            //m_JastrowMat(j,p,1) = -m_JastrowMat(p,j,1);
        }
    }

    m_JastrowGradOld = m_JastrowGrad;

    for (int d = 0; d < m_numberOfDimensions; d++) {
        m_JastrowGrad(p, d) = 0;
        //m_JastrowGrad(p, 1) = 0;

        for (int j=0; j<p; j++) {
            m_JastrowGrad(p, d) += m_JastrowMat(p,j,d);
            //m_JastrowGrad(p, 1) += m_JastrowMat(p,j,1);
        }
        for (int j=p+1; j<m_numberOfParticles; j++) {
            m_JastrowGrad(p, d) += m_JastrowMat(p,j,d);
            //m_JastrowGrad(p, 1) += m_JastrowMat(p,j,1);
        }
        for (int i=0; i<p; i++) {
            m_JastrowGrad(i, d) = m_JastrowGradOld(i,d) - m_JastrowMatOld(i,p,d) + m_JastrowMat(i,p,d);
            //m_JastrowGrad(i, 1) = m_JastrowGradOld(i,1) - m_JastrowMatOld(i,p,1) + m_JastrowMat(i,p,1);
        }
        for (int i=p+1; i<m_numberOfParticles; i++) {
            m_JastrowGrad(i, d) = m_JastrowGradOld(i,d) - m_JastrowMatOld(i,p,d) + m_JastrowMat(i,p,d);
            //m_JastrowGrad(i, 1) = m_JastrowGradOld(i,1) - m_JastrowMatOld(i,p,1) + m_JastrowMat(i,p,1);
        }
    }
}

double ManyElectronsCoefficients::harmonicOscillatorBasis(double x, int nx, int d) {

    double nFac = factorial(nx);

    double n2 = pow(2., nx);
    double pi4 = pow(M_PI, -0.25);
    double omega4 = pow(m_omega, 0.25);
    double constant = omega4*pi4/sqrt(nFac*n2);

    //double xAbs2 = x*x;

    double wavefunc = m_expFactorsDim[d];//exp(-0.5*m_omega*xAbs2);

    double phi = constant*wavefunc*m_system->getHamiltonian()->computeHermitePolynomial(nx, x);

    return phi;


//    double x2 = x*x;
//    double expFactor = exp(-0.5*m_omega*x2);

//    double waveFunction = expFactor*m_system->getHamiltonian()->computeHermitePolynomial(nx, x);

//    return waveFunction;
}

double ManyElectronsCoefficients::harmonicOscillatorBasisDerivative(vec r, vec n, int d) {

    vec constant(m_numberOfDimensions);
    for (int dim = 0; dim < m_numberOfDimensions; dim++) {
        double nFac = factorial(n[dim]);
        double n2 = pow(2., n[dim]);
        double pi4 = pow(M_PI, -0.25);
        double omega4 = pow(m_omega, 0.25);
        constant[dim] = omega4*pi4/sqrt(nFac*n2);
    }

    //double xAbs2 = r[d]*r[d];

    double wavefunc = m_expFactorsDim[d];//exp(-0.5*m_omega*xAbs2);

    double hermitePolynomial = m_system->getHamiltonian()->computeHermitePolynomial(n[d], r[d]);
    double hermitePolynomialDerivative = m_system->getHamiltonian()->computeHermitePolynomialDerivative(n[d], r[d]);

    double phiD = constant[d]*(wavefunc*hermitePolynomialDerivative
                           -m_omega*r[d]*wavefunc*hermitePolynomial);
    for (int dim = 0; dim < m_numberOfDimensions; dim++) {
        if (d != dim) phiD *= constant[dim]*m_expFactorsDim[dim]*m_system->getHamiltonian()->computeHermitePolynomial(n[dim], r[dim]);
    }

    return phiD;
}

double ManyElectronsCoefficients::harmonicOscillatorBasisDoubleDerivative(vec r, vec n, int d) {

    vec constant(m_numberOfDimensions);
    for (int dim = 0; dim < m_numberOfDimensions; dim++) {
        double nFac = factorial(n[dim]);
        double n2 = pow(2., n[dim]);
        double pi4 = pow(M_PI, -0.25);
        double omega4 = pow(m_omega, 0.25);
        constant[dim] = omega4*pi4/sqrt(nFac*n2);
    }

    //double xAbs2 = r[d]*r[d];

    double wavefunc = m_expFactorsDim[d];//exp(-0.5*m_omega*xAbs2);

    double hermitePolynomial = m_system->getHamiltonian()->computeHermitePolynomial(n[d], r[d]);
    double hermitePolynomialDerivative = m_system->getHamiltonian()->computeHermitePolynomialDerivative(n[d], r[d]);
    double hermitePolynomialDoubleDerivative = m_system->getHamiltonian()->computeHermitePolynomialDoubleDerivative(n[d], r[d]);

    double phiDD = constant[d]*(wavefunc*hermitePolynomialDoubleDerivative
                            -m_omega*wavefunc*hermitePolynomial
                            -2*m_omega*r[d]*wavefunc*hermitePolynomialDerivative
                            +m_omega*m_omega*r[d]*r[d]*wavefunc*hermitePolynomial);
    for (int dim = 0; dim < m_numberOfDimensions; dim++) {
        if (d != dim) phiDD *= constant[dim]*m_expFactorsDim[dim]*m_system->getHamiltonian()->computeHermitePolynomial(n[dim], r[dim]);
    }

    return phiDD;
}

//double ManyElectronsCoefficients::computeHermitePolynomial(int nValue, double position) {
//    // Computes Hermite polynomials.
//    double alphaSqrt = 1;//sqrt(m_alpha);
//    double omegaSqrt = sqrt(m_omega);
//    double factor = 2*alphaSqrt*omegaSqrt*position;

//    double HermitePolynomialPP = 0;                 // H_{n-2}
//    double HermitePolynomialP = 1;                  // H_{n-1}
//    double HermitePolynomial = HermitePolynomialP;  // H_n

//    for (int n=1; n <= nValue; n++) {
//        HermitePolynomial = factor*HermitePolynomialP - 2*(n-1)*HermitePolynomialPP;
//        HermitePolynomialPP = HermitePolynomialP;
//        HermitePolynomialP = HermitePolynomial;
//    }

//    return HermitePolynomial;

//}

//double ManyElectronsCoefficients::computeHermitePolynomialDerivative(int nValue, double position) {
//    // Computes Hermite polynomials differentiated w.r.t. position.
//    double alphaSqrt = 1;//sqrt(m_alpha);
//    double omegaSqrt = sqrt(m_omega);
//    double factor1 = 2*alphaSqrt*omegaSqrt;
//    double factor2 = 2*alphaSqrt*omegaSqrt*position;

//    double HPDerivativePP = 0;              // d/dx H_{n-2}
//    double HPDerivativeP = 0;               // d/dx H_{n-1}
//    double HPDerivative = HPDerivativeP;    // d/dx H_n

//    for (int n=1; n <= nValue; n++) {
//        HPDerivative = factor1*computeHermitePolynomial(n-1, position)
//                      +factor2*HPDerivativeP
//                      -2*(n-1)*HPDerivativePP;
//        HPDerivativePP = HPDerivativeP;
//        HPDerivativeP = HPDerivative;
//    }

//    return HPDerivative;

//}

//double ManyElectronsCoefficients::computeHermitePolynomialDoubleDerivative(int nValue, double position) {
//    // Computes Hermite polynomials twice differentiated w.r.t. position.
//    double alphaSqrt = 1;//sqrt(m_alpha);
//    double omegaSqrt = sqrt(m_omega);
//    double factor1 = 4*alphaSqrt*omegaSqrt;
//    double factor2 = 2*alphaSqrt*omegaSqrt*position;

//    double HPDoubleDerivativePP = 0;                    // d/dx d/dx H_{n-2}
//    double HPDoubleDerivativeP = 0;                     // d/dx d/dx H_{n-1}
//    double HPDoubleDerivative = HPDoubleDerivativeP;    // d/dx d/dx H_n

//    for (int n=1; n <= nValue; n++) {
//        HPDoubleDerivative = factor1*computeHermitePolynomialDerivative(n-1, position)
//                            +factor2*HPDoubleDerivativeP
//                            -2*(n-1)*HPDoubleDerivativePP;
//        HPDoubleDerivativePP = HPDoubleDerivativeP;
//        HPDoubleDerivativeP = HPDoubleDerivative;
//    }

//    return HPDoubleDerivative;

//}
