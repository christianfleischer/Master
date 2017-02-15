#include "particle.h"
#include <cassert>

Particle::Particle() {
}

void Particle::setPosition(const std::vector<double> &position) {
    int size = position.size();
    assert(size == m_numberOfDimensions);
    m_position = position;
}

void Particle::adjustPosition(double change, int dimension) {
    m_position.at(dimension) += change;
}

void Particle::setNewPosition(double change, int dimension) {
    m_newPosition = m_position;
    m_newPosition.at(dimension) += change;
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}


SUITE(ParticleTest) {
    class ParticleFixture {
    public:
        Particle* particleTest = new Particle();
        int numberOfDimensions = 2;
        int numberOfParticles = 2;
        void initParticle() {
            std::vector<double> initVector(2);
            initVector[0] = 0.;
            initVector[1] = 0.;

            particleTest->setNumberOfDimensions(numberOfDimensions);
            particleTest->setPosition(initVector);

        }

    };
    TEST_FIXTURE(ParticleFixture, SettingParticlePosition) {
        std::vector<double> testVector(2);
        testVector[0] = 0.;
        testVector[1] = 0.;
        initParticle();
        CHECK_ARRAY_EQUAL(testVector, particleTest->getPosition(), 2);
    }
    TEST_FIXTURE(ParticleFixture, AdjustingParticlePosition) {
        initParticle();
        particleTest->adjustPosition(0.5, 0);

        std::vector<double> testVector(2);
        testVector[0] = 0.5;
        testVector[1] = 0.;
        CHECK_ARRAY_EQUAL(testVector, particleTest->getPosition(), 2);
    }
}
