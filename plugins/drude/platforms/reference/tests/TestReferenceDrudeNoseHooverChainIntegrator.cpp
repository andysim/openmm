/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2016 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**
 * This tests the Reference implementation of DrudeNoseHooverChainIntegrator.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include "openmm/DrudeForce.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/DrudeKernels.h"
#include "openmm/DrudeNoseHooverChainIntegrator.h"
#include "openmm/internal/NoseHooverChainThermostatImpl.h"
#include "openmm/NoseHooverChainThermostat.h"
#include "SimTKOpenMMRealType.h"
#include "SimTKOpenMMUtilities.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerDrudeReferenceKernelFactories();

void testNoseHooverChainPropagation(){
    /* test if the kinetic energy of one degree of freedom converges 
     * to the reference value */
    double velocity;
    const double mass = 1;
    double timestep = 0.001;
    DrudeNoseHooverChainIntegrator dummy_integrator(1,1,1,1,1);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, dummy_integrator, platform);
    NoseHooverChainThermostat thermostat(300, 1);
    NoseHooverChainThermostatImpl thermostat_impl(thermostat);
    ReferenceNoseHooverChainThermostatPropagateKernel(NoseHooverChainThermostatPropagateKernel::Name(), platform);
}

/*
void testHarmonicOscillatorNoDrude() {
    
    const double temperature = 1/0.008314472471220215;
    const double mass = 1.0;
    const double temperatureDrude = 10.0; // TODO:reconsider 

    System system;    
    system.addParticle(mass);
    CustomExternalForce* harmonic_restraint = new CustomExternalForce("x^2+y^2+z^2");
    system.addForce(harmonic_restraint);
    vector<Vec3> positions(1);
    positions[0] = Vec3(0.7, 0.8,0.6);
    vector<Vec3> velocities(1);
    velocities[0] = Vec3(0.0, 0, 0);
    DrudeNoseHooverChainIntegrator integrator(temperature, 1, temperatureDrude, 1, 0.001, 4, 1, 1);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocities(velocities);

    double mean_temperature=0;
    // equilibration
    integrator.step(100000);
    for (size_t i=0; i < 10000; i++){
        integrator.step(100);
        State state = context.getState(State::Energy | State::Velocities | State::Positions);
        const vector<Vec3>& vel = state.getVelocities();
        double kinetic_energy = state.getKineticEnergy();
        double temp = kinetic_energy/(0.5*3*BOLTZ);
        mean_temperature = (i*mean_temperature + temp)/(i+1);
        cout << kinetic_energy << " " << temp << " " << temperature << endl;
    }
    cout << "Mean Temperature: " << mean_temperature;
    ASSERT_USUALLY_EQUAL_TOL(mean_temperature, temperature, 2.0);
    //ASSERT_USUALLY_EQUAL_TOL(3*0.5*BOLTZ*temperatureDrude, keInternal/numSteps, 0.01);
}

void testSinglePair() {
    const double temperature = 300.0;
    const double temperatureDrude = 10.0;
    const double k = ONE_4PI_EPS0*1.5;
    const double charge = 0.1;
    const double alpha = ONE_4PI_EPS0*charge*charge/k;
    const double mass1 = 1.0;
    const double mass2 = 0.1;
    const double totalMass = mass1+mass2;
    const double reducedMass = (mass1*mass2)/(mass1+mass2);
    const double maxDistance = 0.05;
    System system;
    system.addParticle(mass1);
    system.addParticle(mass2);
    DrudeForce* drude = new DrudeForce();
    drude->addParticle(1, 0, -1, -1, -1, charge, alpha, 1, 1);
    system.addForce(drude);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 0, 0);
    DrudeNoseHooverChainIntegrator integ(temperature, 20.0, temperatureDrude, 20.0, 0.003);
    integ.setMaxDrudeDistance(maxDistance);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    context.setPositions(positions);

    // Equilibrate.

    integ.step(1000);
    // Compute the internal and center of mass temperatures.
    double keCM = 0, keInternal = 0;
    int numSteps = 10000;
    for (int i = 0; i < numSteps; i++) {
        integ.step(10);
        State state = context.getState(State::Velocities | State::Positions);
        const vector<Vec3>& vel = state.getVelocities();
        Vec3 velCM = vel[0]*(mass1/totalMass) + vel[1]*(mass2/totalMass);
        keCM += 0.5*totalMass*velCM.dot(velCM);
        Vec3 velInternal = vel[0]-vel[1];
        keInternal += 0.5*reducedMass*velInternal.dot(velInternal);
        Vec3 delta = state.getPositions()[0]-state.getPositions()[1];
        double distance = sqrt(delta.dot(delta));
        ASSERT(distance <= maxDistance*(1+1e-6));
    }
    ASSERT_USUALLY_EQUAL_TOL(3*0.5*BOLTZ*temperature, keCM/numSteps, 0.1);
    ASSERT_USUALLY_EQUAL_TOL(3*0.5*BOLTZ*temperatureDrude, keInternal/numSteps, 0.01);
}

void testWater() {
    // Create a box of SWM4-NDP water molecules.  This involves constraints, virtual sites,
    // and Drude particles.
    
    const int gridSize = 3;
    const int numMolecules = gridSize*gridSize*gridSize;
    const double spacing = 0.6;
    const double boxSize = spacing*(gridSize+1);
    const double temperature = 300.0;
    const double temperatureDrude = 10.0;
    System system;
    NonbondedForce* nonbonded = new NonbondedForce();
    DrudeForce* drude = new DrudeForce();
    system.addForce(nonbonded);
    system.addForce(drude);
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(1.0);
    for (int i = 0; i < numMolecules; i++) {
        int startIndex = system.getNumParticles();
        system.addParticle(15.6); // O
        system.addParticle(0.4);  // D
        system.addParticle(1.0);  // H1
        system.addParticle(1.0);  // H2
        system.addParticle(0.0);  // M
        nonbonded->addParticle(1.71636, 0.318395, 0.21094*4.184);
        nonbonded->addParticle(-1.71636, 1, 0);
        nonbonded->addParticle(0.55733, 1, 0);
        nonbonded->addParticle(0.55733, 1, 0);
        nonbonded->addParticle(-1.11466, 1, 0);
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < j; k++)
                nonbonded->addException(startIndex+j, startIndex+k, 0, 1, 0);
        system.addConstraint(startIndex, startIndex+2, 0.09572);
        system.addConstraint(startIndex, startIndex+3, 0.09572);
        system.addConstraint(startIndex+2, startIndex+3, 0.15139);
        system.setVirtualSite(startIndex+4, new ThreeParticleAverageSite(startIndex, startIndex+2, startIndex+3, 0.786646558, 0.106676721, 0.106676721));
        drude->addParticle(startIndex+1, startIndex, -1, -1, -1, -1.71636, ONE_4PI_EPS0*1.71636*1.71636/(100000*4.184), 1, 1);
    }
    vector<Vec3> positions;
    for (int i = 0; i < gridSize; i++)
        for (int j = 0; j < gridSize; j++)
            for (int k = 0; k < gridSize; k++) {
                Vec3 pos(i*spacing, j*spacing, k*spacing);
                positions.push_back(pos);
                positions.push_back(pos);
                positions.push_back(pos+Vec3(0.09572, 0, 0));
                positions.push_back(pos+Vec3(-0.023999, 0.092663, 0));
                positions.push_back(pos);
            }
    
    // Simulate it and check the temperature.
    
    DrudeNoseHooverChainIntegrator integ(temperature, 50.0, temperatureDrude, 50.0, 0.0005);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    context.setPositions(positions);
    context.applyConstraints(1e-5);
    
    // Equilibrate.
    
    integ.step(500);
    
    // Compute the internal and center of mass temperatures.
    
    double ke = 0;
    int numSteps = 4000;
    for (int i = 0; i < numSteps; i++) {
        integ.step(1);
        ke += context.getState(State::Energy).getKineticEnergy();
    }
    ke /= numSteps;
    int numStandardDof = 3*3*numMolecules-system.getNumConstraints();
    int numDrudeDof = 3*numMolecules;
    int numDof = numStandardDof+numDrudeDof;
    double expectedTemp = (numStandardDof*temperature+numDrudeDof*temperatureDrude)/numDof;
    ASSERT_USUALLY_EQUAL_TOL(expectedTemp, ke/(0.5*numDof*BOLTZ), 0.03);
} 

void testForceEnergyConsistency() {
    // Create a box of polarizable particles.
    
    const int gridSize = 3;
    const int numAtoms = gridSize*gridSize*gridSize;
    const double spacing = 0.6;
    const double boxSize = spacing*(gridSize+1);
    const double temperature = 300.0;
    const double temperatureDrude = 10.0;
    System system;
    vector<Vec3> positions;
    NonbondedForce* nonbonded = new NonbondedForce();
    DrudeForce* drude = new DrudeForce();
    system.addForce(nonbonded);
    system.addForce(drude);
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    nonbonded->setNonbondedMethod(NonbondedForce::PME);
    nonbonded->setCutoffDistance(1.0);
    nonbonded->setUseSwitchingFunction(true);
    nonbonded->setSwitchingDistance(0.9);
    nonbonded->setEwaldErrorTolerance(5e-5);
    for (int i = 0; i < numAtoms; i++) {
        int startIndex = system.getNumParticles();
        system.addParticle(1.0);
        system.addParticle(1.0);
        nonbonded->addParticle(1.0, 0.3, 1.0);
        nonbonded->addParticle(-1.0, 0.3, 1.0);
        nonbonded->addException(startIndex, startIndex+1, 0, 1, 0);
        drude->addParticle(startIndex+1, startIndex, -1, -1, -1, -1.0, 0.001, 1, 1);
    }
    for (int i = 0; i < gridSize; i++)
        for (int j = 0; j < gridSize; j++)
            for (int k = 0; k < gridSize; k++) {
                Vec3 pos(i*spacing, j*spacing, k*spacing);
                positions.push_back(pos);
                positions.push_back(pos);
            }
    
    // Simulate it and check that force and energy remain consistent.
    
    DrudeNoseHooverChainIntegrator integ(temperature, 50.0, temperatureDrude, 50.0, 0.001);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    context.setPositions(positions);
    State prevState;
    for (int i = 0; i < 100; i++) {
        State state = context.getState(State::Energy | State::Forces | State::Positions);
        if (i > 0) {
            double expectedEnergyChange = 0;
            for (int j = 0; j < system.getNumParticles(); j++) {
                Vec3 delta = state.getPositions()[j]-prevState.getPositions()[j];
                expectedEnergyChange -= 0.5*(state.getForces()[j]+prevState.getForces()[j]).dot(delta);
            }
            ASSERT_EQUAL_TOL(expectedEnergyChange, state.getPotentialEnergy()-prevState.getPotentialEnergy(), 0.05);
        }
        prevState = state;
        integ.step(1);
    }
}*/

int main() {
    try {
        //registerDrudeReferenceKernelFactories();

        //testHarmonicOscillatorNoDrude();
        //testSinglePair();
        //testWater();
        //testForceEnergyConsistency();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
