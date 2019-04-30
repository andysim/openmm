#ifndef OPENMM_DRUDENOSEHOOVERCHAININTEGRATOR_H_
#define OPENMM_DRUDENOSEHOOVERCHAININTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2019 Stanford University and the Authors.      *
 * Authors: Andreas Kraemer and Andrew C. Simmonett                           *
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

#include "openmm/Integrator.h"
#include "openmm/Kernel.h"
#include "openmm/internal/windowsExportDrude.h"

namespace OpenMM {

/**
 *
 * Nose Hoover chain velocity Verlet integrator.
 *
 * This Integrator can be used to simulate systems that include Drude particles. If no Drude particles
 * are present in the system, it thermostats all particles in the system using a single reference temperature.
 * If Drude Particles are present, it applies two different Nose-Hoover chain
 * thermostats to different parts of the system.  The first is applied to ordinary particles (ones that
 * are not part of a Drude particle pair), as well as to the center of mass of each Drude particle pair.
 * A second thermostat, typically with a much lower temperature, is applied to the relative internal
 * displacement of each pair.
 *
 * This integrator can optionally set an upper limit on how far any Drude particle is ever allowed to
 * get from its parent particle.  This can sometimes help to improve stability.  The limit is enforced
 * with a hard wall constraint.
 * 
 * This Integrator requires the System to include a DrudeForce, which it uses to identify the Drude
 * particles.
 * TODO: no requirement of Drude force?!
 */

class OPENMM_EXPORT_DRUDE DrudeNoseHooverChainIntegrator : public Integrator {
public:
    /**
     * Create a DrudeNoseHooverChainIntegrator.  See 
     *
     *     G. J. Martyna, M. E. Tuckerman, D. J. Tobias, and Michael L. Klein  "Explicit reversible
     *     integrators for extended systems dynamics", Molecular Physics, 87 1117-1157 (1996)
     *     http://dx.doi.org/10.1080/00268979600100761
     *
     *     and
     *
     *     G. J. Martyna, M. L. Klein, and M. Tuckerman "Nose–Hoover chains: The canonical
     *     ensemble via continuous dynamics",
     *     Journal of Chemical Physics 97, 2635-2643 (1992)
     *     http://dx.doi.org/10.1063/1.463940
     *
     * @param temperature    the temperature of the main heat bath (in Kelvin)
     * @param frictionCoeff  the friction coefficient which couples the system to the main heat bath (in inverse picoseconds)
     * @param drudeTemperature    the temperature of the heat bath applied to internal coordinates of Drude particles (in Kelvin)
     * @param drudeFrictionCoeff  the friction coefficient which couples the system to the heat bath applied to internal coordinates of Drude particles (in inverse picoseconds)
     * @param stepSize       the step size with which to integrator the system (in picoseconds)
     * @param chainLength    the number of beads in the thermostat chain
     * @param numYoshidaSuzuki the number of terms in the Yoshida-Suzuki decomposition (allowed values are 1,3,5): higher value increases both accuracy and expense
     * @param numMTS         the number of terms used in the multi time step decomposition in the integrator: higher value increases both accuracy and expense
     */
    DrudeNoseHooverChainIntegrator(double temperature, double frictionCoeff, double drudeTemperature, double drudeFrictionCoeff, double stepSize, int chainLength=10, int numYoshidaSuzuki=3, int numMTS=3);
    /**
     * Get the temperature of the main heat bath (in Kelvin).
     *
     * @return the temperature of the heat bath, measured in Kelvin
     */
    double getTemperature() const {
        return temperature;
    }
    /**
     * Set the temperature of the main heat bath (in Kelvin).
     *
     * @param temp    the temperature of the heat bath, measured in Kelvin
     */
    void setTemperature(double temp) {
        temperature = temp;
    }
    /**
     * Get the friction coefficient which determines how strongly the system is coupled to
     * the main heat bath (in inverse ps).
     *
     * @return the friction coefficient, measured in 1/ps
     */
    double getFriction() const {
        return friction;
    }
    /**
     * Set the friction coefficient which determines how strongly the system is coupled to
     * the main heat bath (in inverse ps).
     *
     * @param coeff    the friction coefficient, measured in 1/ps
     */
    void setFriction(double coeff) {
        friction = coeff;
    }
    /**
     * Get the temperature of the heat bath applied to internal coordinates of Drude particles (in Kelvin).
     *
     * @return the temperature of the heat bath, measured in Kelvin
     */
    double getDrudeTemperature() const {
        return drudeTemperature;
    }
    /**
     * Set the temperature of the heat bath applied to internal coordinates of Drude particles (in Kelvin).
     *
     * @param temp    the temperature of the heat bath, measured in Kelvin
     */
    void setDrudeTemperature(double temp) {
        drudeTemperature = temp;
    }
    /**
     * Get the friction coefficient which determines how strongly the internal coordinates of Drude particles
     * are coupled to the heat bath (in inverse ps).
     *
     * @return the friction coefficient, measured in 1/ps
     */
    double getDrudeFriction() const {
        return drudeFriction;
    }
    /**
     * Set the friction coefficient which determines how strongly the internal coordinates of Drude particles
     * are coupled to the heat bath (in inverse ps).
     *
     * @param coeff    the friction coefficient, measured in 1/ps
     */
    void setDrudeFriction(double coeff) {
        drudeFriction = coeff;
    }
    /**
     * Get the length of the Nose-Hoover chain of beads (unitless)
     */
    int getChainLength() const {
        return chainLength;
    }
    /**
     * Set the length of the Nose-Hoover chain of beads (unitless)
     *
     * @param length    the number of beads in the Nose-Hoover chain
     */
    void setChainLength(int length) {
        chainLength = length;
    }
    /**
     * Get the number of terms in the Yoshida-Suzuki decomposition
     */
    int getNumYoshidaSuzuki() const {
        return numYoshidaSuzuki;
    }
    /**
     * Set the number of terms in the Yoshida-Suzuki decomposition
     *
     * @param numTerms  the number of terms in the Yoshida-Suzuki decomposition
     */
    void setNumYoshidaSuzuki(int numTerms) {
        numYoshidaSuzuki = numTerms;
    }
    /**
     * Get the number of terms in the multi time step decomposition
     */
    int getNumMTS() const {
        return numMTS;
    }
    /**
     * Set the number of terms in the multi time step decomposition
     *
     * @param numTerms   the number of terms in the multi time step decomposition
     */
    void setNumMTS(int numTerms) {
        numMTS = numTerms;
    }
    /**
     * Get the maximum distance a Drude particle can ever move from its parent particle, measured in nm.  This is implemented
     * with a hard wall constraint.  If this distance is set to 0 (the default), the hard wall constraint is omitted.
     */
    double getMaxDrudeDistance() const;
    /**
     * Set the maximum distance a Drude particle can ever move from its parent particle, measured in nm.  This is implemented
     * with a hard wall constraint.  If this distance is set to 0 (the default), the hard wall constraint is omitted.
     */
    void setMaxDrudeDistance(double distance);
    /**
     * Advance a simulation through time by taking a series of time steps.
     *
     * @param steps   the number of time steps to take
     */
    void step(int steps);
protected:
    /**
     * This will be called by the Context when it is created.  It informs the Integrator
     * of what context it will be integrating, and gives it a chance to do any necessary initialization.
     * It will also get called again if the application calls reinitialize() on the Context.
     */
    void initialize(ContextImpl& context);
    /**
     * This will be called by the Context when it is destroyed to let the Integrator do any necessary
     * cleanup.  It will also get called again if the application calls reinitialize() on the Context.
     */
    void cleanup();
    /**
     * When the user modifies the state, we need to mark that the forces need to be recalculated.
     */
    void stateChanged(State::DataType changed);
    /**
     * Get the names of all Kernels used by this Integrator.
     */
    std::vector<std::string> getKernelNames();
    /**
     * Compute the kinetic energy of the system at the current time.
     */
    double computeKineticEnergy();
private:
    double temperature, friction, drudeTemperature, drudeFriction, maxDrudeDistance;
    int chainLength, numYoshidaSuzuki, numMTS;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_DRUDENOSEHOOVERCHAININTEGRATOR_H_*/
