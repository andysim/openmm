
#ifndef OPENMM_NOSEHOOVERCHAINTHERMOSTATIMPL_H_
#define OPENMM_NOSEHOOVERCHAINTHERMOSTATIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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

#include "openmm/internal/ForceImpl.h"
#include "openmm/NoseHooverChainThermostat.h"
#include "openmm/Kernel.h"
#include <utility>
#include <map>
#include <string>

namespace OpenMM {

/**
 * This is the internal implementation of NoseHooverChainThermostat.
 */

/**
 * NOTE TO SELF:
 *
 * The Impl class is only know to the ContextImpl which triggers
 * updateContextState and calcForcesAndEnergy calls on all forces
 * (calcForcesAndEnergy can also be called on groups of forces). 
 * These calls in turn are triggered by the integrator.
 * 
 * updateContextState is meant to give forces a hook to change state variables
 * and is usually called at the start of each integration step.
 * In the AndersenThermostat, updateContextState triggers the velocity scaling
 * (we cannot do the same thing for NoseHoover, because the propagation has
 * to be called before and after the integration step).
 *
 * Instead, calcForcesAndEnergies will trigger the velocity scaling for NH.
 * 
 * The NoseHooverChainThermostatImpl should be responsible for making the 
 * kernel that propagates kinetic energy through the Nose-Hoover beads.
 *
 * As opposed to our original design decision, the kernel will communicate
 * with the integrator kernel only through the context, rather than
 * passing parameters and return values.
 *
 * All instantaneous data associated with the NH beads should only be 
 * stored in the Context (if for some reason we should revert this decision
 * updateParametersInContext could be called after the NH propagation.)
 *
 * Open questions:
 * - How does NoseHooverChainThermostatImpl::calcForcesAndEnergies know about
 *   the time step? --> has to be a member (again)
 * 
 * The relation between the integrator and the NoseHooverChainThermostat is 
 * as follows:
 * - The constructor of NHCIntegrator will take a system class and inject
 *   the NoseHooverChainThermostat classes as forces with a specific force 
 *   group, and suffixes.
 * - If Drudes are present, one thermostat instance will be created for
 *   each set of particles (Drude and non-Drude), but they can go into
 *   the same force group.
 * -
 *
 */

class NoseHooverChainThermostatImpl : public ForceImpl {
public:
    
    NoseHooverChainThermostatImpl(const NoseHooverChainThermostat& owner, std::string suffix="");
    ~NoseHooverChainThermostatImpl();
    void initialize(ContextImpl& context);
    const NoseHooverChainThermostat& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context, bool& forcesInvalid) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters();
    std::vector<std::string> getKernelNames();
    void updateParametersInContext(ContextImpl& context);
    /**
     * has to be called every time the number of thermostated particles changes (and for initialization)
     */
    void updateQ();

private:
    friend class ReferenceNoseHooverChainThermostatPropagateKernel;
    std::vector<double> Q_, G_, vxi_, xi_;
    const std::string suffix_;
    const NoseHooverChainThermostat& owner;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_NOSEHOOVERCHAINTHERMOSTATIMPL_H_*/
    
