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

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/NoseHooverChainThermostat.h"
#include "openmm/internal/NoseHooverChainThermostatImpl.h"
#include "openmm/kernels.h"
#include <sstream>

using namespace OpenMM;
using std::map;
using std::pair;
using std::vector;
using std::set;
using std::string;
using std::stringstream;

NoseHooverChainThermostatImpl::NoseHooverChainThermostatImpl(const NoseHooverChainThermostat& owner, int numDOFs, std::string suffix) : owner(owner), suffix_(suffix) {
    int chainLength = owner.getChainLength();
    G_.resize(chainLength);
    vxi_.resize(chainLength);
    xi_.resize(chainLength);
    updateQ();
}

NoseHooverChainThermostatImpl::~NoseHooverChainThermostatImpl() {
}

void NoseHooverChainThermostatImpl::initialize(ContextImpl& context) {
// TODO: copy parameters from context
//
//    kernel = context.getPlatform().createKernel(CalcCustomExternalForceKernel::Name(), context);
//
//    // Check for errors in the specification of bonds.
//
//    const System& system = context.getSystem();
//    vector<double> parameters;
//    int numParameters = owner.getNumPerParticleParameters();
//    for (int i = 0; i < owner.getNumParticles(); i++) {
//        int particle;
//        owner.getParticleParameters(i, particle, parameters);
//        if (particle < 0 || particle >= system.getNumParticles()) {
//            stringstream msg;
//            msg << "CustomExternalForce: Illegal particle index: ";
//            msg << particle;
//            throw OpenMMException(msg.str());
//        }
//        if (parameters.size() != numParameters) {
//            stringstream msg;
//            msg << "CustomExternalForce: Wrong number of parameters for particle ";
//            msg << i;
//            throw OpenMMException(msg.str());
//        }
//    }
//    kernel.getAs<CalcCustomExternalForceKernel>().initialize(context.getSystem(), owner);
}

double NoseHooverChainThermostatImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    return 0.0;
}

vector<string> NoseHooverChainThermostatImpl::getKernelNames() {
    vector<string> names;
    //TODO: reconsider empty vector
    //names.push_back(CalcCustomExternalForceKernel::Name());
    return names;
}

map<string, double> NoseHooverChainThermostatImpl::getDefaultParameters() {
    map<string, double> parameters;
    parameters[NoseHooverChainThermostat::ChainLength(suffix_)] = static_cast<double>(owner.getChainLength());
    parameters[NoseHooverChainThermostat::NumDOFs(suffix_)] = static_cast<double>(owner.getNumDOFs());
    parameters[NoseHooverChainThermostat::NumMTS(suffix_)] = static_cast<double>(owner.getNumMTS());
    parameters[NoseHooverChainThermostat::NumYoshidaSuzuki(suffix_)] = static_cast<double>(owner.getNumYoshidaSuzuki());
    parameters[NoseHooverChainThermostat::Temperature(suffix_)] = owner.getTemperature();
    parameters[NoseHooverChainThermostat::CollisionFrequency(suffix_)] = owner.getCollisionFrequency();
    for (int i = 0; i < owner.getChainLength(); i++){
        parameters[NoseHooverChainThermostat::G(suffix_, i)] = G_[i];
        parameters[NoseHooverChainThermostat::Q(suffix_, i)] = Q_[i];
        parameters[NoseHooverChainThermostat::Xi(suffix_, i)] = xi_[i];
        parameters[NoseHooverChainThermostat::Vxi(suffix_, i)] = vxi_[i];
    }
    return parameters;
}

void NoseHooverChainThermostatImpl::updateParametersInContext(ContextImpl& context) {
    // TODO: I'm not so sure about this. There is another method updateContextState() and I am presently not clear about the differences.
    if (static_cast<int>(context.getParameter(NoseHooverChainThermostat::ChainLength(suffix_))) != owner.getChainLength()) {
        throw OpenMMException("Nose-Hoover chain length was changed.");
    }
    if (static_cast<int>(context.getParameter(NoseHooverChainThermostat::NumDOFs(suffix_))) != owner.getNumDOFs()) {
        throw OpenMMException("Nose-Hoover number of degrees of freedom  was changed.");
    }
    context.setParameter(NoseHooverChainThermostat::NumMTS(suffix_), static_cast<double>(owner.getNumMTS()));
    context.setParameter(NoseHooverChainThermostat::NumYoshidaSuzuki(suffix_), static_cast<double>(owner.getNumYoshidaSuzuki()));
    context.setParameter(NoseHooverChainThermostat::Temperature(suffix_), owner.getTemperature());
    context.setParameter(NoseHooverChainThermostat::CollisionFrequency(suffix_), owner.getCollisionFrequency());
    for (int i = 0; i < owner.getChainLength(); i++){
        context.setParameter(NoseHooverChainThermostat::G(suffix_, i), G_[i]);
        context.setParameter(NoseHooverChainThermostat::Q(suffix_, i), Q_[i]);
        context.setParameter(NoseHooverChainThermostat::Xi(suffix_, i), xi_[i]);
        context.setParameter(NoseHooverChainThermostat::Vxi(suffix_, i), vxi_[i]);
    }
    updateQ();
    // We don't want to call systemChanged() here: we don't invalidate forces, positions, velocities, or energies of the particles in this function.
    //context.systemChanged();
}

void NoseHooverChainThermostatImpl::updateQ(){
    double kT = BOLTZ * owner.getTemperature();
    double frequency = owner.getCollisionFrequency();
    int chainLength = owner.getChainLength();
    int numDOFs = owner.getNumDOFs();
    Q_.resize(chainLength, kT / (frequency * frequency));
    Q_[0] *= numDOFs;
}
