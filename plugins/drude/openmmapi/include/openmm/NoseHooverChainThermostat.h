#ifndef NOSEHOOVERCHAINTHERMOSTAT_H_
#define NOSEHOOVERCHAINTHERMOSTAT_H_

#include "openmm/Force.h"
#include "SimTKOpenMMRealType.h"
#include "openmm/OpenMMException.h"

#include <sstream>

namespace OpenMM {

/**
 */
class NoseHooverChainThermostat : public Force {
   public:
   /**
    */
    NoseHooverChainThermostat(double temperature, int numDOFs, double frequency=50.0, int chainLength=10, int numMTS=2, int numYoshidaSuzuki=3) {
        setChainLength(chainLength);
        setNumDOFs(numDOFs);
        setNumYoshidaSuzuki(numYoshidaSuzuki);
        setNumMTS(numMTS);
        setCollisionFrequency(frequency);
        setTemperature(temperature);
    }
    /**
     * Get the temperature of the main heat bath (in Kelvin).
     *
     * @return the temperature of the heat bath, measured in Kelvin
     */
    double getTemperature() const {
        return temperature_;
    }
    /**
     * Set the temperature of the main heat bath (in Kelvin).
     *
     * @param temp    the temperature of the heat bath, measured in Kelvin
     */
    void setTemperature(double temp) {
        temperature_ = temp;
    }
    /**
     * Get the collision frequency which determines how strongly the system is coupled to
     * the main heat bath (in inverse ps).
     *
     * @return the collision frequency, measured in 1/ps
     */
    double getCollisionFrequency() const {
        return frequency_;
    }
    /**
     * Set the collision frequency which determines how strongly the system is coupled to
     * the main heat bath (in inverse ps).
     *
     * @param frequency    the collision frequency, measured in 1/ps
     */
    void setCollisionFrequency(double frequency) {
        frequency_ = frequency;
    }
    /**
     * Get the length of the Nose-Hoover chain of beads (unitless)
     */
    int getChainLength() const {
        return chainLength_;
    }
    /**
     * Get the number of degrees of freedom associated with the particles that this system is acting on 
     */
    int getNumDOFs() const {
        return numDOFs_;
    }
    /**
     * Get the number of terms in the Yoshida-Suzuki decomposition
     */
    int getNumYoshidaSuzuki() const {
        return numYoshidaSuzuki_;
    }
    /**
     * Set the number of terms in the Yoshida-Suzuki decomposition
     *
     * @param numTerms  the number of terms in the Yoshida-Suzuki decomposition
     */
    void setNumYoshidaSuzuki(int numTerms) {
        numYoshidaSuzuki_ = numTerms;
    }
    /**
     * Get the number of terms in the multi time step decomposition
     */
    int getNumMTS() const {
        return numMTS_;
    }
    /**
     * Set the number of terms in the multi time step decomposition
     *
     * @param numTerms   the number of terms in the multi time step decomposition
     */
    void setNumMTS(int numTerms) {
        numMTS_ = numTerms;
    }
    /**
     * set destructor
     */
    virtual ~NoseHooverChainThermostat() {}
    /**
    * This is the name of the parameter which stores the current temperature of the 
    *  heat bath (in Kelvin). 
    */ 
    static const std::string Temperature(std::string suffix) { 
        std::stringstream key;
        key << "NoseHooverTemperature" << suffix; 
        return key.str(); 
    } 
    /** 
    * This is the name of the parameter which store the current collision frequency (in 1/ps). 
    */ 
    static const std::string CollisionFrequency(std::string suffix) {
        std::stringstream key;
        key << "NoseHooverCollisionFrequency" << suffix; 
        return key.str(); 
    } 
    /** 
    * This is the name of the parameter which stores the chain length
    */ 
    static const std::string ChainLength(std::string suffix) {
        std::stringstream key;
        key << "NoseHooverChainLength" << suffix; 
        return key.str(); 
    } 
    /** 
    * This is the name of the parameter which stores the number of degrees of freedom associated with this thermostat 
    */ 
    static const std::string NumDOFs(std::string suffix) {
        std::stringstream key;
        key << "NoseHooverNumDOFs" << suffix; 
        return key.str(); 
    } 
    /** 
    * This is the name of the parameter which stores the number of multi time steps 
    */ 
    static const std::string NumMTS(std::string suffix) {
        std::stringstream key;
        key << "NoseHooverNumMTS" << suffix; 
        return key.str(); 
    } 
    /** 
    * This is the name of the parameter which stores the number of terms in the Yoshida-Suzuki decomposition 
    */ 
    static const std::string NumYoshidaSuzuki(std::string suffix) {
        std::stringstream key;
        key << "NoseHooverNumYoshidaSuzuki" << suffix; 
        return key.str(); 
    } 
    /**
    * This is the name of the parameters which store the forces on Nose-Hoover beads
    */ 
    static const std::string G(std::string suffix, int index) { 
        std::stringstream key;
        key << "NoseHooverG" << suffix << index; 
        return key.str(); 
    } 
    /**
    * This is the name of the parameters which store the masses of Nose-Hoover beads
    */ 
    static const std::string Q(std::string suffix, int index) { 
        std::stringstream key;
        key << "NoseHooverQ" << suffix << index; 
        return key.str(); 
    } 
    /**
    * This is the name of the parameters which store the positions of Nose-Hoover beads
    */ 
    static const std::string Xi(std::string suffix, int index) { 
        std::stringstream key;
        key << "NoseHooverXi" << suffix << index; 
        return key.str(); 
    } 
    /**
    * This is the name of the parameters which store the velocities of Nose-Hoover beads
    */ 
    static const std::string Vxi(std::string suffix, int index) { 
        std::stringstream key;
        key << "NoseHooverVxi" << suffix << index; 
        return key.str(); 
    } 
    /**
     * Compute the weights of the Yoshida-Suzuki decomposition.
     */
    static const std::vector<double> getYoshidaSuzukiWeights(int numYS) {
        switch (numYS) {
            case 1:
                return {1};
            case 3:
                return {0.828981543588751, -0.657963087177502, 0.828981543588751};
            case 5:
                return {0.2967324292201065, 0.2967324292201065, -0.186929716880426, 0.2967324292201065,
                        0.2967324292201065};
            default:
                throw OpenMMException("The number of Yoshida-Suzuki weights must be 1,3, or 5.");
        }
    }
protected:
    virtual ForceImpl* createImpl() const override {
        return NULL;
        //TODO: I guess that this needs to be implemented, in order to provide an instance that lives in a context
    }
    /**
     * Set the number of degrees of freedom associated with the particles that this thermostat is acting on.
     *
     * @param numDOFs the number of translational and rotational degrees of freedom minus the number of constraints.
     */
    void setNumDOFs(int numDOFs) {
        numDOFs_ = numDOFs;
    }
    /**
     * Set the length of the Nose-Hoover chain of beads (unitless)
     *
     * @param length    the number of beads in the Nose-Hoover chain
     */
    void setChainLength(int length) {
        chainLength_ = length;
    }

   private:
    double temperature_, frequency_;
    int chainLength_, numDOFs_, numYoshidaSuzuki_, numMTS_;
};
}

#endif /* NOSEHOOVERCHAINTHERMOSTAT_H_ */
