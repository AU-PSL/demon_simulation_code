/*===- ShieldedCoulombForce.h - libSimulation -=================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef SHIELDEDCOULOMBFORCE_H
#define SHIELDEDCOULOMBFORCE_H

#include "Force.h"
#include "VectorCompatibility.h"

class ShieldedCoulombForce : public Force {
public:
	ShieldedCoulombForce(Cloud * const C, const double shieldingConstant);
	~ShieldedCoulombForce();

	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

private:
	double shielding; // [m^-1]
	
	static const double coulomb;

    void force(const cloud_index currentParticle,
               const double iCharge, const double jCharge,
               const double displacementX, const double displacementY);
    void force(const cloud_index currentParticle, const cloud_index iParticle, 
               const doubleV currentCharge, const doubleV iCharge,
               const doubleV displacementX, const doubleV displacementY);
    void forcer(const cloud_index currentParticle, const cloud_index iParticle, 
                const doubleV currentCharge, const doubleV iCharge,
                const doubleV displacementX, const doubleV displacementY);
    
    static void plusEqualr_pd(double * const a, const doubleV b);
    static void minusEqualr_pd(double * const a, const doubleV b);
    
    void coulombForce1(const cloud_index i, const cloud_index j);
    void block1(const cloud_index ilow, const cloud_index ihigh, 
                const cloud_index jlow, const cloud_index jhigh);
    
    void coulombForce2(const cloud_index i, const cloud_index j);
    void block2(const cloud_index ilow, const cloud_index ihigh, 
                const cloud_index jlow, const cloud_index jhigh);
    
    void coulombForce3(const cloud_index i, const cloud_index j);
    void block3(const cloud_index ilow, const cloud_index ihigh, 
                const cloud_index jlow, const cloud_index jhigh);
    
    void coulombForce4(const cloud_index i, const cloud_index j);
    void block4(const cloud_index ilow, const cloud_index ihigh, 
                const cloud_index jlow, const cloud_index jhigh);
};

#endif // SHIELDEDCOULOMBFORCE_H
