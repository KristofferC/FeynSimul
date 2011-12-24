"""
Coupled harmonic oscillators
"""


class CoupHarmOsc():

    def __init__(self, nbrOfParticles, k = 1.0):
        """Initiates a coupled harmonic oscillator system in 1D.
    
        Contains the potential and virial energy operator for a system with coupled harmonic oscillator.
        Every particle interacts with the others with a harmonic oscillator potential.

        @type nbrOfParticles: integer
        @param nbrOfParticles: The number of particles that are interacting

        @type k: number
        @param k: Strength of interaction between particles
    
        """

	

        # Generate potential
        potential = "0.5f*" + str(k) + "f * ("
        for i in range(0, nbrOfParticles):
            for j in range(0, i):
                extra = '(x' + str(j + 1) + ' - x' + str(i + 1) + ')'
            	potential += extra + '*' + extra + '+'

        self.potential = potential[:-1] + ")"
        self.energyOp = self.potential[5:]
        self.DOF = nbrOfParticles
        self.userCode = ""
