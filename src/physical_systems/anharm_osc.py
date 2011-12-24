"""
Anharmonic oscillator
"""

class AnHarmOsc:

    def __init__(self, k = 1.0, lambd = 1.0):
        """Initiates an anharmonic harmonic oscillator system in 1D.

        Contains the potential and virial energy operator for an anharmonic oscillator 
        of the form V(x) = 1/2 * k * x + lambd * x ^ 4. 

        @type k: number
        @param k: The spring constant of the harmonic ocsillator.
        @type lambd: number
        @param lambd: Strength of distrubance to the harmonic oscillator.
        """
        self.potential = ("0.5f * x1 * x1 + " + str(lambd) + "f * x1 * x1 * x1 *"
        + "x1")
        self.energyOp = "x1 * x1 + 3.0 * " + str(lambd) + "f * x1 * x1 * x1 * x1"
        self.userCode = ""
        self.DOF = 1

#test = AnHarmOsc(1.0);
#print test.potential
#print test.energyOp
