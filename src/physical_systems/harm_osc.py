class HarmOsc:
    
    def __init__(self, k = 1.0,lb=False):
        """Initiates a harmonic oscillator system in 1D.

        Contains the potential and virial energy operator for a harmonic oscillator of
        the form V(x) = 1/2 * k * x ^ 2. 

        @type k: number
        @param k: The spring constant of the harmonic ocsillator.

        """
        self.potential = "0.5f * " + str(k) + "f * x1 * x1"
        if lb:
            self.potential = "0.5f * (1.0f + sqr(eps)/12.0f)*sqr(x1)"
        self.userCode = "float sqr(float x){return x*x;}" 
        self.energyOp = str(k) + "f * x1 * x1"
        if lb:
            self.energyOp = "2.0f * 0.5f*(1.0f + sqr(eps)/4.0f)*sqr(x1)"
        self.DOF = 1
