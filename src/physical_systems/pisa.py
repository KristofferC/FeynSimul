class Pisa:

    def __init__(self):
	"""Initiates a system with two particles interacting in 3D with 
        gaussian potentials.

	The system is written in Jacobi coordinates so only two three degrees
        of freedom. The potential is taken from arXiv:1106.3853v1 
        [physics.atm-clus] and can be found at U{http://arxiv.org/abs/1106.3853.}
        """


        self.userCode = """ inline float sqr(float x){return x * x;}
        		    inline float quad(float x){return x * x * x * x;} """
        
        self.energyOp =  """ native_exp(-0.8604556619274776f * (sqr(x1) + sqr(x2) + sqr(x3)))*(-1.2269999999999996f + 1.0557790971850147f * (sqr(x1) + sqr(x2) + sqr(x3))) """

        self.potential = """-1.2269999999999996f * native_exp(-0.8604556619274776f * (sqr(x1) + sqr(x2) + sqr(x3)))"""
 
        self.DOF = 3

    # groundStateEnergy=-0.001302
