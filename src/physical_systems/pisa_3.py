# -*- coding: utf-8 -*-


class Pisa:

    """ Sum of two gaussian potentials.Three particle system but two
    degrees
    of freedowm with Jacobi coordinates
    """

    def __init__(self,factor=1.0):

        self.userCode = """ inline float sqr(float x){return x*x;}
        					inline float quad(float x){return x*x*x*x;}
        """
        #self.energyOp ="""native_sqrt(2.0f*x1*x1+2.0f*x2*x2+2.0f*x3*x3)"""
        
        self.energyOp = str(factor)+"""* native_exp(-0.8604556619274776f*sqr(x1) - 0.8604556619274776f*sqr(x2) - 0.8604556619274776f*sqr(x3))*(-1.2269999999999996f + 1.0557790971850147f*sqr(x1) + 1.0557790971850147f*sqr(x2) + 1.0557790971850147f*sqr(x3)) """

        self.potential = str(factor)+"""*-1.2269999999999996f*native_exp(-0.8604556619274776f*sqr(x1) - 0.8604556619274776f*sqr(x2) - 0.8604556619274776f*sqr(x3))"""
 
        
        #self.potential = self.LBpotential
        #self.energyOp = self.LBenergyOp      


 
        self.DOF = 3
        self.recomendedBeta=1000 #=1mK
        self.groundStateEnergy=-0.001302
