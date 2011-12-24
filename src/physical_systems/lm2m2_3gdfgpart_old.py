# -*- coding: utf-8 -*-


class Lm2m2_3part:

    """ Sum of two gaussian potentials.Three particle system but two
    degrees
    of freedowm with Jacobi coordinates
    """

    def __init__(self):

        self.userCode = """ inline float sqr(float x){return x*x;}
        		    inline float quad(float x){return x*x*x*x;}
                            #include "lm2m2.h"                          
         """
        self.DOF = 6

        self.energyOp = """(dlm2m2dr(0.34813790744841544*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))*sqrt(sqr(x1) + sqr(x2) + sqr(x3))*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6)*(40.399681885235644*sqr(x1) + 40.399681885235644*sqr(x2) + 40.399681885235644*sqr(x3) + 121.19904565570693*sqr(x4) + 121.19904565570693*sqr(x5) + 121.19904565570693*sqr(x6) + 139.94860326969626*x1*x4 + 139.94860326969626*x2*x5 + 139.94860326969626*x3*x6) + sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6)*((dlm2m2dr(0.6962758148968309*sqrt(sqr(x1) + sqr(x2) + sqr(x3)))*(80.79936377047129*sqr(x1) + 80.79936377047129*sqr(x2) + 80.79936377047129*sqr(x3)) + (232.09010579361728*lm2m2(0.6962758148968309*sqrt(sqr(x1) + sqr(x2) + sqr(x3))) + 232.09010579361728*lm2m2(0.34813790744841544*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6)) + 232.09010579361728*lm2m2(0.34813790744841544*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6)))*sqrt(sqr(x1) + sqr(x2) + sqr(x3)))*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + dlm2m2dr(0.34813790744841544*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6))*sqrt(sqr(x1) + sqr(x2) + sqr(x3))*(40.399681885235644*sqr(x1) + 40.399681885235644*sqr(x2) + 40.399681885235644*sqr(x3) + 121.19904565570693*sqr(x4) + 121.19904565570693*sqr(x5) + 121.19904565570693*sqr(x6) - 139.94860326969626*x1*x4 - 139.94860326969626*x2*x5 - 139.94860326969626*x3*x6)))/(Sqrt(sqr(x1) + sqr(x2) + sqr(x3))*Sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6)*Sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))"""
        self.potential = """232.09010579361728*lm2m2(0.6962758148968309*sqrt(sqr(x1) + sqr(x2) + sqr(x3))) + 232.09010579361728*lm2m2(0.34813790744841544*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6)) + 232.09010579361728*lm2m2(0.34813790744841544*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))"""
