# This file is part of FeynSimul.
#
# FeynSimul is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FeynSimul is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FeynSimul.  If not, see <http://www.gnu.org/licenses/>.

"""
Gaussian interaction between two particles
"""

class Gauss:

 
    def __init__(self, V1, V2, sigma1, sigma2):
        """ Initiates a sum of two 1D gaussian functions potential.
	
        Three particles interact with gaussian potentials in 1D. The system is given in
        Jacobi coordinates so there are only two degrees of freedom.
	
        @type V1: number
        @param V1: Height of first gaussian function.
        @type V2: number
        @param V2: Height of second gaussian function.
        @type sigma1: number
        @param sigma1: Width of first gaussian function.
        @type sigma2: number
        @param sigma2: Width of second gaussian function.
        """


        self.userCode = """
        inline float operatorFunction(float V, float sigma,
                                           float x, float y)
        {
            const float sqrt3 = 1.73205081f;
            const float piapr = 3.1415926535f;
            const float D = 2.0f * sigma * sigma;
            const float D_inv = 1.0f / D;

            float A = 4.0f * x * x;
            float B = (x - sqrt3 * y) * (x - sqrt3 * y);
            float C = (x + sqrt3 * y) * (x + sqrt3 * y);

            return -V * D_inv / (native_sqrt(piapr) * sigma) *
            ((A - D) * native_exp(-A * D_inv) + (B - D) *
            native_exp(-B * D_inv) + (C - D) * native_exp(-C * D_inv));
        }

        inline float potentialFunction(float V, float sigma,
                                            float x, float y)
        {
            const float sqrt3 = 1.73205081f;
            const float piapr = 3.1415926535f;
            const float D_inv = 0.5f/(sigma*sigma);

            float A = 4.0f*x*x;
            float B = (x-sqrt3*y)*(x-sqrt3*y);
            float C = (x+sqrt3*y)*(x+sqrt3*y);

            return V / (native_sqrt(piapr) * sigma) * (native_exp(-A * D_inv) +
            native_exp(- B* D_inv) + native_exp(-C * D_inv));


        }
        """

        self.potential = ("potentialFunction(" + str(V1) + "f, " + str(sigma1)
                          + "f, x1, x2)" + " + " + "potentialFunction("
                          + str(V2) + "f, " + str(sigma2) + "f, x1, x2)")
        self.energyOp = ("operatorFunction(" + str(V1) + "f, " + str(sigma1)
                         + "f, x1, x2)"  + " + " + "operatorFunction("
                         + str(V2) + "f, " + str(sigma2) + "f, x1, x2)")
        self.DOF = 2

