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
