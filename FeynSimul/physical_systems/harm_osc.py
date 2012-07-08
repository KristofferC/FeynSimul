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

class HarmOsc:
    
    def __init__(self, k = 1.0,lb=False):
        """Initiates a harmonic oscillator system in 1D.

        Contains the potential and virial energy operator for a harmonic oscillator of
        the form V(x) = 1/2 * k * x ^ 2. 

        @type k: number
        @param k: The spring constant of the harmonic ocsillator.

        """
        self.potential = "0.5 * " + str(k) + " * x1 * x1"
        if lb:
            self.potential = "0.5 * (1.0 + sqr(eps)/12.0)*sqr(x1)"
        self.userCode = "float sqr(float x){return x*x;}" 
        self.energyOp = str(k) + " * x1 * x1"
        if lb:
            self.energyOp = "2.0 * 0.5*(1.0 + sqr(eps)/4.0)*sqr(x1)"
        self.DOF = 1
