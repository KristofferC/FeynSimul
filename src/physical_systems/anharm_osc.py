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
Anharmonic oscillator
"""

class AnHarmOsc:

    def __init__(self, lambd = 1.0):
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

