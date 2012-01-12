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
