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

# -*- coding: utf-8 -*-


class Lm2m2_4part:

    def __init__(self):
        #lm2m2_try2 was used below, not so much simplification of energyop...
        #kernel beta of 1000 gives real temp of 50 mK
        #lm2m2partsys.nb was used for generating potential and energy operator
        self.DOF = 9
        #source for following three values:http://iopscience.iop.org/0953-4075/35/3/305/pdf/b20305.pdf
        #(they seem to be far from exact)
        self.groundStateEnergy=-0.4#K
        self.meanSquaredRadius=7.15**2.0#Angstrom^2
        self.meanSquaredAtomDist=10.7**2.0#Angstrom^2
        self.xUnit=4.92341e-11 #Meter
        self.potentialUnit=6.90325e-22 #Joule
        self.meanSquaredRadiusOp="""0.25*(sqr(x1) + sqr(x2) + sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + sqr(x7) + sqr(x8) + sqr(x9))"""
        self.meanRadiusOp="""0.072168784*(3.0f*sqrt(sqr(x7) + sqr(x8) + sqr(x9)) + sqrt(8.0f*sqr(x4) + 8.0f*sqr(x5) + 8.0f*sqr(x6) + sqr(x7) + sqr(x8) + sqr(x9) - 5.6568542*x4*x7 - 5.6568542*x5*x8 - 5.6568542*x6*x9) + sqrt(6.0f*sqr(x1) + 6.0f*sqr(x2) + 6.0f*sqr(x3) + 2.0f*sqr(x4) + 2.0f*sqr(x5) + 2.0f*sqr(x6) + sqr(x7) + sqr(x8) + sqr(x9) + 6.9282032*x1*x4 + 6.9282032*x2*x5 + 6.9282032*x3*x6 + 4.8989795*x1*x7 + 2.8284271*x4*x7 + 4.8989795*x2*x8 + 2.8284271*x5*x8 + 4.8989795*x3*x9 + 2.8284271*x6*x9) + sqrt(6.0f*sqr(x1) + 6.0f*sqr(x2) + 6.0f*sqr(x3) + 2.0f*sqr(x4) + 2.0f*sqr(x5) + 2.0f*sqr(x6) + sqr(x7) + sqr(x8) + sqr(x9) + 2.8284271*x4*x7 - 2.0f*x1*(3.4641016*x4 + 2.4494897*x7) + 2.8284271*x5*x8 - 2.0f*x2*(3.4641016*x5 + 2.4494897*x8) + 2.8284271*x6*x9 - 3.4641016*x3*(2.0f*x6 + 1.4142136*x9)))"""
        self.meanSquaredAtomDistOp="""?"""
        self.userCode = """
        inline float sqr(float x){return x*x;}                 
        float lm2m2(float r)
        {
	        //source: physics.nist.gov
	        float k=8.6173324e-5f; //eV/K boltzmann constant, 
	
	        // source: lm2m2-part of table V in JChemPhys_94_8047.pdf
	        // A, beta, and alpha have a star in JChemPhys_94_8047.pdf
	        // but not in SR. A. Aziz and H. H. Chen, J. Chern. Phys. 67, 		// 5719 (1977) were the expression for the base potential is.
	        // c6, c8 and c10 were capital letters in one of the
	        // sources but not in the other.
	        float A=1.89635353e5f;
	        float alpha=10.70203539f;
	        float c6=1.34687065f;
	        float c8=0.41308398f;
	        float c10=0.17060159f;
	        float D=1.4088f;
	        float beta=-1.90740649f;
	        float epsilon_over_k=10.97f;//K
	        float rm=2.9695f; //Angstrom;
	
	        float Aa=0.0026000000f;
            float xa1 = 1.0035359490f;
	        float xa2 = 1.4547903690f;
	
	        float x=r/rm;
	        float x2=x*x;
	        float x6=x2*x2*x2;
	        float x8=x6*x2;
	        float x10=x8*x2;
	
	        //Vb_star=base potential, reduced form
	        //source: SR. A. Aziz and H. H. Chen, J. Chern. Phys. 67, 		// 5719 (1977)
	        float F_times_rational;
	        if(x<D)
	        {
		        if(x>0.0005f/rm)
			        F_times_rational=exp(-(D/x-1.0f)*(D/x-1.0f))*(c6/x6+c8/x8+c10/x10);
		        else
			        F_times_rational=0.0f;
	        }
	        else
		        F_times_rational=c6/x6+c8/x8+c10/x10;
            float Vb_star=A*exp(-alpha*x+beta*x2)-F_times_rational;
	
            //Va_star=addon potential, reduced form
            //source: JChemPhys_94_8047.pdf
            float Va_star;
            if(x<xa1 || x>xa2)
            	Va_star=0.0f;
            else
            {
                float B=2.0f*M_PI/(xa2-xa1);
            	Va_star=Aa*(sin(B*(x-xa1)-M_PI/2.0f)+1.0f);
            }
            
            return epsilon_over_k*k*(Vb_star+Va_star);
        }

        float dlm2m2dr(float r)
        {
	        //source: physics.nist.gov
	        float k=8.6173324e-5f; //eV/K boltzmann constant, 
	
	        // source: lm2m2-part of table V in JChemPhys_94_8047.pdf
	        // A, beta, and alpha have a star in JChemPhys_94_8047.pdf
	        // but not in SR. A. Aziz and H. H. Chen, J. Chern. Phys. 67,
	        // 5719 (1977) were the expression for the base potential is.
	        // c6, c8 and c10 alsa where capital letters in one of the
	        // sources but not in the other.
	        float A=1.89635353e5f;
	        float alpha=10.70203539f;
	        float c6=1.34687065f;
	        float c8=0.41308398f;
	        float c10=0.17060159f;
	        float D=1.4088f;
	        float beta=-1.90740649f;
	        float epsilon_over_k=10.97f;//K
	        float rm=2.9695f; //Angstrom;
	
	        float Aa=0.0026000000f;
            float xa1 = 1.0035359490f;
	        float xa2 = 1.4547903690f;
	
	        float x=r/rm;
	        float x2=x*x;
	        float x6=x2*x2*x2;
	        float x7=x6*x;
	        float x8=x6*x2;
	        float x9=x8*x;
	        float x10=x8*x2;
	        float x11=x10*x;
	
	        //Vb_star=base potential, reduced form
	        //source: SR. A. Aziz and H. H. Chen, J. Chern. Phys. 67, 		// 5719 (1977)
	        float dF_times_rationaldx;
	        if(x<D)
	        {
		        if(x>0.0005f/rm)
			        dF_times_rationaldx=exp(-(D/x-1.0f)*(D/x-1.0f))*(-6.0f*c6/x7-8.0f*c8/x9-10.0f*c10/x11)+
			        exp(-(D/x-1.0f)*(D/x-1.0f))*(-2.0f*(D/x-1.0f)*(-D/x2))*(c6/x6+c8/x8+c10/x10);
		        else
			        dF_times_rationaldx=0.0f;
	        }
	        else
		        dF_times_rationaldx=-6.0f*c6/x7-8.0f*c8/x9-10.0f*c10/x11;
            float dVb_stardx=A*exp(-alpha*x+beta*x2)*(-alpha+beta*2.0f*x)-dF_times_rationaldx;
	
            //Va_star=addon potential, reduced form
            //source: JChemPhys_94_8047.pdf
            float dVa_stardx;
            if(x<xa1 || x>xa2)
            	dVa_stardx=0.0f;
            else
            {
                float B=2.0f*M_PI/(xa2-xa1);
            	dVa_stardx=Aa*(cos(B*(x-xa1)-M_PI/2.0f)*B+1.0f);
            }
            
            return epsilon_over_k*k*(dVb_stardx+dVa_stardx)/rm;
        }
        """
        self.energyOp = """232.09010579361725*lm2m2(0.6962758148968308*sqrt(sqr(x1) + sqr(x2) + sqr(x3))) + 232.09010579361725*lm2m2(0.49234135028973836*sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9))) + 232.09010579361725*lm2m2(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6)) + 232.09010579361725*lm2m2(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6)) + 232.09010579361725*lm2m2(0.4019950291609113*sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9)) + 232.09010579361725*lm2m2(0.20099751458045564*sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9)) + 0.5*(x4*((-80.79936377047126*dlm2m2dr(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6))*(1.7320508075688772*x1 - 3.0f*x4))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + (80.79936377047126*dlm2m2dr(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))*(1.7320508075688772*x1 + 3.0f*x4))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6) + (46.64953442323208*dlm2m2dr(0.49234135028973836*sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)))*(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7))/sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)) + (93.29906884646415*dlm2m2dr(0.4019950291609113*sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9))*(x4 - 1.4142135623730951*x7))/sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9) + (46.649534423232076*dlm2m2dr(0.20099751458045564*sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9))*(1.7320508*x1 + x4 + 2.8284271247461903*x7))/sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9)) + x7*((131.94480851945073*dlm2m2dr(0.49234135028973836*sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)))*(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7))/sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)) - (93.29906884646415*dlm2m2dr(0.4019950291609113*sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9))*(1.4142135623730951*x4 - 2.0f*x7))/sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9) + (93.29906884646415*dlm2m2dr(0.20099751458045564*sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9))*(2.4494897427831783*x1 + 1.4142135623730951*x4 + 4.0f*x7))/sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9)) + x1*((161.59872754094252*dlm2m2dr(0.6962758148968308*sqrt(sqr(x1) + sqr(x2) + sqr(x3)))*x1)/sqrt(sqr(x1) + sqr(x2) + sqr(x3)) + (80.79936377047126*dlm2m2dr(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6))*(x1 - 1.7320508075688772*x4))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + (80.79936377047126*dlm2m2dr(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))*(x1 + 1.7320508*x4))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6) - (80.79936377047125*dlm2m2dr(0.49234135028973836*sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)))*(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7))/sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)) + (46.649534423232076*dlm2m2dr(0.20099751458045564*sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9))*(3.0f*x1 + 1.7320508075688772*x4 + 4.898979485566357*x7))/sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9)) + x5*((-80.79936377047126*dlm2m2dr(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6))*(1.7320508075688772*x2 - 3.0f*x5))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + (80.79936377047126*dlm2m2dr(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))*(1.7320508075688772*x2 + 3.0f*x5))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6) + (46.64953442323208*dlm2m2dr(0.49234135028973836*sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)))*(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8))/sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)) + (93.29906884646415*dlm2m2dr(0.4019950291609113*sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9))*(x5 - 1.4142135623730951*x8))/sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9) + (46.649534423232076*dlm2m2dr(0.20099751458045564*sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9))*(1.7320508*x2 + x5 + 2.8284271247461903*x8))/sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9)) + x8*((131.94480851945073*dlm2m2dr(0.49234135028973836*sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)))*(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8))/sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)) - (93.29906884646415*dlm2m2dr(0.4019950291609113*sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9))*(1.4142135623730951*x5 - 2.0f*x8))/sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9) + (93.29906884646415*dlm2m2dr(0.20099751458045564*sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9))*(2.4494897427831783*x2 + 1.4142135623730951*x5 + 4.0f*x8))/sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9)) + x2*((161.59872754094252*dlm2m2dr(0.6962758148968308*sqrt(sqr(x1) + sqr(x2) + sqr(x3)))*x2)/sqrt(sqr(x1) + sqr(x2) + sqr(x3)) + (80.79936377047126*dlm2m2dr(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6))*(x2 - 1.7320508075688772*x5))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + (80.79936377047126*dlm2m2dr(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))*(x2 + 1.7320508*x5))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6) - (80.79936377047125*dlm2m2dr(0.49234135028973836*sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)))*(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8))/sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)) + (46.649534423232076*dlm2m2dr(0.20099751458045564*sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9))*(3.0f*x2 + 1.7320508075688772*x5 + 4.898979485566357*x8))/sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9)) + x6*((-80.79936377047126*dlm2m2dr(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6))*(1.7320508075688772*x3 - 3.0f*x6))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + (80.79936377047126*dlm2m2dr(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))*(1.7320508075688772*x3 + 3.0f*x6))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6) + (46.64953442323208*dlm2m2dr(0.49234135028973836*sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)))*(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9))/sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)) + (93.29906884646415*dlm2m2dr(0.4019950291609113*sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9))*(x6 - 1.4142135623730951*x9))/sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9) + (46.649534423232076*dlm2m2dr(0.20099751458045564*sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9))*(1.7320508*x3 + x6 + 2.8284271247461903*x9))/sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9)) + x9*((131.94480851945073*dlm2m2dr(0.49234135028973836*sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)))*(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9))/sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)) - (93.29906884646415*dlm2m2dr(0.4019950291609113*sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9))*(1.4142135623730951*x6 - 2.0f*x9))/sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9) + (93.29906884646415*dlm2m2dr(0.20099751458045564*sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9))*(2.4494897427831783*x3 + 1.4142135623730951*x6 + 4.0f*x9))/sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9)) + x3*((161.59872754094252*dlm2m2dr(0.6962758148968308*sqrt(sqr(x1) + sqr(x2) + sqr(x3)))*x3)/sqrt(sqr(x1) + sqr(x2) + sqr(x3)) + (80.79936377047126*dlm2m2dr(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6))*(x3 - 1.7320508075688772*x6))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + (80.79936377047126*dlm2m2dr(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))*(x3 + 1.7320508*x6))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6) - (80.79936377047125*dlm2m2dr(0.49234135028973836*sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)))*(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9))/sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9)) + (46.649534423232076*dlm2m2dr(0.20099751458045564*sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9))*(3.0f*x3 + 1.7320508075688772*x6 + 4.898979485566357*x9))/sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9)))"""
        self.potential = """232.09010579361725*lm2m2(0.6962758148968308*sqrt(sqr(x1) + sqr(x2) + sqr(x3))) + 232.09010579361725*lm2m2(0.49234135028973836*sqrt(sqr(-0.70710678*x1 + 0.40824829*x4 + 1.1547005*x7) + sqr(-0.70710678*x2 + 0.40824829*x5 + 1.1547005*x8) + sqr(-0.70710678*x3 + 0.40824829*x6 + 1.1547005*x9))) + 232.09010579361725*lm2m2(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6)) + 232.09010579361725*lm2m2(0.3481379074484154*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6)) + 232.09010579361725*lm2m2(0.4019950291609113*sqrt(sqr(x4) + sqr(x5) + sqr(x6) + 2.0f*sqr(x7) + 2.0f*sqr(x8) + 2.0f*sqr(x9) - 2.8284271*x4*x7 - 2.8284271*x5*x8 - 2.8284271*x6*x9)) + 232.09010579361725*lm2m2(0.20099751458045564*sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 8.0f*sqr(x7) + 8.0f*sqr(x8) + 8.0f*sqr(x9) + 3.4641016*x3*x6 + 5.6568542*x4*x7 + 3.4641016*x1*(x4 + 2.8284271*x7) + 5.6568542*x5*x8 + 3.4641016*x2*(x5 + 2.8284271*x8) + 9.797959*x3*x9 + 5.6568542*x6*x9))"""
