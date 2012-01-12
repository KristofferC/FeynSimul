# -*- coding: utf-8 -*-


class Lm2m2_3part:
    def __init__(self):
        #kernel beta of 1000 gives real temp of 15 mK
        self.DOF = 6
        self.groundStateEnergy=-0.126#K
        self.xUnit=8.98888e-11 #Meter
        self.potentialUnit=2.07098e-22 #Joule
        self.meanSquaredRadiusOp="""0.33333333333333*(sqr(x1) + sqr(x2) + sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6))"""
        self.meanRadiusOp="""0.13608276*(2.0f*sqrt(sqr(x4) + sqr(x5) + sqr(x6)) + sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + sqrt(3.0f*sqr(x1) + 3.0f*sqr(x2) + 3.0f*sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))"""
        self.meanSquaredAtomDistOp="""sqr(x1) + sqr(x2) + sqr(x3) + sqr(x4) + sqr(x5) + sqr(x6)"""
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
        self.energyOp = """773.633685978724*lm2m2(1.2712199002142859*sqrt(sqr(x1) + sqr(x2) + sqr(x3))) + 773.633685978724*lm2m2(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6)) + 773.633685978724*lm2m2(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6)) + 0.5*(x1*((983.4585370922837*dlm2m2dr(1.2712199002142859*sqrt(sqr(x1) + sqr(x2) + sqr(x3)))*x1)/sqrt(sqr(x1) + sqr(x2) + sqr(x3)) + (491.72926854614184*dlm2m2dr(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6))*(x1 - 1.7320508075688772*x4))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + (491.72926854614184*dlm2m2dr(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))*(x1 + 1.7320508*x4))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6)) + x4*((dlm2m2dr(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6))*(-851.7000766905983*x1 + 1475.1878056384255*x4))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + (dlm2m2dr(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))*(851.7000766905983*x1 + 1475.1878056384255*x4))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6)) + x2*((983.4585370922837*dlm2m2dr(1.2712199002142859*sqrt(sqr(x1) + sqr(x2) + sqr(x3)))*x2)/sqrt(sqr(x1) + sqr(x2) + sqr(x3)) + (491.72926854614184*dlm2m2dr(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6))*(x2 - 1.7320508075688772*x5))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + (491.72926854614184*dlm2m2dr(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))*(x2 + 1.7320508*x5))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6)) + x5*((dlm2m2dr(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6))*(-851.7000766905983*x2 + 1475.1878056384255*x5))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + (dlm2m2dr(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))*(851.7000766905983*x2 + 1475.1878056384255*x5))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6)) + x3*((983.4585370922837*dlm2m2dr(1.2712199002142859*sqrt(sqr(x1) + sqr(x2) + sqr(x3)))*x3)/sqrt(sqr(x1) + sqr(x2) + sqr(x3)) + (491.72926854614184*dlm2m2dr(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6))*(x3 - 1.7320508075688772*x6))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + (491.72926854614184*dlm2m2dr(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))*(x3 + 1.7320508*x6))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6)) + x6*((dlm2m2dr(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6))*(-851.7000766905983*x3 + 1475.1878056384255*x6))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6) + (dlm2m2dr(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))*(851.7000766905983*x3 + 1475.1878056384255*x6))/sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*sqr(x4) + 3.0f*sqr(x5) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6)))"""
        self.potential = """773.633685978724*lm2m2(1.2712199002142859*sqrt(sqr(x1) + sqr(x2) + sqr(x3))) + 773.633685978724*lm2m2(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) - 3.4641016*x1*x4 - 3.4641016*x2*x5 - 3.4641016*x3*x6)) + 773.633685978724*lm2m2(0.6356099501071429*sqrt(sqr(x1) + sqr(x2) + sqr(x3) + 3.0f*(sqr(x4) + sqr(x5)) + 3.0f*sqr(x6) + 3.4641016*x1*x4 + 3.4641016*x2*x5 + 3.4641016*x3*x6))"""
