package distribution;

public class PolyaGamma {
	public double TRUNC;
	public double cutoff;
	
	//for Gaussian random variable
	public double m_dGset;
	public int m_iSet;
	
	public PolyaGamma()
	{
		TRUNC = 0.64;
		cutoff = 1/TRUNC;
		m_iSet = 0;
	}

	// Sample from PG(n, Z) using Devroye-like method.
	// n is a natural number and z is a positive real.
	public double nextPG(int n, double z)
	{
		double dRes = 0;
		for(int i=0; i<n; i++)
		{
			dRes += nextPG1(z);
		}
		return dRes;
	}
	
	// sample from PG(1, z)
	double nextPG1(double zVal)
	{
		double z = Math.abs(zVal) * 0.5;
		double fz = (Math.PI * Math.PI * 0.125 + z * z * 0.5);

		double X = 0;
		int numTrials = 0;
		while ( true ) 
		{
			numTrials ++;
			double dU = Math.random();
			if ( dU < texpon(z) ) {
				X = TRUNC + rexp1() / fz;
			} else {
				X = rtigauss(z);
			}
//			if(Double.isNaN(X))
//				System.out.println("x wrong nan!!!");

			double S = a(0, X);
//			if(Double.isNaN(S))
//				System.out.println("S wrong nan!!!");
			double Y = Math.random() * S;
			int n = 0;
			while ( true ) {
				n ++;
				if ( n % 2 == 1 ) {
	                S = S - a(n, X);
	                if ( Y <= S ) break;
				} else {
	               S = S + a(n,X);
	               if ( Y>S ) break;
	 			}
//				if(n>=1000)
//					break;
				//if ( n % 1000 == 0 ) System.out.println("n is "+n+" (z: "+z+", X: "+X+", s: "+S+", Y: "+Y+")\n");
			};

	       if ( Y <= S ) break;
		  // if ( numTrials % 1000 == 0 )
				//printf("# trials: %d\n", numTrials);
	 	};

		return 0.25 * X;
	}
	
	// rtigauss - sample from truncated Inv-Gauss(1/abs(Z), 1.0) 1_{(0, TRUNC)}.
	double rtigauss(double Z)
	{
		double R = TRUNC;

		Z = Math.abs(Z);
		double mu = 1 / Z;
		double X = R + 1;
		if ( mu > R ) {
			double alpha = 0;
			while ( Math.random() > alpha ) {
				double E1 = rexp1();
				double E2 = rexp1();
				while ( Math.pow(E1,2.0) > 2*E2 / R) {
					E1 = rexp1();
					E2 = rexp1();
				}
				X = R / Math.pow((1 + R*E1), 2.0);
				alpha = Math.exp(-0.5 * Z * Z * X);
			}
		} else {
			while ( X > R || X <= 0 ) {
				//double lambda = 1;
				double Y = Math.pow(rnorm(), 2.0);
				double muY = mu * Y;
				X = mu * (1 + 0.5*muY /*/ lambda*/ - 0.5 /*/ lambda*/ * Math.sqrt(4 * muY * /*lambda **/ (1 + muY)));
				if ( Math.random() > mu / (mu + X) ) {
					X = Math.pow(mu, 2.0) / X;
				}
			}
		}

	    return X;
	}
	
	double texpon(double Z)
	{
	   double x = TRUNC;
	   double fz = (Math.PI*Math.PI*0.125 + Z*Z*0.5);
	   double b = Math.sqrt(1.0 / x) * (x * Z - 1);
	   double a = -1.0 * Math.sqrt(1.0 / x) * (x * Z + 1);
	
	   double x0 = Math.log(fz) + fz * TRUNC;
	   double xb = x0 - Z + pnorm(b, true);
	   double xa = x0 + Z + pnorm(a, true);
	
	   double qdivp = 4 / Math.PI * ( Math.exp(xb) + Math.exp(xa) );
	
	   return (1.0 / (1.0 + qdivp));
	}

	//the cumulative density function for standard normal
	double pnorm(double x, boolean bUseLog)
	{
		final double c0 = 0.2316419;
		final double c1 = 1.330274429;
		final double c2 = 1.821255978;
		final double c3 = 1.781477937;
		final double c4 = 0.356563782;
		final double c5 = 0.319381530;
		final double c6 = 0.398942280401;
		final double negative = (x < 0 ? 1.0 : 0.0);
		final double xPos = (x < 0.0 ? -x : x);
		final double k = 1.0 / ( 1.0 + (c0 * xPos));
		final double y1 = (((((((c1*k-c2)*k)+c3)*k)-c4)*k)+c5)*k;
		final double y2 = 1.0 - (c6*Math.exp(-0.5*xPos*xPos)*y1);
	
		if ( bUseLog ) {
			return Math.log(((1.0-negative)*y2) + (negative*(1.0-y2)));
		} else {
			return ((1.0-negative)*y2) + (negative*(1.0-y2));
		}
	}

	//draw a sample from standard norm distribution.
	double rnorm()
	{
		if ( m_iSet == 0 ) {
			double dRsq = 0;
			double v1, v2;
			do {
				v1 = 2.0 * Math.random() - 1.0;
				v2 = 2.0 * Math.random() - 1.0;
				dRsq = v1 * v1 + v2 * v2;
			} while (dRsq > 1.0 || dRsq < 1e-300);
	
			double dFac = Math.sqrt(-2.0 * Math.log(dRsq) / dRsq);
			m_dGset = v1 * dFac;
			m_iSet = 1;
			return /*dMu + dSigma * */ v2 * dFac;
		} else {
			m_iSet = 0;
			return /*dMu + dSigma **/ m_dGset;
		}
	}

	//draw a sample from an exponential distribution with parameter lambda
	double rexp1(/*double lambda*/)
	{
		double dval = 0;
		while ( dval>=1 || dval <= 0 ) {
			dval = Math.random();
		}
		return (-Math.log( dval )/*/ lambda*/);
	}

	
	//Calculate coefficient n in density of PG(1.0, 0.0), i.e. J* from Devroye.
	double a(int n, double x)
	{
		double dRes = 0;
	
		if ( x>TRUNC )
			dRes = Math.PI * (n+0.5) * Math.exp( - Math.pow((n+0.5)*Math.PI, 2.0) * x * 0.5 );
		else
			dRes = Math.pow((2/Math.PI/x), 1.5) * Math.PI * (n+0.5) * Math.exp( -2* Math.pow((n+0.5), 2.0) / x );
		
//		if(Double.isNaN(dRes))
//			System.out.println("x wrong nan!!!");
	
		return dRes;
	}   

//sample from normal distribution

}
