package distribution;

public class InverseGaussian {
	public double m_dMu;
	public double m_dScale;
	
	//for Gaussian random variable
	public double m_dGset;
	public int m_iSet;
	
	public void reset(double dMu, double dScale)
	{
		m_iSet = 0;
		m_dMu = dMu;
		m_dScale = dScale;
	}
	
	double sample()
	{
		  double v = nextGaussian();   // sample from a normal distribution with a mean of 0 and 1 standard deviation

		  double y = v*v;
	      double x = m_dMu + (m_dMu * m_dMu * y) / (2* m_dScale) - (m_dMu/(2*m_dScale)) * Math.sqrt(4*m_dMu*m_dScale*y + m_dMu*m_dMu*y*y);

		  double test = Math.random();  // sample from a uniform distribution between 0 and 1

		  if (test <= (m_dMu)/(m_dMu + x))
	             return x;
	      else
	             return (m_dMu*m_dMu) / x;
	}

	double nextGaussian()
	{
		//double dMu = 0;
		//double dSigma = 1;

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

}
