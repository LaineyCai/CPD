package distribution;

public class MVGaussian {
	//for Gaussian random variable
	public int m_iSet;
	public double m_dGset;
	
	public MVGaussian()
	{
		m_iSet = 0;
	}
	
	public double[] nextMVGaussian(double[] mean, double[][] precision) 
	{
		double[] res = null;
		double[][] precisionLowerTriangular = Cholesky.cholesky(precision);

		res = nextMVGaussianWithCholesky(mean, precisionLowerTriangular);
		return res;
	}

	public double[] nextMVGaussianWithCholesky(double[] mean, double[][] precisionLowerTriangular) 
	{		
		// Initialize vector z to standard normals
		//  [NB: using the same array for z and x]
		int n = mean.length;
		double[] res = new double[n];
		for (int i = 0; i < n; i++) {
			res[i] = nextGaussian();
		}

		// Now solve trans(L) x = z using back substitution
		double innerProduct;

		for (int i = n-1; i >= 0; i--) {
			innerProduct = 0;
			for (int j = i+1; j < n; j++) {
				// the cholesky decomp got us the precisionLowerTriangular triangular
				//  matrix, but we really want the transpose.
				innerProduct += res[j] * precisionLowerTriangular[j][i];
			}
			res[i] = (res[i] - innerProduct) / precisionLowerTriangular[i][i];
		}

		for (int i = 0; i < n; i++) {
			res[i] += mean[i];
		}
		return res;
	}

	public double nextGaussian()
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
			return v2 * dFac;
		} else {
			m_iSet = 0;
			return m_dGset;
		}	
	}

}
