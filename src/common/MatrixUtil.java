package common;

public class MatrixUtil {
	// irregular array
	public static int[][] getArray() {
		int [][] num={{1,2,3},{4,5},{2}};
		for(int i = 0; i < num.length; i++) {
			for(int j = 0; j < num[i].length; j++)
				System.out.println(num[i][j]);
		}
		return num;
	}
	public static void printArray(int [][] num) {
		//int [][] num={{1,2,3},{4,5},{2}};
		for(int i = 0; i < num.length; i++) {
			for(int j = 0; j < num[i].length; j++)
				System.out.print(num[i][j] + "\t");
			System.out.println();
		}
	}
	public static void printArray(short [][] num) {
		//int [][] num={{1,2,3},{4,5},{2}};
		for(int i = 0; i < num.length; i++) {
			for(int j = 0; j < num[i].length; j++)
				System.out.print(num[i][j] + "\t");
			System.out.println();
		}
	}
	
	public static void printArray(int[] num) {
		for(int i = 0; i < num.length; i++) {
			System.out.print(num[i] + "\t");
		}
		System.out.println();
	}
	
	public static void printArray(long[] num) {
		for(int i = 0; i < num.length; i++) {
			System.out.print(num[i] + "\t");
		}
		System.out.println();
	}
	
	public static void printArray(float[] num) {
		for(int i = 0; i < num.length; i++) {
			System.out.print(num[i] + "\t");
		}
		System.out.println();
	}
	public static void printArray(double[] num) {
		for(int i = 0; i < num.length; i++) {
			System.out.print(num[i] + "\t");
		}
		System.out.println();
	}
	public static void printArray(boolean[][] bs) {
		for(int i = 0; i < bs.length; i++) {
			for(int j = 0; j < bs[i].length; j++) {
				if(bs[i][j])
					System.out.print("1\t");
				else
					System.out.print("0\t");
			}
			System.out.println();
		}
	}
	public static float sumRow(int[][] matrix, int u) {
		float a = 0.0F;
		for(int m = 0; m < matrix[u].length; m++) {
			a += matrix[u][m];
		}
		return a;
	}
	public static float sumColumn(int[][] matrix, int j, int R) {
		float a = 0.0F;
		for(int i=0; i<R; i++)
		{
			a += matrix[i][j];
		}
		return a;
	}
	
	public static float sumRow(int[][][] matrix, int u, int v) {
		float a = 0.0F;
		for(int m = 0; m < matrix[u][v].length; m++) {
			a += matrix[u][v][m];
		}
		return a;
	}
	public static float sum(int[] nW) {
		long a = 0l;
		for(int i = 0; i < nW.length; i++) {
			a += nW[i];
		}
		return a;
	}
	public static int max(int[] flag) {
		int max = flag[0];
		for(int i = 1; i < flag.length; i++) {
			if(flag[i] > max)
				max = flag[i];
		}
		return max;
	}
	public static int min(int[] flag) {
		int min = flag[0];
		for(int i = 1; i < flag.length; i++) {
			if(flag[i] < min)
				min = flag[i];
		}
		return min;
	}
	
	public static double vectorTimes(double[] a, double[] b)
	{
		double res = 0D;
		for(int i=0; i<a.length; i++)
		{
			res += (a[i]*b[i]);
		}
		return res;
	}
	
	//return the jth column of matrix a
	public static double[] getColumn(double[][] a, int j)
	{
		int n = a.length;
		double[] res = new double[n];
		for(int i=0; i<n; i++)
		{
			res[i] = a[i][j];
		}
		return res;
	}
	
	//return matrix a times column vector b, a:N*N, b:N*1
	public static double[] MatrixTimesVector(double[][] a, double[] b)
	{
		int N = b.length;
		double[] res = new double[N];
		for(int i=0;i<N; i++)
		{
			res[i] = 0;
			for(int j=0; j<N; j++)
			{
				res[i] += (a[i][j]*b[j]);
			}
		}
		return res;
	}
	public static int[] getColumn(int[][] a, int j)
	{
		int n = a.length;
		int[] res = new int[n];
		for(int i=0; i<n; i++)
		{
			res[i] = a[i][j];
		}
		return res;
	}
	
	//dot product of the jth column of matrix a and the kth column of matrix a, for each column, the row index start from 0 to d
	public static double ColumnDotProduct(double[][] a, int j, int k, int d)
	{
		double res = 0;
		for(int i=0; i< d; i++)
		{
			res+= (a[i][j]*a[i][k]);
		}
		return res;
	}
	
	// res = ab'
	public static int[][] cov_product(int[] a, int[] b)
	{
		int N = a.length;
		int[][] res = new int[N][N];
		for ( int i=0; i<N; i++ ) {
			int av = a[i];
			for ( int j=0; j<N; j++ ) {
				res[i][j] = av * b[j];
			}
		}
		return res;
	}
	
	// res = ab'
	public static double[][] cov_product(double[] a, double[] b)
	{
		int N = a.length;
		double[][] res = new double[N][N];
		for ( int i=0; i<N; i++ ) {
			double av = a[i];
			for ( int j=0; j<N; j++ ) {
				res[i][j] = av * b[j];
			}
		}
		return res;
	}
	
	//res = diag(ab')
	public static double[] cov_product_diag(double[] a, double[] b)
	{
		int N = a.length;
		double[] res = new double[N];
		for ( int i=0; i<N; i++ ) {
			res[i] = a[i]*b[i];
		}
		return res;
	}
	
	// res = vec(a)*vec(b)'
	public static double[][] cov_product(double[][] a, double[][] b)
	{
		int N = a.length;
		double[][] res = new double[N*N][N*N];
		for ( int i=0; i<N; i++ ) {
			for(int j=0; j<N; j++)
			{
				double av = a[i][j];
				for(int k=0; k<N; k++)
				{
					for(int l=0; l<N; l++)
					{
						res[i*N+j][k*N+l] = av*b[k][l];
					}
				}
			}
		}
		return res;
	}
	
	// res = diag(vec(a)*vec(b)')
	public static double[] cov_product_diag(double[][] a, double[][] b)
	{
		int N = a.length;
		double[] res = new double[N*N];
		for ( int i=0; i<N; i++ ) {
			for(int j=0; j<N; j++)
			{
				res[i*N+j] = a[i][j]*b[i][j];
			}
		}
		return res;
	}
	
//	/* the inverse of a matrix. */
//	double[][] inverse_cholydec(double[][] A, double[][] lowerTriangle)
//	{
//		int n= A.length;
//		double[][] a = new double[n][n];
//		// upper-triangle matrix
//		for ( int i=0; i<n; i++ ) {
//			for ( int j=0; j<n; j++ ) {
//				if ( j < i ) a[i][j] = 0;
//				else a[i][j] = A[i][j];
//			}
//		}
//
//	    if( spdmatrixcholesky(a, n, true) ) {
//			// get cholesky decomposition result.
//			double[] dPtr = NULL;
//			for ( int i=0; i<n; i++ ) {
//				dPtr = lowerTriangle[i];
//				for ( int j=0; j<=i; j++ ) {
//					dPtr[j] = a[j][i];
//				}
//			}
//
//			// inverse
//	        if( spdmatrixcholeskyinverse(a, n, true) ) {
//	        } else {
//				printf("Inverse matrix error!");
//			}
//		} else {
//			printf("Non-PSD matrix!");
//		}
//
//		for ( int i=0; i<n; i++ ) {
//			for ( int j=0; j<n; j++ ) {
//				if ( j < i ) res[i][j] = a(j, i);
//				else res[i][j] = a(i, j);
//			}
//		}
//	}
}