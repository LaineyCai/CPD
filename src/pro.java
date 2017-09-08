import common.FileUtil;

public class pro {
	
	public static void main(String args[])
	{
		for(int train=0; train<=9; train++)
		{
			System.out.println(train);
//			String pathin = "data/dblp/input/train"+train+"/";
			String pathin = "data/toy/input/";
			String ftn = pathin+"tpo.txt";
//			String ftp = pathin+"tp.txt";
			
			
			float[][] tn = FileUtil.readLinesFMatrix(ftn, "\t");
			int M = tn.length, N = tn[0].length;
			float[][] ntn = new float[M][N];
			
			for(int i=0; i<M; i++)
			{
				float s = 0f;
				for(int j=0; j<N; j++)
				{
					s+=tn[i][j];
				}
				for(int j=0; j<N; j++)
				{
					ntn[i][j] = Math.round(tn[i][j]/s*10000f)/10000f;
				}
			}
			
			FileUtil.writeFMatrix(pathin+"tp.txt", ntn);
		}

	}

}
