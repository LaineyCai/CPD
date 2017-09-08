package cpdtwitter;
/*
 * Community profiling - Hongyun Cai & 14/09/2015
 */

import gnu.trove.map.hash.TIntObjectHashMap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import common.ModelUtil;

import distribution.PolyaGamma;
import common.ComUtil;
import common.FileUtil;
import common.MatrixUtil;

import de.bwaldvogel.liblinear.Feature;
import de.bwaldvogel.liblinear.FeatureNode;
import de.bwaldvogel.liblinear.Linear;
import de.bwaldvogel.liblinear.Parameter;
import de.bwaldvogel.liblinear.Problem;
import de.bwaldvogel.liblinear.SolverType;

public class Model {
	// user related parameters
	// public Document Doc;
	public int T; // no of topics
	public int U; // no of users
	public int V; // vocabulary size
	public int M; // no of time stamps
	public int A; //no of communities
	public int Nf = 2;	//no of features for each user
	public int N2f = 2*Nf;
	public int NLu = 3589811;
	public int NLt = 992522;

	// hyperparameters
	public float alpha;
	public float beta;	
	public float rho;


	// model parameters   
	public int niters; // number of Gibbs sampling iteration

	// Estimated/Inferenced parameters
	public float[][] theta; // community-topic distribution, C*T
	public float[][] vPhi; // topic-word distribution, T*V
	public float[][] pai; //user-community distribution, U x C
	public double[][][] eta; //community-community retweet link distribution on topic k,  T x C x C 
	public double[] nu; // individual factor  
	//augmented variable
	public float[] lambda;	// U x U (any two users), NLu
	public float[] delta;	//N x N (any two tweets), NLt

	// Temp variables while sampling
	public int[][] Z; // U x N
	public int[][] C; //U x N
	public int[][] NTW; // T x V, sum is: SNTW[T]
	public int[][] NUC; // sum U x C, sum is: SNUC[U]; c_hat before normalisation
	public int[][] NCT; //C x T, sum is SNCT[C]
	public int[][] popularity;	//T x M, topic popularity at timeslice t
	
	public float[] SNTW;
	public float[] SNUC;
	public float[] SNCT;
	public float[] SNTC;	//for normalize hat_z_k
	public float[] SPTM;	//for normalize popularity
	
	//intermediate variables for training
	public PolyaGamma[] m_uPGsampler;
	public PolyaGamma[] m_tPGsampler;
	
	//variable record the neighbour of each user and tweet   
	TIntObjectHashMap<int[]> userNeighbour; // user, the list of users who follow her, store the index of user, not the real user id, the index is the row number of each user in the file userlist.txt
	TIntObjectHashMap<int[][]> userInNeighbours; //users and the list of users who she follows, also the index of userNeighbour to ref lambda
	
	HashMap<String, List<int[]>> tweetInNeighbours; //the tweet and the list of tweets which retweet it, for each int[], we store the user and the index of the tweet in the user's tweets, also for the second int[], we stire the ubdex to refdelta
	int[] ltRefIx;
	int[] luRefIx;
	int[][] ugroups;
//	double[][] indfea;
	float[][] tpf;
	float[][] tnf;
	int[] nvs;
	float[] cijp;	//positive cij in logistic regression
	float[] cijn;	//negative cij in logistic regression
	float[] popp;	//positive popularity in logistic regression
	float[] popn;	//positive popularity in logistic regression
	HashSet<String> retweet;
 	
	//individual statistics
//	double[][] ind_fea;	//U * number of features; feature need normalize to keep in the same scale
	
	public Model(int no, int comno, float betaV, int[][] myugroups) {
		T = no;
		A = comno;
		beta = betaV;
		alpha = (float) (50.0 / T);
		rho = (float) (50.0 / A);
		ugroups = myugroups;
	}
	
	/*protected boolean init(User[] users, int vocabSize, int tsNo,TIntObjectHashMap<int[]> myuserNeighbour, int[] myltRefIx,int[] myluRefIx, HashSet<String> mynoretweet, double[][] myindfea) {
		U = users.length;
		V = vocabSize;
		M = tsNo;
		// assign topics and communities randomly
		Z = new int[U][];
		C = new int[U][];		
		for (int i = 0; i < U; i++) {
			Z[i] = new int[users[i].getDocWords().length];
			C[i] = new int[Z[i].length];
			for (int j = 0; j < Z[i].length; j++) {
 				Z[i][j] = (int) Math.floor(Math.random() * T);
				if (Z[i][j] < 0)
					Z[i][j] = 0;
				if (Z[i][j] > T - 1)
					Z[i][j] = (int) (T - 1);
				C[i][j] = (int) Math.floor(Math.random() * A);
				if (C[i][j] < 0)
					C[i][j] = 0;
				if (C[i][j] > A - 1)
					C[i][j] = (int) (A - 1);				
			}
		}
		m_uPGsampler = new PolyaGamma[NLu];
		for(int i=0; i<NLu; i++)
		{
			m_uPGsampler[i] = new PolyaGamma();
		}
		m_tPGsampler = new PolyaGamma[NLt];
		for(int i=0; i<NLt; i++)
		{
			m_tPGsampler[i] = new PolyaGamma();
		}
		userNeighbour = myuserNeighbour;
		Graphs g = new Graphs();
		userInNeighbours = ModelUtil.generateuserInNeighbourT(userNeighbour,U);
		tweetInNeighbours = g.generatetweetInNeighbour(users);
		ltRefIx = myltRefIx;
		luRefIx = myluRefIx;
		retweet = mynoretweet;
		cleanTempPrmts();
		computeTempPrmts(users, Z,C);
		computeSum(users, T);
		indfea = myindfea;

		return true;
	}*/
	

	/**
	 * initialize the model
	 */
	protected boolean init(User[] users, int vocabSize, int tsNo,TIntObjectHashMap<int[]> myuserNeighbour, int[] myltRefIx,int[] myluRefIx, HashSet<String> mynoretweet, float[][] mytpf, float[][] mytnf, int[] mynvs) {
		U = users.length;
		V = vocabSize;
		M = tsNo;
		// assign topics and communities randomly
		Z = new int[U][];
		C = new int[U][];		
		for (int i = 0; i < U; i++) {
			Z[i] = new int[users[i].getDocWords().length];
			C[i] = new int[Z[i].length];
			for (int j = 0; j < Z[i].length; j++) {
 				Z[i][j] = (int) Math.floor(Math.random() * T);
				if (Z[i][j] < 0)
					Z[i][j] = 0;
				if (Z[i][j] > T - 1)
					Z[i][j] = (int) (T - 1);
				C[i][j] = (int) Math.floor(Math.random() * A);
				if (C[i][j] < 0)
					C[i][j] = 0;
				if (C[i][j] > A - 1)
					C[i][j] = (int) (A - 1);				
			}
		}
		m_uPGsampler = new PolyaGamma[NLu];
		for(int i=0; i<NLu; i++)
		{
			m_uPGsampler[i] = new PolyaGamma();
		}
		m_tPGsampler = new PolyaGamma[NLt];
		for(int i=0; i<NLt; i++)
		{
			m_tPGsampler[i] = new PolyaGamma();
		}
		userNeighbour = myuserNeighbour;
		Graphs g = new Graphs();
		userInNeighbours = ModelUtil.generateuserInNeighbourT(userNeighbour,U);
		tweetInNeighbours = g.generatetweetInNeighbour(users);
		ltRefIx = myltRefIx;
		luRefIx = myluRefIx;
		retweet = mynoretweet;
		cleanTempPrmts();
		computeTempPrmts(users, Z,C);
		computeSum(users, T);
//		indfea = myindfea;
		tpf = mytpf;
		tnf = mytnf;
		NLt = tpf.length;
		nvs = mynvs;
		cijp = new float[NLt];
		cijn = new float[NLt];
		popp = new float[NLt];
		popn = new float[NLt];
		return true;
	}
	
	/**
	 * initialize the model from files
	 */
	protected boolean init_frIter(User[] users, int vocabSize, int tsNo,TIntObjectHashMap<int[]> myuserNeighbour, String iterfile, int[] myltRefIx, int[] myluRefIx,  HashSet<String> mynoretweet) {
		U = users.length;
		V = vocabSize;
		M = tsNo;
		// assign topics and communities randomly
		Z = new int[U][];
		C = new int[U][];
		for (int i = 0; i < U; i++) {
			Z[i] = new int[users[i].getDocWords().length];
			C[i] = new int[Z[i].length];
		}
		
		FileUtil.readMatrix(iterfile+"Z", Z);
		FileUtil.readMatrix(iterfile+"C", C);
		m_uPGsampler = new PolyaGamma[NLu];
		for(int i=0; i<NLu; i++)
		{
			m_uPGsampler[i] = new PolyaGamma();
		}
		m_tPGsampler = new PolyaGamma[NLt];
		for(int i=0; i<NLt; i++)
		{
			m_tPGsampler[i] = new PolyaGamma();
		}
		userNeighbour = myuserNeighbour;
		Graphs g = new Graphs();
		userInNeighbours = ModelUtil.generateuserInNeighbourT(userNeighbour,U);
		tweetInNeighbours = g.generatetweetInNeighbour(users);
		ltRefIx = myltRefIx;
		luRefIx = myluRefIx;
		retweet = mynoretweet;
		cleanTempPrmts();
		//computeTempPrmts(users, Z,C);
		FileUtil.readMatrix(iterfile+"nuc", NUC);
		FileUtil.readMatrix(iterfile+"nct", NCT);
		FileUtil.readMatrix(iterfile+"ntw", NTW);
		FileUtil.readMatrix(iterfile+"popu", popularity);
		computeSum(users, T);
		
		FileUtil.readMatrix(iterfile+"theta", theta);
		FileUtil.readMatrix(iterfile+"vPhi", vPhi);
		FileUtil.readMatrix(iterfile+"nu", nu);
		FileUtil.readMatrix(iterfile+"pai", pai);
		FileUtil.readMatrix(iterfile+"eta", eta);
		FileUtil.readMatrix(iterfile+"lambda", lambda);
		FileUtil.readMatrix(iterfile+"delta", delta);
		
		return true;
	}

	//clean value for variational parameters
	public void cleanTempPrmts() {
		// initial parameters 
		NTW = new int[T][];
		vPhi = new float[T][];
		popularity = new int[T][];
		eta = new double[T][][];
		for (int t = 0; t < T; t++) {
			NTW[t] = new int[V];
			vPhi[t] = new float[V];
			for (int i = 0; i < V; i++) {
				NTW[t][i] = 0;
				vPhi[t][i] = 0.0f;
			}
			popularity[t] = new int[M];
			for (int i = 0; i < M; i++) {
				popularity[t][i] = 0;
			}
			eta[t] = new double[A][];
			for(int i=0; i<A;i++)
			{
				eta[t][i] = new double[A];
				for(int j=0; j<A; j++)
					eta[t][i][j] = 1.0f;
			}			
		}

		NCT = new int[A][];
		theta = new float[A][];
		for (int i = 0; i < A; i++) {
			NCT[i] = new int[T];
			theta[i] = new float[T];			
			for (int t = 0; t < T; t++) {
				NCT[i][t] = 0;
				theta[i][t] = 0.0f;
			}
		}
		
		NUC = new int[U][];
		pai = new float[U][];
		for (int i = 0; i < U; i++) {
			NUC[i] = new int[A];
			pai[i] = new float[A];
			for (int t = 0; t < A; t++) {
				NUC[i][t] = 0;
				pai[i][t] = 0.0f;
			}
		}
		nu = new double[N2f+2];
		for(int i=0; i<N2f+2; i++)
		{
			nu[i] = 1.0d;
		}
		lambda = new float[NLu];
		for(int i=0; i<NLu;i++)
		{
			lambda[i] = 1.0f;
		}
		delta = new float[NLt];
		for(int i=0; i<NLt;i++)
		{
			delta[i] = 1.0f;
		}
	}
	
	//clean value for variational parameters
	public void cleanCounter() {
		// initial parameters NW[] NWT[][] NIT[][] NY[] NUT[][]
		NTW = new int[T][];
		popularity = new int[T][];
		for (int t = 0; t < T; t++) {
			NTW[t] = new int[V];
			for (int i = 0; i < V; i++) {
				NTW[t][i] = 0;
			}
			popularity[t] = new int[M];
			for (int i = 0; i < M; i++) {
				popularity[t][i] = 0;
			}
		}

		NCT = new int[A][];
		for (int i = 0; i < A; i++) {
			NCT[i] = new int[T];
			for (int t = 0; t < T; t++) {
				NCT[i][t] = 0;
			}
		}
		
		NUC = new int[U][];
		for (int i = 0; i < U; i++) {
			NUC[i] = new int[A];
			for (int t = 0; t < A; t++) {
				NUC[i][t] = 0;
			}
		}
	}
	
//	public void printnoretweet(ArrayList<User> users, String outputDir)
//	{
//		ArrayList<String> retweet = new ArrayList<String>();
//		for (int u = 0; u < U; u++) {
//			for (int n = 0; n < users.get(u).getDocWords().length; n++) {
//				if(!getNoretweets(u, n, users))
//					retweet.add(u+"\t"+n);
//			}			
//		}
//		FileUtil.writeLines(outputDir+"retweet.txt", retweet);
//	}
	
	public void SampleCT(User[] users, double arho, double talpha, double vbeta, int[] UIdx)
	{
//		for (int u = 0; u < U; u++) {
		for(int u: UIdx){
			double[] gVal = ComputeNeighborLhoodG(users, u, luRefIx[u]);
			int tRefIx = ltRefIx[u];
			for (int n = 0; n < users[u].getDocWords().length; n++) {			
				SampleCommunity( u, n, users, tRefIx, arho, talpha, gVal);
				SampleTopic(users[u].getDocWords()[n], users[u]
						.getDocTimeStamp()[n], u, n, users, tRefIx,talpha, vbeta);
//				if(users.get(u).getRetweets()!=null)
//				{
					if(users[u].getRetweets().containsKey(n))
					{
						tRefIx++;
					}
//				}
			}
		}
	}
	
	public void inference(int iteration, User[] users,
			int saveStep, int saveTimes, String outputDir, int initer) {
		if (iteration < saveStep * saveTimes) {
			System.err.println("iteration should be at least: " + saveStep * saveTimes);
			System.exit(0);
		}
		niters = iteration;
		double arho = A * rho, talpha = T * alpha, vbeta = V * beta;
		//multithreading
		int Nthread = ugroups.length;
		ThreadUser[] tt = new ThreadUser[Nthread];
		for(int n=0; n<Nthread; n++)
		{
			tt[n] = new ThreadUser(this, users, arho, talpha, vbeta, ugroups[n]); 
		}
		ExecutorService executor = Executors.newFixedThreadPool(Nthread);
		
		for (int i = initer+0; i < iteration; i++) {
			System.out.println("iteration " + i);
			long begintime = System.nanoTime();			
			drawEta(users);
			drawLambda(users);
			drawDelta(users);			
			
			for(int n=0; n<Nthread; n++)
			{
				tt[n].UpdateData(this);
				executor.submit(tt[n]);
			}
								
			boolean finished = false;
			while (! finished) {
				
				try {
					Thread.sleep(10);
				} catch (InterruptedException e) {
					
				}
				
				finished = true;
				
				// Are all the threads done?
				for (int nt = 0; nt < Nthread; nt++) {				
					finished = finished && tt[nt].isFinish();
				}
			}
			
			//update the shared variables 

			cleanTempPrmts();
			computeTempPrmts(users, Z,C);
			computeSum(users, T);
			
			//update cij and pop
			assigncijpop(users);
			drawNu();
			
			if (i >= (iteration - 1 - (saveStep * (saveTimes - 1))))
			{
				if ((iteration - i - 1) % saveStep == 0) {
					System.out.println("Saving the model at " + (i + 1) + "-th iteration");
//					uRefIx = 0; 
//					tRefIx =0;
					drawEta(users);
					drawLambda(users);
					drawDelta(users);
					for (int u = 0; u < U; u++) {
						//approximation, here we do not consider gVal for each single tweet, NUC[u][c]-- was need before the gVal calculation
						double[] gVal = ComputeNeighborLhoodG(users, u, luRefIx[u]);
						int tRefIx = ltRefIx[u];
						for (int n = 0; n < users[u].getDocWords().length; n++) {
							//sample p(c_un=c)
							SampleCommunity_Final(u, n, users, tRefIx,arho, talpha,gVal);
							//sample p(z_un=z)
							SampleTopic_Final(users[u].getDocWords()[n], users[u]
									.getDocTimeStamp()[n], u, n, users, tRefIx,talpha, vbeta);
							if(users[u].getRetweets().containsKey(n))
							{
								tRefIx++;
							}
						}
//						tRefIx+=(users.get(u).getRetweetcount());
					}
					
					//update cij and pop
					assigncijpop(users);
					drawNu();
					
/*					List<int[]> instances = new ArrayList<int[]>();
					instances = getInstance("data/toy2/TMinput/");
					LearnNu(instances, users);*/
					try {
						computeModelParameter();
						saveModel(outputDir,i + 1);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			long endtime = System.nanoTime();
			System.out.println("used "+ (endtime-begintime)+" nano seconds.");			
		}
		executor.shutdownNow();
	}

	public void computeModelParameter() {
		System.out.println("computing model parameters...");

		for (int t = 0; t < T; t++) {
			for (int w = 0; w < V; w++)
				vPhi[t][w] = (float) ((NTW[t][w] + beta) / (SNTW[t] + V * beta));
		}

		for (int a = 0; a < A; a++) {
			for (int t = 0; t < T; t++) {
				theta[a][t] = (float) ((NCT[a][t] + alpha) / (SNCT[a] + T * alpha));
			}
		}
		
		for(int i=0; i<U;i++) {
			for(int a=0; a<A;a++) {
				pai[i][a] = (float) ((NUC[i][a] + rho) / (SNUC[i] + A*rho));
			}
		}		
		
		System.out.println("model parameters are computed");
	}

	private void computeSum(User[] users, int T) {
		SNUC = new float[users.length];
		for (int i = 0; i < users.length; i++) {
			SNUC[i] = MatrixUtil.sumRow(NUC, i);
		}
		SNTW = new float[T];
		for (int t = 0; t < T; t++) {
			SNTW[t] = MatrixUtil.sumRow(NTW, t);
		}
		SNCT = new float[A];
		for(int a=0; a<A; a++)
		{
			SNCT[a]=MatrixUtil.sumRow(NCT, a);
		}
		SNTC = new float[T];
		for(int t=0; t<T; t++)
		{
			SNTC[t] = MatrixUtil.sumColumn(NCT, t, A);
		}
		SPTM = new float[T];
		for(int t=0; t<T; t++)
		{
			SPTM[t] = MatrixUtil.sumRow(popularity,t);
		}
	}

	private void computeTempPrmts(User[] users, int[][] newZ, int[][] newC) {
		for (int i = 0; i < U; i++) {
			for (int j = 0; j < users[i].getDocWords().length; j++) {
				NUC[i][newC[i][j]]++;
				NCT[newC[i][j]][newZ[i][j]]++;
				for (int k = 0; k < users[i].getDocWords()[j].length; k++)
					NTW[newZ[i][j]][users[i].getDocWords()[j][k]]++;
				popularity[Z[i][j]][users[i].getDocTimeStamp()[j]]++;
			}
		}
	}

	public void getResFromLastIteration(User[] users) {
		System.out.println("getting results from last interation...");
		//cleanTempPrmts(users);
		cleanCounter();
		computeTempPrmts(users, Z,C);
	}


	//likelihood of friendship link when sample community/topic
	private double[] ComputeNeighborLhoodG(User[] users, int userInd, int RefIx)
	{
		double[] dLhoodVal = new double[A];
		double sumval = 0D;
		
		for(int a=0; a<A; a++)
		{
			dLhoodVal[a] = 0D;
		}

		int[] neighbours = null;
		
		if(userNeighbour.get(userInd)!=null)
		{
			neighbours = userNeighbour.get(userInd);	//neighbours of user userInd, the one who she follows
			for(int i=0; i< neighbours.length; i++)
			{
				int uNeighbor = neighbours[i];
				//pai_u * pai_v
				//double tweight = uDiscFun(NUC[userInd], NUC[uNeighbor], users.get(userInd).getDocWords().length, users.get(uNeighbor).getDocWords().length);  // discriminant function value between pair <d, nNeighbor>
				double unitlength = (double)(users[userInd].getDocWords().length* users[uNeighbor].getDocWords().length);
				double tweight = uDiscFun(NUC[userInd], NUC[uNeighbor])/unitlength;  // discriminant function value between pair <d, nNeighbor>
				for(int a=0; a<A; a++)
				{
					double fweight = tweight+ (NUC[uNeighbor][a]/unitlength);	//the paiu*paiv after assign community a to user u
					dLhoodVal[a] += (NUC[uNeighbor][a]/unitlength - lambda[RefIx+i]*(fweight*fweight-tweight*tweight));
				}
			}
		}		
		
		if(userInNeighbours.get(userInd)!=null)
		{
			int[][] inneighbours = null;
			//for all the neighbours of user userInd, the one who follow her
			inneighbours = userInNeighbours.get(userInd);
			for(int i=0; i< inneighbours.length; i++)
			{
				int uNeighbor = inneighbours[i][0];
				//pai_u * pai_v
				//double tweight = uDiscFun(NUC[userInd], NUC[uNeighbor], users.get(userInd).getDocWords().length, users.get(uNeighbor).getDocWords().length);  // discriminant function value between pair <d, nNeighbor>
				double unitlength = (double)(users[userInd].getDocWords().length* users[uNeighbor].getDocWords().length);
				double tweight = uDiscFun(NUC[userInd], NUC[uNeighbor])/ unitlength;  // discriminant function value between pair <d, nNeighbor>
				for(int a=0; a<A; a++)
				{
					double fweight = tweight + (NUC[uNeighbor][a]/unitlength);
					dLhoodVal[a] += (NUC[uNeighbor][a]/unitlength - lambda[inneighbours[i][1]]*(fweight*fweight-tweight*tweight));
				}
			}	
		}
			

		for(int a=0; a<A; a++)
		{
			dLhoodVal[a] = Math.exp(0.5*dLhoodVal[a]);	
			sumval += dLhoodVal[a];
		}
			

		for(int a=0; a<A; a++)
		{
			dLhoodVal[a] = dLhoodVal[a]/sumval;
			if(dLhoodVal[a]==0 || Double.isNaN(dLhoodVal[a]) ||Double.isInfinite(dLhoodVal[a]))
			{
				for(int i=0; i<A; i++)
				{
					dLhoodVal[i] = 1;
				}
				return dLhoodVal;
			}

		}			

		return dLhoodVal;
	}
	
	//likelihood of diffusion link when sample community
	private double[] ComputeNeighborLhoodE_C(User[] users, int userInd, int tweetInd, int RefIx)
	{
		double[] dLhoodVal = new double[A];
		TIntObjectHashMap<int[]>  retweetinfo = users[userInd].getRetweets();
		//for all tweets which retweet current tweet
		List<int[]> intweetneighbours = new ArrayList<int[]>();
		int[] upair = new int[]{userInd, tweetInd};
		intweetneighbours = tweetInNeighbours.get(ComUtil.generateKey(upair));
		
		for(int a=0; a<A; a++)
		{
			dLhoodVal[a] = 0D;
		}
		double sumval = 0D;
			

		//for the tweet where this tweet retweet, if any
		if(retweetinfo!=null)
		{
			if(retweetinfo.containsKey(tweetInd))
			{
				int[] p = retweetinfo.get(tweetInd);
				int v = p[0];	//the user where the tweet is from
				int z = Z[userInd][tweetInd];	//the topic label of this retweet
//					double[] f= new double[2*Nf];
//					//f = concatIndFea(userInd,v);
//					f=IndFea(v);
				//calculate sum_c sum_c' (c_hat{u,c}z_hat{c,k}eta{c,c',k}c_hat{v,c'}z_hat{c',k}+popularity{k,t}+nu*f)
				float popu =  (float)popularity[z][users[userInd].getDocTimeStamp()[tweetInd]]/SPTM[z];
				//double weight = tDiscFun(NUC[userInd],NUC[v],users.get(userInd).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu);
				int[] z_hat_k = new int[A];
				z_hat_k = MatrixUtil.getColumn(NCT, z);
				double tweight = tDiscFun(NUC[userInd],NUC[v],z,z_hat_k);
				double unitlength = (double) (users[userInd].getDocWords().length*users[v].getDocWords().length*SNTC[z]*SNTC[z]);
				tweight = (tweight/unitlength)*nu[0] + popu*nu[1];
				for(int nn=0; nn<N2f; nn++)
				{
					tweight+=(nu[nn+2]*tpf[RefIx][nn]);
				}

				for(int a=0; a<A; a++)
				{
					double difweight = 0D;	//difference weight, full-current(i.e., remove current tweet)
					for (int cp=0; cp<A; cp++)
					{
						difweight += eta[z][a][cp]*NUC[v][cp]*z_hat_k[a]*z_hat_k[cp];
					}
					difweight /= unitlength;
					double fweight = tweight + difweight;
					dLhoodVal[a] += (difweight- delta[RefIx]*(fweight*fweight-tweight*tweight));
				}
				RefIx++;

			}
		}
		
		//for all tweets which retweet current tweet
		if(intweetneighbours!=null)
		{
			for(int i=0; i<intweetneighbours.size();i++)
			{
				int v = intweetneighbours.get(i)[0];
				int vtIx = intweetneighbours.get(i)[1];
				int z = Z[v][vtIx];
//					double[] f= new double[2*Nf];
//					//f = concatIndFea(v,userInd);
//					f=IndFea(userInd);
				float popu =  (float)popularity[z][users[v].getDocTimeStamp()[vtIx]]/SPTM[z];
				//double weight = tDiscFun(NUC[userInd],NUC[v],users.get(userInd).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu);
				//dLhoodVal[a] += (weight - delta[intweetneighbours.get(i)[2]]*weight*weight);
				int[] z_hat_k = new int[A];
				z_hat_k = MatrixUtil.getColumn(NCT, z);
				double tweight = tDiscFun(NUC[userInd],NUC[v],z,z_hat_k);
				double unitlength = (double) (users[userInd].getDocWords().length*users[v].getDocWords().length*SNTC[z]*SNTC[z]);

				tweight = (tweight/unitlength)*nu[0] + popu*nu[1];
				for(int nn=0; nn<N2f; nn++)
				{
					tweight+=(nu[nn+2]*tpf[intweetneighbours.get(i)[2]][nn]);
				}

				for(int a=0; a<A; a++)
				{
					double difweight = 0D;	//difference weight, full-current(i.e., remove current tweet)
					for (int cp=0; cp<A; cp++)
					{
						difweight += eta[z][cp][a]*NUC[v][cp]*z_hat_k[a]*z_hat_k[cp];
					}
					difweight /= unitlength;
					double fweight = tweight + difweight;
					dLhoodVal[a] += (difweight- delta[intweetneighbours.get(i)[2]]*(fweight*fweight-tweight*tweight));
				}					

			}
		}


		
		for(int a=0; a<A; a++)
		{
			dLhoodVal[a] = Math.exp((0.5*dLhoodVal[a]));
			sumval+=dLhoodVal[a];
		}
		
		for(int a=0; a<A; a++)
		{
			dLhoodVal[a] = dLhoodVal[a]/sumval;
			if(dLhoodVal[a]==0 || Double.isNaN(dLhoodVal[a]) ||Double.isInfinite(dLhoodVal[a]))
			{
				for(int i=0; i<A; i++)
				{
					dLhoodVal[i] = 1;
				}
				return dLhoodVal;
			}
		}	

		
		return dLhoodVal;
	}
	
	
	//likelihood of diffusion link when sample topic
	private double[] ComputeNeighborLhoodE_T(User[] users, int userInd, int tweetInd, int RefIx)
	{
		double[] dLhoodVal = new double[T];
		TIntObjectHashMap<int[]>  retweetinfo = users[userInd].getRetweets();
		//for all tweets which retweet current tweet
		List<int[]> intweetneighbours = new ArrayList<int[]>();
		int[] upair = new int[]{userInd, tweetInd};
		intweetneighbours = tweetInNeighbours.get(ComUtil.generateKey(upair));
		
		for(int a=0; a<T; a++)
		{
			dLhoodVal[a] = 0D;
		}
		double sumval = 0D;
			

		//for the tweet where this tweet retweet, if any
		if(retweetinfo!=null)
		{
			if(retweetinfo.containsKey(tweetInd))
			{
				int[] p = retweetinfo.get(tweetInd);
				int v = p[0];	//the user where the tweet is from
				int z = Z[userInd][tweetInd];	//the topic label of this retweet
				int Cu = C[userInd][tweetInd];
				int Cv = C[v][p[1]];
//					double[] f= new double[2*Nf];
//					//f = concatIndFea(userInd,v);
//					f=IndFea(v);
				//calculate sum_c sum_c' (c_hat{u,c}z_hat{c,k}eta{c,c',k}c_hat{v,c'}z_hat{c',k}+popularity{k,t}+nu*f)
				//current tweet's popularity
				float popu =  (float)popularity[z][users[userInd].getDocTimeStamp()[tweetInd]]/SPTM[z];
				//double weight = tDiscFun(NUC[userInd],NUC[v],users.get(userInd).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu);
				int[] z_hat_k = new int[A];
				z_hat_k = MatrixUtil.getColumn(NCT, z);
				double tweight = tDiscFun(NUC[userInd],NUC[v],z,z_hat_k);
				double unitlength = (double) (users[userInd].getDocWords().length*users[v].getDocWords().length*SNTC[z]*SNTC[z]);
				tweight = (tweight/unitlength)*nu[0] + popu*nu[1];
				for(int nn=0; nn<N2f; nn++)
				{
					tweight+=(nu[nn+2]*tpf[RefIx][nn]);
				}

				for(int t=0; t<T; t++)
				{
					double difweight = 0D;	//difference weight, full-current(i.e., remove current tweet)
					
					if(Cu==Cv)
					{
						difweight = eta[z][Cu][Cv]*NUC[v][Cv]*NUC[userInd][Cu]+2*(eta[z][Cu][Cv]*NUC[v][Cv]*NUC[userInd][Cu]*z_hat_k[Cu]);
					}
					else
					{
						difweight = eta[z][Cu][Cv]*z_hat_k[Cv]*NUC[v][Cv]*NUC[userInd][Cu];
					}
					difweight/=unitlength;
					double fweight = 0D;
					if(t==z)
					{
						fweight = tweight + difweight -popu + (float)(popularity[z][users[userInd].getDocTimeStamp()[tweetInd]]+1)/(SPTM[z]+1);
					}
					else
					{
						fweight = tweight + difweight;
					}
					dLhoodVal[t] += (difweight- delta[RefIx]*(fweight*fweight-tweight*tweight));						
				}
				RefIx++;

			}
		}
		
		//for all tweets which retweet current tweet
		if(intweetneighbours!=null)
		{
			for(int i=0; i<intweetneighbours.size();i++)
			{
				int v = intweetneighbours.get(i)[0];
				int vtIx = intweetneighbours.get(i)[1];
				int z = Z[v][vtIx];
				int Cu = C[userInd][tweetInd];
				int Cv = C[v][vtIx];
//					double[] f= new double[2*Nf];
//					//f = concatIndFea(v,userInd);
//					f=IndFea(userInd);
				//the tweet who retweet current tweet, hence it is the popularity at that moment, have nothing to do currently
				float popu =  (float)popularity[z][users[v].getDocTimeStamp()[vtIx]]/SPTM[z];
				//double weight = tDiscFun(NUC[userInd],NUC[v],users.get(userInd).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu);
				//dLhoodVal[a] += (weight - delta[intweetneighbours.get(i)[2]]*weight*weight);
				int[] z_hat_k = new int[A];
				z_hat_k = MatrixUtil.getColumn(NCT, z);
				double tweight = tDiscFun(NUC[userInd],NUC[v],z,z_hat_k);
				double unitlength = (double) (users[userInd].getDocWords().length*users[v].getDocWords().length*SNTC[z]*SNTC[z]);
				tweight = (tweight/unitlength)*nu[0] + popu*nu[1];
				for(int nn=0; nn<N2f; nn++)
				{
					tweight+=(nu[nn+2]*tpf[intweetneighbours.get(i)[2]][nn]);
				}

				for(int t=0; t<T; t++)
				{
					double difweight = 0D;
					if(Cu==Cv)
					{
						difweight = eta[z][Cv][Cu]*NUC[v][Cv]*NUC[userInd][Cu]+2*(eta[z][Cv][Cu]*NUC[v][Cv]*NUC[userInd][Cu]*z_hat_k[Cu]);
					}
					else
					{
						difweight = eta[z][Cv][Cu]*z_hat_k[Cv]*NUC[v][Cv]*NUC[userInd][Cu];
					}
					difweight/=unitlength;
					double fweight = tweight + difweight;
					
					dLhoodVal[t] += (difweight- delta[intweetneighbours.get(i)[2]]*(fweight*fweight-tweight*tweight));						
				}				
			}
		}


		
		for(int a=0; a<T; a++)
		{
			dLhoodVal[a] = Math.exp((0.5*dLhoodVal[a]));
			sumval+=dLhoodVal[a];
		}
		
		for(int a=0; a<T; a++)
		{
			dLhoodVal[a] = dLhoodVal[a]/sumval;
			if(dLhoodVal[a]==0 || Double.isNaN(dLhoodVal[a]) ||Double.isInfinite(dLhoodVal[a]))
			{
				for(int i=0; i<T; i++)
				{
					dLhoodVal[i] = 1;
				}
				return dLhoodVal;
			}
		}	

		
		return dLhoodVal;
	}

	
	//to record the tweets which are not retweets and haven't been retweeted by others
	private boolean getNoretweets(int u, int n, User[] users)
	{
		boolean res = false;
		boolean hasretweet = users[u].getRetweets().containsKey(n);
		//for all tweets which retweet current tweet
		List<int[]> intweetneighbours = new ArrayList<int[]>();
		int[] upair = new int[]{u, n};
		intweetneighbours = tweetInNeighbours.get(ComUtil.generateKey(upair));
		

		if(!hasretweet&&intweetneighbours==null)
		{
			res = true;
		}
		return res;
	}
	
	public void printretweet(User[] users, String outputDir)
	{
		ArrayList<String> retweet = new ArrayList<String>();
		for (int u = 0; u < U; u++) {
			for (int n = 0; n < users[u].getDocWords().length; n++) {
				if(!getNoretweets(u, n, users))
					retweet.add(u+"\t"+n);
			}			
		}
		FileUtil.writeLines(outputDir, retweet);
	}
	
	private boolean SampleCommunity(int u, int n, User[] users, int tRefIx,double arho, double talpha, double[] gVal) {
		
		int topic = Z[u][n];
		int community = C[u][n];

		NUC[u][community]--;
		SNUC[u]--;
		NCT[community][topic]--;
		SNCT[community]--;
		SNTC[topic]--;

		// get p(Z_{u,n} = z|Z_c, W, Y, I)
		double[] pt = new double[A];
		double NUCsumRowU = SNUC[u];
		boolean noret = false;
		double[] eVal = null;
		if(!retweet.contains(u+"\t"+n))
		{
			noret = true;
		}
		else
		{			
			eVal = ComputeNeighborLhoodE_C(users, u, n, tRefIx);
		}		

		for (int a = 0; a < A; a++) {
			double p1 = (double) (NUC[u][a] + rho) / (NUCsumRowU + arho);			
			double p2 = (double) (NCT[a][topic] + alpha) / (SNCT[a] + talpha);
			if(noret)
				pt[a] = p1 * p2 * gVal[a];
			else
				pt[a] = p1 * p2 * gVal[a] * eVal[a];
		}
		//System.out.print("Sample community ");
		// cummulate multinomial parameters
		int sample = ComUtil.sample(pt, A, gVal, eVal);
//		assert (sample >= 0 && sample < A) : "sample community value error:" + sample+ pt.toString();

		C[u][n] = sample;
		community = sample;

		NUC[u][community]++;
		SNUC[u]++;
		NCT[community][topic]++;
		SNCT[community]++;
		SNTC[topic]++;
		return true;
	}
	
	private boolean SampleTopic(int[] words, int ts, int u, int n, User[] users, int tRefIx, double talpha, double vbeta) {
		
		int topic = Z[u][n];
		int community = C[u][n];
		int timestamp = users[u].getDocTimeStamp()[n];
		// get words and their count in [u,n]
		ArrayList<Integer> tempUniqueWords = new ArrayList<Integer>();
		ArrayList<Integer> tempCounts = new ArrayList<Integer>();
		uniqe(words, tempUniqueWords, tempCounts);

		NCT[community][topic]--;
		SNCT[community]--;
		SNTC[topic]--;
		for (int w1 = 0; w1 < tempUniqueWords.size(); w1++) {
			NTW[topic][tempUniqueWords.get(w1)] -= tempCounts.get(w1);
			SNTW[topic] -= tempCounts.get(w1);
		}
		popularity[topic][timestamp] --;
		SPTM[topic]--;
		
		// get p(Z_{u,n} = z|Z_c, W, Y, I)
		double[] pt = new double[T];
		double NCTsumRowC = SNCT[community];
		boolean noret = false;
		double[] eVal = null;
		if(!retweet.contains(u+"\t"+n))
		{
			noret = true;
		}
		else
		{			
			eVal = ComputeNeighborLhoodE_T(users, u, n, tRefIx);
		}	
		
		for (int i = 0; i < T; i++) {
			int wcount = 0;
			double p1 = (double) (NCT[community][i] + alpha) / (NCTsumRowC + talpha);
			double p2 = 1.0D;
			for (int w = 0; w < tempUniqueWords.size(); w++) {
				int tempvalue = NTW[i][tempUniqueWords.get(w)];
				// double sumRow = MatrixUtil.sumRow(NTW, i);
				double NTWsumRowT = SNTW[i];
				// checkEqual(sumRow, MatrixUtil.sumRow(NTW, i), "NTW");
				for (int numC = 0; numC < tempCounts.get(w); numC++) {
					p2 = p2 * ((double) (tempvalue + beta + numC) / ((double) NTWsumRowT
									+ vbeta + wcount));
					wcount++;
				}
			}
			if(p2==0)
				p2=1;
			if(noret)
				pt[i] = p1 * p2;
			else
				pt[i] = p1 * p2 * eVal[i];

		}

		// cummulate multinomial parameters
		int sample = ComUtil.sample(pt, T);
//		assert (sample >= 0 && sample < T) : "sample topic value error:" + sample+pt.toString();

		Z[u][n] = sample;
		topic = sample;

		// update NTW[T][W](y=1) NTI[T][M] NUT[U][T] in {u,n}
		NCT[community][topic]++;
		SNCT[community]++;
		for (int w1 = 0; w1 < tempUniqueWords.size(); w1++) {
			NTW[topic][tempUniqueWords.get(w1)] += tempCounts.get(w1);
			SNTW[topic] += tempCounts.get(w1);
		}
		popularity[topic][timestamp] ++;
		SPTM[topic]++;
		SNTC[topic]++;
		tempUniqueWords.clear();
		tempCounts.clear();		
		
		return true;
	}
	
	//after sampling both community and topic, we record the new cij and pop
	public void assigncijpop(User[] users)
	{		
		int count = 0;
		for(int u=0; u<U; u++)
		{
			TIntObjectHashMap<int[]>  retweetinfo = users[u].getRetweets();
			if(retweetinfo!=null)
			{
				for(int utIx: retweetinfo.keys())
				{
					//positive cij and pop
					int[] p = retweetinfo.get(utIx);
					int v = p[0];	//the user where the tweet is from
					int z = Z[u][utIx];
					int[] z_hat_k = new int[A];
					z_hat_k = MatrixUtil.getColumn(NCT, z);
					double tweight = tDiscFun(NUC[u],NUC[v],z,z_hat_k);
					double unitlength = (double) (users[u].getDocWords().length*users[v].getDocWords().length*SNTC[z]*SNTC[z]);
					cijp[count] = (float)(tweight/unitlength);
					popp[count] = (float)popularity[z][users[u].getDocTimeStamp()[utIx]]/SPTM[z];
					
					//negative cij and pop
					v = nvs[count];
					z = Z[v][0];
					z_hat_k = MatrixUtil.getColumn(NCT, z);
					tweight = tDiscFun(NUC[u],NUC[v],z,z_hat_k);
					unitlength = (double) (users[u].getDocWords().length*users[v].getDocWords().length*SNTC[z]*SNTC[z]);
					cijn[count] = (float)(tweight/unitlength);
					popn[count] = (float)popularity[z][users[u].getDocTimeStamp()[utIx]]/SPTM[z];
					
					count++;
				}
			}

		}		
	}
	
	private boolean SampleCommunity_Final(int u, int n, User[] users, int tRefIx, double arho, double talpha, double[] gVal) {
		int topic = Z[u][n];
		int community = C[u][n];
		
		NUC[u][community]--;
		SNUC[u]--;
		NCT[community][topic]--;
		SNCT[community]--;
		SNTC[topic] --;
		// get p(Z_{u,n} = z|Z_c, W, Y, I)
		double[] pt = new double[A];
		double NUCsumRowU = SNUC[u];
		boolean noret = false;
		double[] eVal = null;
		if(!retweet.contains(u+"\t"+n))
		{
			noret = true;
		}
		else
		{			
			eVal = ComputeNeighborLhoodE_C(users, u, n, tRefIx);
		}

		
		for (int a = 0; a < A; a++) {
			double p1 = (double) (NUC[u][a] + rho) / (NUCsumRowU + arho);
			double p2 = (double) (NCT[a][topic] + alpha) / (SNCT[a]+ talpha);
			if(noret)
				pt[a] = p1 * p2 * gVal[a];
			else
				pt[a] = p1 * p2 * gVal[a] * eVal[a];
		}
		// cummulate multinomial parameters
		//int sample = ComUtil.sample(pt, A);
		int sample = ComUtil.findMax(pt, A);
//		assert (sample >= 0 && sample < A) : "sample community value error:" + sample;

		C[u][n] = sample;
		community = sample;
		NUC[u][community]++;
		SNUC[u]++;
		NCT[community][topic]++;
		SNCT[community]++;
		SNTC[topic]++;
		return true;
	}
	
	private boolean SampleTopic_Final(int[] words, int ts, int u, int n, User[] users, int tRefIx,double talpha, double vbeta) {
		int topic = Z[u][n];
		int community = C[u][n];
		int timestamp = users[u].getDocTimeStamp()[n];
		// get words and their count in [u,n]
		ArrayList<Integer> tempUniqueWords = new ArrayList<Integer>();
		ArrayList<Integer> tempCounts = new ArrayList<Integer>();
		uniqe(words, tempUniqueWords, tempCounts);

		NCT[community][topic]--;
		SNCT[community]--;
		for (int w1 = 0; w1 < tempUniqueWords.size(); w1++) {
			NTW[topic][tempUniqueWords.get(w1)] -= tempCounts.get(w1);
			SNTW[topic] -= tempCounts.get(w1);
		}
		popularity[topic][timestamp] --;
		SPTM[topic]--;
		SNTC[topic]--;
		
		// get p(Z_{u,n} = z|Z_c, W, Y, I)
		double[] pt = new double[T];
		double NCTsumRowC = SNCT[community];
		boolean noret = false;
		double[] eVal = null;
		if(!retweet.contains(u+"\t"+n))
		{
			noret = true;
		}
		else
		{			
			eVal = ComputeNeighborLhoodE_T(users, u, n, tRefIx);
		}	
		for (int i = 0; i < T; i++) {
			int wcount = 0;
			double p1 = (double) (NCT[community][i] + alpha) / (NCTsumRowC + talpha);
			double p2 = 1.0D;
			for (int w = 0; w < tempUniqueWords.size(); w++) {
				int tempvalue = NTW[i][tempUniqueWords.get(w)];
				double NTWsumRowT = SNTW[i];
				for (int numC = 0; numC < tempCounts.get(w); numC++) {
					p2 = p2 * ((double) (tempvalue + beta + numC) / ((double) NTWsumRowT
									+ vbeta + wcount));
					wcount++;
				}
			}
			if(p2==0)
				p2=1;
			if(noret)
				pt[i] = p1 * p2;
			else
				pt[i] = p1 * p2 * eVal[i];
		}

		// cummulate multinomial parameters
		//int sample = ComUtil.sample(pt, T);
		int sample = ComUtil.findMax(pt,T);
//		assert (sample >= 0 && sample < T) : "sample topic value error:" + sample;

		Z[u][n] = sample;
		topic = sample;

		// update NTW[T][W](y=1) NTI[T][M] NUT[U][T] in {u,n}
		NCT[community][topic]++;
		SNCT[community]++;
		for (int w1 = 0; w1 < tempUniqueWords.size(); w1++) {
			NTW[topic][tempUniqueWords.get(w1)] += tempCounts.get(w1);
			SNTW[topic] += tempCounts.get(w1);
		}
		popularity[topic][timestamp] ++;
		SPTM[topic]++;
		SNTC[topic]++;
		tempUniqueWords.clear();
		tempCounts.clear();
		return true;
	}
	
	//lambda for user-user following friendship, assign lambda value for each user-user pair
	private void drawLambda(User[] users)
	{
		double uDiscFuncVal;
		int uPIx = 0;
		for ( int u = 0; u < U; u++) {
			if(userNeighbour.get(u)!=null)
			{
				for ( int v=0; v<userNeighbour.get(u).length; v++ ) {
					int uNeighbor = userNeighbour.get(u)[v];
					//pai_u * pai_v
					uDiscFuncVal = uDiscFun(NUC[u], NUC[uNeighbor], users[u].getDocWords().length, users[uNeighbor].getDocWords().length);  // discriminant function value between pair <d, nNeighbor>
					lambda[uPIx] = (float) m_uPGsampler[uPIx].nextPG(1, uDiscFuncVal);
					uPIx ++;
				}
			}
		}
	}
	
	//delta for retweet, assign delta for each retweet
	private void drawDelta(User[] users)
	{
		double tDiscFuncVal;
		int tPIx = 0;
		for(int u=0; u<U; u++)
		{
			TIntObjectHashMap<int[]>  retweetinfo = users[u].getRetweets();
			if(retweetinfo!=null)
			{
				for(int utIx: retweetinfo.keys())
				{
					int[] p = retweetinfo.get(utIx);
					int v = p[0];	//the user where the tweet is from
					int z = Z[u][utIx];	//the topic label of this retweet
					
//					double[] f= new double[2*Nf];
//					//f = concatIndFea(u,v);
//					f=IndFea(v);
					float popu =  (float)popularity[z][users[u].getDocTimeStamp()[utIx]]/SPTM[z];
					
					int[] z_hat_k = new int[A];
					z_hat_k = MatrixUtil.getColumn(NCT, z);
					tDiscFuncVal = tDiscFun(NUC[u],NUC[v],z,z_hat_k);
					double unitlength = (double) (users[u].getDocWords().length*users[v].getDocWords().length*SNTC[z]*SNTC[z]);
					tDiscFuncVal = (tDiscFuncVal/unitlength)*nu[0] + popu*nu[1];
					for(int nn=0; nn<N2f; nn++)
					{
						tDiscFuncVal+=(nu[nn+2]*tpf[tPIx][nn]);
					}					
					
					if(!Double.isInfinite(tDiscFuncVal)&&!Double.isNaN(tDiscFuncVal))
					{
						delta[tPIx] = (float) m_tPGsampler[tPIx].nextPG(1, tDiscFuncVal);
					}
					tPIx ++;
				}
			}
		}
	}
	
	//generate fixed user feature matrix for all retweets
/*	public void generateUserFeature(User[] users, String output)
	{
		Random rand = new Random();
		HashMap<Integer, HashSet<Integer>> userretweets = new HashMap<Integer, HashSet<Integer>>();
		int[] userretweetscount = new int[U];
		int cnlt= 0;
		for(int u=0; u<U; u++)
		{
			TIntObjectHashMap<int[]>  retweetinfo = users[u].getRetweets();
			if(retweetinfo!=null)
			{
				for(int utIx: retweetinfo.keys())
				{
					cnlt++;
				}
			}
		}
		double[][] tp = new double[cnlt][N2f];	//positive instances
		double[][] tn = new double[cnlt][N2f];	//negative instances
		int count = 0;
		for(int u=0; u<U; u++)
		{
			int indi = 0;
			TIntObjectHashMap<int[]>  retweetinfo = users[u].getRetweets();
			if(retweetinfo!=null)
			{
				HashSet<Integer> r = new HashSet<Integer>();
				for(int utIx: retweetinfo.keys())
				{
					int[] p = retweetinfo.get(utIx);
					int v = p[0];	//the user where the tweet is from
					for (int i=0; i<Nf; i++)	//concat u and v's feature
					{
						tp[count][i] = indfea[u][i];
						tp[count][Nf+i] = indfea[v][i];
					}
					r.add(v);
					count++;
					indi++;
				}
				userretweets.put(u, r);
			}
			userretweetscount[u] = indi;
		}
		
		//sample negative instances from v who never retweet u
		//record the randomly sampled negative users, later to calculate cij
		ArrayList<Integer> nvs = new ArrayList<Integer>();
		count = 0;
		for(int u=0; u<U; u++)
		{
			int indi = userretweetscount[u];
			if(indi>0)
			{
				HashSet<Integer> vs = new HashSet<Integer>();
				vs = userretweets.get(u);
				int indic = 0;
				while(indic<indi)	//same amount of negarive instances as positive instances
				{
					int rv = rand.nextInt(U);
					if(!vs.contains(rv))	//rv never retweet u
					{
						for (int i=0; i<Nf; i++)	//concat u and rv's feature
						{
							tn[count][i] = indfea[u][i];
							tn[count][Nf+i] = indfea[rv][i];
						}
						count++;
						indic++;
						nvs.add(rv);
					}
				}				
			}
		}
		
		FileUtil.writeDMatrix(output+"tp.txt", tp);
		FileUtil.writeDMatrix(output+"tn.txt", tn);
		FileUtil.writeLines(output+"negavs.txt", nvs);
		
	}
	*/
	
	//nu for parameters of the logistic regression, nu*f+cif*popu
	private void drawNu()
	{		
		int M = 2*NLt, N = N2f+2;
		Feature[][] featureMatrix = new Feature[M][N];
        for(int i=0; i<NLt; i++)
        {
        	featureMatrix[i][0] = new FeatureNode(1,cijp[i]);
        	featureMatrix[i][1] = new FeatureNode(2,popp[i]);
        	for(int j=0; j<N2f; j++)
        	{
        		featureMatrix[i][j+2]=new FeatureNode(j+3,tpf[i][j]);
        	}
        	featureMatrix[i+NLt][0] = new FeatureNode(1,cijn[i]);
        	featureMatrix[i+NLt][1] = new FeatureNode(2,popn[i]);
        	for(int j=0; j<N2f; j++)
        	{
        		featureMatrix[i+NLt][j+2]=new FeatureNode(j+3,tnf[i][j]);
        	}
        }
        
        //loading target value
        double[] targetValue = new double[M];
        for(int i=0; i<NLt; i++)
        {
        	targetValue[i] = 1;	//positive instances
        	targetValue[NLt+i]=-1;	//negative instances
        }
        
        Problem problem = new Problem();
        problem.l = M; // number of training examples
        problem.n = N; // number of features
        problem.x = featureMatrix; // feature nodes
        problem.y = targetValue; // target values

        SolverType solver = SolverType.L1R_LR; // -s 0
        double C = 1.0;    // cost of constraints violation
        double eps = 0.01; // stopping criteria
//        int maxIter = 1500;           
      
        Parameter parameter = new Parameter(solver, C, eps);
        de.bwaldvogel.liblinear.Model lrmodel = Linear.train(problem, parameter);      
        nu = lrmodel.getFeatureWeights();
	}
	
	//delta for retweet, assign delta for each retweet
	private void drawEta(User[] users)
	{

		for (int t = 0; t < T; t++) {
			for(int i=0; i<A;i++)
			{
				for(int j=0; j<A; j++)
					eta[t][i][j] = 0;
			}			
		}
		for(int u=0; u<U; u++)
		{
			TIntObjectHashMap<int[]>  retweetinfo = users[u].getRetweets();
			if(retweetinfo!=null)
			{
				for(int utIx: retweetinfo.keys())
				{
					int[] p = retweetinfo.get(utIx);
					int v = p[0];	//the user where the tweet is from
					int vtIx = p[1];
					int z = Z[u][utIx];	//the topic label of this retweet
					eta[z][C[u][utIx]][C[v][vtIx]] ++;
				}
			}
		}
		
		//normalize to make eta<1
		double maxeta = 0;
		for (int t = 0; t < T; t++) {
			for(int i=0; i<A;i++)
			{
				for(int j=0; j<A; j++)
					
				{
					if(eta[t][i][j]>maxeta)
						maxeta = eta[t][i][j];
				}
			}			
		}
		
		for (int t = 0; t < T; t++) {
			for(int i=0; i<A;i++)
			{
				for(int j=0; j<A; j++)
					eta[t][i][j] = eta[t][i][j]/maxeta;
			}			
		}
	}
	
//	//concat the feature for user u and user v
//	private double[] concatIndFea(int u, int v)
//	{
//		double[] f= new double[2*Nf];
//		for(int i=0; i<Nf; i++)
//		{
//			f[i] = ind_fea[u][i]*paranu;
//			f[Nf+i] = ind_fea[v][i]*paranu;
//		}
//		return f;
//	}
	
//	//use only the statistics of users being retweeted
//	private double[] IndFea(int u)
//	{
//		double[] f= new double[Nf];
//		for(int i=0; i<Nf; i++)
//		{
//			f[i] = ind_fea[u][i]*paranu;
//		}
//		return f;
//	}
	
	
	//calculate c_hat u * c_hat v
	private double uDiscFun(int[] c_hat_u, int[] c_hat_v, int N_u, int N_v)
	{
		double sum = 0D;
		for(int i=0; i< c_hat_u.length; i++)
		{
			sum += (double) c_hat_u[i]*c_hat_v[i];
		}
		sum /= (double)(N_u*N_v);
		return sum;
	}
	
	//calculate c_hat u * c_hat v (nuc(u,c)++), this is to speed up ComputeNeighborLhoodG
	private double uDiscFun(int[] c_hat_u, int[] c_hat_v)
	{
		double sum = 0D;
		for(int i=0; i< c_hat_u.length; i++)
		{
			sum += (double) c_hat_u[i]*c_hat_v[i];
		}
		return sum;
	}
	
	//calculate sum_c sum_c' (c_hat{u,c}z_hat{c,k}eta{c,c',k}c_hat{v,c'}z_hat{c',k}+popularity{k,t}+nu*f)
	private double tDiscFun(int[] c_hat_u, int[] c_hat_v, int N_u, int N_v, float N_k, int topick, int[] z_hat_k, float popularity/*, double[] nu, double[] f*/)
	{
		double res = 0D;
		for(int c=0; c<A; c++)
		{
			for(int cp=0; cp<A; cp++)
			{
				res+=(double)(c_hat_u[c]*z_hat_k[c]*eta[topick][c][cp]*c_hat_v[cp]*z_hat_k[cp]);
			}
		}
		res /= (double)(N_u*N_v*N_k*N_k);
		res += popularity;
		return res;
	}
	
	//calculate sum_c sum_c' (c_hat{u,c}z_hat{c,k}eta{c,c',k}c_hat{v,c'}z_hat{c',k}+popularity{k,t}+nu*f), to speed up computeNeighbourhoodEC
	private double tDiscFun(int[] c_hat_u, int[] c_hat_v, int topick,int[] z_hat_k)
	{
		double res = 0D;
		for(int c=0; c<A; c++)
		{
			for(int cp=0; cp<A; cp++)
			{
				res+=(double)(c_hat_u[c]*eta[topick][c][cp]*c_hat_v[cp]*z_hat_k[c]*z_hat_k[cp]);
			}
		}
		return res;
	}
	
	private double tDiscFun(int[] c_hat_u, int[] c_hat_v, int N_u, int N_v, float N_k, int topick, int[] z_hat_k, float popularity, double[] nu, double[] f)
	{
		double res = 0D;
		for(int c=0; c<A; c++)
		{
			for(int cp=0; cp<A; cp++)
			{
				res+=(double)(c_hat_u[c]*z_hat_k[c]*eta[topick][c][cp]*c_hat_v[cp]*z_hat_k[cp]);
			}
		}
		res /= (double)(N_u*N_v*N_k*N_k);
		res += (popularity*+MatrixUtil.vectorTimes(nu, f));			//need to normalize to make the three factors in the same scale
		return res;
	}
		

	public static void uniqe(int[] words, 
			ArrayList<Integer> tempUniqueWords, ArrayList<Integer> tempCounts) {
		for (int i = 0; i < words.length; i++) {
			if (tempUniqueWords.contains(words[i])) {
				int index = tempUniqueWords.indexOf(words[i]);
				tempCounts.set(index, tempCounts.get(index) + 1);
			} else {
				tempUniqueWords.add(words[i]);
				tempCounts.add(1);
			}

		}
	}
	
	public double classify(int z, int u, int v, int t, double[] f, ArrayList<User> users)
	{
		double res = 0;
		float popu =  (float)popularity[z][t]/SPTM[z];
		double prob = tDiscFun(NUC[u],NUC[v],users.get(u).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu,nu,f);
		
		res = 1/(1+Math.exp(-prob));
		return res;
	}
	
	List<int[]> getInstance(String dir)
	{
		List<int[]> instances = new ArrayList<int[]>();
		ArrayList<String> Pinstance = new ArrayList<String>();
		ArrayList<String> Ninstance = new ArrayList<String>();
		FileUtil.readLines(dir+"PInstance.txt", Pinstance);
		FileUtil.readLines(dir+"NInstance.txt", Ninstance);
		for (String ins: Pinstance)
		{
			String[] terms = ins.split(",");
			int[] in = new int[5];
			for(int i=0; i<4; i++)
			{
				in[i] = Integer.parseInt(terms[i]);
			}
			in[4] = 1;
			instances.add(in);
		}
		for (String ins: Ninstance)
		{
			String[] terms = ins.split(",");
			int[] in = new int[5];
			for(int i=0; i<4; i++)
			{
				in[i] = Integer.parseInt(terms[i]);
			}
			in[4] = 0;
			instances.add(in);
		}
		
		return instances;
	}
//	//logistic regression to learn \nu
//	public void LearnNu(List<int[]> instances, ArrayList<User> users)
//	{
//		double rate = 0.00035;
//		for (int n = 0; n < 2000; n++)
//		{
//			double lik = 0.0;
//			for(int i=0; i<instances.size(); i++)
//			{
//				int[] ins = instances.get(i);
//				double[] f= new double[2*Nf];
//				//f=IndFea(ins[2]);
//				f = concatIndFea(ins[0],ins[2]);
//				double predicted = classify(Z[ins[0]][ins[1]],ins[0],ins[2],ins[3], f,users);
//				int label = ins[4];
//				
//				for (int j=0; j<nu.length; j++) {
//                    nu[j] = nu[j] + rate * (label - predicted) * f[j];
//                }
//				 lik += label * Math.log(classify(Z[ins[0]][ins[1]],ins[0],ins[2],ins[3], f,users)) + (1-label) * Math.log(1- classify(Z[ins[0]][ins[1]],ins[0],ins[2],ins[3],f, users));
//			}
//           // System.out.println("iteration: " + n + " " + Arrays.toString(nu) + " mle: " + lik);
//		}
//	}
	
//	//predict if the retweet will happen
//	public double ProbRetweet(int z, int u, int v, int t, ArrayList<User> users)
//	{
//		double res = 0;
//		double[] f= new double[2*Nf];
//		f = concatIndFea(u,v);
//		//f=IndFea(v);
//		float popu =  (float)popularity[z][t]/SPTM[z];
//		double prob = tDiscFun(NUC[u],NUC[v],users.get(u).getDocWords().length, users.get(v).getDocWords().length,/*SNTC[z],*/z,MatrixUtil.getColumn(NCT, z),popu,nu,f);
//		
//		res = 1/(1+Math.exp(-prob));
//		return res;
//	}
	
//	public void PrintRetweetPro(ArrayList<User> users)
//	{
//		for (int u = 0; u < U; u++) {	
//			for (int v=0; v<U;v++) {
//				if(u!=v)
//				{
//					for(int z=0; z<T;z++)
//					{
//						//for(int t=0; t<M;t++)
//						//{
//							System.out.println(/*"The probability of user "+(u+1)+" retweet user "+(v+1)+" on topic "+z+" is "+*/ProbRetweet(z, u, v, 1, users));
//						//}
//					}
//					System.out.println();
//				}
//			}
//		}
//	}

	/**
	 * output Model paramters
	 * */
	public void outputModelRes() {
		// output Z
		System.out.println("Z[u][n]: ");
		MatrixUtil.printArray(Z);
	}

	public void outTaggedDoc(ArrayList<User> users,
			/*ArrayList<String> uniWordMap,*/ String outputDir) {
		ArrayList<String> datalines = new ArrayList<String>();
		for (int i = 0; i < users.size(); i++) {
			for (int j = 0; j < users.get(i).getDocWords().length; j++) {
				String tmpline = "Community "+C[i][j]+", Topic " + Z[i][j] + ": ";
				tmpline += users.get(i).getDocTimeStamp()[j] +" ";
				for (int k1 = 0; k1 < users.get(i).getDocWords()[j].length; k1++) {
						tmpline += /*uniWordMap.get(*/
								users.get(i).getDocWords()[j][k1]/*)*/ +" ";
				}
				datalines.add(tmpline);
				System.out.println(tmpline);
			}
			System.out.println("");
			FileUtil.writeLines(outputDir + users.get(i).getId(),
					datalines);
			datalines.clear();
		}
	}

	void saveModelRes(String string) throws Exception {
		BufferedWriter writer = null;
		writer = new BufferedWriter(new FileWriter(new File(string)));
		writer.write("Z[u][n]: \n");
		for (int i = 0; i < Z.length; i++) {
			for (int j = 0; j < Z[i].length; j++)
				writer.write(Z[i][j] + "\t");
			writer.write("\n");
		}
		writer.flush();
		writer.close();
	}

	public boolean saveModel(String output, int iter) throws Exception {
		output = output+"iter"+iter+"/";
		(new File(output)).mkdirs();
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.theta")));
		ModelComFunc.writeData(theta, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.pai")));
		ModelComFunc.writeData(pai, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.eta")));
		ModelComFunc.writeData(eta, writer);
		writer.close();

		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.vPhi")));
		ModelComFunc.writeData(vPhi, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.nu")));
		ModelComFunc.writeData(nu, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.lambda")));
		ModelComFunc.writeData(lambda, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.delta")));
		ModelComFunc.writeData(delta, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.nuc")));
		ModelComFunc.writeData(NUC, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.nct")));
		ModelComFunc.writeData(NCT, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.ntw")));
		ModelComFunc.writeData(NTW, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.popu")));
		ModelComFunc.writeData(popularity, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.Z")));
		ModelComFunc.writeData(Z, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.C")));
		ModelComFunc.writeData(C, writer);
		writer.close();

		return true;
	}
	
	public boolean saveModel(String output/*, ArrayList<String> uniWordMap*/) throws Exception {
		ArrayList<Integer> rankList = new ArrayList<Integer>();
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.theta")));
		ModelComFunc.writeData(theta, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.pai")));
		ModelComFunc.writeData(pai, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.eta")));
		ModelComFunc.writeData(eta, writer);
		writer.close();

		/*writer = new BufferedWriter(new FileWriter(new File(output
				+ "model-topic-words.txt")));
		for (int t = 0; t < vPhi.length; t++) {
			ComUtil.getTop(vPhi[t], rankList, 20);
			writer.write("Topic " + t + "\n");
			ModelComFunc.writeData(vPhi[t], uniWordMap, rankList, writer, "\t");
			rankList.clear();
		}
		writer.close();*/
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.vPhi")));
		ModelComFunc.writeData(vPhi, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.nu")));
		ModelComFunc.writeData(nu, writer);
		writer.close();

		return true;
	}


	/*public void output(User[] users, ArrayList<String> uniWords,
			ArrayList<String> uniItems, String outputDir) {
		ArrayList<String> datalines = new ArrayList<String>();
		for (int i = 0; i < users.length; i++) {
			for (int j = 0; j < users[i].getDocWords().length; j++) {
				String tmpline = "";
				for (int k1 = 0; k1 < users.get(i).getDocWords()[j].length; k1++) {
					tmpline += uniWords.get(users.get(i).getDocWords()[j][k1])
							+ " ";
				}
				tmpline += uniItems.get(users.get(i).getDocTimeStamp()[j]);
				datalines.add(tmpline);
			}
			FileUtil.writeLines(outputDir + users.get(i).getId() + ".txt",
					datalines);
			datalines.clear();
		}
	}*/
	
	class ThreadUser extends Thread{
		private Model m;
		double arho;
		double talpha;
		double vbeta;
		User[] users;
		int[] lu;
		boolean isfinish;
		
		ThreadUser(Model m, User[] users, double arho, double talpha, double vbeta, int[] lu){
			this.m = m;
			this.talpha = talpha;
			this.vbeta = vbeta;
			this.arho = arho;
			this.users = users;
			this.lu = lu;
			this.isfinish = false;
		}
		
		public void UpdateData(Model m)
		{
			this.m = m;
			this.isfinish = false;
		}
		public boolean isFinish()
		{
			return isfinish;
		}
		public void run()
		{
			m.SampleCT(users, arho, talpha, vbeta,lu);
			this.isfinish = true;
		}
	}

}
