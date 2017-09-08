package cpdtoy;
/*
 * 2016-05-16 ** Clean version parallel, divide dataset by documents, for tweets dataset
 * For community profiling, sigmoid link function for user-user and tweet-tweet, polya-gamma data agumented & gibbs sampling
 */
import gnu.trove.map.hash.TIntObjectHashMap;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import common.FileUtil;
import common.ModelUtil;
import common.Stopwords;

public class CommunityLDA {
	
	public static void main(String args[]) throws IOException {

		String dataset = "toy";
		String foldername = "input/";
		String indir = "data/"+dataset+"/"+foldername;		
		String outdir = "data/"+dataset;
		String tempOutDir = outdir + "/output/";			//output folder
		String tRefIxpath = indir +"tRefIx.txt";		//tRefIx list
		String uRefIxpath = indir +"uRefIx.txt";		//tRefIx list
		String retweetpath = indir+"retweet.txt";
		
		(new File(tempOutDir)).mkdirs();

		int wordMapSize = 9;
		int tsMapSize = 7;		

		float beta = 0.1f;		
		int iteration = 200;
		int saveStep = 200;
		int saveTimes = 1;
		
		int A = Integer.parseInt(args[0]);
		int T = 3;
		
		tempOutDir += ("C"+A+ "/");
		(new File(tempOutDir)).mkdirs();
		
		//read user papers from file
		// the array to store all users
		ArrayList<String> bodies = new ArrayList<String>();
		FileUtil.readLines(indir+"alltweets.txt", bodies);
		int U = 8;
		User[] users = new User[U];
		for (int i = 0; i < U; i++) {
//				System.out.println("Reading user "+i);
				User user = new User(i, bodies.get(i));
				users[i] = user;
			}
		bodies.clear();	
		
		//user-user graph
		TIntObjectHashMap<int[]> userNeighbour = new TIntObjectHashMap<int[]>();
		userNeighbour = FileUtil.readTIntObjectHashMapILI(indir+"uunetwork.txt", ",");	
		

		int[][] lgroups = ModelUtil.getLgroups(indir);
//		double[][] indfea = FileUtil.readLinesDMatrix(indir+"ufor.txt", ",");
		//load serial positive and negative instances feature
		float[][] tpf = FileUtil.readLinesFMatrix(indir+"tp.txt", "\t");
		float[][] tnf = FileUtil.readLinesFMatrix(indir+"tn.txt", "\t");
		int NLt = tnf.length;
		int[] nvs = FileUtil.readLinesI(indir+"negavs.txt", NLt);
		
		//read user list from file
		
		HashSet<String> retweet = new HashSet<String>();
		FileUtil.readLines(retweetpath, retweet);
		
		int[] luRefIx = FileUtil.readLinesI(uRefIxpath, U);
		int[] ltRefIx = FileUtil.readLinesI(tRefIxpath, U);

		String outputDir = tempOutDir;
		FileUtil.mkdir(new File(outputDir));

		Model model = new Model(T, A, beta, lgroups);
		int initer = 0;
//		model.init(users, wordMapSize, tsMapSize, userNeighbour, ltRefIx,luRefIx,retweet,indfea);	
		model.init(users, wordMapSize, tsMapSize, userNeighbour, ltRefIx,luRefIx,retweet,tpf, tnf,nvs);		
		model.inference(iteration, users, saveStep, saveTimes, outputDir,initer);
//		model.generateUserFeature(users, dindir);
		System.out.println("done");
		
	}

}
