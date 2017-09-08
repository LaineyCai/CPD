package cpdtwitter2;
/*
 * 2016-05-16 ** Clean version parallel, divide dataset by documents, for tweets dataset
 * For community profiling, sigmoid link function for user-user and tweet-tweet, polya-gamma data agumented & gibbs sampling
 * all clarified, but skip nu*f
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

		String dataset = /*"toy2"*/"twitter";
		String foldername = "input/";
		String indir = "data/"+dataset+"/"+foldername;		
		String outdir = "data/"+dataset;
		int train = Integer.parseInt(args[0]);
		String dindir = indir +"train"+train+"/";
		String tempOutDir = outdir + "/output/train"+train+"/";			//output folder
		String tRefIxpath = dindir +"tRefIx"+train+".txt";		//tRefIx list
		String uRefIxpath = dindir +"uRefIx"+train+".txt";		//tRefIx list
		String retweetpath = dindir+"retweet"+train+".txt";
		
		(new File(tempOutDir)).mkdirs();

		int wordMapSize = 2316020;
		int tsMapSize = 1633;		

		float beta = 0.1f;		
		int iteration = 1000;
		int saveStep = 1000;
		int saveTimes = 1;
		
		int A = Integer.parseInt(args[1]);
		int T = A;
		
		tempOutDir += ("C"+A+ "/");
		(new File(tempOutDir)).mkdirs();
		
		//read user papers from file
		// the array to store all users
		ArrayList<String> bodies = new ArrayList<String>();
		FileUtil.readLines(dindir+"alltweets_train"+train+".txt", bodies);
		int U = 137325;
		User[] users = new User[U];
		for (int i = 0; i < U; i++) {
//				System.out.println("Reading user "+i);
				User user = new User(i, bodies.get(i));
				users[i] = user;
			}
		bodies.clear();	
		
		//user-user graph
		TIntObjectHashMap<int[]> userNeighbour = new TIntObjectHashMap<int[]>();
		userNeighbour = FileUtil.readTIntObjectHashMapILI(dindir+"uunetwork_train"+train+".txt", ",");	
		

		int[][] lgroups = ModelUtil.getLgroups(indir);
		
		//read user list from file
		
		HashSet<String> retweet = new HashSet<String>();
		FileUtil.readLines(retweetpath, retweet);
		
		int[] luRefIx = FileUtil.readLinesI(uRefIxpath, U);
		int[] ltRefIx = FileUtil.readLinesI(tRefIxpath, U);

		String outputDir = tempOutDir;
		FileUtil.mkdir(new File(outputDir));

		Model model = new Model(T, A, beta, lgroups);
		int initer = 0;
		model.init(users, wordMapSize, tsMapSize, userNeighbour, ltRefIx,luRefIx,retweet);
		model.inference(iteration, users, saveStep, saveTimes, outputDir,initer);

		System.out.println("done");
		
	}

}
