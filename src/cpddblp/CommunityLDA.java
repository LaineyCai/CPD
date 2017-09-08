package cpddblp;
/*
 * 2015-12-14 ** Clean version parallel, divide dataset by documents
 * For community profiling, sigmoid link function for user-user and tweet-tweet, polya-gamma data agumented & gibbs sampling
 */
import gnu.trove.map.hash.TIntObjectHashMap;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import common.FileUtil;
import common.ModelUtil;




public class CommunityLDA {

	public static void main(String args[]) throws IOException {

		String dataset = /*"toy2"*/"dblp";
		String foldername = "input/";
		String indir = "data/"+dataset+"/"+foldername;		
		String outdir = "data/"+dataset;
		int train = Integer.parseInt(args[0]);
		String dindir = indir +"train"+train+"/";
		String tempOutDir = outdir + "/output/train"+train+"/";			//output folder
		String tRefIxpath = dindir +"tRefIx"+train+".txt";		//tRefIx list
		String uRefIxpath = dindir +"uRefIx"+train+".txt";		//tRefIx list
		String noretweetpath = dindir+"noretweet"+train+".txt";
		
		(new File(tempOutDir)).mkdirs();

		int wordMapSize = 330334;
		int tsMapSize = 77;		

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
		FileUtil.readLines(dindir+"papers_train"+train+".txt", bodies);
		int U = 916907;
		User[] users = new User[U];
		for (int i = 0; i < U; i++) {
				System.out.println("Reading user "+i);
				User user = new User(i, bodies.get(i));
				users[i] = user;
			}
		bodies.clear();		
		
		int[][] lgroups = ModelUtil.getLgroups(indir);
//		double[][] indfea = FileUtil.readLinesDMatrix(indir+"ufor.txt", ",");
		
		float[][] tpf = FileUtil.readLinesFMatrix(dindir+"tp.txt", "\t");
		float[][] tnf = FileUtil.readLinesFMatrix(dindir+"tn.txt", "\t");
		int NLt = tnf.length;
		int[] nvs = FileUtil.readLinesI(dindir+"negavs.txt", NLt);
		
		//read user list from file
		
		HashSet<String> noretweet = new HashSet<String>();
		FileUtil.readLines(noretweetpath, noretweet);
		
		int[] luRefIx = FileUtil.readLinesI(uRefIxpath, U);
		int[] ltRefIx = FileUtil.readLinesI(tRefIxpath, U);
		
		//user-user graph
		TIntObjectHashMap<int[]> userNeighbour = new TIntObjectHashMap<int[]>();
		userNeighbour = FileUtil.readTIntObjectHashMapILI(dindir+"uunetwork_train"+train+".txt", ",");

		String outputDir = tempOutDir;
		FileUtil.mkdir(new File(outputDir));

		Model model = new Model(T, A, beta, lgroups);
		int initer = 0;
//		model.init(users, wordMapSize, tsMapSize, userNeighbour, ltRefIx,luRefIx,noretweet, indfea);
		model.init(users, wordMapSize, tsMapSize, userNeighbour, ltRefIx,luRefIx,noretweet, tpf, tnf,nvs);
		model.inference(iteration, users, saveStep, saveTimes, outputDir,initer);
//		model.generateUserFeature(users, dindir);
		System.out.println("done");
		
	}

}
