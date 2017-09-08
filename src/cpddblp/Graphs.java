package cpddblp;

import gnu.trove.map.hash.TIntObjectHashMap;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import common.ComUtil;
import common.FileUtil;

public class Graphs {
	
	//return the number of positive links to decare lambda
	public HashMap<Integer, List<Integer>> readUserGraph(String UserNetworkFileName, String userListFileName)
	{
		HashMap<Integer, List<Integer>> userNeighbour = new HashMap<Integer, List<Integer>>();	
		
		ArrayList<String> usernetworks = new ArrayList<String>();
		FileUtil.readLines(UserNetworkFileName, usernetworks);
		
		ArrayList<String> users = new ArrayList<String>();
		FileUtil.readLines(userListFileName, users);
		
		for(String pair: usernetworks)
		{
			String[] userpair = pair.split(",");
			int userid1 = users.indexOf(userpair[0]);
			int userid2 = users.indexOf(userpair[1]);
			
			if(userNeighbour.containsKey(userid1))
			{
				userNeighbour.get(userid1).add(userid2);	
			}
			else
			{
				List<Integer> followers = new ArrayList<Integer>();
				followers.add(userid2);
				userNeighbour.put(userid1, followers);
			}
		}
		
		return userNeighbour;
	}
	
	public HashMap<Integer, List<Integer>> writeUserGraph(String UserNetworkFileName, String userListFileName)
	{
		HashMap<Integer, List<Integer>> userNeighbour = new HashMap<Integer, List<Integer>>();	
		
		ArrayList<String> usernetworks = new ArrayList<String>();
		FileUtil.readLines(UserNetworkFileName, usernetworks);
		
		ArrayList<String> users = new ArrayList<String>();
		FileUtil.readLines(userListFileName, users);
		
		for(String pair: usernetworks)
		{
			String[] userpair = pair.split(",");
			int userid1 = users.indexOf(userpair[0]);
			int userid2 = users.indexOf(userpair[1]);
			
			if(userNeighbour.containsKey(userid1))
			{
				userNeighbour.get(userid1).add(userid2);	
			}
			else
			{
				List<Integer> followers = new ArrayList<Integer>();
				followers.add(userid2);
				userNeighbour.put(userid1, followers);
			}
		}
		
		return userNeighbour;
	}
	
	public HashMap<Integer, List<int[]>> generateuserInNeighbour(HashMap<Integer, List<Integer>> userNeighbour, int U)
	{
		HashMap<Integer, List<int[]>> userInNeighbours =  new HashMap<Integer, List<int[]>>();
		
		int uPIx = 0;
		for ( int u = 0; u < U; u++) {

			if(userNeighbour.get(u)!=null)
			{
				for ( int v=0; v<userNeighbour.get(u).size(); v++ ) {
					int uNeighbor = userNeighbour.get(u).get(v);
					int[] pair = new int[]{u,uPIx};
					
					if(userInNeighbours.containsKey(uNeighbor))
					{					
						userInNeighbours.get(uNeighbor).add(pair);
					}
					else
					{
						List<int[]> followees = new ArrayList<int[]>();
						followees.add(pair);
						userInNeighbours.put(uNeighbor, followees);
					}
					uPIx ++;
				}
			}
		}
		return userInNeighbours;
	}
	
	public HashMap<String, List<int[]>> generatetweetInNeighbour(ArrayList<User> users)
	{
		HashMap<String, List<int[]>> tweetInNeighbours = new HashMap<String, List<int[]>>();
		int tPIx = 0;
		int U = users.size();
		for(int u=0; u<U; u++)
		{
			TIntObjectHashMap<int[][]>  retweetinfo = users.get(u).getRetweets();
			if(retweetinfo!=null)
			{
				for(int utIx: retweetinfo.keys())
				{
					int[][] p = retweetinfo.get(utIx);
					for(int i=0; i<p.length; i++)
					{
						int v =p[i][0];	//the user where the tweet is from
						int vtIx = p[i][1];	//the tweet index of user v
						
						int[] vpair = new int[]{v, vtIx};
						int[] upair = new int[]{u, utIx, tPIx};
						String svpair = ComUtil.generateKey(vpair);
						if(tweetInNeighbours.containsKey(svpair))
						{
							tweetInNeighbours.get(svpair).add(upair);
						}
						else
						{
							List<int[]> beretweeted = new ArrayList<int[]>();
							beretweeted.add(upair);
							tweetInNeighbours.put(svpair, beretweeted);
						}
						
						tPIx ++;
					}
				}
			}
			
		}
		return tweetInNeighbours;
	}
	
	public HashMap<String, List<int[]>> generatetweetInNeighbour(User[] users)
	{
		HashMap<String, List<int[]>> tweetInNeighbours = new HashMap<String, List<int[]>>();
		int tPIx = 0;
		int U = users.length;
		for(int u=0; u<U; u++)
		{
			TIntObjectHashMap<int[][]>  retweetinfo = users[u].getRetweets();
			if(retweetinfo!=null)
			{
				for(int utIx: retweetinfo.keys())
				{
					int[][] p = retweetinfo.get(utIx);
					for(int i=0; i<p.length; i++)
					{
						int v =p[i][0];	//the user where the tweet is from
						int vtIx = p[i][1];	//the tweet index of user v
						
						int[] vpair = new int[]{v, vtIx};
						int[] upair = new int[]{u, utIx, tPIx};
						String svpair = ComUtil.generateKey(vpair);
						if(tweetInNeighbours.containsKey(svpair))
						{
							tweetInNeighbours.get(svpair).add(upair);
						}
						else
						{
							List<int[]> beretweeted = new ArrayList<int[]>();
							beretweeted.add(upair);
							tweetInNeighbours.put(svpair, beretweeted);
						}
						
						tPIx ++;
					}
				}
			}
			
		}
		return tweetInNeighbours;
	}
	
	public void generatePositiveRetweet(ArrayList<User> users, String OutPath)
	{
		ArrayList<String> pinstance = new ArrayList<String>();
		ArrayList<String> ppair = new ArrayList<String>();
		int U = users.size();
		for(int u=0; u<U; u++)
		{
			TIntObjectHashMap<int[][]>  retweetinfo = users.get(u).getRetweets();
			if(retweetinfo!=null)
			{
				for(int utIx: retweetinfo.keys())
				{	
					int[][] p = retweetinfo.get(utIx);
					for(int i=0; i<p.length; i++)
					{
						int v = p[i][0];	//the user where the tweet is from
						int t = users.get(u).getDocTimeStamp()[utIx];
						String ins = u+"\t"+utIx+"\t" + v+"\t"+ t;
						String par = u+"\t" + v;
						pinstance.add(ins);
						ppair.add(par);
					}
				}
			}			
		}
		
		FileUtil.writeLines(OutPath+"PInstance.txt", pinstance);
		FileUtil.writeLines(OutPath+"PPair.txt", ppair);		
	}
	
	//use #follower, #tweets
	public double[][] readIndividualFeature_back(String featureFileName, ArrayList<String> users, int N)
	{
		int U = users.size();
		double[][] indfea = new double[U][N];
		ArrayList<String> features = new ArrayList<String>();
		FileUtil.readLines(featureFileName, features);
		double max_c1 = 0, max_c2=0;
		for(String f: features)
		{
			String[] terms = f.split("\t");
			/////!!!!!!!!! need to refine how to define popularity and activity of a user
			//popularity, follower/followee
			int u = users.indexOf(terms[0]);
/*			indfea[u][0] = Double.parseDouble(terms[3])/Double.parseDouble(terms[2]);
			
			//activity
			indfea[u][1] = Double.parseDouble(terms[4])/Double.parseDouble(terms[2]);*/
			indfea[u][0] = Double.parseDouble(terms[3]);
			indfea[u][1] = Double.parseDouble(terms[4]);
			if(indfea[u][0]>max_c1)
				max_c1 = indfea[u][0];
			if(indfea[u][1]>max_c2)
				max_c2 = indfea[u][1];
		}
		for(int i=0; i<indfea.length; i++)
		{
			indfea[i][0]/=max_c1;
			indfea[i][1]/=max_c2;
		}
		return indfea;
	}
	//read formated features
		public double[][] readIndFea(String featureFileName, int U, int N)
		{
			double[][] indfea = new double[U][N];
			BufferedReader reader = null;

			int i=0;
			try {
				reader = new BufferedReader(new FileReader(new File(featureFileName)));

				String line = null;
				while ((line = reader.readLine()) != null) {
					String[] terms = line.split(",");
					for(int k=0; k<N; k++)
					{
						indfea[i][k] = Double.parseDouble(terms[k]);
					}
					i++;
				}

			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			} finally {
				if (reader != null) {
					try {
						reader.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}
			return indfea;
		}
	//user follower/friends
	public int[][] readIndividualFeature(String featureFileName, ArrayList<String> users, int N, String outfile)
	{
		int U = users.size();
		int[][] indfea = new int[U][N];
		ArrayList<String> features = new ArrayList<String>();
		FileUtil.readLines(featureFileName, features);
		//double max_c1 = 0;
		int i=0;
		for(String f: features)
		{
			String[] terms = f.split("\t");
			/////!!!!!!!!! need to refine how to define popularity and activity of a user
			//popularity, follower/followee
			int u = users.indexOf(terms[0]);
			indfea[u][0] = Integer.parseInt(terms[2]);
			indfea[u][1] = Integer.parseInt(terms[3]);
			System.out.println(i++);
			//activity
//			indfea[u][1] = Double.parseDouble(terms[4])/Double.parseDouble(terms[2]);
//			indfea[u][0] = Double.parseDouble(terms[3]);
//			indfea[u][1] = Double.parseDouble(terms[4]);
//			if(indfea[u][0]>max_c1)
//				max_c1 = indfea[u][0];
//			if(indfea[u][1]>max_c2)
//				max_c2 = indfea[u][1];
		}
//		for(int i=0; i<indfea.length; i++)
//		{
//			indfea[i][0]/=max_c1;
//			System.out.print(indfea[i][0]+" ");
////			indfea[i][1]/=max_c2;
//		}
//		System.out.println();
		FileUtil.writeLines(outfile, indfea);
		return indfea;
	}
	
	public double[][] readIndividualFeature_back2(String featureFileName, ArrayList<String> users, int N)
	{
		int U = users.size();
		double[][] indfea = new double[U][N];
		ArrayList<String> features = new ArrayList<String>();
		FileUtil.readLines(featureFileName, features);
		double max_c1 = 0;
		for(String f: features)
		{
			String[] terms = f.split("\t");
			/////!!!!!!!!! need to refine how to define popularity and activity of a user
			//popularity, follower/followee
			int u = users.indexOf(terms[0]);
			indfea[u][0] = Double.parseDouble(terms[3])/Double.parseDouble(terms[2]);
			
			//activity
//			indfea[u][1] = Double.parseDouble(terms[4])/Double.parseDouble(terms[2]);
//			indfea[u][0] = Double.parseDouble(terms[3]);
//			indfea[u][1] = Double.parseDouble(terms[4]);
			if(indfea[u][0]>max_c1)
				max_c1 = indfea[u][0];
//			if(indfea[u][1]>max_c2)
//				max_c2 = indfea[u][1];
		}
		for(int i=0; i<indfea.length; i++)
		{
			indfea[i][0]/=max_c1;
			System.out.print(indfea[i][0]+" ");
//			indfea[i][1]/=max_c2;
		}
		System.out.println();
		return indfea;
	}
	
}
