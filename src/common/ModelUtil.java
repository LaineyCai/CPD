package common;

import gnu.trove.map.hash.TIntObjectHashMap;

import java.util.ArrayList;
import java.util.List;


public class ModelUtil {
	public static int[][] getLgroups(String path)
	{
		ArrayList<String> lines = new ArrayList<String>();
		FileUtil.readLines(path+"ugroups.txt", lines);
		
		int N = lines.size();
		int[][] lgroups = new int[N][];
		
		for(int i=0; i<N; i++)
		{
			String[] users = lines.get(i).split(",");
			lgroups[i] = new int[users.length];
			for(int j=0; j<users.length; j++)
			{
				lgroups[i][j] = Integer.parseInt(users[j]);
			}
		}
		
		return lgroups;
	}
	
	public static int[][] getLgroups(String path, String file)
	{
		ArrayList<String> lines = new ArrayList<String>();
		FileUtil.readLines(path+file, lines);
		
		int N = lines.size();
		int[][] lgroups = new int[N][];
		
		for(int i=0; i<N; i++)
		{
			String[] users = lines.get(i).split(",");
			lgroups[i] = new int[users.length];
			for(int j=0; j<users.length; j++)
			{
				lgroups[i][j] = Integer.parseInt(users[j]);
			}
		}
		
		return lgroups;
	}
	
	public static TIntObjectHashMap<int[][]> generateuserInNeighbourT(TIntObjectHashMap<int[]> userNeighbour, int U)
	{
		TIntObjectHashMap<List<int[]>> tempuserInNeighbours =  new TIntObjectHashMap<List<int[]>>();
		TIntObjectHashMap<int[][]> userInNeighbours =  new TIntObjectHashMap<int[][]>();
		
		int uPIx = 0;
		for ( int u = 0; u < U; u++) {
			if(userNeighbour.get(u)!=null)
			{
				for ( int v=0; v<userNeighbour.get(u).length; v++ ) {
					int uNeighbor = userNeighbour.get(u).length;
					int[] pair = new int[]{u,uPIx};
					
					if(tempuserInNeighbours.containsKey(uNeighbor))
					{					
						tempuserInNeighbours.get(uNeighbor).add(pair);
					}
					else
					{
						List<int[]> followees = new ArrayList<int[]>();
						followees.add(pair);
						tempuserInNeighbours.put(uNeighbor, followees);
					}
					uPIx ++;
				}
			}
		}
		
		for(int u: tempuserInNeighbours.keys())
		{
			List<int[]> lpairs = tempuserInNeighbours.get(u);
			int[][] pairs = new int[lpairs.size()][2];
			for(int i=0; i<lpairs.size(); i++)
			{
				pairs[i] = lpairs.get(i);
			}
			userInNeighbours.put(u, pairs);
		}
		
		return userInNeighbours;
	}

}
