package common;

import gnu.trove.map.hash.TIntObjectHashMap;

import java.util.*;
import java.io.*;

import common.FileUtil;

public class FileUtil {
	
	/*public static void WriteMap(String file, HashMap<Integer, List<Integer>> m)
	{
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));

			Set<?> s = hashMap.entrySet();
			Iterator<?> it = s.iterator();
			while (it.hasNext()) {
				Map.Entry m = (Map.Entry) it.next();
				writer.write(m.getKey() + "\t" + m.getValue() + "\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}*/
	
	public static int[][] readLines(String file, int N) {
		int[][] lines = new int[N][];
		BufferedReader reader = null;

		int count = 0;
		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				int M = terms.length;
				lines[count] = new int[M];
				for(int m=0; m<M; m++)
				{
					lines[count][m] = Integer.parseInt(terms[m]);
				}
				count++;
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
		return lines;

	}
	
	public static HashMap<String, List<String>> readHashMapSLS(String file) {
		BufferedReader reader = null;
		String line = "";
		HashMap<String, List<String>> m =new HashMap<String, List<String>>();
		try {
			reader = new BufferedReader(new FileReader(new File(file)));
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				List<String> l = new ArrayList<String>();
				for(String t: terms[1].split(","))
				{
					l.add(t);
				}
				m.put(terms[0],l);
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
		return m;
	}
	
	public static void writeHashMapSS(String file, HashMap<?, ?> hashMap) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));

			Set<?> s = hashMap.entrySet();
			Iterator<?> it = s.iterator();
			while (it.hasNext()) {
				Map.Entry m = (Map.Entry) it.next();
				writer.write(m.getKey() + "," + m.getValue() + "\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void writeHashMapSS(String file, HashMap<?, ?> hashMap, String delim) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));

			Set<?> s = hashMap.entrySet();
			Iterator<?> it = s.iterator();
			while (it.hasNext()) {
				Map.Entry m = (Map.Entry) it.next();
				writer.write(m.getKey() + delim + m.getValue() + "\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void writeHashMapTIntHI(String file, TIntObjectHashMap<HashSet<Integer>> hm, String delim) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));
			
			for(int k: hm.keys())
			{
				String line = k+"\t";
				for(int u: hm.get(k))
				{
					line+=(u+delim);
				}
				writer.write(line.substring(0, line.length()-delim.length()) + "\n");
			}			
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static HashMap<String, String> readHashMapSS1(String file) {
		BufferedReader reader = null;
		String line = "";
		HashMap<String, String> m =new HashMap<String, String>();
		try {
			reader = new BufferedReader(new FileReader(new File(file)));
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split(",");
				m.put(terms[0],terms[1]);
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
		return m;
	}
	
	public static HashMap<String, Double> readHashMapSD(String file) {
		BufferedReader reader = null;
		String line = "";
		HashMap<String, Double> m =new HashMap<String, Double>();
		try {
			reader = new BufferedReader(new FileReader(new File(file)));
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				m.put(terms[0],Double.parseDouble(terms[1]));
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
		return m;
	}
	
	public static TIntObjectHashMap<int[]> readTIntObjectHashMapILI(String file, String delim) {
		BufferedReader reader = null;
		String line = "";
		TIntObjectHashMap<int[]> m =new TIntObjectHashMap<int[]>();
		try {
			reader = new BufferedReader(new FileReader(new File(file)));
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				String[] ls = terms[1].split(delim);
				int[] als = new int[ls.length];
				for(int i=0; i<ls.length; i++)
				{
					als[i] = Integer.parseInt(ls[i]);
				}
				m.put(Integer.parseInt(terms[0]), als);
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
		return m;
	}
	
	public static int[] readLinesI(String file, int Cap) {
		BufferedReader reader = null;
		int[] lines = new int[Cap];
		int i = 0;
		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null && i<Cap) {
				lines[i++] = Integer.parseInt(line);
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
		return lines;
	}

	public static void readLines(String file, ArrayList<String> lines) {
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				lines.add(line);
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

	}
	
	public static int[][] readLinesAAINT(String file) {
		ArrayList<String> lines = new ArrayList<String>();
		readLines(file, lines);
		int[][] res = new int[lines.size()][];
		for(int i=0; i<lines.size(); i++)
		{
			String[] terms = lines.get(i).split(",");
			res[i] = new int[terms.length];
			for(int j=0; j<terms.length; j++)
			{
				res[i][j] = Integer.parseInt(terms[j]);
			}
		}
		return res;
	}
	
	public static double[][] readLinesDMatrix(String file, String sep) {
		ArrayList<String> lines = new ArrayList<String>();
		readLines(file, lines);
		double[][] res = new double[lines.size()][];
		for(int i=0; i<lines.size(); i++)
		{
			String[] terms = lines.get(i).split(sep);
			res[i] = new double[terms.length];
			for(int j=0; j<terms.length; j++)
			{
				res[i][j] = Double.parseDouble(terms[j]);
			}
		}
		return res;
	}
	
	public static float[][] readLinesFMatrix(String file, String sep) {
		ArrayList<String> lines = new ArrayList<String>();
		readLines(file, lines);
		float[][] res = new float[lines.size()][];
		int M = lines.size();
		int N = lines.get(0).split(sep).length;
		for(int i=0; i<M; i++)
		{
			String[] terms = lines.get(i).split(sep);
			res[i] = new float[N];
			for(int j=0; j<N; j++)
			{
				res[i][j] = Float.parseFloat(terms[j]);
			}
		}
		return res;
	}
	
	public static ArrayList<HashSet<String>> readLinesAHS(String file) {
		ArrayList<String> lines = new ArrayList<String>();
		readLines(file, lines);
		ArrayList<HashSet<String>> res = new ArrayList<HashSet<String>>();
		for(String line: lines)
		{
			String[] terms = line.split(",");
			HashSet<String> hs = new HashSet<String>();
			for(String term: terms)
			{
				hs.add(term);
			}
			res.add(hs);
		}
		return res;
	}
	
	public static int[][] readUserAction(String userActionFile, int U)
	{
		int[][] userActions = new int[U][];
		BufferedReader reader = null;

		int i=0;
		try {
			reader = new BufferedReader(new FileReader(new File(userActionFile)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split(",");
				int N = terms.length;
				userActions[i] = new int[N];
				
				for(int k=0; k<N; k++)
				{
					userActions[i][k] = Integer.parseInt(terms[k]);
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
		return userActions;
	}
	
	public static void readLines2(String file, ArrayList<ArrayList<Integer>> lines) {
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				ArrayList<Integer> lu = new ArrayList<Integer>();
				for(String u: line.split(","))
				{
					lu.add(Integer.parseInt(u));
				}
				lines.add(lu);
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

	}
	
	public static void readLines2(String file, ArrayList<ArrayList<Integer>> lines, String delim) {
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				ArrayList<Integer> lu = new ArrayList<Integer>();
				for(String u: line.split(delim))
				{
					lu.add(Integer.parseInt(u));
				}
				lines.add(lu);
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

	}
	
	public static void readLines4(String file, ArrayList<ArrayList<String>> lines) {
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				ArrayList<String> lu = new ArrayList<String>();
				for(String u: line.split(","))
				{
					lu.add(u);
				}
				lines.add(lu);
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

	}
	
	public static void readLines3(String file, ArrayList<int[]> arint) {
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				int[] a = new int[terms.length];
				for(int i=0; i<terms.length; i++)
				{
					a[i] = Integer.parseInt(terms[i]);
				}
				arint.add(a);
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

	}
	
	public static void readLines1(String file, ArrayList<Integer> lines) {
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				lines.add(Integer.parseInt(line));
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

	}
	
	public static void readLines(String file, HashSet<String> lines) {
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				lines.add(line);
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

	}
	
	public static void readLinesI(String file, HashSet<Integer> lines) {
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				lines.add(Integer.parseInt(line));
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

	}
	
	public static void readTweetsBody(String file, ArrayList<String> lines) {
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				if(terms[0].equals("O"))
					lines.add(terms[2]);
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

	}
	
	public static boolean readTweetsBodyAll(String file, ArrayList<String> lines) {
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				if(terms[0].equals("O"))
					lines.add(terms[2]);
				else
					lines.add("not original");
			}
			return true;

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return false;
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
				return true;
			}
			else
			{
				return false;
			}
		}

	}
	public static void readUsers(String file, ArrayList<String> lines) {
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				lines.add(line.split("\t")[0]);
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

	}
	
	public static String readLines(String file) {
		BufferedReader reader = null;
		String line = "";
		try {
			reader = new BufferedReader(new FileReader(new File(file)));

		line = reader.readLine();


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
		return line;
	}
	
	public static void readContents(String file, ArrayList<String> lines) {
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				System.out.println(file+"\t"+line);
				lines.add(terms[2].replaceFirst("RT  ", ""));
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

	}
	
		public static void writeLines(String file, HashMap<?, ?> hashMap) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));

			Set<?> s = hashMap.entrySet();
			Iterator<?> it = s.iterator();
			while (it.hasNext()) {
				Map.Entry m = (Map.Entry) it.next();
				writer.write(m.getKey() + "\t" + m.getValue() + "\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void writeHashMap(String file, HashMap<?, ?> hashMap) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));

			Set<?> s = hashMap.entrySet();
			Iterator<?> it = s.iterator();
			while (it.hasNext()) {
				Map.Entry m = (Map.Entry) it.next();
				String v = "";
				for(Integer va: (List<Integer>)m.getValue())
				{
					v+=(va+",");
				}
				v = v.substring(0,v.length()-1);
				writer.write(m.getKey() + "\t" + v + "\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void writeHashMapSLS(String file, HashMap<?, ?> hashMap) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));

			Set<?> s = hashMap.entrySet();
			Iterator<?> it = s.iterator();
			while (it.hasNext()) {
				Map.Entry m = (Map.Entry) it.next();
				String v = "";
				for(String va: (List<String>)m.getValue())
				{
					v+=(va+",");
				}
				v = v.substring(0,v.length()-1);
				writer.write(m.getKey() + "\t" + v + "\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void writeHashMapSSS(String file, HashMap<?, ?> hashMap) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));

			Set<?> s = hashMap.entrySet();
			Iterator<?> it = s.iterator();
			while (it.hasNext()) {
				Map.Entry m = (Map.Entry) it.next();
				String v = "";
				for(String va: (HashSet<String>)m.getValue())
				{
					v+=(va+",");
				}
				v = v.substring(0,v.length()-1);
				writer.write(m.getKey() + "\t" + v + "\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static HashMap<Integer, List<Integer>> readHashMap(String file) {
		BufferedReader reader = null;
		String line = "";
		HashMap<Integer, List<Integer>> m =new HashMap<Integer, List<Integer>>();
		try {
			reader = new BufferedReader(new FileReader(new File(file)));
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				List<Integer> l = new ArrayList<Integer>();
				for(String t: terms[1].split(","))
				{
					l.add(Integer.parseInt(t));
				}
				m.put(Integer.parseInt(terms[0]),l);
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
		return m;
	}
	
		
	public static void writeHashMapSI(String file, HashMap<?, ?> hashMap) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));

			Set<?> s = hashMap.entrySet();
			Iterator<?> it = s.iterator();
			while (it.hasNext()) {
				Map.Entry m = (Map.Entry) it.next();
				writer.write(m.getKey() + "\t" + m.getValue() + "\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static HashMap<String, String> readHashMapSS(String file) {
		BufferedReader reader = null;
		HashMap<String, String> m =new HashMap<String, String>();
		String line = "";
		try {
			reader = new BufferedReader(new FileReader(new File(file)));
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				m.put(terms[0], terms[1]);
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
		return m;
	}
	
	public static HashMap<String, String> readReverseHashMapSS(String file) {
		BufferedReader reader = null;
		HashMap<String, String> m =new HashMap<String, String>();
		String line = "";
		try {
			reader = new BufferedReader(new FileReader(new File(file)));
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				m.put(terms[1], terms[0]);
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
		return m;
	}
	
	public static HashSet<String> readHMValueSet(String file) {
		BufferedReader reader = null;
		HashSet<String> m = new HashSet<String>();
		String line = "";
		try {
			reader = new BufferedReader(new FileReader(new File(file)));
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				m.add(terms[1]);
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
		return m;
	}
	
	public static HashMap<String, Integer> readHashMapSI(String file) {
		BufferedReader reader = null;
		HashMap<String, Integer> m =new HashMap<String, Integer>();
		String line = "";
		try {
			reader = new BufferedReader(new FileReader(new File(file)));
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				m.put(terms[0], Integer.parseInt(terms[1]));
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
		return m;
	}
	
	public static void writeLines(String file, int[][] a) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));
			for(int i=0; i<a.length; i++)
			{
				writer.write(a[i][0] + "\t" + a[i][1] + "\n");
			}				

		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void writeLines(String file, double[] a) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));
			for(int i=0; i<a.length; i++)
			{
				writer.write(a[i]+ "\n");
			}				

		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void writeLines(String file, int[] a) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));
			for(int i=0; i<a.length; i++)
			{
				writer.write(a[i]+ "\n");
			}				

		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void writeLines(String file, long[] a) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));
			for(int i=0; i<a.length; i++)
			{
				writer.write(a[i]+ "\n");
			}				

		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void writeLines(String file, String[] a) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));
			for(int i=0; i<a.length; i++)
			{
				writer.write(a[i]+ "\n");
			}				

		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void writeLines(String file, String line) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));
			writer.write(line + "\n");

		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void appendLines(String file, String line) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file),true));

			writer.append(line+"\n");
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	public static void writeLines(String file, ArrayList<?> counts) {
		BufferedWriter writer = null;

		try {

			writer = new BufferedWriter(new FileWriter(new File(file)));

			for (int i = 0; i < counts.size(); i++) {
				writer.write(counts.get(i) + "\n");
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

	}
	
	public static void writeLines(String file, HashSet<String> counts) {
		BufferedWriter writer = null;

		try {

			writer = new BufferedWriter(new FileWriter(new File(file)));

			for (String term: counts) {
				writer.write(term + "\n");
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

	}
	
	public static void writeLines1(String file, ArrayList<?> counts) {
		BufferedWriter writer = null;
		String line = "";
		try {

			writer = new BufferedWriter(new FileWriter(new File(file)));

			for (int i = 0; i < counts.size(); i++) {
				line += (counts.get(i) + ",");				
			}
			line = line.substring(0, line.length()-1);
			writer.write(line);

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

	}
	public static void writeLLI(String file, ArrayList<List<Integer>> lli) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));

			for (int i=0; i<lli.size(); i++)
			{
				List<Integer> lu = lli.get(i);
				String v = lu.toString();
				v = v.substring(1,v.length()-1);
				writer.write( v + "\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void writeLHS(String file, ArrayList<HashSet<String>> lli) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));

			for (int i=0; i<lli.size(); i++)
			{
				HashSet<String> lu = lli.get(i);
				String v = lu.toString();
				v = v.substring(1,v.length()-1);
				writer.write( v + "\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void writeLHI(String file, ArrayList<HashSet<Integer>> lli) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));

			for (int i=0; i<lli.size(); i++)
			{
				HashSet<Integer> lu = lli.get(i);
				String v = lu.toString();
				v = v.substring(1,v.length()-1);
				writer.write( v + "\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void writeMatrix(String file, int[][] m)
	{
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file),true));
			for (int i = 0; i < m.length; i++) {
				for(int j=0; j<m[i].length;j++)
				{
					writer.write(m[i][j]+"\t");
				}
				writer.write("\n");
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}		
	}
	
	public static void writeDMatrix(String file, double[][] m)
	{
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file),true));
			for (int i = 0; i < m.length; i++) {
				for(int j=0; j<m[i].length;j++)
				{
					writer.write(m[i][j]+"\t");
				}
				writer.write("\n");
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}		
	}
	
	public static void writeFMatrix(String file, float[][] m)
	{
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file),false));
			for (int i = 0; i < m.length; i++) {
				for(int j=0; j<m[i].length;j++)
				{
					writer.write(m[i][j]+"\t");
				}
				writer.write("\n");
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}		
	}
	
	public static void writeMatrix(String file, double[] m)
	{
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file),true));
			for (int i = 0; i < m.length; i++) {
				writer.write(m[i]+"\n");

			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}		
	}
	
	public static void writeMatrix(String file, float[][] m)
	{
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file),false));
			for (int i = 0; i < m.length; i++) {
				for(int j=0; j<m[i].length;j++)
				{
					writer.write(m[i][j]+"\t");
				}
				writer.write("\n");
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}		
	}
	
	public static void writeMatrix(String file, double[][] m)
	{
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file),false));
			for (int i = 0; i < m.length; i++) {
				for(int j=0; j<m[i].length;j++)
				{
					writer.write(m[i][j]+"\t");
				}
				writer.write("\n");
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}		
	}
	
	public static void appendMatrix(String file, double[][] m)
	{
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file),true));
			for (int i = 0; i < m.length; i++) {
				for(int j=0; j<m[i].length;j++)
				{
					writer.write(m[i][j]+"\t");
				}
				writer.write("\n");
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}		
	}
	
	public static void readMatrix(String file, boolean[][]m)
	{
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));
			String line = null;
			int linenb = 0;
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				for(int i=0; i<terms.length;i++)
				{
					m[linenb][i] = Boolean.parseBoolean(terms[i]);
				}
				linenb++;
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
	}
	
	public static void readMatrix(String file, double[][]m)
	{
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));
			String line = null;
			int linenb = 0;
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				for(int i=0; i<terms.length;i++)
				{
					m[linenb][i] = Double.parseDouble(terms[i]);
				}
				linenb++;
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
	}
	
	public static void readMatrix(String file, double[] m)
	{
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));
			String line = null;
			int linenb = 0;
			while ((line = reader.readLine()) != null) {
				m[linenb] = Double.parseDouble(line);
				linenb++;
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
	}
	
	public static void readMatrix(String file, float[] m)
	{
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));
			String line = null;
			int linenb = 0;
			while ((line = reader.readLine()) != null) {
				m[linenb] = (float) Double.parseDouble(line);
				linenb++;
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
	}
	
	public static void readMatrix1(String file, float[] m)
	{
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));
			String line = reader.readLine();
			String[] terms = line.split(";");
			for(int i=0; i<terms.length; i++)
			{
				m[i] = Float.parseFloat(terms[i]) ;
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
	}
	
	public static void readDBLPBody(String file, HashMap<Integer, ArrayList<String>> body, ArrayList<String> userlist)
	{
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				String newline = "";
				for(int i=1; i<terms.length; i++)
				{
					newline +=(terms[i]+"\t");
				}
				newline = newline.substring(0,newline.length()-1);
				int id = userlist.indexOf(terms[0]);
				if(body.containsKey(id))
				{					
					body.get(id).add(newline);
				}
				else
				{
					ArrayList<String> lb = new ArrayList<String>();
					lb.add(newline);
					body.put(id, lb);
				}
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
	}
	
	//change all the author id to the index in the userlist, and write all papers belonging to one author to one line for speed up the reading process
	public static void rewriteDBLPBody(String path, String file,  ArrayList<String> userlist, String fileout)
	{
		HashMap<Integer, ArrayList<String>> body = new HashMap<Integer, ArrayList<String>>();
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(path+file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				String newline = "";
				for(int i=1; i<=2; i++)
				{
					newline +=(terms[i]+"\t");
				}
				for(int i=3; i<terms.length; i+=2)
				{
					newline +=(userlist.indexOf(terms[i])+"\t"+terms[i+1]+"\t");
				}
				newline = newline.substring(0,newline.length()-1);
				int id = userlist.indexOf(terms[0]);
				if(body.containsKey(id))
				{					
					body.get(id).add(newline);
				}
				else
				{
					ArrayList<String> lb = new ArrayList<String>();
					lb.add(newline);
					body.put(id, lb);
				}
			}
			
			ArrayList<String> alluserpapers = new ArrayList<String>();
			for(int i=0; i<userlist.size(); i++)
			{
				System.out.println(i);
				String uline = "";
				for(String s: body.get(i))
				{
					uline +=(s+";");
				}
				uline = uline.substring(0, uline.length()-1);
				alluserpapers.add(uline);
			}
			FileUtil.writeLines(path+fileout, alluserpapers);

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
	}
	
	public static void writeMatrix(String file, boolean[][] m)
	{
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));
			for (int i = 0; i < m.length; i++) {
				for(int j=0; j<m[i].length;j++)
				{
					writer.write(m[i][j]+"\t");
				}
				writer.write("\n");
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}		
	}
	
	public static void readMatrix(String file, int[][]m)
	{
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));
			String line = null;
			int linenb = 0;
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				m[linenb] = new int[terms.length];
				for(int i=0; i<terms.length;i++)
				{
					m[linenb][i] = Integer.parseInt(terms[i]);
				}
				linenb++;
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
	}
	
	
	
	public static void readMatrix(String file, float[][]m)
	{
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));
			String line = null;
			int linenb = 0;
			int N= m.length;
			while ((line = reader.readLine()) != null && linenb<N) {
				String[] terms = line.split("\t");
				for(int i=0; i<terms.length;i++)
				{
					m[linenb][i] = (float) Double.parseDouble(terms[i]);
				}
				linenb++;
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
	}
	
	public static void readMatrix(String file, float[][]m, String delim)
	{
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));
			String line = null;
			int linenb = 0;
			int N= m.length;
			while ((line = reader.readLine()) != null && linenb<N) {
				String[] terms = line.split(delim);
				for(int i=0; i<terms.length;i++)
				{
					m[linenb][i] = (float) Double.parseDouble(terms[i]);
				}
				linenb++;
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
	}
	public static void writeMatrix(String file, boolean[][][] m)
	{
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));
			for (int i = 0; i < m.length; i++) {
				for(int j=0; j<m[i].length;j++)
				{
					for(int k=0; k<m[i][j].length; k++)
					{
						writer.write(m[i][j][k]+"\t");
					}
					writer.write("\n");
				}
				writer.write("\n");
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}		
	}
	
	public static void writeMatrix(String file, double[][][] m)
	{
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file),false));
			for (int i = 0; i < m.length; i++) {
				for(int j=0; j<m[i].length;j++)
				{
					for(int k=0; k<m[i][j].length; k++)
					{
						writer.write(m[i][j][k]+"\t");
					}
					writer.write("\n");
				}
				writer.write("\n");
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}		
	}
	
	public static void readMatrix(String file, boolean[][][] m)
	{
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));
			String line = null;
			int linei = 0, linej=0;
			while ((line = reader.readLine()) != null) {
				if(line.equals(""))
				{
					linei++;
					linej=0;
				}
				else
				{
					String[] terms = line.split("\t");
					for(int i=0; i<terms.length;i++)
					{
						m[linei][linej][i] = Boolean.parseBoolean(terms[i]);
					}
					linej++;
				}				
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
	}
	
	public static void writeMatrix(String file, int[][][] m)
	{
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));
			for (int i = 0; i < m.length; i++) {
				for(int j=0; j<m[i].length;j++)
				{
					for(int k=0; k<m[i][j].length; k++)
					{
						writer.write(m[i][j][k]+"\t");
					}
					writer.write("\n");
				}
				writer.write("\n");
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}		
	}
	
	public static void readMatrix(String file, int[][][] m)
	{
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));
			String line = null;
			int linei = 0, linej=0;
			while ((line = reader.readLine()) != null) {
				if(line.equals("null"))
				{}
				else if(line.equals(""))
				{
					linei++;
					linej=0;
				}
				else
				{
					String[] terms = line.split("\t");
					for(int i=0; i<terms.length;i++)
					{
						m[linei][linej][i] = Integer.parseInt(terms[i]);
					}
					linej++;
				}				
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
	}
	
	public static void readMatrix(String file, float[][][] m)
	{
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));
			String line = null;
			int linei = 0, linej=0;
			while ((line = reader.readLine()) != null) {
				if(line.equals(""))
				{
					linei++;
					linej=0;
				}
				else
				{
					String[] terms = line.split("\t");
					for(int i=0; i<terms.length;i++)
					{
						m[linei][linej][i] = (float) Double.parseDouble(terms[i]);
					}
					linej++;
				}				
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
	}
	
	public static void readMatrix(String file, double[][][] m)
	{
		BufferedReader reader = null;

		try {

			reader = new BufferedReader(new FileReader(new File(file)));
			String line = null;
			int linei = 0, linej=0;
			while ((line = reader.readLine()) != null) {
				if(line.equals(""))
				{
					linei++;
					linej=0;
				}
				else
				{
					String[] terms = line.split("\t");
					for(int i=0; i<terms.length;i++)
					{
						m[linei][linej][i] = Double.parseDouble(terms[i]);
					}
					linej++;
				}				
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
	}
	
	public static void writeTwoArrays(String file, ArrayList<?> a1, ArrayList<?> a2) {
		BufferedWriter writer = null;

		try {

			writer = new BufferedWriter(new FileWriter(new File(file)));

			for (int i = 0; i < a1.size(); i++) {
				writer.write(a1.get(i) +"\t"+a2.get(i)+ "\n");
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

	}
	
	public static void writeLinesAppend(String file, ArrayList<?> counts) {
		BufferedWriter writer = null;

		try {

			writer = new BufferedWriter(new FileWriter(new File(file),true));

			for (int i = 0; i < counts.size(); i++) {
				writer.write(counts.get(i) + "\n");
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

	}

	public static void writeLines(String file, ArrayList<String> uniWordMap,
			ArrayList<Integer> uniWordMapCounts) {
		BufferedWriter writer = null;

		try {

			writer = new BufferedWriter(new FileWriter(new File(file)));

			for (int i = 0; i < uniWordMap.size()
					|| i < uniWordMapCounts.size(); i++) {
				writer.write(uniWordMap.get(i) + "\t" + uniWordMapCounts.get(i)
						+ "\n");
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	@SuppressWarnings("unchecked")
	public static void writeLinesSorted(String file, ArrayList<?> uniWordMap,
			ArrayList<?> uniWordMapCounts, int flag) {
		// flag = 0 decreasing order otherwise increasing
		HashMap map = new HashMap();
		if (uniWordMap.size() != uniWordMapCounts.size()) {
			System.err.println("Array sizes are not equal!!! Function returned.");
		} else {
			for (int i = 0; i < uniWordMap.size(); i++) {
				map.put(uniWordMap.get(i), uniWordMapCounts.get(i));
			}
			map = (HashMap<String, Integer>) ComUtil.sortByValue(map, flag);
			writeLines(file, map);
			map.clear();
		}
	}

	public static void tokenize(String line, ArrayList<String> tokens) {
		StringTokenizer strTok = new StringTokenizer(line);
		while (strTok.hasMoreTokens()) {
			String token = strTok.nextToken();
			tokens.add(token);
		}
	}

	public static void print(ArrayList<?> tokens) {
		for (int i = 0; i < tokens.size(); i++) {
			System.out.print(tokens.get(i) + " ");
		}
		System.out.print("\n");
	}

	// HashMap Operations
	public static void printHash(HashMap<String, Integer> hashMap) {
		Set<?> s = hashMap.entrySet();
		Iterator<?> it = s.iterator();
		while (it.hasNext()) {
			Map.Entry m = (Map.Entry) it.next();
			System.out.println(m.getKey() + "\t" + m.getValue());
		}
	}

	public static ArrayList<String> getHashMap(HashMap<?, ?> hm) {
		ArrayList<String> a = new ArrayList<String>();
		Set<?> s = hm.entrySet();
		Iterator<?> it = s.iterator();
		while (it.hasNext()) {
			Map.Entry m = (Map.Entry) it.next();
			a.add(m.getKey() + "\t" + m.getValue());
		}
		return a;
	}

	public static String getKeysFromValue(HashMap<Integer, String> hm,
			String value) {
		Set<?> s = hm.entrySet();
		// Move next key and value of HashMap by iterator
		Iterator<?> it = s.iterator();
		while (it.hasNext()) {
			// key=value separator this by Map.Entry to get key and value
			Map.Entry m = (Map.Entry) it.next();
			if (m.getValue().equals(value))
				return m.getKey() + "";
		}
		System.err.println("Error, can't find the data in Hashmap!");
		return null;
	}

	public static void readHash(String type_map, HashMap<String, String> typeMap) {

		ArrayList<String> types = new ArrayList<String>();
		ArrayList<String> tokens = new ArrayList<String>();

		if (type_map != null) {
			readLines(type_map, types);
			for (int i = 0; i < types.size(); i++) {
				if (!types.get(i).isEmpty()) {
					FileUtil.tokenize(types.get(i), tokens);
					if (tokens.size() != 0) {
						if (tokens.size() != 2) {
							for (int j = 0; j < tokens.size(); j++) {
								System.out.print(tokens.get(j) + " ");
							}
							System.err
									.println(type_map
											+ " Error ! Not two elements in one line !");
							return;
						}
						if (!typeMap.containsKey(tokens.get(0)))
							typeMap.put(tokens.get(0), tokens.get(1));
						else {
							System.out.println(tokens.get(0) + " "
									+ tokens.get(1));
							System.err.println(type_map
									+ " Error ! Same type in first column !");
							return;
						}
					}
					tokens.clear();
				}
			}
		}
	}

	public static void readHash2(String type_map,
			HashMap<String, Integer> hashMap) {

		ArrayList<String> types = new ArrayList<String>();
		ArrayList<String> tokens = new ArrayList<String>();

		if (type_map != null) {
			readLines(type_map, types);
			for (int i = 0; i < types.size(); i++) {
				if (!types.get(i).isEmpty()) {
					FileUtil.tokenize(types.get(i), tokens);
					if (tokens.size() != 0) {
						if (tokens.size() != 2) {
							for (int j = 0; j < tokens.size(); j++) {
								System.out.print(tokens.get(j) + " ");
							}
							System.err
									.println(type_map
											+ " Error ! Not two elements in one line !");
							return;
						}
						if (!hashMap.containsKey(tokens.get(0)))
							hashMap.put(tokens.get(0),
									new Integer(tokens.get(1)));
						else {
							System.out.println(tokens.get(0) + " "
									+ tokens.get(1));
							System.err.println(type_map
									+ " Error ! Same type in first column !");
							return;
						}
					}
					tokens.clear();
				}
			}
		}
	}

	public static void readHash3(String type_map,
			HashMap<String, Integer> hashMap) {

		ArrayList<String> types = new ArrayList<String>();
		ArrayList<String> tokens = new ArrayList<String>();

		if (type_map != null) {
			readLines(type_map, types);
			for (int i = 0; i < types.size(); i++) {
				if (!types.get(i).isEmpty()) {
					FileUtil.tokenize(types.get(i), tokens);
					if (tokens.size() != 0) {
						if (tokens.size() < 2) {
							for (int j = 0; j < tokens.size(); j++) {
								System.out.print(tokens.get(j) + " ");
							}
							System.err
									.println(type_map
											+ " Error ! Not two elements in one line !");
							return;
						}
						String key = tokens.get(0);
						String value = tokens.get(tokens.size() - 1);
						for (int no = 1; no < tokens.size() - 1; no++) {
							key += " " + tokens.get(no);
						}
						if (!hashMap.containsKey(key))
							hashMap.put(key, new Integer(value));
						else {
							System.out.println(key + " " + value);
							System.err.println(type_map
									+ " Error ! Same type in first column !");
							return;
						}
					}
					tokens.clear();
				}
			}
		}
	}

	/**
	 * Create a directory by calling mkdir();
	 * 
	 * @param dirFile
	 */
	public static void mkdir(File dirFile) {
		try {
			// File dirFile = new File(mkdirName);
			boolean bFile = dirFile.exists();
			if (bFile == true) {
				System.err.println("The folder exists.");
			} else {
				System.err
						.println("The folder do not exist,now trying to create a one...");
				bFile = dirFile.mkdir();
				if (bFile == true) {
					System.out.println("Create successfully!");
				} else {
					System.err
							.println("Disable to make the folder,please check the disk is full or not.");
				}
			}
		} catch (Exception err) {
			System.err.println("ELS - Chart : unexpected error");
			err.printStackTrace();
		}
	}

	/**
	 * 
	 * @param path
	 * @return
	 */
	static public boolean deleteDirectory(File path) {
		if (path.exists()) {
			File[] files = path.listFiles();
			for (int i = 0; i < files.length; i++) {
				if (files[i].isDirectory()) {
					deleteDirectory(files[i]);
				} else {
					files[i].delete();
				}
			}
		}
		return (path.delete());
	}

	/**
	 * List files in a given directory
	 * 
	 * */
	static public String[] listFiles(String inputdir) {
		File dir = new File(inputdir);

		String[] children = dir.list();
		if (children == null) {
			// Either dir does not exist or is not a directory
		} else {
			for (int i = 0; i < children.length; i++) {
				// Get filename of file or directory
				String filename = children[i];
			}
		}

		return children;
	}

	/**
	 * List files in a given directory
	 * 
	 * */
	static public String[] listFilteredFiles(String inputdir,
			final String filterCondition) {
		File dir = new File(inputdir);

		String[] children = dir.list();
		// It is also possible to filter the list of returned files.
		// This example does not return any files that start with `.'.
		FilenameFilter filter = new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.endsWith(filterCondition);
			}
		};
		children = dir.list(filter);

		return children;
	}

	/**
	 * List files recursively in a given directory
	 * 
	 * */
	static public void listFilesR() {
		File dir = new File("directoryName");

		String[] children = dir.list();

		// The list of files can also be retrieved as File objects
		File[] files = dir.listFiles();

		// This filter only returns directories
		FileFilter fileFilter = new FileFilter() {
			public boolean accept(File file) {
				return file.isDirectory();
			}
		};
		files = dir.listFiles(fileFilter);

	}

	/**
	 * Frequently used functions
	 * */
	static public int count(String a, String contains) {
		int i = 0;
		int count = 0;
		while (a.contains(contains)) {
			i = a.indexOf(contains);
			a = a.substring(0, i)
					+ a.substring(i + contains.length(), a.length());
			count++;
		}
		return count;
	}

	public static void print(String[] files) {

		for (int i = 0; i < files.length; i++) {
			System.out.print(files[i] + " ");
		}
		System.out.print("\n");
	}

	public static void print(int[] c1) {
		for (int i = 0; i < c1.length; i++) {
			System.out.print(c1[i] + " ");
		}
		System.out.println();
	}

	public static void test() {
		String a = "fdsfdsaf";
		a += "\nfdsaf fd fd";
		a += "\nfd sf fd fd\n";
		System.out.println(a);
		a = a.replaceAll("\n+", " ");
		System.out.println(a);
		System.exit(0);
	}

	public static void readHash(String type_map, HashMap<String, String> typeMap, boolean flag) {

		ArrayList<String> types = new ArrayList<String>();
		ArrayList<String> tokens = new ArrayList<String>();
		
		if(type_map != null) {
			readLines(type_map, types);
			for (int i = 0; i < types.size(); i++) {
				if (!types.get(i).isEmpty()) {
					FileUtil.tokenize(types.get(i), tokens);
					if(tokens.size() != 0) {
						if (tokens.size() != 2) {
							for(int j = 0; j < tokens.size(); j++) {
								System.out.print(tokens.get(j)+" ");
							}
							System.err
									.println(type_map + " Error ! Not two elements in one line !");
							return;
						}
						String tokens0 = "";
						String tokens1 = "";
						if(flag) {
							tokens0 = tokens.get(0).trim();
							tokens1 = tokens.get(1).trim();
						} else {
							tokens0 = tokens.get(1).trim();
							tokens1 = tokens.get(0).trim();
						}
						if (!typeMap.containsKey(tokens0))
							typeMap.put(tokens0, tokens1);
						else {
							System.err.println(tokens0 + " " + tokens1);
							System.err
									.println(type_map + " Ignore this one ! Same type in first column !");
						}
					}
					tokens.clear();
				}
			}
		}
	}
	
	public static void writeHashMapS(String file, HashMap<?, ?> hashMap) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));

			Set<?> s = hashMap.entrySet();
			Iterator<?> it = s.iterator();
			while (it.hasNext()) {
				Map.Entry m = (Map.Entry) it.next();
				String v = "";
				for(String va: (List<String>)m.getValue())
				{
					v+=(va+"\t");
				}
				v = v.substring(0,v.length()-1);
				writer.write(m.getKey() + "," + v + "\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static HashMap<String, List<String>> readHashMap1(String file) {
		BufferedReader reader = null;
		String line = "";
		HashMap<String, List<String>> m =new HashMap<String, List<String>>();
		try {
			reader = new BufferedReader(new FileReader(new File(file)));
			while ((line = reader.readLine()) != null) {
				String[] terms = line.split("\t");
				List<String> l = new ArrayList<String>();
				for(String t: terms[1].split(","))
				{
					l.add(t);
				}
				m.put(terms[0],l);
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
		return m;
	}
	
	public static String readLine(String file, int linenb) {
		BufferedReader reader = null;
		String res = "";
		int count = 0;
		try {

			reader = new BufferedReader(new FileReader(new File(file)));

			String line = null;
			while ((line = reader.readLine()) != null) {
				if(count ==linenb)
				{
					res = line;
					break;
				}
				count ++;
				
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
		return res;
	}
	
	public static void writeLines1(String file, int[][] a) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(new File(file)));
			for(int i=0; i<a.length; i++)
			{
				String line = "";
				for(int j=0; j<a[i].length; j++)
				{
					line += (a[i][j]+"\t");
				}
				line = line.substring(0, line.length()-1);
				writer.write(line+"\n");
			}				

		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
}