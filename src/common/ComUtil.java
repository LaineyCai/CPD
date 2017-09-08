package common;
/**
 * 	This class contains common operations for data structure like array, arrayList, HashMap etc.
 * 
 * */

import java.io.File;
import java.io.FileFilter;
import java.io.FilenameFilter;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

public class ComUtil {
	/**
	 * String Operations
	 * */
	public static void tokenize(String line, ArrayList<String> tokens) {
		StringTokenizer strTok = new StringTokenizer(line);
		while (strTok.hasMoreTokens()) {
			String token = filter(strTok.nextToken());
			// for toy
			//if(token.length()>2)
				tokens.add(token);
		}
	}
	
    /** Create a new string by eliminating non-alphanumeric chars */
    public static String filter(String s) {
            // Create a string buffer
            StringBuffer strBuf = new StringBuffer();

            // Examine each char in the string to skip alphanumeric char
            for (int i = 0; i < s.length(); i++) {
                    if (Character.isLetterOrDigit(s.charAt(i))) {
                    strBuf.append(s.charAt(i));
                    }
            }

            // Return a new filtered string
            return strBuf.toString();
    }

	/**
	 * Print
	 * */
	public static void print(ArrayList tokens) {
		for (int i = 0; i < tokens.size(); i++) {
			System.out.print(tokens.get(i) + " ");
		}
		System.out.print("\n");
	}

	public static void print(String[] files) {

		for (int i = 0; i < files.length; i++) {
			System.out.print(files[i] + " ");
		}
		System.out.print("\n");
	}

	/**
	 * HashMap Operations
	 * */
	public static void printHash(HashMap<String, Integer> hashMap) {
		System.out.println("Print HashMap");
		Set s = hashMap.entrySet();
		Iterator it = s.iterator();
		while (it.hasNext()) {
			Map.Entry m = (Map.Entry) it.next();
			System.out.println(m.getKey() + "\t" + m.getValue());
		}
	}

	public static ArrayList<String> getHashMap(HashMap<String, String> hm) {
		ArrayList<String> a = new ArrayList<String>();
		Set s = hm.entrySet();
		Iterator it = s.iterator();
		while (it.hasNext()) {
			Map.Entry m = (Map.Entry) it.next();
			a.add(m.getKey() + "\t" + m.getValue());
		}
		return a;
	}

	public static ArrayList<String> getHashMap2(HashMap<String, Integer> hm) {
		ArrayList<String> a = new ArrayList<String>();
		Set s = hm.entrySet();
		Iterator it = s.iterator();
		while (it.hasNext()) {
			Map.Entry m = (Map.Entry) it.next();
			a.add(m.getKey() + "\t" + m.getValue());
		}
		return a;
	}

	public static String getKeysFromValue(HashMap<Integer, String> hm,
			String value) {
		Set s = hm.entrySet();
		// Move next key and value of HashMap by iterator
		Iterator it = s.iterator();
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
			FileUtil.readLines(type_map, types);
			for (int i = 0; i < types.size(); i++) {
				if (!types.get(i).isEmpty()) {
					ComUtil.tokenize(types.get(i), tokens);
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
			FileUtil.readLines(type_map, types);
			for (int i = 0; i < types.size(); i++) {
				if (!types.get(i).isEmpty()) {
					ComUtil.tokenize(types.get(i), tokens);
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
			HashMap<String, Double> hashMap) {

		ArrayList<String> types = new ArrayList<String>();
		ArrayList<String> tokens = new ArrayList<String>();

		if (type_map != null) {
			FileUtil.readLines(type_map, types);
			for (int i = 0; i < types.size(); i++) {
				if (!types.get(i).isEmpty()) {
					ComUtil.tokenize(types.get(i), tokens);
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
									new Double(tokens.get(1)));
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

	public static double readHashTopValue(HashMap<String, Integer> scores, int k) {
		List list = new LinkedList(scores.entrySet());
		int count = 0;
		int value = 0;
		double res = 0;
		for (Iterator it = list.iterator(); count < k && it.hasNext();) {
			Map.Entry entry = (Map.Entry) it.next();
			value = (Integer) entry.getValue();
			res += (double) value * Math.log(2) / Math.log(count + 2);
			// res += (Integer) entry.getValue();
			count++;
		}
		return res;
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

	@SuppressWarnings("unchecked")
	public static HashMap<?,?> sortByValue(HashMap<?,?> map, final int flag) {
		// flag = 0 decreasing order otherwise increasing
		List list = new LinkedList(map.entrySet());
		Collections.sort(list, new Comparator() {
			public int compare(Object o1, Object o2) {
				if(flag == 0 )
					return ((Comparable) ((Map.Entry) (o2)).getValue())
						.compareTo(((Map.Entry) (o1)).getValue());
				else
					return ((Comparable) ((Map.Entry) (o1)).getValue())
					.compareTo(((Map.Entry) (o2)).getValue());
			}
		});

		HashMap result = new LinkedHashMap();
		for (Iterator it = list.iterator(); it.hasNext();) {
			Map.Entry entry = (Map.Entry) it.next();
			result.put(entry.getKey(), entry.getValue());
		}
		return result;
	}

	public static double getSumValue(HashMap<String, Double> map) {
		Double count = 0.0D;
		List list = new LinkedList(map.entrySet());
		for (Iterator it = list.iterator(); it.hasNext();) {
			Map.Entry entry = (Map.Entry) it.next();
			count += map.get(entry.getKey());
		}
		return count;
	}

	public static int getFrequentElement(int[] bcp) {
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		ArrayList<Integer> count = new ArrayList<Integer>();
		ArrayList<Integer> uniId = new ArrayList<Integer>();
		int id = 0;

		for (int col = 0; col < bcp.length; col++) {
			// System.out.print(bcp[col] + "\t");
			int no = 0;
			if (!map.containsKey(bcp[col])) {
				map.put(bcp[col], id++);
				count.add(1);
				uniId.add(bcp[col]);
			} else {
				no = map.get(bcp[col]);
				count.set(no, count.get(no) + 1);
			}
		}

		int maximum = Integer.MIN_VALUE;
		int maxId = Integer.MIN_VALUE;
		for (int i = 0; i < count.size(); i++) {
			// System.out.print(uniId.get(i) + ":" + count.get(i) + ",\t");
			if (maximum < count.get(i)) {
				maximum = count.get(i);
				maxId = uniId.get(i);
			}
		}
		// System.out.println();

		map.clear();
		uniId.clear();
		count.clear();
		return maxId;
	}

	public static void getFrequentElement(int[][] bcp, int[] res, char flag) {
		if (flag == 'r') {
			for (int row = 0; row < bcp.length; row++) {
				res[row] = getFrequentElement(bcp[row]);
			}
		} else {
			int colL = bcp[0].length;
			int[] column = new int[bcp.length];
			for (int col = 0; col < colL; col++) {
				for (int row = 0; row < bcp.length; row++) {
					column[row] = bcp[row][col];
				}
				res[col] = getFrequentElement(column);
			}
		}
	}

	public static short getFrequentElement(short[] bcp) {
		HashMap<Short, Short> map = new HashMap<Short, Short>();
		ArrayList<Short> count = new ArrayList<Short>();
		ArrayList<Short> uniId = new ArrayList<Short>();
		short id = 0;

		for (short col = 0; col < bcp.length; col++) {
			// System.out.print(bcp[col] + "\t");
			short no = 0;
			if (!map.containsKey(bcp[col])) {
				map.put(bcp[col], id++);
				count.add((short) 1);
				uniId.add(bcp[col]);
			} else {
				no = map.get(bcp[col]);
				count.set(no, (short) (count.get(no) + 1));
			}
		}

		short maximum = Short.MIN_VALUE;
		short maxId = Short.MIN_VALUE;
		for (int i = 0; i < count.size(); i++) {
			// System.out.print(uniId.get(i) + ":" + count.get(i) + ",\t");
			if (maximum < count.get(i)) {
				maximum = count.get(i);
				maxId = uniId.get(i);
			}
		}
		// System.out.println();

		map.clear();
		uniId.clear();
		count.clear();
		return maxId;
	}

	public static boolean getFrequentElementBinary(int[] sample) {
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		ArrayList<Integer> count = new ArrayList<Integer>();
		ArrayList<Integer> uniId = new ArrayList<Integer>();
		int id = 0;

		for (int col = 0; col < sample.length; col++) {
			// System.out.print(bcp[col] + "\t");
			int no = 0;
			if (!map.containsKey(sample[col])) {
				map.put(sample[col], id++);
				count.add(1);
				uniId.add(sample[col]);
			} else {
				no = map.get(sample[col]);
				count.set(no, count.get(no) + 1);
			}
		}

		int maximum = Integer.MIN_VALUE;
		int maxId = Integer.MIN_VALUE;
		for (int i = 0; i < count.size(); i++) {
			// System.out.print(uniId.get(i) + ":" + count.get(i) + ",\t");
			if (maximum < count.get(i)) {
				maximum = count.get(i);
				maxId = uniId.get(i);
			}
		}
		// System.out.println();

		map.clear();
		uniId.clear();
		count.clear();
		if(maxId == 1)
			return true;
		else
			return false;
	}
	
	public static int[] CountElmt(ArrayList<Integer> newScores1,
			ArrayList<Integer> scores) {
		int a[] = new int[scores.size()];
		for (int i = 0; i < scores.size(); i++) {
			a[i] = 0;
		}
		for (int i = 0; i < newScores1.size(); i++) {
			int value = newScores1.get(i);
			int pos = scores.indexOf(value);
			a[pos]++;
		}
		return a;
	}

	public static int countCommElmts(ArrayList<Integer> newScores1,
			ArrayList<Integer> newScores2) {
		int count = 0;
		for (int i = 0; i < newScores1.size(); i++) {
			if (newScores1.get(i) == newScores2.get(i))
				count++;
		}
		return count;
	}

	public static void uniqe(int[] words, ArrayList<Integer> tempUniqueWords,
			ArrayList<Integer> tempCounts) {
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

	public static void uniqe(ArrayList<Integer> items) {
		// add elements to al, including duplicates
		HashSet<Integer> hs = new HashSet<Integer>();
		hs.addAll(items);
		items.clear();
		items.addAll(hs);
	}

	public static void getTop(float[] array, ArrayList<Integer> rankList, int i) {
		int index = 0;
		int count = 0;
		HashSet<Integer> scanned = new HashSet<Integer>();
		float max = Float.MIN_VALUE;
		for (int m = 0; m < i && m < array.length; m++) {
			max = Float.MIN_VALUE;
			for (int no = 0; no < array.length; no++) {
				if (array[no] >= max && !scanned.contains(no)) {
					index = no;
					max = array[no];
				}
			}
			scanned.add(index);
			rankList.add(index);
			//System.out.println(m + "\t" + index);
		}
	}
	
	public static void getTop(float[] array, int[] rankList, int i) {
		int index = 0;
		int count = 0;
		HashSet<Integer> scanned = new HashSet<Integer>();
		//int[] rankedList = new int[(i>array.length?array.length:i)];
		int[] rankedList = new int[i];	//faster, but has to gurantee i<=array.length
		float max = Float.MIN_VALUE;
		for (int m = 0; m < i && m < array.length; m++) {
			max = Float.MIN_VALUE;
			for (int no = 0; no < array.length; no++) {
				if (array[no] >= max && !scanned.contains(no)) {
					index = no;
					max = array[no];
				}
			}
			scanned.add(index);
			rankList[m] = index;
			//System.out.println(m + "\t" + index);
		}
	}
	public static void getTop1(float[] array, int[] rankList, int i) {
		int index = 0;
		int count = 0;
		HashSet<Integer> scanned = new HashSet<Integer>();
		//int[] rankedList = new int[(i>array.length?array.length:i)];
		int[] rankedList = new int[i];	//faster, but has to gurantee i<=array.length
		float max = Float.MIN_VALUE;
		int pindex=-1;
		for (int m = 0; m < i && m < array.length; m++) {
			max = Float.MIN_VALUE;
			for (int no = 0; no < array.length; no++) {
				if (array[no] >= max && !scanned.contains(no)) {
					index = no;
					max = array[no];
				}
			}
			if(index!=pindex)
			{
				scanned.add(index);
				rankList[m] = index;
				pindex = index;
			}
			else
			{
				int rindex = (int) Math.floor(Math.random() * (array.length));
				scanned.add(rindex);
				rankList[m] = rindex;
			}
			//System.out.println(m + "\t" + index);
		}
	}
	
	public static int[] getTop(float[] array, float value) {
		int index = 0;
		int count = 0;
		HashSet<Integer> scanned = new HashSet<Integer>();
		ArrayList<Integer> rl = new ArrayList<Integer>();
		//int[] rankedList = new int[(i>array.length?array.length:i)];
//		int[] rankedList = new int[i];	//faster, but has to gurantee i<=array.length
		float max = Float.MIN_VALUE;
		for (int m = 0; m < array.length; m++) {
			max = Float.MIN_VALUE;
			for (int no = 0; no < array.length; no++) {
				if (array[no] >= max && !scanned.contains(no)) {
					index = no;
					max = array[no];
				}
			}
			scanned.add(index);
			if(array[index]>=value)
			{
				rl.add(index);
			}
			else
			{
				break;
			}
			
		}
		int[] rankList = new int[rl.size()];
//		System.out.println(rl.size());
		for(int i=0; i<rl.size(); i++)
		{
			rankList[i] = rl.get(i);
		}	
		return rankList;
	}

	public static int sample(float[] p, int T) {
		float[] pt = new float[T];
		pt[0] = p[0];
		for (int i = 1; i < T; i++) {
			pt[i] = p[i] + pt[i - 1];
		}

		// scaled sample because of unnormalized p[]
		double rouletter = (double) (Math.random() * pt[T - 1]);
		short sample = 0;
		for (sample = 0; sample < T; sample++) {
			if (pt[sample] > rouletter)
				break;
		}
		return sample;
	}

	public static int sample(double[] p, int T) {
		double[] pt = new double[T];
		pt[0] = p[0];
		for (int i = 1; i < T; i++) {
			pt[i] = p[i] + pt[i - 1];
		}

		// scaled sample because of unnormalized p[]
		double rouletter = (double) (Math.random() * pt[T - 1]);
		short sample = 0;
		for (sample = 0; sample < T; sample++) {
			if (pt[sample] > rouletter)
				break;
		}
//		if(sample <0 || sample>=T)
//			System.out.println("wrong!!!"+rouletter);
		return sample;
	}
	
	public static int sample(double[] p, int T, double[] g, double[] e) {
		double[] pt = new double[T];
		pt[0] = p[0];
		for (int i = 1; i < T; i++) {
			pt[i] = p[i] + pt[i - 1];
		}

		// scaled sample because of unnormalized p[]
		double rouletter = (double) (Math.random() * pt[T - 1]);
		short sample = 0;
		for (sample = 0; sample < T; sample++) {
			if (pt[sample] > rouletter)
				break;
		}
//		if(sample <0 || sample>=T)
//		{
//			for(int i=0; i<T; i++)
//			{
//				System.out.print(p[i]+", "+g[i]+", "+e[i]+"\t");
//			}
//			System.out.println();
//		}
//			System.out.println(rouletter+"\t"+p.toString()+"\t"+g.toString()+"\t"+e.toString());
		return sample;
	}
	
	
	
	public static int findMax(double[] p, int T) {
		double maxvalue = -1;
		int maxind = -1;
		for(int i=0; i<T; i++)
		{
			if(p[i]>maxvalue)
			{
				maxvalue = p[i];
				maxind = i;
			}
		}
		return maxind;
	}
	
	public static String getMonthNumber(String month)
	{
		String monthnumber = "";
		switch(month)
		{
			case "Jan":
				monthnumber = "01";
				break;
			case "Feb":
				monthnumber = "02";
				break;
			case "Mar":
				monthnumber = "03";
				break;
			case "Apr":
				monthnumber = "04";
				break;
			case "May":
				monthnumber = "05";
				break;
			case "Jun":
				monthnumber = "06";
				break;
			case "Jul":
				monthnumber = "07";
				break;
			case "Aug":
				monthnumber = "08";
				break;
			case "Sep":
				monthnumber = "09";
				break;
			case "Oct":
				monthnumber = "10";
				break;
			case "Nov":
				monthnumber = "11";
				break;
			case "Dec":
				monthnumber = "12";
				break;
		}
		return monthnumber;
	}
	
	//dtstring example: "Thu Jan 01 02:25:38 CDT 2011"; use hour as unit
	public static int getTimeStamp(String dtstring)
	{
		int timestamp = 0;
		String[] terms = dtstring.split(" ");
		String year = terms[5], month = terms[1], date = terms[2], time = terms[3];
		String monthnumber = getMonthNumber(month);
		String dt = year+"-"+monthnumber+"-"+date+" "+time;
		
		String firstdate = "2006-07-01 00:00:00";
		SimpleDateFormat sf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		Long diff = null, h = null, d = null;
		try {
			diff = sf.parse(dt).getTime()-sf.parse(firstdate).getTime();
			h = diff/1000/60/60%24;	//hour
			d = diff/1000/60/60/24;	//day
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		timestamp = (int) (d*24+h+1);
		return timestamp;
	}
	
	//use day as unit
	public static int getDayStamp(String dtstring)
	{
		int timestamp = 0;
		String[] terms = dtstring.split(" ");
		String year = terms[5], month = terms[1], date = terms[2], time = terms[3];
		String monthnumber = getMonthNumber(month);
		String dt = year+"-"+monthnumber+"-"+date+" "+time;
		
		String firstdate = "2006-07-01 00:00:00";
		SimpleDateFormat sf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		Long diff = null, h = null, d = null;
		try {
			diff = sf.parse(dt).getTime()-sf.parse(firstdate).getTime();
//			h = diff/1000/60/60%24;	//hour
			d = diff/1000/60/60/24;	//day
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		timestamp = (int) (long) d;
		return timestamp;
	}
	
	//generate key for hashmap
	public static String generateKey(int[] key)
	{
		String res = "";
		for(int item: key)
		{
			res+=(item+"\t");
		}
		return res;
	}
}
