package cpddblp;
import gnu.trove.map.hash.TIntObjectHashMap;


/**
 * A Document object represents a user, the document has been preprocessed, all words format as wordid and time format as timestampid
 */

public class User {
	
	// the ID or title of this document
	private int id;

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	private int [][] docWords;	//tweets 
	//use hashmap instead of two-d array for quick search when sample z, we need to find the retweet info for each tweet, need to be efficient
	TIntObjectHashMap<int[][]> retweets;
	int retweetcount;
	//cite is different from retweet, one tweet can only be a retweet from at most one tweet, but one paper can cite multiple papers
	
	private int [] docTimestamp;
	
	// Init. Document -> get DataLines from Document
	public User(int id, String bodys) {
//		int Maxusernb = 1000;
		this.id = id;
		
		String[] body = bodys.split(";");
		int totallinenb = body.length;
		
		this.docWords = new int[totallinenb][];
		this.docTimestamp = new int[totallinenb];
		retweets = new TIntObjectHashMap<int[][]> ();
		
		retweetcount = 0;
		
		for(int lineNb = 0; lineNb<totallinenb; lineNb++)
		{
			String line = body[lineNb];
			String[] terms = line.split("\t");
			docTimestamp[lineNb] = Integer.parseInt(terms[0]);
			String title = terms[1];
			String[] words = title.split(" ");
			docWords[lineNb] = new int[words.length];
			for(int w = 0; w < words.length; w++)
				docWords[lineNb][w] = Integer.parseInt(words[w]);
			int[][] pair = new int[(terms.length-2)/2][2];
			if(pair.length>0)
			{
				for(int i=2; i<terms.length-1; i+=2)
				{	//add one index, indicate the ith citation 
					pair[(i-2)/2][0]= Integer.parseInt(terms[i]);
					pair[(i-2)/2][1] = Integer.parseInt(terms[i+1]);
					retweetcount++;
				}	
	
				retweets.put(lineNb, pair);
			}
		}
		
}
	

	public int getRetweetcount() {
		return retweetcount;
	}

	public void setRetweetcount(int retweetcount) {
		this.retweetcount = retweetcount;
	}

	public int[][] getDocWords() {
		return docWords;
	}

	public void setDocWords(int[][] docWords) {
		this.docWords = docWords;
	}

	public int[] getDocTimeStamp() {
		return docTimestamp;
	}

	public void setDocTimeStamp(int[] docts) {
		this.docTimestamp = docts;
	}

	public TIntObjectHashMap<int[][]> getRetweets() {
		return retweets;
	}

	public void setRetweets(TIntObjectHashMap<int[][]> retweets) {
		this.retweets = retweets;
	}

	public int[] getDocTimestamp() {
		return docTimestamp;
	}

	public void setDocTimestamp(int[] docTimestamp) {
		this.docTimestamp = docTimestamp;
	}


	
}