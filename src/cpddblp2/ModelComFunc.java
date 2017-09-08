package cpddblp2;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.util.ArrayList;


public class ModelComFunc {
	public static void writeData(float[] array, ArrayList<String> strings,
			ArrayList<Integer> rankList, BufferedWriter writer, String prefix) {
		PrintWriter writer2 = new PrintWriter(writer);
		for (int row = 0; row < rankList.size(); row++) {
			writer2.printf("%s\t%s\t%f\n", prefix, strings.get(rankList.get(row)),
					array[rankList.get(row)]);
//			writer2.printf(prefix + "\t",
//					strings.get(rankList.get(row)) + "\t" + array[rankList.get(row)] + "\n");
		}
	}

	public static void writeData(float[] vPhiB2, BufferedWriter writer) {
		PrintWriter writer2 = new PrintWriter(writer);
		for (int row = 0; row < vPhiB2.length; row++) {
//			writer2.printf("\t" + vPhiB2[row] + "\n");
			writer2.printf("\t%f\n", vPhiB2[row]);
		}
	}
	
	public static void writeData(double[] vPhiB2, BufferedWriter writer) {
		PrintWriter writer2 = new PrintWriter(writer);
		for (int row = 0; row < vPhiB2.length; row++) {
//			writer2.printf("\t" + vPhiB2[row] + "\n");
			writer2.printf("\t%f\n", vPhiB2[row]);
		}
	}
	
//	public static void writeData(ArrayList<Integer>[] cNP2, BufferedWriter writer) {
//		PrintWriter writer2 = new PrintWriter(writer);
//		writer2 = new PrintWriter(writer);
//		for (int i = 0; i < cNP2.length; i++) {
//			// writer2.printf("%d-th topic:\n", i);
//			for (int j = 0; j < cNP2[i].size(); j++) {
//				// writer2.printf("%s,\t", Doc.getNps().get(cNP2[i].get(j)));
//			}
//			writer2.print("\n\n");
//		}
//	}

	public static void writeData(int[][] phi2, PrintWriter writer2) {
		for (int row = 0; row < phi2.length; row++) {
			// writer2.printf("%d", row);
			for (int col = 0; col < phi2[row].length; col++) {
				writer2.printf("%d\t", phi2[row][col]);
//				writer2.printf(phi2[row][col] + "\t");
			}
			writer2.print("\n");
		}
	}
	
	public static void writeData(float[][] array, BufferedWriter writer) {
		PrintWriter writer2 = new PrintWriter(writer);
		for (int row = 0; row < array.length; row++) {
			// writer2.printf("%d\t", row);
			for (int col = 0; col < array[row].length; col++) {
				writer2.printf("%f\t", array[row][col]);
//				writer2.printf(array[row][col] + "\t");
			}
			writer2.print("\n");
		}
	}
	
	public static void writeData(int[][] array, BufferedWriter writer) {
		PrintWriter writer2 = new PrintWriter(writer);
		for (int row = 0; row < array.length; row++) {
			// writer2.printf("%d\t", row);
			for (int col = 0; col < array[row].length; col++) {
				writer2.printf("%d\t", array[row][col]);
//				writer2.printf(array[row][col] + "\t");
			}
			writer2.print("\n");
		}
	}
	
	public static void writeData(float[][][] array, BufferedWriter writer) {
		PrintWriter writer2 = new PrintWriter(writer);
		for (int row = 0; row < array.length; row++) {
			for (int col = 0; col < array[row].length; col++) {
				for(int i=0; i<array[row][col].length;i++) {
					writer2.printf("%f\t", array[row][col][i]);
				}
				writer2.print("\n");
			}
			writer2.print("\n");
		}
	}
	
	public static void writeData(double[][][] array, BufferedWriter writer) {
		PrintWriter writer2 = new PrintWriter(writer);
		for (int row = 0; row < array.length; row++) {
			for (int col = 0; col < array[row].length; col++) {
				for(int i=0; i<array[row][col].length;i++) {
					writer2.printf("%f\t", array[row][col][i]);
				}
				writer2.print("\n");
			}
			writer2.print("\n");
		}
	}
	
	public static void writeData(double[][] array, BufferedWriter writer) {
		PrintWriter writer2 = new PrintWriter(writer);
		for (int row = 0; row < array.length; row++) {
			// writer2.printf("%d\t", row);
			for (int col = 0; col < array[row].length; col++) {
				writer2.printf("%f\t", array[row][col]);
//				writer2.printf(array[row][col] + "\t");
			}
			writer2.print("\n");
		}
	}

	public static void writeData(double[][] vph2, PrintWriter writer2) {
		for (int row = 0; row < vph2.length; row++) {
			// writer2.printf("%d", row);
			for (int col = 0; col < vph2[row].length; col++) {
				writer2.printf("\t%d", vph2[row][col]);
//				writer2.printf("\t" + vph2[row][col]);
			}
			writer2.print("\n");
		}
	}

	public static void writeData(double[] phi2, PrintWriter writer2) {
		for (int row = 0; row < phi2.length; row++) {
			writer2.printf("\t%d", phi2[row]);
//			writer2.printf("\t" + phi2[row]);
		}
	}
}
