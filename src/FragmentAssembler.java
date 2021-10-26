import java.util.*;

/**
* Main class of FragmentAssembler package. 
*
* @author 	HUYLENBROECK Florent
*/
class FragmentAssembler{

	private static Collection collection;

	public static void main(String[] args){
		String path_in = "";
		String path_out = "";
		String path_out_ic = "";

		if(args.length==5){
			path_in = args[0];
			if(args[1].equals("-out")){
				path_out=args[2];
			}
			else{
				System.out.println("Error while parsing command. Expected first flag to be \"-out\". Exiting.");
				System.exit(1);
			}
			if(args[3].equals("-out-ic")){
				path_out_ic=args[4];
			}
			else{
				System.out.println("Error while parsing command. Expected second flag to be \"-out-ic\". Exiting.");
				System.exit(1);
			}
			System.out.println("FragmentAssembler - HUYLENBROECK Florent - Group 6B");
			System.out.println("\t1. Opening file \""+path_in+"\".");
			collection = new Collection(FastaIO.openFasta(path_in));
			String collection_n = path_in.replaceAll("[^0-9S]", "");
			System.out.println("\t2. Generating overlap graph.");
			int[][] overlap_graph = getOverlapGraph();
			System.out.println("\t3. Finding a greedy hamiltonian path amongst overlap graph.");
			int[] path = greedyHamiltonianPath(overlap_graph);
			System.out.println("\t4. Aligning fragments.");
			AlignmentLinkedList alignment = new AlignmentLinkedList(collection, path);
			System.out.println("\t5. Building consensus contig using majority vote.");
			String contig = alignment.getContig();
			String contig_ic = invertAndComplement(contig);
			System.out.println("\t6. Saving contig to \""+path_out+"\".");
			FastaIO.writeFasta(path_out, contig, collection_n);
			System.out.println("\t7. Inverting and complementing contig.");
			System.out.println("\t8. Saving inverted and complemented contig to \""+path_out_ic+"\".");
			FastaIO.writeFasta(path_out_ic, contig_ic, collection_n);
			System.out.println("Done.");
		}
		else{
			System.out.println("Error while parsing command. Incorrect number of argument found.");
			System.out.println("Please format you command as follows :");
			System.out.println("java -jar FragmentAssembler.jar <file.fasta> -out <out.fasta> -out-ic <out_ic.fasta>");
			System.exit(1);
		}		
	}

	/**
	* Inverts and complements a contig.
	* A <-> T
	* C <-> G
	*
	* @param contig 	String, the contig to invert and complement.
	* @return 			String, the inverted and complemented contig.
	*/
	private static String invertAndComplement(String contig){
		String ic = "";
		for(int i=contig.length()-1; i>=0; i--){
			switch(contig.charAt(i)){
				case 'a' : ic+='t'; break;
				case 'c' : ic+='g'; break;
				case 'g' : ic+='c'; break;
				case 't' : ic+='a'; break;
				default :;
			}
		}
		return ic;
	}

	/**
	* Builds the overlap graph for the collection, as an adjacency matrix. 
	* Vertices are the pairs of indexes and edges are the value in the array at each pair of index.
	* The matrix being symetrical, only the upper part of the matrix is computed, and mirrored to fill the rest of the matrix.
	* Also, the diagonal is filled with zeroes and ignored during the computation.
	*
	* @return 	int[][], adjacency matrix of the overlap graph.
	*/
	private static int[][] getOverlapGraph(){

		int length = collection.length();

		int[][] graph = new int[length][length];

		for(int i =0; i<length; i++){
			for(int j=i; j<length; j++){
				if(i!=j){
					int[] tmp_score = collection.getFragment(i).semiGlobalAlignmentScore(collection.getFragment(j));
					graph[i][j]=tmp_score[0];
					graph[j][i]=tmp_score[1];
				}
			}
		}

		return graph;
	}

	/**
	* Finds a hamiltonian path in a graph given it's adjacency matrix using greedy heuristic.
	* Algorithm is described at slide 25-26 of the project's presentation slides.
	*
	* @param graph 	int[][], adjacency matrix of the graph.
	* @return 		int[][], an array containing the selected vertices, if the form [f, g]
	*/
	private static int[] greedyHamiltonianPath(int[][] graph){

		int length = collection.length();

		byte[] in = new byte[length], out = new byte[length]; 
		ArrayList<int[]> sets = new ArrayList<int[]>();
		ArrayList<int[]> vertices = new ArrayList<int[]>();

		for(int i=0; i<length; i++){
			sets.add(new int[] {i});

			for(int j=0; j<length; j++){
				if(i!=j){
					vertices.add(new int[] {graph[i][j], i, j});
				}
			}
		}

		vertices.sort(Comparator.comparing(a -> -a[0]));

		int[][] greedy_hamiltionian_path_vertice = new int[length-1][2];
		int greedy_index = 0;

		for(int[] vertex:vertices){
			int f = vertex[1]; 
			int g = vertex[2]; 
			if(in[g]==0 && out[f]==0){
				int[] setF = findSet(sets, f); 
				int[] setG = findSet(sets, g); 
				if(!setF.equals(setG)){
					greedy_hamiltionian_path_vertice[greedy_index] = new int[] {f, g};
					greedy_index++;
					in[g]=1;
					out[f]=1;
					if(union(sets, setF, setG)==1){
						break; 
					}
				}
			}
		}

		// translate vertice list in a path
		int[] greedy_hamiltionian_path = new int[length];


		// find starting point (the one that has no entry in 'in' but has one in 'out')
		for(int i=0; i<length; i++){
			if(in[i]==0){
				greedy_hamiltionian_path[0]=i;
				break;
			}
		}

		// then fill vector 
		for(int i=1; i<length; i++){
			for(int[] vertex : greedy_hamiltionian_path_vertice){
				if(vertex[0]==greedy_hamiltionian_path[i-1]){
					greedy_hamiltionian_path[i]=vertex[1];
					break;
				}
			}
		}


		return greedy_hamiltionian_path;
	}

	/**
	* Finds an element in a set of sets and returns the set containing that element.
	*
	* @param sets 	ArrayList<int[]> the set of sets
	* @param elem 	int, the element to find
	* @return 		int[], set containing the element
	*/
	private static int[] findSet(ArrayList<int[]> sets, int elem){
		int[] ret = {};
		for(int[] set:sets){
			for(int i=0; i<set.length; i++){
				if (set[i]==elem){
					ret=set;
				}
			}
		}
		return ret;
	}

	/**
	* Merges two sets into one.
	*
	* @param sets 	ArrayList<int[]> the set inside which the sets are merged
	* @param set1 	int[], the first set to merge
	* @param set2 	int[], the second set to merge
	* @return 		int, cardianlity of the set containing the merged sets
	*/
	private static int union(ArrayList<int[]> sets, int[] set1, int[] set2){
		int new_length = set1.length+set2.length;
		int[] new_set = new int[new_length];

		int fill_index=0;
		for(int i=0; i<set1.length; i++){
			new_set[fill_index]=set1[i];
			fill_index++;
		}
		for(int i=0; i<set2.length; i++){
			new_set[fill_index]=set2[i];
			fill_index++;
		}

		sets.remove(set1);
		sets.remove(set2);

		sets.add(new_set);

		return sets.size();
	}
}