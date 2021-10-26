import java.util.*;

/**
* Class representing a LinkedList used to build the contig.
*
* @author HUYLENBROECK Florent
*/
class AlignmentLinkedList{

	/**
	* Sentinel node to begin the list.
	*/
	private AlignmentNode head;

	/**
	* Sentinel node to end the list.
	*/
	private AlignmentNode tail;

	/**
	* Fragments are stored in reverse order, since they are aligned backwards using the alignment matrix.
	* 
	* @param collection 	Collection containing the fragments to align.
	* @param path 			int[], a greedy hamiltonian path amongst the collection semi-global alignment scores.
	*/
	public AlignmentLinkedList(Collection collection, int[] path){
		int owner_f = path[0];
		Fragment f = collection.getFragment(owner_f);

		head = new AlignmentNode((byte)0, -1);
		tail = new AlignmentNode((byte)0, -1);
		
		AlignmentNode current = head;
		for(int i=f.length()-1; i>=0; i--){
			AlignmentNode new_node = new AlignmentNode(f.bitAt(i), owner_f);
			current.setNext(new_node);
			current = new_node;
		}
		current.setNext(tail);

		for(int i=1; i<path.length; i++){
			Fragment g = collection.getFragment(path[i]);
			align(f, g, owner_f, path[i]);
			f=g;
			owner_f=path[i];
		}
	}

	/**
	* Aligns the fragment G to the fragment F.
	* 
	* @param f 			Fragment, the fragment to align upon.
	* @param g 			Fragment to be aligned.
	* @param owner_f 	int used to mark which nodes belong to F.
	* @param owner_g 	int used to mark which nodes belong to G.
	*/
	public void align(Fragment f, Fragment g, int owner_f, int owner_g){

		int[][] a = f.semiGlobalAlignmentMatrix(g);

		int tmp_max=-2*a.length*a[0].length;
		int index_f=0;
		int index_g=0;

		// Finding entry point in alignment matrix
		for(int j=0; j<a[0].length; j++){
			if(tmp_max<=a[a.length-1][j]){
				tmp_max=a[a.length-1][j];
				index_f=a.length-1;
				index_g=j;
			}
		}

		AlignmentNode current = head;

		// Treating nucleids (at the end) of G unmatched to nucleids of F
		for(int i=0; i<a[0].length-index_g-1; i++){
			AlignmentNode unmatched = new AlignmentNode(g.bitAt(g.length()-i-1), owner_g);
			unmatched.setNext(current.getNext());
			current.setNext(unmatched);
			current = unmatched;
		}

		while(index_f>0 && index_g>0){

			// Find which move gave the best score.
			int left = a[index_f][index_g-1];
			int leftup = a[index_f-1][index_g-1];
			int up = a[index_f-1][index_g];

			int max = Math.max(left, Math.max(up, leftup));

			if(max==leftup){
				// Go to the next node of F and adds G's data to it.
				current=current.getNext(owner_f);
				current.addData(g.bitAt(index_g-1), owner_g);
				index_f--;
				index_g--;
			}
			else if(max==left){
				// Create a new node for G's data and insert it before next node of F.
				AlignmentNode new_g = new AlignmentNode(g.bitAt(index_g-1), owner_g);
				new_g.setNext(current.getNext(owner_f)); // MIGHT LOSE NODES HERE ?
				current.setNext(new_g);
				current=new_g;
				index_g--;
			}
			else if(max==up){
				// Find next F's node.
				current=current.getNext(owner_f); 
				index_f--;
			}
		}
		while(index_g>0){
			// When we reached the end of F without reaching the end of G
			if(current.getNext().equals(tail)){
				// If tail is next, insert G's data.
				AlignmentNode new_g = new AlignmentNode(g.bitAt(index_g-1), owner_g);
				current.setNext(new_g);
				new_g.setNext(tail);
				current=new_g;
			}
			else{
				// If next is not tail, add data to next.
				current=current.getNext();
				current.addData(g.bitAt(index_g-1), owner_g);
			}
			index_g--;
		}
	}

	/**
	* toString override, for printing purpose.
	*
	* @return 	String representing the list.
	*/
	public String toString(){
		String ret = "[HEAD]\n";
		AlignmentNode current = head.getNext();
		while(!current.equals(tail)){
			ret+=current.toString()+"\n";
			current=current.getNext();
		}
		return ret+"[TAIL]";
	}

	/**
	* Computes the consensus contig.
	*
	* @return 	String, the consensus contig.
	*/
	public String getContig(){
		String gitnoc = "";
		AlignmentNode current = head.getNext();
		while(!current.equals(tail)){
			gitnoc+=current.consensus();
			current=current.getNext();
		}
		String contig = "";
		for(int i=gitnoc.length()-1; i>=0; i--){
			contig+=gitnoc.charAt(i);
		}
		return contig;
	}

	/**
	* Class that represents a node in the AlignmentLinkedList.
	*
	* @author 	HUYLENBROECK Florent
	*/
	private class AlignmentNode{

		/**
		* Next node in the list.
		*/
		private AlignmentNode next;

		/**
		* array of length 4. One entry matches to the number of occurence of a nucleid in the node.
		*/
		private int[] data;

		/**
		* Lists the fragments that use this node to store data.
		*/
		private ArrayList<Integer> owners;

		/**
		* @param data 	byte, the data to initialize the node with.
		* @param owner 	int, the first fragment to uses this node.
		*/
		public AlignmentNode(byte data, int owner){

			this.data = new int[4];
			owners = new ArrayList<Integer>();
			this.addData(data, owner);
		}

		/**
		* Getter for the next node.
		*
		* @return 	AlignmentNode following this one in the list.
		*/
		public AlignmentNode getNext(){
			return next;
		}

		/**
		* Getter for the next node that is used by a given fragment.
		*
		* @param owner 	int, the fragment that uses the next node.
		* @return 		AlignmentNode following this one in the list.
		*/
		public AlignmentNode getNext(int owner){
			AlignmentNode current = next;
			while(!current.hasOwner(owner)){
				current=current.getNext();
			}
			return current;
		}

		/**
		* Setter for next node.
		*
		* @param node 	AlignmentNode, next node.
		*/
		public void setNext(AlignmentNode node){
			this.next=node;
		}

		/**
		* Adds data to the node. Increments data array and adds owner to owners list if the owner is not already using this fragment.
		*
		* @param data 	byte, data to add.
		* @param owner 	int, the fragment that adds the new data to the node.
		*/
		public void addData(byte data, int owner){
			if(!owners.contains(owner)){
				this.data[(int)data]++;
				owners.add(owner);
			}
		}

		/**
		* Tells if a node is used by a certain fragment.
		*
		* @param owner 	int, the potential owner.
		* @return 		boolean that tells whether if the fragment uses the node or not.
		*/
		private boolean hasOwner(int owner){
			return (owners.contains(owner));
		}

		/**
		* toString override, for printing purpose.
		*
		* @return 	String representing the node.
		*/
		public String toString(){
			String ret="[";
			ret+="(A:"+data[0];
			ret+=")(C:"+data[1];
			ret+=")(G:"+data[2];
			ret+=")(T:"+data[3];
			ret+=")-owners:"+owners.toString();
			return ret+"]";
		}

		/**
		* Goes trough the node's data to figure out the consensus nucleid by majority vote. If tied : a>c>t>g.
		*
		* @return 	char, consensus nucleid by majority vote.
		*/
		public char consensus(){
			int winner=0;
			for(int i=1; i<4; i++){
				winner = data[i] > data[winner] ? i : winner;
			}
			return Fragment.twoBitsToChar((byte)winner);
		}
	}
}