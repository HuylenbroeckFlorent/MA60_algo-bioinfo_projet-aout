/**
* Class that describes a sequence of nucleids.
* A fragment stores a sequence of nucleids as an array of bytes. Each bytes represents at most 4 nucleids (every byte except the last one holds strictly 4).
* One byte being 8 bit, each nucleid is mapped to a two-bit representation.
* A : 00
* C : 01
* G : 10
* T : 11
* The last byte has its leftmost bits being nucleids, bounded by the length of the sequence. It won't read further into the byte than the length allows it, 
* thus no extra A (00) will be read.
*
* @author 	HUYLENBROECK Florent
*/
class Fragment{

	private int length;
	private byte[] fragment;

	/**
	* @param seq 	String that describes the sequence of nucleides. Usually read from a .fasta file.
	*/
	public Fragment(String seq){

		// Figures out fragment's length.
		length = seq.length();
		fragment = new byte[(int)Math.ceil(length/4.0)];

		// Adds data to fragment, 4 nucleids at a time.
		byte subfragment = (byte)0; 
		for(int i=0; i<fragment.length; i++){
			for(int j=0; j<4; j++){

				// If there are still nucleids to read, reads the next one.
				if(i*4+j<length){
					subfragment += charToTwoBits(seq.charAt(i*4+j));
				}

				// Shifts the bits in the current subfragment if 3> nucleids processed in current subfragment.
				if(j<3){
					subfragment = (byte)(subfragment<<2);
				}
			}

			// Add the subfragment to fragment list and resets current subfragment.
			fragment[i]+=subfragment;
			subfragment=(byte)0;
		}
	}

	/**
	* Getter for the fragment's length value.
	*
	* @return 	int, length of the fragment.
	*/
	public int length(){
		return length;
	}

	/**
	* Maps a nucleid's character representation to it's byte value.
	*
	* @param c 	char representing the nucleid. Can be a, c, t, g.
	* @return 	byte, the value given to the input nucleid.
	*/
	public static byte charToTwoBits(char c){
		byte ret = 0;
		switch(c){
			case 'a' : ret = 0; break;
			case 'c' : ret = 1; break;
			case 'g' : ret = 2; break;
			case 't' : ret = 3; break;
			default : ret = -1;
		}
		return (byte)ret;
	}

	/**
	* Maps a nucleid's byte value to it's character representation.
	*
	* @param b 	byte, the nucleid byte value.
	* @return 	char representing the nucleid.
	*/
	public static char twoBitsToChar(byte b){
		char ret = ' ';
		switch(b){
			case 0 : ret = 'a'; break;
			case 1 : ret = 'c'; break;
			case 2 : ret = 'g'; break;
			case 3 : ret = 't'; break;
			default : ret = ' ';
		}
		return ret;
	}

	/**
	* toString override, for printing purpose.
	*
	* @return 	String representing the fragment.
	*/
	public String toString(){
		String ret="";
		for(int i=0; i<length; i++){
			ret+=nucleidAt(i);
		}
		return ret;
	}

	/**
	* Getter for the raw data stored in the fragment object.
	*
	* @return 	byte[], array of byte each holding at most 4 nucleids.
	*/
	public byte[] getFragment(){
		return fragment;
	}

	/**
	* Gives the value of the two-bits representation of a certain nucleid within the fragment.
	*
	* @param index 	int, the index of the nucleid within the sequence.
	* @return 		byte that has it's two rightmost bits being the two-bits representation of the nucleid and other bits set to zero.
	*/
	public byte bitAt(int index){
		byte ret = (byte)0;
		if(index<length){

			// Finds which subfragment contains target nucleid.
			int i = index/4;
			byte subfragment = fragment[i];

			// Finds index of nucleid in subfragment.
			int j = index%4;


			// Shifts the subfragment to have the target nucleid as first nucleid then converts it to char.
			for(int k=0; k<3-j; k++){
				subfragment = (byte)(subfragment>>2);
			}
			ret = (byte)(subfragment&3);
		}
		return ret;
	}

	/**
	* Gives the character representation of a certain nucleid within the fragment.
	*
	* @param index 	int, the index of the nucleid within the sequence.
	* @return 		char representation of the nucleid.
	*/
	public char nucleidAt(int index){
		if(index<length)
			return twoBitsToChar(bitAt(index));
		else
			return 'x';
	}

	/**
	* Computes the semiglobal alignment score of the fragment object with another fragment. This algorithm is optimized to only store one row at a time instead of 
	* the whole matrix. 
	*
	* @param f2 	Fragment to align with the fragment object.
	* @return 		int, semiglobal alignment score.
	*/
	public int[] semiGlobalAlignmentScore(Fragment f2){

		int n = f2.length(), gap_score = -2, mismatch_score = -1, match_score = 1, treshold = 3;

		int[] a = new int[n+1];

		for(int j = 0; j<=n; j++){
			a[j]= 0; 
		}

		int temp_max_g_f = length*n*gap_score;

		for(int i = 1; i<=length; i++){
			int old = a[0];
			a[0] = 0; 
			for(int j = 1; j<=n; j++){
				int temp = a[j];
				int p = (this.bitAt(i-1) == f2.bitAt(j-1)) ? match_score:mismatch_score; // -1 to shift the entries in the sims tab
				a[j] =  Math.max(a[j]+gap_score, Math.max(old+p, a[j-1]+gap_score));
				old = temp;
			}
			temp_max_g_f = Math.max(temp_max_g_f, a[n]);
		}

		int temp_max_f_g = length*n*gap_score;

		for(int i=1; i<a.length; i++){
			temp_max_f_g = Math.max(temp_max_f_g, a[i]);
		}

		return new int[] {temp_max_f_g, temp_max_g_f};
	}

	public int[][] semiGlobalAlignmentMatrix(Fragment f2){

		int n = f2.length(), gap_score = -2, mismatch_score = -1, match_score = 1, treshold = 3;

		int[][] a = new int[length+1][n+1];

		for(int i = 0; i<=length; i++){
			a[i][0] = 0;
		}

		for(int j = 0; j<=n; j++){
			a[0][j]= 0;
		}

		for(int i = 1; i<=length; i++){
			for(int j = 1; j<=n; j++){

				int p = (this.bitAt(i-1) == f2.bitAt(j-1)) ? match_score:mismatch_score;
				a[i][j] =  Math.max(a[i-1][j]+gap_score, Math.max(a[i-1][j-1]+p, a[i][j-1]+gap_score));
			}
		}

		return a;
	}
}
