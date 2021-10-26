import java.io.*;
import java.util.*;

/**
* Input/output class for .fasta file.
*
* @author 	HUYLENBROECK Florent
*/
class FastaIO{
	/**
	* Opens a .fasta file and imports the fragments that it describes.
	*
	* @param path 	String, path to the .fasta file.
	* @return 		Arraylist<String> that contains all the retreived fragments.
	*/
	public static ArrayList<String> openFasta(String path){
		ArrayList<String> fragments = new ArrayList<String>();
		try{
			BufferedReader fastaReader = new BufferedReader(new FileReader(path));
			String line;
			String fragment = "";
			while((line = fastaReader.readLine()) != null){
				if(line.length()>0){
					if (line.charAt(0)=='>'){
						if (fragment!=""){
							fragments.add(fragment.toLowerCase());
							fragment="";
						}
					}
					else{
						fragment+=line;
					}
				}
			}
			if (fragment!="")
				fragments.add(fragment.toLowerCase());
		} catch(Exception e) {
			e.printStackTrace();
		}

		return fragments;
	}

	/**
	* Writes a .fasta file given an input sequence
	*
	* @param path 				String, path to the .fasta file.
	* @param sequence 			String, the sequence to be written.
	* @param collection_number 	String, the collection number from which the fragments have been read.
	*/
	public static void writeFasta(String path, String sequence, String collection_number){
		String fasta = "> Groupe-6B Collection "+collection_number+" longueur "+sequence.length()+"\n";
		for(int i=1; i<=sequence.length(); i++){
			fasta += sequence.charAt(i-1);
			if(i%80==0){
				fasta+="\n";
			}
		}
		try{
			FileWriter fastaWriter = new FileWriter(path);
			fastaWriter.write(fasta);
			fastaWriter.close();
		} catch(Exception e) {
			e.printStackTrace();
		}

		
	}
}