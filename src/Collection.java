import java.util.*;

/**
* Class that describes the collection object. A collection is a set of fragments.
*
* @author 	HUYLENBROECK Florent
*/
class Collection{

	private Fragment[] collection;
	private int length;

	/**
	* @param fragments 	ArrayList<String>. One String in the list is one fragment to be put in the collection.
	*/
	public Collection(ArrayList<String> fragments){

		length = fragments.size();
		collection = new Fragment[length];

		for(int i=0; i<fragments.size(); i++){
			collection[i] = new Fragment(fragments.get(i));
		}
	}

	/**
	* Getter for the collection's length value, being the number of fragment it holds.
	*
	* @return 	int, the number of fragment in the collection.
	*/
	public int length(){
		return length;
	}

	/**
	* toString override, for printing purpose.
	*
	* @return 	String representing the collection.
	*/
	public String toString(){
		String ret = "";
		for(int i=0; i<length; i++){
			ret+=">fragment "+i+"\n"+collection[i].toString()+"\n";
		}
		return ret;
	}

	/**
	* Gets a certain fragment in the collection.
	*
	* @param index 	int, the index of the fragment in the colelction.
	* @return 		Fragment at input index.
	*/
	public Fragment getFragment(int index){
		if(index<length){
			return collection[index];
		}
		else{
			return collection[0];
		}
	}
}