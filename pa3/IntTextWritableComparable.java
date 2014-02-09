import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.WritableComparable;

public class IntTextWritableComparable implements WritableComparable<IntTextWritableComparable> {

   private IntWritable firstInt;
	private Text secondText;

	public IntTextWritableComparable() {
	   
	   firstInt = new IntWritable();
	   secondText = new Text();
	}
	
	public IntTextWritableComparable(IntTextWritableComparable that) {
      
      firstInt = new IntWritable(that.getFirstInt().get());
      secondText = new Text(that.getSecondText().toString());
   }
	
	public IntTextWritableComparable(IntWritable firstInt, Text secondText) {
	   this.firstInt = firstInt;
	   this.secondText = secondText;
	}

	public IntWritable getFirstInt() {
		return firstInt;
	}

	public Text getSecondText() {
		return secondText;
	}

	
	//	 inherited from interface org.apache.hadoop.io.Writable
	@Override
	public void readFields(DataInput in) throws IOException {
	   firstInt.readFields(in);
		secondText.readFields(in);
	}

	@Override
	public void write(DataOutput out) throws IOException {
		secondText.write(out);
		firstInt.write(out);
	}

	//	 inherited from interface java.lang.Comparable
	@Override
	public int compareTo(IntTextWritableComparable that) {
	   //this-that
	   //mapreduce default is smaller first, but we want large score and large title first, so minus sign
		int cmp = -firstInt.compareTo(that.firstInt);
		if (cmp == 0) {
			return -secondText.compareTo(that.secondText);
		}
		return cmp;
	}

//	@Override
//	public boolean equals(IntTextWritableComparable obj) {
//
//		if (obj instanceof IntTextWritableComparable) {
//		   IntTextWritableComparable that = (IntTextWritableComparable) obj;
//			return (secondText.equals(that.secondText) && firstInt.equals(that.firstInt));
//		}
//
//		return false;
//	}

	@Override
	public int hashCode() {
		return secondText.hashCode() + 13 * firstInt.hashCode();
	}

	@Override
	public String toString() {
	   return firstInt.toString()+"\t"+secondText.toString();
	}
}
