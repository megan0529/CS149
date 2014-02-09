import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;



public class Ngram {
  
   
   public static void main(String[] args) throws Exception {
      //    hadoop jar ngram.jar Ngram 4 query1.txt /wikipedia/8gb output
      //                               0 1          2              3
      Configuration conf = new Configuration();
      
      /*local mode for debugging*/
      conf.set("mapreduce.jobtracker.address", "local");
      conf.set("fs.defaultFS", "file:///");
      /*local mode for debugging*/
      
      //delete the old output:
      FileSystem fs = new Path(args[3]).getFileSystem(conf);
      fs.delete(new Path(args[3]), true);

      conf.setInt("n", Integer.parseInt(args[0]));
      conf.setInt("k", 20);
      conf.set("queryFile", args[1]);
      conf.set("xmlStart", "<page>");
      conf.set("xmlEnd", "</page>");

      //setup job: 
      Job job = new Job(conf, "Ngram");
      job.setJarByClass(Ngram.class);

      job.setMapOutputKeyClass(IntTextWritableComparable.class);// <(cnt, title), null>
      job.setMapOutputValueClass(NullWritable.class); //consistant with Mapper<1,2,3,4>types
      job.setOutputKeyClass(IntWritable.class);//cnt
      job.setOutputValueClass(Text.class);//title

      job.setMapperClass(NgramMapper.class);
      job.setReducerClass(NgramReducer.class);

      job.setInputFormatClass(XmlInputFormat.class);//self-implemented, extends TextInputFormat
      job.setOutputFormatClass(TextOutputFormat.class);

      FileInputFormat.addInputPath(job, new Path(args[2]));
      FileOutputFormat.setOutputPath(job, new Path(args[3]));

      job.waitForCompletion(true);
   }

  
   
   
   // Mapper class
   public static class NgramMapper extends
         Mapper<LongWritable, Text, IntTextWritableComparable, NullWritable> {

      private Set<String> queryNgramSet = new HashSet<String>();// 3-gram: word1_word2_word3
      private int n;
      
      @Override
      protected void setup(Context context) 
               throws IOException, InterruptedException{
         System.out.println("\n\n");
         super.setup(context);
         Configuration conf = context.getConfiguration();
         n = conf.getInt("n", 0);
         Path path = new Path(conf.get("queryFile"));//new Path("query1.txt")
         FileSystem fs = path.getFileSystem(conf);
         BufferedReader reader = new BufferedReader(new InputStreamReader(fs.open(path), "UTF-8"));
         
         StringBuilder sb = new StringBuilder();
         String line;
         while ((line = reader.readLine())!=null) {
            sb.append(line+"\n");
            System.out.println(line);
         }
         reader.close();
         System.out.println("\n\n");
         
      }

      @Override
      protected void map(LongWritable key, Text value, Context context) 
               throws IOException, InterruptedException {

      }

      protected void cleanup(Context context) 
               throws IOException, InterruptedException {
      
      }
   }
   
   // Reducer Class
   public static class NgramReducer extends
         Reducer<IntTextWritableComparable, NullWritable, IntWritable, Text> {

      @Override
      protected void setup(Context context) {

      }

      @Override
      public void reduce(IntTextWritableComparable key, Iterable<NullWritable> values, Context context) 
               throws IOException, InterruptedException {
      
      }

      @Override
      protected void cleanup(Context context) 
              throws IOException, InterruptedException {
        
      }
   } 

   
}
