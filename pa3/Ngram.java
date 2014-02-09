import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Collection;
import java.util.Deque;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;
import java.util.LinkedHashSet;
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
  
   private static Set<String> queryNgramSet = new LinkedHashSet<String>();// 3-gram: word1_word2_word3
   private static int n;
   
   public static void main(String[] args) throws Exception {
      //    hadoop jar ngram.jar Ngram 4 query1.txt /wikipedia/8gb output
      //                               0 1          2              3
      Configuration conf = new Configuration();
      
//      /*local mode for debugging*/
conf.set("mapreduce.jobtracker.address", "local");
conf.set("fs.defaultFS", "file:///");
//      /*local mode for debugging*/
      
      //delete the old output:
      FileSystem fs = new Path(args[3]).getFileSystem(conf);
      fs.delete(new Path(args[3]), true);

      conf.setInt("n", Integer.parseInt(args[0]));
      conf.setInt("k", 20);
      conf.set("queryFile", args[1]);
      conf.set("xmlStart", "<page>");
      conf.set("xmlEnd", "</page>");
      setupQuery(conf);
      //setup job: 
      Job job = new Job(conf, "Ngram");
      job.setJarByClass(Ngram.class);

      job.setMapOutputKeyClass(IntTextWritableComparable.class);// <(cnt, title), null>
      job.setMapOutputValueClass(IntWritable.class); //consistant with Mapper<1,2,3,4>types
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
  
   public static void setupQuery(Configuration conf) throws IOException, InterruptedException
   {
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
      
      Tokenizer tokenizer = new Tokenizer(sb.toString());
      Deque<String> ngramList = new LinkedList<String>();
      while(tokenizer.hasNext()) {
         String token = tokenizer.next();
         ngramList.addLast(token);
         if(ngramList.size() > n) {
            ngramList.removeFirst();
         }
         if (ngramList.size() == n) {
            StringBuilder sb2 = new StringBuilder();
            for(String s : ngramList) {
               sb2.append(s);
               sb2.append(",");
            }
            
            queryNgramSet.add(sb2.toString());
         }
      }
      
      //print ngram
      Iterator<String> iter = queryNgramSet.iterator();
      int i=0;
      while(iter.hasNext()) {
         i++;
         System.out.println("Ngram_"+i+"is:"+iter.next());
      }
      System.out.println("\n\n");  
      //print end
   }
   
   // Mapper class
   public static class NgramMapper extends
         Mapper<LongWritable, Text, IntTextWritableComparable, IntWritable> {

      private TreeSet<IntTextWritableComparable> resultToEmit = new TreeSet<IntTextWritableComparable>();
      
      @Override
      protected void setup(Context context) 
               throws IOException, InterruptedException{
         System.out.println("\n");
         super.setup(context);
      }

      @Override
      protected void map(LongWritable key, Text value, Context context) 
               throws IOException, InterruptedException {
         String page = value.toString();
         String titleString=null;
//         System.out.println(page);
         Pattern pattern = Pattern.compile("<title>(.*)</title>");
         Matcher matcher = pattern.matcher(page);
         if ((matcher.find()) && (matcher.groupCount() > 0)) {
            titleString = matcher.group(1);
//            System.out.println("::title:"+titleString);
            page = page.substring(matcher.end());
         }
         if (titleString ==null) 
            return;
         
         Tokenizer tokenizer = new Tokenizer(page);
         Deque<String> ngram = new LinkedList<String>();
         int cnt = 0;
         while (tokenizer.hasNext()){
            String tokenString = tokenizer.next();
            ngram.addLast(tokenString);
            while (ngram.size() > n) {
               ngram.removeFirst();
            }
            if(ngram.size()==n ) {
               StringBuilder sb = new StringBuilder();
               Iterator<String> i = ngram.iterator();
               while(i.hasNext()) {
                  sb.append(i.next());
                  sb.append(",");
               }
               if (queryNgramSet.contains(sb.toString())) {
                  cnt +=1;
               }
            }
         }
         if(cnt > 0) {
            IntTextWritableComparable pair = new IntTextWritableComparable(new IntWritable(cnt), new Text(titleString));
            System.out.println("==>emit:"+pair.toString());
//            context.write(pair, NullWritable.get());
            
            resultToEmit.add(pair);
            while(resultToEmit.size()>20) {
               resultToEmit.pollFirst();
            }
         }
         
      }
      
      @Override
      protected void cleanup(Context context) 
               throws IOException, InterruptedException {
         for (IntTextWritableComparable pair : resultToEmit) {
            System.out.println(pair.getFirstInt().toString()+", "+pair.getSecondText().toString());
            context.write(new IntTextWritableComparable(pair), new IntWritable(1));
         }
      
         
      }
   
   }
      // Reducer Class
   public static class NgramReducer extends
         Reducer<IntTextWritableComparable, IntWritable, IntWritable, Text> {

      public static int writeCnt = 0;

      @Override
      public void reduce(IntTextWritableComparable key, Iterable<IntWritable> values, Context context) 
               throws IOException, InterruptedException {
         System.out.println("-----receive:"+key.toString());
         while(writeCnt < 20) {
            IntWritable count = key.getFirstInt();
            Text  title = key.getSecondText();
            context.write(count, title);
            writeCnt += 1;
         }
      
      }
      
   }
   
   
}
