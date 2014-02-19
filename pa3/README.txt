Song Han  (songhan)
Artem Vasilye (tema8)


[Preprocessing]
We split the input file by <page> and </page> by the class XmlInputFormat which
extends TextInputFormat. Then we wrote a method setupQuery to construct a list
of Ngrams in the query document. Since the query document is pretty small, these
Ngrams could be easily fit into memory.


[Mapper]
The input of the mapper is the <LongWritable id, Text document> pairs. Only the
document text is useful. Then we scan through the document to match each Ngram
to the query list. If there's a hit, we increase the score by one. When the
whole document has been scanned and the score is not zero, we write out the
output: <LongWritable score, Text titleString>

[Reducer]
We are using the decreasing comparator by setting:
job.setSortComparatorClass(LongWritable.DecreasingComparator.class);
So after shuffling the reducer is guaranteed to receive the score in reverse
order. We need to further sort the titleString, which is done within the reducer.
The output of the reducer is <LongWritable score, Text titleString>.
One alternative could be creating a user defined class IntText which implements
WritableComparable and use the value to key inversion and secondary sort.



[Scalability]
Our implementation performs a single pass of Map-Reduce, during which each Mapper
builds NGram from the input query file, so the time it takes is proportional to
the input size q and since all mappers are doing it only once and in parallel,
this time is independent from the number of processor. The input pages from
Wikipedia are divided across P Mappers by Hadoop. Assuming that all P Mappers
are running at the same time, it will take O(n/P) time to process the entire set.
Each Mapper will produce a list of up to 20 best, so we'll get P such lists to
be reduced. If we have P processors, the reduction of P lists can be done using
a tree structure ( like in merge sort ) in O( log P), because the height of
binary tree is log P. All this stages are done one after another, so the resulting
time for the algorithm is the sum of 3 stages or: O(q + n/P + logP)

[Extra credit]
We implemented the first extra credit to output top 20 matches
