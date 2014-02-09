#!/bin/sh

HADOOP_HOME="/usr/local/hadoop-1.2.1"

# Clean up the directory
find . -name '*.class' -print0 | xargs -0 rm -f
mkdir -p class_dir

# Compile the program
find . -name '*.java' -and -not -name '.*' -print0 | xargs -0 javac -cp "${HADOOP_HOME}/hadoop-core-1.2.1.jar" -d class_dir

jar -cvf ngram.jar -C class_dir .

hadoop jar ngram.jar Ngram 4 query1.txt ./wikipedia/1gb output
cat output/part-*
