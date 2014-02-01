#!/bin/sh

# Clean up the directory
find . -name '*.class' -print0 | xargs -0 rm -f

# Compile the program
find . -name '*.java' -and -not -name '.*' -print0 | xargs -0 javac -cp deuceAgent-1.3.0.jar

# Run the program
echo 'Checking correctness...'
java -javaagent:deuceAgent-1.3.0.jar -cp . Driver --test=correctness --key-range=1000000 --read-percent=5
java -javaagent:deuceAgent-1.3.0.jar -cp . Driver --test=correctness --key-range=10000 --read-percent=95
java -javaagent:deuceAgent-1.3.0.jar -cp . Driver --test=correctness --key-range=1000 --read-percent=99

echo
echo 'Checking performance...'
java -javaagent:deuceAgent-1.3.0.jar -cp . Driver --test=performance --key-range=1000000 --read-percent=5
java -javaagent:deuceAgent-1.3.0.jar -cp . Driver --test=performance --key-range=10000 --read-percent=95
java -javaagent:deuceAgent-1.3.0.jar -cp . Driver --test=performance --key-range=1000 --read-percent=99
