#!/bin/bash

rm -rf *.class
find . -name '*.java' -print0 | xargs -0 javac
java -cp . ChatServer
