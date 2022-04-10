#!/bin/zsh



for (( i=0; i<100; i=i+1 )); do
     #echo $i
     ./demo.out $(( $i + 1 )) $(( $i * 200 ))
done 

