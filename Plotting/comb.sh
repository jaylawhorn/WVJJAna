#!/bin/bash

echo "combine -M Significance wjj_XX.txt -t -1 --expectSignal=1" > results_XX.txt
combine -M Significance wjj_XX.txt -t -1 --expectSignal=1 >> results_XX.txt
echo "combine -M Significance wv_XX.txt -t -1 --expectSignal=1" >> results_XX.txt
combine -M Significance wv_XX.txt -t -1 --expectSignal=1 >> results_XX.txt

echo "combine -M Significance zjj_XX.txt -t -1 --expectSignal=1" >> results_XX.txt
combine -M Significance zjj_XX.txt -t -1 --expectSignal=1 >> results_XX.txt
echo "combine -M Significance zv_XX.txt -t -1 --expectSignal=1" >> results_XX.txt
combine -M Significance zv_XX.txt -t -1 --expectSignal=1 >> results_XX.txt

combineCards.py mWjj=wjj_XX.txt mWV=wv_XX.txt > comb_W_XX.txt 
echo "combine -M Significance comb_W_XX.txt -t -1 --expectSignal=1" >> results_XX.txt
combine -M Significance comb_W_XX.txt -t -1 --expectSignal=1 >> results_XX.txt

combineCards.py mZjj=zjj_XX.txt mZV=zv_XX.txt > comb_Z_XX.txt 
echo "combine -M Significance comb_Z_XX.txt -t -1 --expectSignal=1" >> results_XX.txt
combine -M Significance comb_Z_XX.txt -t -1 --expectSignal=1 >> results_XX.txt

combineCards.py mWjj=wjj_XX.txt mWV=wv_XX.txt mZjj=zjj_XX.txt mZV=zv_XX.txt > comb_all_XX.txt
combine -M Significance comb_all_XX.txt -t -1 --expectSignal=1 >> results_XX.txt

