#!/bin/bash

echo "combine -M Significance wjj_XX.txt -t -1 --expectSignal=1" > results_XX.txt
combine -M Significance wjj_XX.txt -t -1 --expectSignal=1 >> results_XX.txt
echo "combine -M Significance wv_XX.txt -t -1 --expectSignal=1" >> results_XX.txt
combine -M Significance wv_XX.txt -t -1 --expectSignal=1 >> results_XX.txt

echo "combine -M Significance zjj_XX.txt -t -1 --expectSignal=1" >> results_XX.txt
combine -M Significance zjj_XX.txt -t -1 --expectSignal=1 >> results_XX.txt
echo "combine -M Significance zv_XX.txt -t -1 --expectSignal=1" >> results_XX.txt
combine -M Significance zv_XX.txt -t -1 --expectSignal=1 >> results_XX.txt

#combineCards.py mWjj=wjj_2017.txt mWV=wv_2017.txt > comb_W_2017.txt 
#echo "combine -M Significance comb_W_2017.txt -t -1 --expectSignal=1" >> results_2017.txt
#combine -M Significance comb_W_2017.txt -t -1 --expectSignal=1 >> results_2017.txt
#
#combineCards.py mZjj=zjj_2017.txt mZV=zv_2017.txt > comb_Z_2017.txt 
#echo "combine -M Significance comb_Z_2017.txt -t -1 --expectSignal=1" >> results_2017.txt
#combine -M Significance comb_Z_2017.txt -t -1 --expectSignal=1 >> results_2017.txt
#
#combineCards.py mWjj=wjj_2017.txt mWV=wv_2017.txt mZjj=zjj_2017.txt mZV=zv_2017.txt > comb_all_2017.txt
#combine -M Significance comb_all_2017.txt -t -1 --expectSignal=1 >> results_2017.txt

