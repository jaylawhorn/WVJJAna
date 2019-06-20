#!/bin/bash

#mVV fit
echo "combine -M Significance wjj_0.txt -t -1 --expectSignal=1" > result_0.txt
combine -M Significance wjj_0.txt -t -1 --expectSignal=1 >> result_0.txt
echo "combine -M Significance wv_0.txt -t -1 --expectSignal=1" >> result_0.txt
combine -M Significance wv_0.txt -t -1 --expectSignal=1 >> result_0.txt
combineCards.py mWjj=wjj_0.txt mWV=wv_0.txt > comb_0.txt 
echo "combine -M Significance comb_0.txt -t -1 --expectSignal=1" >> result_0.txt
combine -M Significance comb_0.txt -t -1 --expectSignal=1 >> result_0.txt

#mVV fit with dEta and mVBF categories
echo "combine -M Significance wjj_1.txt -t -1 --expectSignal=1" > result_1.txt
combine -M Significance wjj_1.txt -t -1 --expectSignal=1 >> result_1.txt
echo "combine -M Significance wv_1.txt -t -1 --expectSignal=1" >> result_1.txt
combine -M Significance wv_1.txt -t -1 --expectSignal=1 >> result_1.txt
combineCards.py mWjj=wjj_1.txt mWV=wv_1.txt > comb_1.txt 
echo "combine -M Significance comb_1.txt -t -1 --expectSignal=1" >> result_1.txt
combine -M Significance comb_1.txt -t -1 --expectSignal=1 >> result_1.txt