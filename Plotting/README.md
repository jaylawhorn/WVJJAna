First step is to `hadd` relevant files and update the text files, for example `files2016.txt`. The utility `RemoveDuplicateEvents.C` removes duplicated events from the combination of multiple datasets.

Next, `DoDrawing.C` defines the histogram contents and binning, and further selection requirements. (Should do same job as `SignalExtraction/MakeCards.C`, which is a bit dangerous.) 

Finally, `MakePlots.C` handles actually making the plots and their style, etc.

