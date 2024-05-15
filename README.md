# CRISPR_Screen_QA-QC
This project will determine correlation coefficients and slopes of regression lines between a pre-analyzed, publicly available DrugZ or MaGeCK output and a freshly generated counterpart.
Forewarning: All files selected MUST be csv files. Maybe future updates will make this file type more malleable. Only time will tell :)

The input is chosen through a GUI. The input options will require a csv file containing a list of authors in a column called "First_Last" with the name of each concatenated together by an underscore. For example, "Heer_Bindra." Then it will require you to select a directory that has folders named after each author pair, another folder denoted by screen identifiers such as Screen#_Drug_Cellline (i.e. 1_Cisplatin_U2OS) which will then have files within the directory named the same but ending in "_pDZ", "_gDZ", "_pMGK", and/or "_gMKG". Lastly, you will choose an output directory to write the output csv file with the calculated values.
