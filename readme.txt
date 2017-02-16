
### cORF_predictio_pipeline ####

usage: cORF_predictio_pipeline.py [input data, options] 

input data:

-circ , -ci, -c      
 A file with the circles locations in the format as in example file : circ_example_dm3.txt
    
-annotation, -a, -A
Annotation file in the format as in example file : annotation_example_dm3.txt
 
-system,-sys,-Sys
full path to fasta file of the genome

-output, -o, -O
output file name.

optional:

--force_exact
Use only circRNA with junction on annotaed exons boundaries

--collapse_frames
report only one cORF for each circRNA. Choose the longest ORF.



#### Output #####

<output file name>.GARY conatain the list of cORF found and in addition:
information about the predicted sequence, strat and stop positions and annotation data.

