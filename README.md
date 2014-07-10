This project contains two approaches for SV classification: rule-based and machine-learning based.

They work well with paired result. ML approach can support unpaired result but with relatively low accuracy. Normally, I use paired result for training and paired+unpaired result for testing. (Since there are few "WOFS" in paired result, some entries of this class should be copied from testing file into training file).
They support SVs within single chromosome. Cross-chromosomal events are not supported. 

The data used are generated by the scripts in the folder of simulate_script. 




Rule-based:
svdetect.py is the main file for rule-based classification.
Usage: ./svdetect.py Socrates_paired_result.txt corresponding_breakpoints.txt

The comparison is done between the result of svdetect.py and the breakpoints.txt which is generated during data simulating.
Check index.html for the detail of lastest comparison.

example:
Mos-MacBook-Pro-2:SVdetect HELE$ python ./svdetect.py ./simulated_data/set1/results_Socrates_paired_simulated_chr12_1_simseq_s_long_sc_l25_q5_m5_i95.txt ./simulated_data/set1/simulated_chr12_1.fa_breakpoints.txt





ML-approach:
script1_training_data.sh and script1_test_data.sh are for generating training data and test data respectively. 
script1_training_data.sh uses paired Socrates result only. 
script1_test_data.sh uses both paired and unpired result.
If wanna used paired result for testing, just use script1_training_data.sh to generate the test data. There is no difference between these two scripts expect script1_test_data.sh utilizes unpaired result of Scorates.

The output of these scripts is .ARFF file which can be used by WEKA to do classification. (I used WEKA GUI to do the classification. No code for ML classification within this project)

example:
Mos-MacBook-Pro-2:SVdetect HELE$ ./script1_training_data.sh ./simulated_data/set1/ ./ML_data/set55
 
