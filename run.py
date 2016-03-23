from xml.dom import minidom
import sys
import os
import time

# numInstances and numThreads would be required for cluster mode
doc = minidom.parse("./config.xml")
numInstances = doc.getElementsByTagName("numInstances")[0].firstChild.data
numThreads = doc.getElementsByTagName("numThreads")[0].firstChild.data
	
def run():

	if len(sys.argv)==1 :  
		cmdStr = "/home/bohaozhang/spark/bin/spark-submit " + \
		"--class \"DNASeqAnalyzer\" --master local[*] --driver-memory 32g target/scala-2.10/dnaseqanalyzer_2.10-1.0.jar local"
	elif (len(sys.argv)==2) & (sys.argv[1]=="cluster") :
		cmdStr = "/home/bohaozhang/spark/bin/spark-submit " + \
		"--class \"DNASeqAnalyzer\" --master yarn-cluster --driver-memory 14g --executor-memory 14g --num-executors "+numInstances+" --executor-cores "+numThreads+" target/scala-2.10/dnaseqanalyzer_2.10-1.0.jar cluster"
	else:  
		cmdStr = "/home/bohaozhang/spark/bin/spark-submit " + \
		"--class \"DNASeqAnalyzer\" --master local[*] --driver-memory 128g target/scala-2.10/dnaseqanalyzer_2.10-1.0.jar local"
	
	print cmdStr
	os.system(cmdStr)
	
start_time = time.time()

run()
	
time_in_secs = int(time.time() - start_time)
mins = time_in_secs / 60
secs = time_in_secs % 60

print "{Time taken = " + str(mins) + " mins " + str(secs) + " secs}"
