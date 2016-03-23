import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkConf
import sys.process._

import java.io._
import java.nio.file.{Paths, Files}
import java.net._
import java.util.Calendar

import scala.sys.process.Process
import scala.io.Source
import scala.collection.JavaConversions._
import scala.util.Sorting._

import net.sf.samtools._

import tudelft.utils.ChromosomeRange
import tudelft.utils.DictParser
import tudelft.utils.Configuration
import tudelft.utils.SAMRecordIterator

import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.spark.storage.StorageLevel._
import org.apache.spark.HashPartitioner

import collection.mutable.HashMap

object DNASeqAnalyzer 
{
final val MemString = "-Xmx14336m" 
final val RefFileName = "ucsc.hg19.fasta"
final val SnpFileName = "dbsnp_138.hg19.vcf"
final val ExomeFileName = "gcat_set_025.bed"
//////////////////////////////////////////////////////////////////////////////
def bwaRun (x: String, config: Configuration) : 
	Array[(Int, SAMRecord)] = 
{
	val refFolder = config.getRefFolder
	
	// Create the command string (bwa mem...)and then execute it using the Scala's process package. More help about 
	//	Scala's process package can be found at http://www.scala-lang.org/api/current/index.html#scala.sys.process.package. 
	
	// bwa mem refFolder/RefFileName -p -t numOfThreads fastqChunk > outFileName
	(Process(config.getToolsFolder++"bwa mem "++refFolder++RefFileName++" -p -t "++config.getNumThreads++" "++config.getInputFolder++x) #> new File(config.getTmpFolder++x)).!
	
	val bwaKeyValues = new BWAKeyValues(config.getTmpFolder++x)
	bwaKeyValues.parseSam()
	val kvPairs: Array[(Int, SAMRecord)] = bwaKeyValues.getKeyValuePairs()
	
	// Delete the temporary files
	Process("rm "++config.getTmpFolder++x).!
	
	return kvPairs
}
	 
def writeToBAM(fileName: String, samRecordsSorted: Array[SAMRecord], config: Configuration) : ChromosomeRange = 
{
	val header = new SAMFileHeader()
	header.setSequenceDictionary(config.getDict())
	val outHeader = header.clone()
	outHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
	val factory = new SAMFileWriterFactory();
	val writer = factory.makeBAMWriter(outHeader, true, new File(fileName));
	
	val r = new ChromosomeRange()
	val input = new SAMRecordIterator(samRecordsSorted, header, r)
	while(input.hasNext()) 
	{
		val sam = input.next()
		writer.addAlignment(sam);
	}
	writer.close();
	
	return r
}

def sort(xs: Array[SAMRecord]):Array[SAMRecord] = 
{
	if(xs.length <= 1)
		xs;
    else
	{
		val pivot = xs(xs.length/2);
		Array.concat(sort(xs filter (a=>a.compare(pivot)<0)), xs filter (a=>a.compare(pivot)==0), sort(xs filter (a=>a.compare(pivot)>0)))
    }
}

def sort(xs: Array[(Integer, String)]):Array[(Integer, String)] = 
{
	if(xs.length <= 1)
		xs;
    else
	{
		val pivot = xs(xs.length/2);
		Array.concat(sort(xs filter (a=>a._1<pivot._1)), xs filter (a=>a._1==pivot._1), sort(xs filter (a=>a._1>pivot._1)))
    }
}

def loadBalancer (chunk: Array[(Int,Int)], numPartition: Int) :
	(Array[(Int,Int)]) = //best deviation and the distribution to achieve the best result
{
	var deviation = Array(0)
	var distribution = chunk
	val avg = chunk.map(a=>a._2).sum/numPartition
	
	for (i<-2 to numPartition)
		deviation = deviation ++ Array(0)
	
	for (i<-0 to chunk.length-1)
	{
		distribution = distribution.updated(chunk.apply(i)._1-1,(distribution.apply(chunk.apply(i)._1-1)._1,deviation.indexOf(deviation.min)))
		deviation = deviation.updated(deviation.indexOf(deviation.min), deviation.apply(deviation.indexOf(deviation.min)) + chunk.apply(i)._2)
	}
	return distribution
}

def variantCall (chrRegion: Int, samRecordsSorted: Array[SAMRecord], config: Configuration) : 
	Array[(Integer, (Integer, String))] = 
{	
	val tmpFolder = config.getTmpFolder
	val toolsFolder = config.getToolsFolder
	val refFolder = config.getRefFolder
	val numOfThreads = config.getNumThreads
	
	val X=chrRegion.toString
	
	// Following is shown how each tool is called. Replace the X in regionX with the chromosome region number (chrRegion). 
	// 	You would have to create the command strings (for running jar files) and then execute them using the Scala's process package. More 
	// 	help about Scala's process package can be found at http://www.scala-lang.org/api/current/index.html#scala.sys.process.package.
	//	Note that MemString here is -Xmx14336m, and already defined as a constant variable above, and so are reference files' names.
	
	// SAM records should be sorted by this point
	// val chrRange = writeToBAM(tmpFolder/regionX-p1.bam, samRecordsSorted, config)
	val chrRange = writeToBAM(tmpFolder++"region"++X++"-p1.bam", samRecordsSorted, config)
	
	// Picard preprocessing
	//	java MemString -jar toolsFolder/CleanSam.jar INPUT=tmpFolder/regionX-p1.bam OUTPUT=tmpFolder/regionX-p2.bam
	Process("java "++MemString++" -jar "++toolsFolder++"CleanSam.jar INPUT="++tmpFolder++"region"++X++"-p1.bam OUTPUT="++tmpFolder++"region"++X++"-p2.bam").!
	//	java MemString -jar toolsFolder/MarkDuplicates.jar INPUT=tmpFolder/regionX-p2.bam OUTPUT=tmpFolder/regionX-p3.bam 
	//		METRICS_FILE=tmpFolder/regionX-p3-metrics.txt
	Process("java "++MemString++" -jar "++toolsFolder++"MarkDuplicates.jar INPUT="++tmpFolder++"region"++X++"-p2.bam OUTPUT="++tmpFolder++"region"++X++"-p3.bam METRICS_FILE="++tmpFolder++"region"++X++"-p3-metrics.txt").!
	//	java MemString -jar toolsFolder/AddOrReplaceReadGroups.jar INPUT=tmpFolder/regionX-p3.bam OUTPUT=tmpFolder/regionX.bam 
	//		RGID=GROUP1 RGLB=LIB1 RGPL=ILLUMINA RGPU=UNIT1 RGSM=SAMPLE1
	Process("java "++MemString++" -jar "++toolsFolder++"AddOrReplaceReadGroups.jar INPUT="++tmpFolder++"region"++X++"-p3.bam OUTPUT="++tmpFolder++"region"++X++".bam RGID=GROUP1 RGLB=LIB1 RGPL=ILLUMINA RGPU=UNIT1 RGSM=SAMPLE1").!
	// 	java MemString -jar toolsFolder/BuildBamIndex.jar INPUT=tmpFolder/regionX.bam
	Process("java "++MemString++" -jar "++toolsFolder++"BuildBamIndex.jar INPUT="++tmpFolder++"region"++X++".bam").!
	//	delete tmpFolder/regionX-p1.bam, tmpFolder/regionX-p2.bam, tmpFolder/regionX-p3.bam and tmpFolder/regionX-p3-metrics.txt
	Process("rm "++tmpFolder++"region"++X++"-p1.bam "++tmpFolder++"region"++X++"-p2.bam "++tmpFolder++"region"++X++"-p3.bam "++tmpFolder++"region"++X++"-p3-metrics.txt").!
	
	// Make region file 
	//	val tmpBed = new File(tmpFolder/tmpX.bed)
	val tmpBed = new File(tmpFolder++"tmp"++X++".bed")
	//	chrRange.writeToBedRegionFile(tmpBed.getAbsolutePath())
	chrRange.writeToBedRegionFile(tmpBed.getAbsolutePath())
	//	toolsFolder/bedtools intersect -a refFolder/ExomeFileName -b tmpFolder/tmpX.bed -header > tmpFolder/bedX.bed
	(Process(toolsFolder++"bedtools intersect -a "++refFolder++ExomeFileName++" -b "++tmpFolder++"tmp"++X+".bed -header") #> new File(tmpFolder++"bed"++X++".bed")).!
	//	delete tmpFolder/tmpX.bed
	Process("rm "++tmpFolder++"tmp"++X++".bed").!
	
	// Indel Realignment 
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt numOfThreads -R refFolder/RefFileName 
	//		-I tmpFolder/regionX.bam -o tmpFolder/regionX.intervals -L tmpFolder/bedX.bed
	Process("java "++MemString++" -jar "++toolsFolder++"GenomeAnalysisTK.jar -T RealignerTargetCreator -nt "++numOfThreads++" -R "++refFolder++RefFileName++" -I "++tmpFolder++"region"++X++".bam -o "++tmpFolder++"region"++X++".intervals -L "++tmpFolder++"bed"++X++".bed").!
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T IndelRealigner -R refFolder/RefFileName -I tmpFolder/regionX.bam 
	//		-targetIntervals tmpFolder/regionX.intervals -o tmpFolder/regionX-2.bam -L tmpFolder/bedX.bed
	Process("java "++MemString++" -jar "++toolsFolder++"GenomeAnalysisTK.jar -T IndelRealigner -R "++refFolder++RefFileName++" -I "++tmpFolder++"region"++X++".bam -targetIntervals "++tmpFolder++"region"++X++".intervals -o "++tmpFolder++"region"++X++"-2.bam -L "++tmpFolder++"bed"++X++".bed").!
	//	delete tmpFolder/regionX.bam, tmpFolder/regionX.bai, tmpFolder/regionX.intervals
	Process("rm "++tmpFolder++"region"++X++".bam "++tmpFolder++"region"++X++".bai "++tmpFolder++"region"++X++".intervals").!
	
	// Base quality recalibration 
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T BaseRecalibrator -nct numOfThreads -R refFolder/RefFileName -I 
	//		tmpFolder/regionX-2.bam -o tmpFolder/regionX.table -L tmpFolder/bedX.bed --disable_auto_index_creation_and_locking_when_reading_rods 
	//		-knownSites refFolder/SnpFileName
	Process("java "++MemString++" -jar "++toolsFolder++"GenomeAnalysisTK.jar -T BaseRecalibrator -nct "++numOfThreads++" -R "++refFolder++RefFileName++" -I "++tmpFolder++"region"++X++"-2.bam -o "++tmpFolder++"region"++X++".table -L "++tmpFolder++"/bed"++X++".bed --disable_auto_index_creation_and_locking_when_reading_rods -knownSites "++refFolder++SnpFileName).!
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T PrintReads -R refFolder/RefFileName -I 
	//		tmpFolder/regionX-2.bam -o tmpFolder/regionX-3.bam -BSQR tmpFolder/regionX.table -L tmpFolder/bedX.bed 
	Process("java "++MemString++" -jar "++toolsFolder++"GenomeAnalysisTK.jar -T PrintReads -R "++refFolder++RefFileName++" -I "++tmpFolder++"region"++X++"-2.bam -o "++tmpFolder++"region"++X++"-3.bam -BQSR "++tmpFolder++"region"++X++".table -L "++tmpFolder++"bed"++X++".bed").!
	// delete tmpFolder/regionX-2.bam, tmpFolder/regionX-2.bai, tmpFolder/regionX.table
	Process("rm "++tmpFolder++"region"++X++"-2.bam "++tmpFolder++"region"++X++"-2.bai "++tmpFolder++"region"++X++".table").!
	
	// Haplotype -> Uses the region bed file
	// java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T HaplotypeCaller -nct numOfThreads -R refFolder/RefFileName -I 
	//		tmpFolder/regionX-3.bam -o tmpFolder/regionX.vcf  -stand_call_conf 30.0 -stand_emit_conf 30.0 -L tmpFolder/bedX.bed 
	//		--no_cmdline_in_header --disable_auto_index_creation_and_locking_when_reading_rods
	Process("java "++MemString++" -jar "++toolsFolder++"GenomeAnalysisTK.jar -T HaplotypeCaller -nct "++numOfThreads++" -R "++refFolder++RefFileName++" -I "++tmpFolder++"region"++X++"-3.bam -o "++tmpFolder++"region"++X++".vcf  -stand_call_conf 30.0 -stand_emit_conf 30.0 -L "++tmpFolder++"bed"++X++".bed --no_cmdline_in_header --disable_auto_index_creation_and_locking_when_reading_rods").!
	// delete tmpFolder/regionX-3.bam, tmpFolder/regionX-3.bai, tmpFolder/bedX.bed
	Process("rm "++tmpFolder++"region"++X++"-3.bam "++tmpFolder++"region"++X++"-3.bai "++tmpFolder++"bed"++X++".bed").!
	
	// return the content of the vcf file produced by the haplotype caller.
	val vcfContent = Source.fromFile(tmpFolder++"region"++X++".vcf").getLines.toArray.filter(a=>a.apply(0)!='#').map(a=>(Int.box(if(a.split('\t').apply(0)!="chrX") a.split('\t').apply(0).drop(3).toInt else 23),(Int.box(a.split('\t').apply(1).toInt),a)))
	//	Return those in the form of <Chromsome number, <Chromosome Position, line>>
	return vcfContent
}

def main(args: Array[String]) 
{
	val config = new Configuration()
	config.initialize()
		 
	val conf = new SparkConf().setAppName("DNASeqAnalyzer")
	conf.set("spark.local.dir", "/data/bohaozhang/tmp")
	// For local mode, include the following two lines
	if(args(0)=="local")
	{
		conf.setMaster("local[" + config.getNumInstances() + "]")
		conf.set("spark.cores.max", config.getNumInstances())
		println("Running in local mode!")
	}
	
	// For cluster mode, include the following commented line
	if(args(0)=="cluster")
	{
		conf.set("spark.shuffle.blockTransferService", "nio") 
		println("Running in cluster mode!")
	}
	
	val sc = new SparkContext(conf)
	
	// Rest of the code goes here
	val input=new File(config.getInputFolder);
	val fileList=input.listFiles.filter(_.isFile).map(a=>a.getName);
	
	val tempDir = new File(config.getTmpFolder);
		if(!tempDir.isDirectory)
			tempDir.mkdir;

	val bwaResults = sc.parallelize(fileList).flatMap(a=>bwaRun(a,config))
	
	val bwaResultsNum = bwaResults.groupByKey.map(a=>(a._1,a._2.toArray.length)).collect.sortBy(a=>0-a._2)
	
	val balancedDistribution = loadBalancer(bwaResultsNum, config.getNumInstances().toInt).sortBy(a=>a._1)

	val balancedBWAResults = bwaResults.map(a=>(balancedDistribution.apply(a._1-1)._2,a._2)).groupByKey
	
	val variantResult = (sort(balancedBWAResults.flatMap(a=>variantCall(a._1,sort(a._2.toArray),config)).groupByKey.map(a=>(a._1,sort(a._2.toArray).map(a=>a._2).mkString("\n"))).collect)).map(a=>a._2).mkString("\n")
	
	//Finally, output process goes here
	val outputDir = new File(config.getOutputFolder);
		if(!outputDir.isDirectory)
			outputDir.mkdir;
			
	val finalResult = new PrintWriter(new File(config.getOutputFolder++"FinalResult.vcf"))
	finalResult.println("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ls")
	finalResult.print  (variantResult)
	finalResult.close
	
}
//////////////////////////////////////////////////////////////////////////////
} // End of Class definition
