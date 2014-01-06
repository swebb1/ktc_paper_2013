
java -classpath /opt/software/Trimmomatic/trimmomatic-0.22.jar org.usadellab.trimmomatic.TrimmomaticSE -phred33 -threads 14 -trimlog $2 $1 $3.fastq ILLUMINACLIP:/homes/swebb/tool_files/adapters_illumina_genomic_w_revcomps.fasta:2:20:10 SLIDINGWINDOW:5:20 HEADCROP:4 MINLEN:20

