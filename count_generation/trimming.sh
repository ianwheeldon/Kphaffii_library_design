#Demultiplexing samples using differernt barcodes
#barcodes used in this example: GTAGT
cutadapt -g ^NTAGT --discard-untrimmed --error-rate 0.2 --times 1 example/raw_sequencing_example.fastq -o example/Results/sample1.fastq;

#Removing the 5' sequences from reads:
cutadapt -g GTAGACCGGGGTTCAATTCCCCGTCGCGGAGC --discard-untrimmed --overlap 18 --error-rate 0.2 --times 1 example/Results/sample1.fastq -o example/Results/sample1_5ptrimmed.fastq;

#Removing the 3' sequences from the reads
cutadapt -a GTTTTAGAGCTAGAAATAGC  --discard-untrimmed --overlap 9 --error-rate 0.2 --times 1 example/Results/sample1_5ptrimmed.fastq -o example/Results/sample1_5ptrimmed_3ptrimmed.fastq;

