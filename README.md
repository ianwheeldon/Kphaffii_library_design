# Kphaffii_library_design
1.0 Publication

Functional genomic screening in Komagataella phaffii enabled by high-activity CRISPR-Cas9 library

2.0 Prerequisite

2.1 For library design:

python 2.7.14

scipy 1.2.2

gffutils 0.10.1

mySQL-python 1.2.3

numpy 1.16.6

pandas 0.24.2

scikit-learn 0.18.1

scipy 1.2.2


2.2 For generating sgRNA counts from raw sequencing reads: 

cutadapt 4.2

python 3.9.18


3.0 Running the code

3.1 Running CHOPCHOP v3

The first step in designing the CRISPR-Cas9 sgRNA library using our codes is to be able to run CHOPCHOP v3 on the command line (https://doi.org/10.1093/nar/gkz365). All the necessary requirements and step-by-step installation guides for CHOPCHOP v3 are explained by its authors here: https://bitbucket.org/valenlab/chopchop/src/master/ . After installation, you should be able to run the examples mentioned on their website properly.

To use our source code to design a genome-wide sgRNA library, clone this repository in the chopchop/ directory: 

cd chopchop/

git clone https://github.com/ianwheeldon/Kphaffii_library_design.git/


3.2 Run example on Komagataella phaffii GS115 genomic data

First, copy the files in the Library_design/ directory in the chopchop/ directory:  

cp -r Kphaffii_library_design/Library_design/* .

cp -r example/* ./

Run this command:

./chopchop.py -T 1 -M NGG --maxMismatches 3 -g 20 -G GS115 -o Results -Target CP014715.1:5164-5200 --scoringMethod ALL

You should be able to have this same result:

Rank   Target sequence Genomic location        Strand  GC content (%)  Self-complementarity    MM0     MM1     MM2     MM3     XU_2015 DOENCH_2014     DOENCH_2016     MORENO_MATEOS_2015        CHARI_2015      G_20    ALKAN_2018      ZHANG_2019

1       GAAGCTTTCAGGAGACCTCTTGG CP014715.1:5149 +       50      0       0       0       0       0       0.40    0.13    44.35   0.56    0.00    0.00    24.87   37.09

2       TCACAAAATCGGACTCCAAGAGG CP014715.1:5164 -       45      0       0       0       0       0       0.59    0.65    62.70   0.55    0.00    1.00    23.84   58.97

3       GAATGATGAATTCACAAAATCGG CP014715.1:5175 -       25      0       0       0       0       0       0.47    0.14    45.68   0.61    0.00    0.00    13.66   17.49


3.3 Designing the genome-wide CRISPR-Cas9 sgRNA library

Running the Library_design.py code would start designing the library similar to BEST_LIBRARY.csv and CHOPCHOP_Total.csv files in the example directory:

python Library_design.py


3.4 Generating sgRNA read counts from raw reads 

Navigate to the count_generation directory: cd Kphaffii_library_design/count_generation/

Firstly, Next Generation Sequencing reads should be trimmed from barcodes and 3’ and 5’ redundant sequences:

cutadapt -g ^NTAGT --discard-untrimmed --error-rate 0.2 --times 1 example/raw_sequencing_example.fastq -o example/Results/sample1.fastq;

cutadapt -g GTAGACCGGGGTTCAATTCCCCGTCGCGGAGC --discard-untrimmed --overlap 18 --error-rate 0.2 --times 1 example/Results/sample1.fastq -o example/Results/sample1_5ptrimmed.fastq;

cutadapt -a GTTTTAGAGCTAGAAATAGC  --discard-untrimmed --overlap 9 --error-rate 0.2 --times 1 example/Results/sample1_5ptrimmed.fastq -o example/Results/sample1_5ptrimmed_3ptrimmed.fastq;

Next, run NEMandIEM.py to generate a sample1_ExactMatched_reads.csv file with the number of exact matched counts of the sgRNAs in the library as well as an additional sample1_unmatched_reads.csv file of the reads that were not matched to the library.

python NEMandIEM.py

The final step in generating the counts is to use Bowtie to align the reads that were not exact matched to the sgRNA library. This data will be stored in the Final_sgRNA_count.csv:
python sgRNA_count.py


