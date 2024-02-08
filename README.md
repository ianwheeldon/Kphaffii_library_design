<h1>Kphaffii_library_design</h1>

<h2>1.0 Publication</h2>

Functional genomic screening in Komagataella phaffii enabled by high-activity CRISPR-Cas9 library

<h2>2.0 Prerequisite</h2>

<h3>2.1 For library design:</h3>

python 2.7.14

scipy 1.2.2

gffutils 0.10.1

mySQL-python 1.2.3

numpy 1.16.6

pandas 0.24.2

scikit-learn 0.18.1

scipy 1.2.2


<h3>2.2 For generating sgRNA counts from raw sequencing reads:</h3>

cutadapt 4.2

python 3.9.18


<h2>3.0 Running the code</h2>

<h3>3.1 Running CHOPCHOP v3</h3>

The first step in designing the CRISPR-Cas9 sgRNA library using our codes is to be able to run CHOPCHOP v3 on the command line (https://doi.org/10.1093/nar/gkz365). All the necessary requirements and step-by-step installation guides for CHOPCHOP v3 are explained by its authors here: https://bitbucket.org/valenlab/chopchop/src/master/ . After installation, you should be able to run the examples mentioned on their website properly.

To use our source code to design a genome-wide sgRNA library, clone this repository in the <code>chopchop/</code> directory: 

<code>cd chopchop/</code>

<code>git clone https://github.com/ianwheeldon/Kphaffii_library_design.git/</code>


<h3>3.2 Run example on Komagataella phaffii GS115 genomic data</h3>

First, copy the files in the <code>Library_design/</code> directory in the <code>chopchop/</code> directory:  

<code>cp -r Kphaffii_library_design/Library_design/* .</code>

<code>cp -r example/* ./</code>

Run this command:

<code>./chopchop.py -T 1 -M NGG --maxMismatches 3 -g 20 -G GS115 -o Results -Target CP014715.1:5164-5200 --scoringMethod ALL</code>

You should be able to have this same result:

<code>Rank   Target sequence Genomic location        Strand  GC content (%)  Self-complementarity    MM0     MM1     MM2     MM3     XU_2015 DOENCH_2014     DOENCH_2016     MORENO_MATEOS_2015        CHARI_2015      G_20    ALKAN_2018      ZHANG_2019</code>

<code>1       GAAGCTTTCAGGAGACCTCTTGG CP014715.1:5149 +       50      0       0       0       0       0       0.40    0.13    44.35   0.56    0.00    0.00    24.87   37.09</code>

<code>2       TCACAAAATCGGACTCCAAGAGG CP014715.1:5164 -       45      0       0       0       0       0       0.59    0.65    62.70   0.55    0.00    1.00    23.84   58.97</code>

<code>3       GAATGATGAATTCACAAAATCGG CP014715.1:5175 -       25      0       0       0       0       0       0.47    0.14    45.68   0.61    0.00    0.00    13.66   17.49</code>


<h3>3.3 Designing the genome-wide CRISPR-Cas9 sgRNA library</h3>

Running the <code>Library_design.py</code> code would start designing the library similar to <code>BEST_LIBRARY.csv</code> and <code>CHOPCHOP_Total.csv</code> files in the example directory:

<code>python Library_design.py</code>


<h2>4.0 Generating sgRNA read counts from raw reads</h2>

Navigate to the count_generation directory: 
<code>cd Kphaffii_library_design/count_generation/</code>

Firstly, Next Generation Sequencing reads should be trimmed from barcodes and 3’ and 5’ redundant sequences:

<code>cutadapt -g ^NTAGT --discard-untrimmed --error-rate 0.2 --times 1 example/raw_sequencing_example.fastq.gz -o example/Results/sample1.fastq.gz;</code>

<code>cutadapt -g GTAGACCGGGGTTCAATTCCCCGTCGCGGAGC --discard-untrimmed --overlap 18 --error-rate 0.2 --times 1 example/Results/sample1.fastq.gz -o example/Results/sample1_5ptrimmed.fastq.gz;</code>

<code>cutadapt -a GTTTTAGAGCTAGAAATAGC  --discard-untrimmed --overlap 9 --error-rate 0.2 --times 1 example/Results/sample1_5ptrimmed.fastq.gz -o example/Results/sample1_5ptrimmed_3ptrimmed.fastq.gz;</code>

Next, run NEMandIEM.py to generate a <code>sample1_ExactMatched_reads.csv</code> file with the number of exact matched counts of the sgRNAs in the library as well as an additional <code>sample1_unmatched_reads.csv</code> file of the reads that were not matched to the library.

<code>python NEMandIEM.py</code>

The final step in generating the counts is to use Bowtie to align the reads that were not exact matched to the sgRNA library. This data will be stored in the <code>Final_sgRNA_count.csv</code>:
<code>python sgRNA_count.py</code>



