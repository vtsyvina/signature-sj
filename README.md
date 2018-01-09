# SignatureSJ
SignatureSJ is a tool for estimation of genetic relatedness between members of heterogeneous populations of closely related genomic variants. As input it takes fasta file and creates an edges list of all related sequence pairs as output.

It works for two metrics:
- Edit Distance (Levenshtein Distance)
- Hamming Distance

It has two input types:
- Single sample - tool will find all related sequence pairs within one sample
- Multi-sample - tool will find all related sequence pairs between all given samples such that one sequence will be from one sample and second will be from another sample (tool goes through all possible sample pairs)

# How to run
## How to Run
1) Copy repository to your local machine

``` git clone git@github.com:vyacheslav-tsivina/sequence-comparison.git```

2) Build jar file with Maven
  ```mvn clean install```
  It will create sequence-comparison.jar in the project root.

**OR**

Download jar from <a href="https://drive.google.com/open?id=1NsSTtiq3Fi-69ONS3cip6Um5LKw-hVlx">here</a> (I always update it with the latest version)

## Parameters
There are several available parameters:
- ``-m`` mandatory parameter to specify a method that you want to run. There are four possible values:
  - 'edit-single' - find all related pairs of sequences in a single sample for edit distance
  - 'hamming-single' - find all related pairs of sequences in a single sample for Hamming distance
  - 'edit-multi' - find all related pairs of sequences between all pairs of given samples for edit distance
  - 'hamming-multi' - find all related pairs of sequences between all pairs of given samples for Hamming distance
- ``-in`` the input path. If not specified default ``cleaned_independent_264/AMC_P01_1b.fas`` file will be used.
  It can be relative as well as the absolute path. For single-sample it can be either file or folder with files. If it's a folder, the tool will read all files in the given folder and concatenate them into one sample. For multi-sample version one should give a folder and tool will consider each file as a separated sample.
 - ``k`` a threshold for related sequences, so in output will be only sequences (S, Q) such that d(S, Q) <= k. 10 is a default value
 - ``l`` a length for mers that are used to create signatures. 11 is the default. Not used for Hamming distance (entropy based signature length is used there)
 - ``-outDir`` an output directory. output/ is a default value.
 - ``-threads`` number of threads for parallel  execution. By default number of available cores will be used.
 
 ## Input
 As input the program takes simple fasta format as follows:
 ```
  >read1
  CACCTACAGCAGCCCTAGTGGTATCGCAGTTACTCCGGATCCCACAAGCTGTCGTGGATATGGTGGCGGG
  
  >read2
  CACCTACGACAGCCCTAGTGGTGTCGCAGTTACTCCGGATCCCACAAGCCGTCGTGGATATGGTGGCGGG
  
  ...
  ```
  For single-sample if a folder is specified as input, the program will concatenate all files from the folder in one sample.
  
  ## Output
  The Output is a simple text file that has sample name(file name)+method+output.txt as a name (e.g., ``db6-signature-hamming-output.txt``).
  It contains a set of number pairs, where each number represents a sequence number in the input file. It's done in that way since even in that way output may have a size of several Gb and sequence names from input file will significantely increase this volume.
  If a folder is specified as input for single sample method, the tool will read files in alphabetical order to escape ambiguity.
  
  ## Any questions
For any questions, please, contact: vyacheslav.tsivina@gmail.com
