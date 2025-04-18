# Advanced_Bioinformatics_Assignment
1.1What does  ./../.. stand for ?
A.  Current directory
B.  Up one directory
C.  Up two directories
D.  None of Above

C

1.2What does cd / mean in UNIX? Please explain what the cd command does.

cd/ means we change directory to the root directory.
cd command helps us change directories in UNIX.

1.3What command would you use to get help about the command cp? (please provide an example command)

man cp

1.4 What does the command pwd do?

The command pwd prints the current working directory.

1.5 How do you display a listing of file details such as date, size, and access permissions in a given directory? (please provide an example command)

ls -l

1.6 How do you print on the terminal the first 15 lines of all files ending by .txt? (please provide an example command)

head -n 15 *.txt

1.7 How do you rename a file from new to old? (please provide an example command)

mv new old

1.8 How do you display the contents of a file myfile.txt? (please provide an example command)

less myfile.txt

1.9 How do you create a new directory called flower? (please provide an example command)

mkdir flower

1.10 How do you change the current directory to /usr/local/bin? (please provide an example command)

cd /usr/local/bin

1.11 How can you display a list of all files in the current directory, including the hidden files? (please provide an example command)

ls -a

1.12 What command do you have to use to go to the parent directory? (please provide an example command)

cd ..

1.13 Which command would you use to create a sub-directory in your home directory?  (please provide an example)

mkdir ~/directory

1.14 Which command would you use to list the first lines in a text file? (please provide an example)

head -n 1 textfile

1.15 Which command will display the last lines of the text file file1? (please provide an example)

tail -n 1 file1

1.16 Which command is used to extract a column from a text file? (please provide an example)

cut -f 1 textfile

1.17 How do you copy an entire directory structure? E.g. from Project to Project.backup (please provide an example)

cp -r Project Project.backup

1.18 How would you search for the string Hypertension at the end of the line in a file called diseases.txt? (please provide an example)

grep -w Hypertension diseases.txt | tail -n 1

1.19 How do you see hidden files in your home directory? (please provide an example)

ls -d ~/.*

1.20 How do you run a job that will continue running even if you are logged out? (please provide an example)

nohup command &

2.6 Using an alternative tool

Provide below the new commands used to run the alternative tool and comment on your choice of options and how and if using this tool would affect the results.

#Generate the index files

sudo apt install bowtie2

cd ~/7BBG2016/data/reference

bowtie2-build -f hg19.fa hg19_index

#Replace the aligner with Bowtie2

bowtie2 -x hg19_index -1 ~/7BBG2016/data/trimmed_fastq/NGS0001_trimmed_R_1P -2 ~/7BBG2016/data/trimmed_fastq/NGS0001_trimmed_R_2P -S ~/7BBG2016/data/aligned_data/NGS0001.sam

Bowtie2 is more suitable for aligning short sequences, with a relatively high operating speed.

The results of aligning using BWA Mem may be more accurate.
