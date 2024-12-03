# Tutorial to Generating the Reconciled Notung Gene and Species Phylogenetic Tree along with the Midpoint Rooted Maximum Likelihood Gene Tree!
### Final Repository for the gene PHKA2
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Before we start, it is absolutely critical to mention the most importance of ensuring that you are in the correct directory in Bash when you want to create any changes. To see the directory you are located in, use the command "pwd" to list the absolute path to your current directory. If you need to change directory, use the command "cd" and add the directory of interest! 





## 1. Introduction to the PHKA2 gene and acquiring its data

In lab 3, we started out by downloading the query protein PHKA2 from the National Center for Biotechnology Information (NCBI). 

To achieve this we used the following code:
```
ncbi-acc-download -F fasta -m protein "NP_000283.1"
```

Below is each part of the code explained:
* ncbi-acc-download: command-line tool that downloads sequence data from the NCBI by taking its accession number as input and will retrive the corresponding sequence. 
* -F fasta: specifies the file format for the output 
* -m protein: specifies the type of molecule (eg. protein or nucleotide) 
* "NP_000283.1": specifies the accession of the desired protein (in our case, PHKA2)

Following this, we utilized the following code to perform a BLAST search of the query protein to find potential homologs of PHKA2
```
blastp -db ../allprotein.fas -query NP_000283.1.fa -outfmt 0 -max_hsps 1 -out globins.blastp.typical.out
```
Again, below is each part of the code explained:
* blastp: specifies the tool (Basic Local Alignment Search Tool for Proteins) that will compare our query sequence to a protein database to find regions of similarity/homology. 
* -db ../allprotein.fas: refers to the protein database file located in the parent directory (.. means "go up one directory level"). 
* -query NP_000283.1.fa: specify the input query file to perform the blastp search from. 
* -outfmt 0: specifies the output format
* -max_hsps 1: will limit the number of high-scoring pair alignments, which are regions of local alignment that have the highest similarity score, found to 1.
* -out globins.blastp.typical.out: specifies the output file name where the results will be saved.
  
Next, we used the following code to create different output of the same analysis that is a bit more detailed and readable.
```
blastp -db ../allprotein.fas -query NP_000283.1.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out myprotein.blastp.detail.out
```
The function of the new parts added to the previous the code is explained below: 
* -outfmt: flag that specifies the output format
* "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle": specifies the new format of the output which is organized in a row-column format which will specify the values of each. Each word is its own column (shown in chronological order: sequence ID, percentage of identical matches, alignment length, mismatch number, gap number, e-value, alignment bit score, and full title of the matching sequence).
* -out startingprotein.blastp.detail.out: specifies the output file name where the results will be saved.

After this, we used the following command to filter the output file to have only high-scoring matches in our analyses (those with e-values less than 1e-30).
```
awk '{if ($6< 1e-30)print $1 }' myprotein.blastp.detail.out > myprotein.blastp.detail.filtered.out
```
The awk program will select particular records in a file and perform operations on them.





## 2. Generating the Multiple Sequence Alignment for PHKA2

When generating the alignment, the allproteins.fas file that contains all the proteomes for our study was already provided to us. We will use this file in order to obtain sequences in this file to perform multiple sequence alignments. The command to obtain these sequences is as follows: 
```
seqkit grep --pattern-file ~/lab03-$MYGIT/projectprotein/protein.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/myprotein/myprotein.homologs.fas
```
This command uses seqkit (toolkit that is useful in FASTA/Q file manipulation).
Also it is important to note that **~/lab03-$MYGIT/** and **~/lab04-$MYGIT/** followed by **projectprotein** or **myprotein** and refers to the directory and folder that the desired files are in. 
In addition, the output *myproteins.homologs.fas* contains the unaligned sequences which will be critical throughout for the generation of **both** figures.


Following this command, we performed the multiple sequence alignment by using the muscle program. Essentially, the multiple sequence alignment will align all the sequences with each other along their entire lengths and produce an alignment that predicts the regions of homology and of insertions or deletions with respect to the other sequences. 
```
muscle -align ~/lab04-$MYGIT/myprotein/myprotein.homologs.fas -output ~/lab04-$MYGIT/myprotein/myprotein.homologs.al.fas
```

To remove any identically labeled sequences, that may represent multiple alleles, the following command was used
```
sed 's/ //g' ~/lab04-$MYGIT/myprotein/myprotein.homologs.al.fas | seqkit grep - v -r -p "dupelabel" > ~/lab05-$MYGIT/myprotein/myprotein.homologsf.al.fas
```
This command will essentially remove sequences with identical labels that also contain the word "duplabel" that was added by the muscle program.





## 3. Generating the Optimal Phyogenetic Tree of PHKA2: the Reconciled Notung Gene and Species Phylogenetic Tree.

First  we created the species tree of the species we were using in the study. The organisms used in our study is as follows (taken from the Lab 05 GitHub):

| Organism                        | Organism Abbreviation          |
|---------------------------------|--------------------------------|
| *Carcharodon carcharias*        | C.carcharias                   |
| *Chelonia mydas*                | C.mydas                        |
| *Danio rerio*                   | D.rerio                        |
| *Equus caballus*                | E.caballus                     |
| *Felis catus*                   | F.catus                        |
| *Gallus gallus*                 | G.gallus                       |
| *Gasterosteus aculeatus*        | G.aculeatus                    |
| *Homo sapiens*                  | H.sapiens                      |
| *Salmo salar*                   | S.salar                        |
| *Sphaerodactylus townsendi*     | S.townsendi                    |
| *Xenopus laevis*                | X.laevis                       |

The species tree was created using the following command 
```
echo "((((((F.catus,E.caballus)Laurasiatheria,H.sapiens)Boreoeutheria,(S.townsendi,(C.mydas,G.gallus)Archelosauria)Sauria)Amniota,X.laevis)Tetrapoda,(D.rerio,(S.salar,G.aculeatus)Euteleosteomorpha)Clupeocephala)Euteleostomi,C.carcharias)Gnathostomata;"  > ~/lab05-$MYGIT/species.tre
```
This tree will be useful later on for the reconcilation of the species and gene trees.



Now, we use IQtree to find the maximum likelihood *unrooted* tree estimate by using the following command
```
iqtree -s ~/lab05-$MYGIT/myprotein/myprotein.homologsf.al.fas -bb 1000 -nt 2
```
This command will use the alignment we just generated and will calculate the optimal amino acid frequency and substition model, subsequently performing a tree search. This tree search will estimate the branch length and find the most optimal tree.


In order to root our current unrooted tree, we utilized midpoint rooting in order to specify the oldest divergence event. The root of this midpoint rooted tree will be in the middle of the longest branch on the tree. The command used for midpoint rooting our unrooted tree is as follows:
```
gotree reroot midpoint -i ~/lab05-$MYGIT/myprotein/myprotein.homologsf.al.fas.tre efile -o ~/lab05-$MYGIT/myprotein/myprotein.homologsf.al.mid.treefile
```
The gotree program used in the code will manipulate phylogenetic trees. 


In order to visualize the tree into a graphical display, the following command created an svg image of the tree
```
nw_order -c n ~/lab05-$MYGIT/myprotein/myprotein.homologsf.al.mid.treefile | nw_d isplay -w 1000 -b 'opacity:0' -s > ~/lab05-$MYGIT/myprotein/myprotein.homologsf.al.mid.treefile.svg - 
```

In our study, we converted this svg image into a pdf for easier viewing using the following command:
```
convert ~/lab05-$MYGIT/myprotein/myprotein.homologsf.al.mid.treefile.svg ~/lab0 5-$MYGIT/myprotein/myprotein.homologsf.al.mid.treefile.pdf
```
Shown here: [myprotein.homologsf.al.mid.treefile.pdf](https://github.com/user-attachments/files/17986885/myprotein.homologsf.al.mid.treefile.pdf)


Now that we acquired both the species and gene trees, the next step was to reconcile the trees using the software package called Notung. Notung will do this by estimating the duplication events and loss events. 
The command used to perform the reconciliation is as follows: 
```
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05 -$MYGIT/species.tre -g ~/lab06-$MYGIT/PHK2A/PHK2A.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/PHK2A/
```
This command will take both the species tree we had previously generated and reconcile it to the gene tree we just created. 


Following the reconciliation, we created a RecPhyloXML object, which allows for the representation of all the relationships in the gene and species reconciled tree. The RecPhyloXML object was created using the following command:
```
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/PHK2A/PHK2A.homologsf.al.mid.treefile.rec.ntg --include.species
```

Using this RecPhyloXML object, we used thirdkind to interpret and format the reconciled tree into a reconciliation graphic of a gene within species tree.
This was done by using the following command:
```
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/PHK2A/PHK2A.homologsf.al.mid.tr eefile.rec.ntg.xml -o ~/lab06-$MYGIT/PHK2A/PHK2A.homologsf.al.mid.treefile.rec.svg
```

Finally, we converted this svg image generated of the gene-within-species reconciled tree into a pdf using the following command:
```
convert -density 150 ~/lab06-$MYGIT/PHK2A/PHK2A.homologsf.al.mid.tree file.rec.svg ~/lab06-$MYGIT/PHK2A/PHK2A.homologsf.al.mid.treefile.rec.pdf
```
This command gave us our Figure 1 image as our output: [PHK2A.homologsf.al.mid.treefile.rec.pdf](https://github.com/user-attachments/files/17987512/PHK2A.homologsf.al.mid.treefile.rec.pdf)




## 4. Generating the Midpoint Rooted Maximum Likelihood Gene Tree of PHKA2:

In this section, we were able to generate the second figure of our study: the midpoint tooted maximum likelihood gene tree of our gene of interest PHKA2.
To start, we used RPS-BLAST in order to identify the Pfam domains within the PHKA2 sequence. The Pfam database is used to predict the domains and critical sites by analyzing the functionality of proteins and classifying them into their respective families. 

In order to run the rps-blast to identify these Pfam domains, the following command was utilized:
```
rpsblast -query ~/lab08-$MYGIT/PHK2A/PHK2A.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/PHK2A/PHK2A.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
```
In this command it is important to note that the unaligned sequence generated in the beginning of Part 2, *myprotein.homologs.fas* was be used to run this analysis. Also, in this command, there was an e-value of 0.0000000001 was specified as an output measure for the Pfam domains.


Next, we plotted the predicted Pfam domains on our final midpoint rooted gene tree *myprotein.homologsf.al.mid.treefile* by using R script to plot the predicted the pfam domains and next to their respective proteins on the tree. The command to run the code is below:
```
Rscript --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/PHK2A/myprotein.homologsf.al.mid.treefile ~/lab08-$MYGIT/PHK2A/PHK2A.rps-blast.out ~/lab08-$MYGIT/PHK2A/PHK2A.tree.rps.pdf 
```
This code uses R script. The functions of the code is listed below:
* --vanilla: tells R not save or restore a previous R setting.
* the script that Dr. Rest wrote called plotTreeAndDomains.r
* the threee files (tree file, rps-blast output and output filename, all in order).

The output generated is shown here and was used the second figure of our study: [PHK2A.tree.rps.pdf](https://github.com/user-attachments/files/17987463/PHK2A.tree.rps.pdf)




In order to look at our predicted Pfam domains in more detail, we used the following code to output the information in a nice spreadhseet format:
```
mlr --inidx --ifs "\t" --opprint cat ~/lab08-$MYGIT/PHK2A/PHK2A.rps-blast.out | tail -n +2 | less -S
```


To further analyze our data, we used the following commands to elucidate more regarding the predicted Pfam domains.

To count the most commonly found Pfam domain, we utilized the following code:
```
cut -f 6 ~/lab08-$MYGIT/PHK2A/PHK2A.rps-blast.out | sort | uniq -c
```

To find the longest annotated protein domain, we used the following code:
```
awk '{a=$4-$3;print $1,'\t',a;}' ~/lab08-$MYGIT/PHK2A/PHK2A.r ps-blast.out | sort -k2nr
```

Finally, to find the domain annotation with the best evalue, the following code was used:
```
cut -f 1,5 -d '\t' ~/lab08-$MYGIT/globins/globins.rps-blast.out
```


## 5. Conclusion
It is important to note that some codes may have been cut out as they were simply not as relevant in the replication of our study (for example, using less to see an output, providing the specific commands when moving between folders and directories, and pushing the git to the repository). 
