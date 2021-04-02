### SingleCell RnaSeq Analysis
_________________

`xenomeCreation.sh` create the genome reference based on 10X tutorial that will be used by Cell Ranger. 

`align10x.sh` run Cell Ranger alignment for the several samples. 

`single.sh` remove reads that mapped on mouse and create clean fasta with reads specific to human.

In bash `triggers dir`, there is bash scripts that call (trigger) other bash scripts...