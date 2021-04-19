### SingleCell RnaSeq Analysis
_________________

You will find bash and R scripts in two distincts directories.  

1. `xenomeCreation.sh` create the genome reference based on 10X tutorial that will be used by Cell Ranger. 

2. `align10x.sh` run Cell Ranger alignment for the several samples. 

3. `filterMouseReads.sh` remove reads that mapped on mouse and create clean fasta with reads specific to human.

4. `qualityControl.R` make plots to check quality.

5. `filterCells.R` filter single cells (based on previous generated figures) Seurat Object and save it to rds file on disk. Thresholds are stored in files containted in the data directory.


5. `DE-speudoBulk.R` retrieve individual rds files and do speudo-bulk differential expression analysis between two conditions.


In `bash/triggers`, there is bash scripts that call (trigger) other bash scripts...


* `trigger_filterMouseReads.sh` call filterMouseReads.sh for several samples.

* `trigger_qualityControl.sh` calls sequentially `qualityControl.R`, `filterCells.R` and `DE-speudoBulk.R` for several samples/conditions.
