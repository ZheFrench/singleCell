### SingleCell RnaSeq Analysis
_________________

You will find bash and R scripts in two distincts directories.  

1. `xenomeCreation.sh` create the genome reference based on 10X tutorial that will be used by Cell Ranger. 

2. `align10x.sh` run Cell Ranger alignment for the several samples. 

3. `filterMouseReads.sh` remove reads that mapped on mouse and create clean fasta with reads specific to human.

4. `qualityControl.R` make plots to check quality.

5. `filterCells.R` filter single cells (based on previous generated figures) Seurat Object and save it to *_Clean.rds file on disk. Thresholds are stored in files containted in the data directory.

6. `DE-speudoBulk.R` retrieve individual *_Clean.rds files from previous step and do speudo-bulk differential expression analysis between two conditions.

7. `clusterCells.R` from the rds object containing filtered seurat object for each condition , you apply seurat normalisation and save in an *_integrated.rds file... (* can be 2 or several conditions)

8. `plotclusterCells.R` load the *_integrated.rds file..

In `bash/triggers`, there is bash scripts that call (trigger) other bash scripts...


* `trigger_filterMouseReads.sh` call filterMouseReads.sh for several samples.

* `trigger_qualityControl.sh` calls sequentially `qualityControl.R`, `filterCells.R` and `DE-speudoBulk.R` for several samples/conditions.
