# DATA
### VCF FILES
-[1KG 30x-grch38](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38)<br>
-[phased VCFs](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)<br>
-[sample information](https://www.internationalgenome.org/api/beta/sample/_search/igsr-1000%20genomes%2030x%20on%20grch38.tsv.tsv)<br>
-[pedigree information](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt)<br>
-[PED-SIM, simulated](https://github.com/williamslab/ped-sim)
- following instructions on the README.md, all necessary dependencies should be installed with mamba.<br>
- ped-sim Makefile needs a simple change to include mamba environment include directory. copy of this edit is [here](https://github.com/kristen-schneider/precision-medicine/tree/main/notes/Makefile).<br>
### MAP FILES
-[HapMap genetic maps in cM units from Beagle](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip)<br>
-[ilash analyzer](https://github.com/roohy/ilash_analyzer/blob/master/interpolate_maps.py) to help create map files.<br>
