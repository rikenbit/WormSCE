rule all:
	input:
		'data/paper1.RData'

rule paper1:
	output:
		'data/paper1.RData'
	benchmark:
		'benchmarks/paper1.txt'
	log:
		'logs/paper1.log'
	shell:
		'src/paper1.sh >& {log}'
		
	#'/home/yamaken/software/R-4.0.0/bin/Rscript'
	#'/home/yamaken/software/R-4.0.0/bin/R'
	#download.file("http://waterston.gs.washington.edu/sci_RNA_seq_gene_count_data/Cao_et_al_2017_vignette.RData", "Cao_et_al_2017_vignette.RData")