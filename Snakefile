rule all:
	input:
		'output/sce_paper1.RData',
		'output/sce_paper2.RData',
		'output/sce_paper3.RData'

rule paper1:
	output:
		'output/sce_paper1.RData'
	resources:
		mem_gb=200
	benchmark:
		'benchmarks/paper1.txt'
	log:
		'logs/paper1.log'
	shell:
		'src/paper1.sh >& {log}'

rule paper2:
	output:
		'output/sce_paper2.RData'
	resources:
		mem_gb=200
	benchmark:
		'benchmarks/paper2.txt'
	log:
		'logs/paper2.log'
	shell:
		'src/paper2.sh >& {log}'

rule paper3:
	output:
		'output/sce_paper3.RData'
	resources:
		mem_gb=200
	benchmark:
		'benchmarks/paper3.txt'
	log:
		'logs/paper3.log'
	shell:
		'src/paper3.sh >& {log}'