PAPERS = ["paper1_Neuron", "paper2", "paper3"]
rule all:
	input:
		expand('output/sce_{p}.RData',
			p=PAPERS)
rule paper:
	output:
		'output/sce_{p}.RData'
	benchmark:
		'benchmarks/{p}.txt'
	conda:
		'envs/myenv.yaml'
	resources:
		mem_gb=200
	log:
		'logs/{p}.log'
	shell:
		'src/{wildcards.p}.sh >& {log}'