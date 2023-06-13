version 1.0

task prepareFiles {
	input{
		File ref_load
		File bed
		File bim
		File fam
		Double? overlap
		Int mem_gb = 8
	}

	Float disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB"))) * 1.5	#hoping this works?
	String basename = basename(bed, ".bed")
	
	command <<<
		#get a list of variant names in common between the two, save to extract.txt
    # modify after we know the format of the loadings, just need to pull out variant name
		awk 'FNR==NR{a[$2]; next}{if($2 in a){print $2}}' ~{ref_load} ~{bim} > extract.txt

		#subset bed with --extract extract.txt
		/plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} --extract extract.txt --make-bed --keep-allele-order --out ~{basename}_pcaReady

		#extract variants in-common variants from ref_loadings
    # change format once we know the loadings
		paste -d'\t' <(cut -f2 ~{ref_bim}) ~{P} > tmp
		head tmp
		head ~{ref_bim}
		head ~{P}
		awk 'FNR==NR{a[$1]; next}{if($1 in a){print $0}}' extract.txt tmp | cut -f2- > loadings_pcaReady.txt
    
    #check for overlap, if overlap is less than threshold, stop, default overlap threshold is 0.95
    loadings_count=$(wc -l < ${ref_load})
    new_loadings_count=$(wc -l < loadings_pcaReady.txt)
    #doing this calculation because I'm lousy at bash...
    prop=$(awk -v old=${loadings_count} -v new=${new_loadings_count} '{prop=new/old; print prop}' )
    printf "Variant overlap is ${prop} of original.\n"
    $https://support.terra.bio/hc/en-us/articles/360037484851-Variable-Types-in-WDL#:~:text=When%20working%20with%20optional%20variables%20in%20your%20command%2C,The%20syntax%20for%20that%20is%3A%20%24%20%7Bdefault%3D%22value%22%20variableName%7D
    myoverlap=${default=0.95 overlap}
    exit_code=$(awk -v prop=${prop} -v default_threshold=${myoverlap} '{myexit=0; if(prop < default_threshold){myexit=1}; print myexit }' )
    if [${exit_code} -gt 0]; then
    	printf "SNP overlap ${prop} is lower than
    	exit ${exit_code};
    fi
    
    
	>>>
	
	output {
		File snps_to_keep="extract.txt"
		File subset_bed="~{basename}_pcaReady.bed"
		File subset_bim="~{basename}_pcaReady.bim"
		File subset_fam="~{basename}_pcaReady.fam"
		File subset_log="~{basename}_pcaReady.log"
		File subset_loadings="loadings_pcaReady.txt"
	}
	
	runtime {
    	docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    	#disks: "local-disk " + disk_size + " HDD"
    	memory: mem_gb + " GB"
	}
}

task summary {
	input{
		File P
		File snps
		Int mem_gb = 4
	}
	
	#String my_k = ceil(awk NR==1{print NF} P) #maybe do a check on the file to make sure all lines are the same to make more robust?
	
	command <<<
		my_k=$(head -n 1 ~{P} | awk '{print NF}')
		printf "\nProjection will be run with k=${my_k} clusters.\n"
		printf "${my_k}" > tmp.txt
		snp_count=$(wc -l < ~{snps})
		printf "\nProjection uses ${snp_count} snps.\n"
	>>>
	
	output {
		String k = read_string("tmp.txt")
	}
	
	runtime{
		docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    	#disks: "local-disk " + disk_size + " HDD"
    	memory: mem_gb + " GB"
	}
}

task run_pca_projected {
	input {
    	File bed
    	File bim
    	File fam
    	File loadings
  	}

	Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB")))
	String basename = basename(bed, ".bed")
	#ln --symbolic ${P} ${basename}.${k}.P.in

	command <<<
		command="/plink2 --bfile ${basename} --score loadings --out ${basename}_pca
		printf "${command}\n"
		${command}
	>>>

	runtime {
		docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		#disks: "local-disk " + disk_size + " HDD"
		memory: mem_gb + " GB"
		cpu: n_cpus
  	}

	output {
		#check output file name from --score in plink2
		File pca_projection = "~{basename}_pca.score"
		File projection_log = "~{basename}_pca.log"
	}
}

workflow pca_projection {
	input{
    File ref_loadings
		File ref_bim
		File bed
		File bim
		File fam
		String? mem_gb
		Int? seed # https://wdl-docs.readthedocs.io/en/latest/WDL/different_parameters/
		Int? n_cpus
		Boolean? cv
	}
	
	call prepareFiles {
		input:
      ref_load = ref_loadings,
			bed = bed, 
			bim = bim, 
			fam = fam
	}
	
	call summary{
		input:
			P = P,
			snps = admixReady.snps_to_keep
	}
	
	call run_admixture_projected {
		input:
			bed = admixReady.subset_bed,
    			bim = admixReady.subset_bim,
    			fam = admixReady.subset_fam,
    			P = admixReady.subset_P,
    			k = summary.k,
    			seed = seed,
    			cv = cv,
    			mem_gb = mem_gb,
    			n_cpus = n_cpus
	}
	
	meta {
    author: "Jonathan Shortt"
    email: "jonathan.shortt@cuanschutz.edu"
    description: "## run_projection_admixture\n This workflow is used to project a genetic test dataset (in plink format, i.e., .bed/.bim/.fam) into clusters (\"ancestral populations\") using ADMIXTURE. First, the cluster file (.P produced by ADMIXTURE) and the test dataset are both subset to contain the same set of variants (Note: this workflow assumes that variants from both the .P and test dataset have been previously harmonized such that variants follow the same naming convention, alleles at each site are ordered identically, and variants are sorted). Thenthe test dataset is projected into the clusters determined by the .P."
	}

}
