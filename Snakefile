
configfile: "config.yml"
seeds = [i+10 for i in range(config["nseeds"])]

#TODO: add rule generate cna trees 
ruleorder:  simulate > generatesinglecells

rule all:
    # noinspection PyInterpreter
    input:
        expand("input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/cells.p0",
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust']
        ),

rule simulate:
    input: "cnatrees.txt"
    output:
        tree =  "input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/tree.tsv",
        prop=  "input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/proportions.tsv",
        genotypes ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/node.tsv",

    params:
        cp_thresh = 0.05,
        purity = 0.99,
        sample = 1,
        alpha= 0.001,
        truncalSegs = 2,
        simout_dir = "/Users/annahart/CLionProjects/pharming/simulation_study/input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}",
        scout_dir = "/Users/annahart/CLionProjects/pharming/simulation_study/input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}"
    log:
        std ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/run.log", 
        err ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/err.log" 
    shell:
        "clonesim/build/simulate -r -S {input} "
        " -purity {params.purity} -minProp {params.cp_thresh} "
        " -kk {params.truncalSegs} -f "
        "-s {wildcards.s}  -l {wildcards.mclust} "
        "-k {wildcards.nsegs} -n {wildcards.snvs} -m {params.sample} "
        "-output_file_dir {params.simout_dir}  > {log.std} 2> {log.err}  "


rule generatesinglecells:
    input:
        tree =  "input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/tree.tsv",
        prop=  "input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/proportions.tsv",
        genotypes ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/node.tsv",
    output:
        phi ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/cellAssignments.p0",
        sparse ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/sparse.p0",
        copy_profiles ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/cells.p0",
    params:
        cp_thresh = 0.05,
        purity = 0.99,
        sample = 1,
        alpha= 0.001,
        truncalSegs = 2,
        simout_dir = "input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}",
        scout_dir = "input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}",
    log:
        std ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/run.log", 
        err ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/err.log" 
    shell:
     "clonesim/build/generatesinglecells -num_cells {wildcards.cells} -read_depth {wildcards.cov} "
        "-alpha_fp {params.alpha} -out_dir {params.scout_dir} -in_dir {params.simout_dir} -k {wildcards.nsegs}"
        " -m {params.sample} > {log.std} 2> {log.err} "
    
    


"""
rule preprocess:
    input: 
        tree =  "input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/tree.tsv",
        genotypes ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/node.tsv",
        phi ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/cellAssignments.p0",
        sparse ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/sparse.p0",
        copy_profiles ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/cells.p0",
    output:
        png =  "input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/tree.png",
        data ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/data.pickle",
        ct ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/clonal_tree.pickle",
    log:
        std ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/preprocess.log", 
        err ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/preprocess.err.log" 
    shell:
        "python ../src/preprocess.py -f {input.sparse} -c {input.copy_profiles} "
        "-t {input.tree} --phi {input.phi} -g {input.genotypes} "
        "-D {output.data} -T {output.ct} --draw {output.png} "
        " > {log.std} 2> {log.err} "


rule sens_analysis: 
    input:
        data ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/data.pickle",
        ct ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/clonal_tree.pickle",
    params:
        nreps = 10,
        cell_rates = "10:110:10",
        mut_rates = "10:110:10",
    output: "obj_comp/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}_obj.csv"
    shell:
        "nice python ../src/sens_analysis.py -D {input.data}  "
        "-T {input.ct} -c {params.cell_rates} -m {params.mut_rates} "
        "-s {wildcards.s} -r {params.nreps} "
        "-o {output} "
"""






        






