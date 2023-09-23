configfile: "config.yml"
seeds = [i+10 for i in range(config["nseeds"])]

#TODO: add rule generate cna trees
#ruleorder:  simulate > preprocess

rule all:
    # noinspection PyInterpreter
    input:
        expand("/Users/annahart/CLionProjects/clonesim/simulation/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/",
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust']
        ),

rule simulate:
    input: "/Users/annahart/CLionProjects/clonesim/build/output/cnatrees.txt"
    output:
        tree =  "/Users/annahart/CLionProjects/clonesim/simulation/input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/tree.txt",
        prop=  "/Users/annahart/CLionProjects/clonesim/simulation/input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/prop.csv",
        genotypes ="/Users/annahart/CLionProjects/clonesim/simulation/input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/genotypes.txt",
        phi ="/Users/annahart/CLionProjects/clonesim/simulation/input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/cellAssignments.p0",
        sparse ="/Users/annahart/CLionProjects/clonesim/simulation/input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/sparse.p0",
        copy_profiles ="/Users/annahart/CLionProjects/clonesim/simulation/input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/cells.p0",
    params:
        cp_thresh = 0.05,
        purity = 0.99,
        sample = 1,
        alpha= 0.001,
        truncalSegs = 2,
        simout_dir = "simout/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}",
        scout_dir = "scout/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}"
    log:
        std ="/Users/annahart/CLionProjects/clonesim/simulation/input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/run.log",
        err ="/Users/annahart/CLionProjects/clonesim/simulation/input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/err.log"
    shell:
        "clonesim/build/simulate -r -S {input} "
        " -purity {params.purity} -minProp {params.cp_thresh} "
        " -kk {params.truncalSegs} "
        "-s {wildcards.s}  -l {wildcards.mclust} "
        "-k {wildcards.nsegs} -n {wildcards.snvs} -m {params.sample}  "
        "-out_dir {params.simout_dir}  > {log.std} 2> {log.err} "
        "clonesim/build/generatesinglecells -num_cells {wildcards.cells} -read_depth {wildcards.cov} "
        "-alpha_fp {params.alpha} -out_dir {params.scout_dir} -in_dir {params.simout_dir} -k {wildcards.nsegs}"
        " -m {params.sample} > {log.std} 2> {log.err} "
















