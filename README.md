# Signature SNVs for FEAST

SignatureSNVs is a method to generate signature SNVs for input into FEAST for source tracking ([Shenhav et al. 2019](https://github.com/cozygene/FEAST)). The signature SNVs are selected from SNV output produced by running the metagenomic pipeline MIDAS ([Nayfach et al. 2017](https://github.com/snayfach/MIDAS)).

Source tracking is a broad term for methods that can estimate the percentage of a microbiome of interest that derives from differences potential sources. A sample of an infant's gut microbiome may be a _sink_ of interest 

Two key terms in understanding source tracking and our approach for signature selection are _sink_ and _source_. A _sink_ is a sample that you are interested in investigating the sources of, such as the gut microbiome of an infant. You may be interested in investigating how much the mother, the crib and the dog contribute to this infant's microbiome. You therefore collect samples from all these potential sources. Once you have whole metagenomic shotgun sequencing on these samples, you are ready to begin analysis to find signature SNVs. 

A signature SNV is a SNV that has a higher probability of coming from one source over other sources or only the sink. The output of counts of the alternative and reference allele go into FEAST for source tracking. FEAST can be run as usual with this input. 

The general workflow is as follows:

```mermaid
graph LR;
A[Metagenomic shotgun data] -->B(MIDAS)
    B --> C{Signature-SNVs}
    C --> D{FEAST}
<!-- 
seq=>start: Metagenomic shotgun data
mid=>midas: MIDAS
sig=>SignatureSNVs: Signature-SNVs
fea=>FEAST

seq->mid->sig->fea -->
```




![alt text](https://github.com/garudlab/signature-snvs/tree/main/misc/Official_Abstract.png?raw=true)


## Table of Contents

1. [Tutorial](#tutorial)
1. [Quick Start](#quickstart)
1. [Example Input Files](#exampleinput)
1. [Optional Pre-processing of signature SNV files](#preprocessing)
1. [Optional Post-processing of signature SNV files](#postprocessing)
1. [Example run of FEAST](#feast)


## <a name="tutorial"> Tutorial </a>

1. In your documents folder, git clone this repo to get the example directory

    ```
    git clone https://github.com/garudlab/Signature-SNVs.git
    ```
1. [Optional] Go into directory SNV-FEAST and Create a virtual environment
    Note: to avoid any dependency conflicts, we recommend installing this in a virtual environment
    ```
    python3 -m pip install --user virtualenv"
    python3 -m virtualenv signature_snvs_env
    source ./snv_feast_env/bin/activate
    ```

2. Install SNV-FEAST with pip. (takes about 1 min)
    This also automatically updates any dependencies. 

    ```
    python3 -m pip install Signature-SNVs==0.0.8
    ```

2. [Option 1] Run code as a module. Start up python in command line interface, then import and run module with example1

    In terminal:
    ```
    python3
    ```

    In interactive python:
    ```
    from signature_snvs import signature_snvs
    signature_snvs.signature_snvs_per_species(species="Bacteroides_uniformis_57318", min_reads=5, start_index=1, end_index=200, config_file_path="configs/config.yaml")
    ```

4. [Option 2] Run on command line with example1

    ```
    python <site-packages_directory>/snv_feast/snv_feast_cli.py --species Bacteroides_uniformis_57318 --min_reads 5 --start_index 1 --end_index 200 --config_file_path configs/config.yaml
    ```

    For example, my site-packages directory is `./lib/python3.9/site-packages/`





## <a name="quickstart"> Quick Start </a>

1. Install SNV-FEAST with pip. (takes about 1 min)
    This also automatically updates any dependencies. 

    ```
    python3 -m pip install SNV-FEAST==0.0.8
    ```

2. Set up your directories and config.

    **Required input**
    'example_template' shows how the directory and config should be set up

    There can be a single directory containing the following:

    1. **sink_source.csv**  a comma-delimited file with the sink source configuration. It should have the accession numbers  in the first column for each sink of interest, and in the following columns, the accession numbers for the sources for each sink (example [sink_source.csv](#sinksource))
    2. **midas_output/snps** MIDAS output with a subdirectory called 'snps/' with it's own subdirectory for each species. In side each species subsubdirectory are two bzipped files 'snps_depth.txt.bz2' and 'snps_ref_freq.txt.bz2' output from MIDAS snps and MIDAS merge_snps step
    3. **config.yaml** YAML indicating the directory with input files, the snp directory, and the output directory for the signature SNVs (example [config.yaml](#config))

3. Determine input arguments:
    * **species** : the species you want to get signature SNVs for
    * **min_reads** : minimum reads required at a site to determine signature SNVs. Recommend 10 if sufficiently high coverage sample, otherwise 5 reads.
    * **start_index** : is there a specific region you want to check for signature SNVs? This number is the row in the midas output merged snps file for snp_depth and snps_ref_freq. If you want to check the whole file, provide 0.
    * **end_index** :  end index for the region of interest. This number is the row in the midas output merged snps file for snp_depth and snps_ref_freq. If you want to check the whole file, provide the length of the file, or some high number (e.g. 10000000), or determine the length of the file from [here](#preprocessing)
    * **config_file_path** : the path where the config.yaml is located

3. [Option 1] Import and run module

    ```
    from signature_snvs import signature_snvs 
    signature_snvs.signature_snvs_per_species(species="Bacteroides_uniformis_57318", min_reads=5, start_index=1, end_index=200, config_file_path="configs/config.yaml")
    ```

4. [Option 2] Run on command line

    ```
    python <site-packages_directory>/snv_feast/snv_feast_cli.py --species Bacteroides_uniformis_57318 --min_reads 5 --start_index 1 --end_index 200 --config_file_path configs/config.yaml
    ```

    For example, my site-packages directory is `./lib/python3.9/site-packages/`



## <a name="exampleinput"> Example Inputs </a>

### <a name="sinksource"> Example sink_source.csv file </a>
This is a comma-delimited file where each row of this table represents a single source tracking experiment. The first cell in each row is the accession number for the sink sample (matching the accession number in the MIDAS output). The second and onward cells in each row should be the accession numbers for the sources for each sink


| family_id	| B	| M	| M1 | M2 | M3|
|-----------|---|---|---|---|---|
|Experiment1	|ERR00001	|ERR00010	|ERR00020	|ERR00030	|ERR00040|
|Experiment2	|ERR00002	|ERR00020	|ERR00030	|ERR00040	|ERR00010|
|Experiment3	|ERR00003	|ERR00030	|ERR00040		|ERR00010	|ERR00020|
|Experiment4	|ERR00004	|ERR00040		|ERR00010	|ERR00020	|ERR00030|

In the toy example_1


| family_id	| B	| M	| M1 | M2 | M3|
|-----------|---|---|---|---|---|
|Test1	|Baby1	|Mother1	|Mother2	|Mother3	|Mother4|
|Test2	|Baby2	|Mother2	|Mother3	|Mother4	|Mother1|
|Test3	|Baby3	|Mother3	|Mother4	|Mother1	|Mother2|
|Test4	|Baby4	|Mother4	|Mother1	|Mother2	|Mother3|


### <a name="config"> Example config.yaml </a>

Example config.yaml:

```
input_dir: '/Users/leahbriscoe/Documents/FEASTX/snv-feast/example_1/'
snp_dir: '/Users/leahbriscoe/Documents/FEASTX/snv-feast/example_1/midas_output/snps/'
output_dir: '/Users/leahbriscoe/Documents/FEASTX/snv-feast/example_1/signature_snvs/'
```
## <a name="preprocessing"> Optional Pre-Processing </a>
Files for post-processing of the signature SNV output are [here](https://github.com/garudlab/SNV-FEAST/tree/main/preprocessing).
We had determined a window size of 200,000 bp was helpful for analysis. First we determined the length of the species files

To be run inside your data directory (e.g. example_1)
1. Step 1: Get length of all snps_depth file given a list of species in the midas_output/snps directory
2. Step 2: Generate a file with all the lengths of species files


## <a name="postprocessing"> Optional Post-Processing </a>

Files for post-processing of the signature SNV output are [here](https://github.com/garudlab/SNV-FEAST/tree/main/postprocessing).

1. Step 1: Merge signature SNV files across windows per species per source tracking experiment
2. Step 2: Merge signature SNV files across speciess per source tracking experiment


## <a name="feast"> R Example of running FEAST  </a>

Here is an example FEAST command for R.
```
FEAST(C = snv_count_matrix, metadata =metadata, different_sources_flag = 0, dir_path =input_dir,
                                   outfile="demo",COVERAGE =coverage_min)
```
where coverage min is the minimum total reads per sample of all samples







