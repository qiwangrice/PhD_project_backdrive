# **[Control Theory Metagenome] Pipeline**

Find driver species from real metagenomic data and simulate fecal microbial transplantation FMT process
Script and results for the publilcation are in the data branch. 

## Installation 

```
pip install micom
pip install pulp
```

**Dependencies:** [Python 3.x](https://www.python.org/download/releases/3.0/), database of the genome-scale metabolic reconstruction of human gut microbes [AGORA database](https://github.com/VirtualMetabolicHuman/AGORA),  a python package for metabolic modeling of microbial communities [MICOM](https://github.com/micom-dev/micom), a python package of an LP modeler [PuLP](https://pypi.org/project/PuLP/)

## **Usage**

### **Step1** Ecological Networks Inference 

Infer ecological networks from metagenomic taxonomic classification results 

```
usage: thesis.py interaction [-h] [-m MEDIUM] [-d MODEL] [-p PERCENTAGE]
                             [-f FLAG] [-o OUTPUT]
                             input_file

positional arguments:
  input_file            Input file of a list of taxonomic classification files

optional arguments:
  -h, --help            show this help message and exit
  -m MEDIUM, --medium MEDIUM
                        Medium CSV file, default medium.csv
  -d MODEL, --model MODEL
                        Metabolic model database, default dbs
  -p PERCENTAGE, --percentage PERCENTAGE
                        Percentage of species removed, default 0.1
  -f FLAG, --flag FLAG  Calculate growth rate, default True
  -o OUTPUT, --output OUTPUT
                        output folder, default output_interaction
```

**Example**

### **Step2** Driver Species Identification

Identify driver species from a multilayer ecological network 

```
usage: thesis.py drivers [-h] [-s STRENGTH] [-p PREFIX] [-o OUTPUT]
                         input_folder

positional arguments:
  input_folder          Input Folder of Bacteria Interaction Networks

optional arguments:
  -h, --help            show this help message and exit
  -s STRENGTH, --strength STRENGTH
                        Threshold of Interaction Strength, default 0.2
  -p PREFIX, --prefix PREFIX
                        Output file prefix, default driver_nodes
  -o OUTPUT, --output OUTPUT
                        Output file folder, default output_driver
```

### **Step3** Simulate Fecal Metagenomic Transplantation (FMT) Process

Simulate the FMT process following the Generalized Lotka-Volterra (GLV) model

#### a) FMT donor samples 

Add input donor sample directly to a given disease sample

```
usage: thesis.py fmt_all [-h] [-m MEDIUM] [-d MODEL] [-p PERCENTAGE]
                         [-s STRENGTH] [-o OUTPUT]
                         input_file

positional arguments:
  input_file            Input Disease and Donor Samples File Address

optional arguments:
  -h, --help            show this help message and exit
  -m MEDIUM, --medium MEDIUM
                        Medium CSV file, default medium.csv
  -d MODEL, --model MODEL
                        Metabolic model database, default dbs
  -p PERCENTAGE, --percentage PERCENTAGE
                        Percentage of species removed, default 0.1
  -s STRENGTH, --strength STRENGTH
                        Threshold of Interaction Strength, default 0.2
  -o OUTPUT, --output OUTPUT
                        Output file folder, default fmt_output
```

#### b) FMT driver species 

Add equal amounts of driver species to the disease sample

```
usage: thesis.py fmt_driver [-h] [-i DRIVER] [-a AMOUNT] [-m MEDIUM]
                            [-d MODEL] [-p PERCENTAGE] [-s STRENGTH]
                            [-o OUTPUT]
                            input_file

positional arguments:
  input_file            Input Disease Samples File Address

optional arguments:
  -h, --help            show this help message and exit
  -i DRIVER, --driver DRIVER
                        Input Driver Species
  -a AMOUNT, --amount AMOUNT
                        Input amount of driver species, default 4g
  -m MEDIUM, --medium MEDIUM
                        Medium CSV file, default medium.csv
  -d MODEL, --model MODEL
                        Metabolic model database, default dbs
  -p PERCENTAGE, --percentage PERCENTAGE
                        Percentage of species removed, default 0.1
  -s STRENGTH, --strength STRENGTH
                        Threshold of Interaction Strength, default 0.2
  -o OUTPUT, --output OUTPUT
                        Output file folder, default fmt_driver_output
```

#### c) FMT only

Given ecological networks of after-FMT samples, fmt_only only simulate species abundance changes following the GLV model

```
usage: thesis.py fmt_only [-h] [-s STRENGTH] [-p PREFIX] [-o OUTPUT]
                          input-file

positional arguments:
  input-file            Input file prefix

optional arguments:
  -h, --help            show this help message and exit
  -s STRENGTH, --strength STRENGTH
                        Threshold of Interaction Strength, default 0.2
  -p PREFIX, --prefix PREFIX
                        Output file prefix
  -o OUTPUT, --output OUTPUT
                        Output file prefix
```

## **Description**
**Pipeline**
![](image/pipeline.png)
