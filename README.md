# **[Control Theory Metagenome] Pipeline**

Find driver species from real metagenomic data and simulate fecal microbial transplantation (FMT) process
Data analysis script and final results for the publilcation are in the data branch. 

## Installation 

```
pip install micom
pip install pulp
```

**Dependencies:** [Python 3.x](https://www.python.org/download/releases/3.0/), database of the genome-scale metabolic reconstruction of human gut microbes [AGORA database](https://github.com/VirtualMetabolicHuman/AGORA),  a python package for metabolic modeling of microbial communities [MICOM](https://github.com/micom-dev/micom), a python package of an LP modeler [PuLP](https://pypi.org/project/PuLP/)

## **Usage**

Bacdrive pipeline contains four modules: 
1. ecological network inferences
2. driver species identification
3. FMT process simulation: a) donor sample transplantation b) driver species transplantation

```
usage: bacdrive.py [-h] {interaction,driver,fmt_donor,fmt_driver,fmt_only} ...

positional arguments:
  {interaction,driver,fmt_donor,fmt_driver,fmt_only}
                        sub-command help
    interaction         Bacterial interaction inference using MICOM
    driver              Driver nodes detection using MDSM
    fmt_donor           After-FMT community construction and simulation
                        following the GLV model
    fmt_driver          Afte-driver species transplantation (ADT) community
                        consturction and simulation following the GLV model
    fmt_only            After-FMT or ADT simulation following the GLV model

optional arguments:
  -h, --help            show this help message and exit
```

### **Step1** Ecological Networks Inference 

Infer ecological networks from metagenomic taxonomic classification results 

```
usage: bacdrive.py interaction [-h] [-m MEDIUM] [-d MODEL] [-p PERCENTAGE]
                               [-f FLAG] [-o OUTPUT]
                               input_file

positional arguments:
  input_file            Input file of a list of taxonomic classification file
                        addresses

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

```
python bacdrive.py interaction example/example_donor_input.txt -o example/donor_interaction
```

### **Step2** Driver Species Identification

Identify driver species from a multilayer ecological network 

```
usage: bacdrive.py driver [-h] [-s STRENGTH] [-p PREFIX] [-o OUTPUT]
                          input_folder

positional arguments:
  input_folder          Input folder of bacteria interaction networks

optional arguments:
  -h, --help            show this help message and exit
  -s STRENGTH, --strength STRENGTH
                        Threshold of interaction strength, default 0.2
  -p PREFIX, --prefix PREFIX
                        Output file prefix, default driver_nodes
  -o OUTPUT, --output OUTPUT
                        Output file folder, default output_driver
```
**Example**
```
python bacdrive.py driver example/donor_interaction -o example/donor_drivers
```


### **Step3** Simulate Fecal Metagenomic Transplantation (FMT) Process

Simulate the FMT process following the Generalized Lotka-Volterra (GLV) model

#### a) FMT donor samples 

Add input donor sample directly to a given disease sample

```
usage: bacdrive.py fmt_donor [-h] [-m MEDIUM] [-d MODEL] [-p PERCENTAGE]
                             [-s STRENGTH] [-o OUTPUT]
                             input_file

positional arguments:
  input_file            Input disease and donor sample file addresses

optional arguments:
  -h, --help            show this help message and exit
  -m MEDIUM, --medium MEDIUM
                        Medium CSV file, default medium.csv
  -d MODEL, --model MODEL
                        Metabolic model database, default dbs
  -p PERCENTAGE, --percentage PERCENTAGE
                        Percentage of species removed, default 0.1
  -s STRENGTH, --strength STRENGTH
                        Threshold of interaction strength, default 0.2
  -o OUTPUT, --output OUTPUT
                        Output file folder, default fmt_output
```

#### b) FMT driver species 

Add equal amounts of driver species to the disease sample

```
usage: bacdrive.py fmt_driver [-h] [-i DRIVER] [-a AMOUNT] [-m MEDIUM]
                              [-d MODEL] [-p PERCENTAGE] [-s STRENGTH]
                              [-o OUTPUT]
                              input_file

positional arguments:
  input_file            Input a list of disease sample file addresses

optional arguments:
  -h, --help            show this help message and exit
  -i DRIVER, --driver DRIVER
                        Input Driver Species
  -a AMOUNT, --amount AMOUNT
                        Input amount of driver species, default 40g
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

**Example**
```
python bacdrive.py fmt_driver example/example_disease_input.txt -i example/donor_drivers/driver_nodes.3layer.str02.txt -o example/fmt_driver_output
```

#### c) FMT only

Given ecological networks of after-FMT or ADT samples, fmt_only only simulate species abundance changes following the GLV model

```
usage: bacdrive.py fmt_only [-h] [-s STRENGTH] [-p PREFIX] [-o OUTPUT]
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
