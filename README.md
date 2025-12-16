# Gene2Struct

**Gene2Struct** is an advanced bioinformatics toolkit specifically designed for the investigation of **functional genes in non-model organisms**.  
It provides an integrated workflow that bridges genomic data mining, evolutionary inference, and structural–functional characterization, encompassing:  

- **Accurate recovery of coding sequences (CDS)** from shallow sequencing data to enable robust **phylogenetic reconstruction**  
- **Comprehensive evolutionary analyses**, including **conservation profiling** and **positive selection detection** at both the **site** and **branch** levels  
- **Protein structure prediction** coupled with large-scale **molecular docking simulations**  

By uniting methodologies from **phylogenomics**, **structural biology**, and **functional genomics**, Gene2Struct establishes a versatile platform for elucidating the molecular mechanisms of gene function in diverse non-model lineages.  



## Installation

### 1. Clone the repository

To begin, clone the repository from its GitHub source and navigate into the project directory.

```bash
git clone https://github.com/scu-shiyi/Gene2Struct.git
cd Gene2Struct
```

### 2. Install MGLTools (Optional)

We recommend that users manually download the MGLTools package and place it under the utils/ directory:
- Download  [MGLTools](https://ccsb.scripps.edu/mgltools/downloads/) .
- Place the extracted folder under the `Gene2Struct/gene2struct/utils/` directory.

### 3. Configure the environment

Set up the Conda environment using the provided YAML file to ensure all required dependencies are installed.

```bash
conda env create -f environment.yaml
conda activate gene2struct
pip install .
```

---

## Usage

Coming soon... (examples of phylogenetic analysis, conservation test, structure prediction, and docking will be provided)

---

## Citation

If you use **Gene2Struct** in your research, please cite this repository in your publications.
