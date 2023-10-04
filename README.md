# Cosbin
This is the repository for our application: **Cosbin: Cosine score based iterative normalization of biologically diverse samples**

## Installation
1. Clone or download this GitHub repository
2. Install all the packages required in `Dependencies.R`

## Repository Introduction
- `Cosbin_functions.R` houses all the `Cosbin` functions.
- Check `toy_exmaple.R` and `Cosbin toy example.xlsx` to see how `Cosbin` works step by step.

- Full experiment workflow:
  - `Generate_idealistic_simulation_data_1.R` (or any of your data) 
  - If you are using your own data, you'll need to calculate the average of each group as the input of `Cosbin`
  - Data cleaning (e.g. `data_cleaning()`) & Initial normalization 
  - Apply `cosbin()` function to the data
  - `evaluation.R` 
  - Apply `cosbin_convert()` to get the final results

- Application to real benchmark data:
`Example.R`

- Additional experiments for comparing Cosbin with CSS and Qsmooth:
`Comparison_Experiment.R`


## Paper Overview
### Motivation: 
Data normalization is essential to ensure accurate inference and comparability of gene expressions across samples or conditions. Ideally, gene expressions should be rescaled based on consistently expressed reference genes. However, for normalizing biologically diverse samples, most commonly used reference genes have exhibited striking expression variability, and distribution-based approaches can be problematic when the magnitudes of differentially expressed genes are significantly asymmetric. 
### Results: 
We report an efficient and accurate data-driven method - Cosine score based iterative normalization (Cosbin), to normalize biologically diverse samples. Based on the Cosine scores of cross-group expression patterns, the Cosbin pipeline iteratively eliminates asymmetrically and differentially expressed genes, and accordingly identifies consistently expressed genes and calculates normalization factors. We demonstrate the superior performance and enhanced utility of Cosbin compared with peer methods using both simulation and real multi-omics expression datasets. Implemented in open-source R scripts and specifically designed to address normalization bias due to asymmetric differential expression, the Cosbin tool complements not replaces the existing methods and will allow biologists to detect subtle yet important molecular signals among phenotypic groups.
### illustration:
<img width="935" alt="Screenshot 2023-10-04 at 2 40 34 AM" src="https://github.com/MinjieSh/Cosbin/assets/42553263/10dafeb9-e986-434d-88b0-6d884eeac31f">

![Screen Recording 2023-10-04 at 2 34 58 AM](https://github.com/MinjieSh/Cosbin/assets/42553263/32ffa903-a5b9-4fff-a3e1-662028415567)

