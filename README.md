# Cosbin
This is the repository of our article: [TBD]

## Installation
1. Clone or download this GitHub repository
2. Install all the packages required in `Dependencies.R`

## Repository Introduction
- `Cosbin_functions.R` houses all the `Cosbin` functions.
- Check `toy_exmaple.R` and `Cosbin toy example.xlsx` to see how `Cosbin` works step by step.

- Full experiment workflow:
  - `Generate_idealistic_simulation_data.R` (or any of your data) 
  - Calculate the average of each group as the input of `Cosbin`
  - Apply `cosbin()` function to the data
  - `evaluation.R` 
  - Apply `cosbin_convert()` to get the final results

- Application to real data (workflow):
`Example.R`


## Paper Overview
![Cosbin workflow](https://user-images.githubusercontent.com/42553263/159149958-4fff031f-64b1-491b-9405-4a057ed166e3.png)
