# The Art of Thinking in Time: A Cognitive Modeling Approach to Musical Synchrony in Jazz Trios
**Marc Barcelos** (202302260)

> *All code edits were completed before the official exam deadline: **June 2, 2025**.*

---

## Master's Thesis – Cognitive Science MA (Spring 2025)

This repository contains all data and code associated with the thesis:  
**_The Art of Thinking in Time: A Cognitive Modeling Approach to Musical Synchrony in Jazz Trios_**

---

## Repository Structure

- **`/Data/`**  
  Contains all relevant datasets folders:
  - Original raw data  
  - Filtered and preprocessed datasets  
  - Concatenated files ready for analysis  
  - Descriptive visualizations

- **`/Analysis/`**  
  Contains all relevant analysis files:
  - **`DescriptiveStats_JazzTrios.ipynb`** included all of the data preprocessing prior to model fitting and data analysis
  - **`ThesisADAM_bounded.m`** was the MATLAB file through which both versions of the ADAM were fit to the data
  - **`/SMSToolbox_Piano/`**  contains all of the necessary MATLAB scripts to run **`ThesisADAM_bounded.m`**
  - **`/Statistical_Analyses/`** contains all of the R files used to undergo the statistical analyses outlined in the paper, primarily linear mixed effects models and data visualizations
    
- **`/Outputs/`**  `**
  Contains all the outputs produced in order to produce the thesis:
  - **`ThesisADAM_boundedresults.xlsx`** contains the the ADAM related and descriptive values produced by **`ThesisADAM_bounded.m`**
  - **`/Script #/`** for all 6 R scripts used in this study, containing these folders in all of them:
    - **`/Script #/Plots/`** includes all of the visualizations produced within a specific R  script
    - **`/Script #/Results/`** includes all of the results from the statistical analyses performed within a specific R  script
---

Feel free to reach out for questions or replication assistance.
