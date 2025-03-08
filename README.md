An R implementation of the Multi State Survival Supervised Clustering (MS<sup>3</sup>C) algorithm for cancer subtyping based on longitudinal data.

The code is arranged in the following way:  
**Data**: Helpers function to load and preprocess the dataset (the dataset is not published)  
**Core**: Functions that implement the MS<sup>3</sup>C algorithm  
**Models**: Definitions of the MSM employed in the study  
**Validation**: Implementation of different metrics for the evaluation of the goodness of the clustering  
  
**Main.r**: Entry point for the full pipeline  
**Grid search.r**: Code for hyperparameters tuning
