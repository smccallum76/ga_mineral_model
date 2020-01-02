# ga_mineral_model
Code for defining mineral model input matrix using a GA and application of tuned matrix for mineral predictions. 
This code consists of run files and function files.  
The function files are called by the run files.
The '01..' file is used to tune the elemental fractions using a genetic algorithm.  
After the '01..' file has defined the tuned matrix, it is imported by the '02..' run file and used to make mineral predictions.  
Some modification to the import data must be made depending on the users file structure and naming convention.  
