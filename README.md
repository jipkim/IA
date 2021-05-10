# Impact of Imperfect Foresight on the Optimal DER Deployment, Remuneration and Policy

This code was developed under Julia v1.6/JuMP0.21 by Jip Kim in 2021.
The following packages must be installed:

  - JuMP
  - Distributions
  - LinearAlgebra
  - KNITRO
 
To run the code, execute CCG.jl, or include() it from a Julia prompt.

Testsystem is given as follows in data folder:
  - Node.csv: Transmission network data
  - Line.csv: Transmission line data
  - Generator.csv: Generation data
  - We also specify additional input data in the beginning of MPini.jl/MP.jl/SP.jl
