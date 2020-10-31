module basicDFT

using LinearAlgebra, Printf, CSV, DataFrames
export LinearAlgebra, Printf, CSV, DataFrames
include("numerov.jl")
export Numerov
include("kohn_sham.jl")
export solve_KS

# add missing directories in current folder
workdir = "./"
any(x->x=="Data", readdir(workdir)) || mkdir(string(workdir,"Data"))
any(x->x=="Plots", readdir(workdir)) || mkdir(string(workdir,"Plots"))



end
