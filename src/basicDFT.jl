module basicDFT

using LinearAlgebra, Printf, DataFrames, CSV
include("numerov.jl")

# add missing directories in current folder
workdir = "./"
any(x->x=="Data", readdir(workdir)) || mkdir(workdir+"Data")
any(x->x=="Plots", readdir(workdir)) || mkdir(workdir+"Plots")



end
