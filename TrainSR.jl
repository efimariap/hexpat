######################################################
# Julia file for symbolic regression

# Filippos Sofos
# Department of Physics
#University of Thessaly, GR
######################################################


import Pkg
Pkg.add("DataFrames")
Pkg.add("SymbolicRegression")
Pkg.add("CSV")
Pkg.add("Random")
Pkg.add("LossFunctions")

using DataFrames
using SymbolicRegression
using CSV
using Random
using LossFunctions

# create a function for partitioning the data
function partitionTrainTest(data, at = 0.8)
    n = nrow(data)
    idx = shuffle(1:n)
    train_idx = view(idx, 1:floor(Int, at*n))
    test_idx = view(idx, (floor(Int, at*n)+1):n)
    data[train_idx,:], data[test_idx,:]
end

#Load data from path
Df=CSV.File("....Datapath") |> DataFrame

######################################################
Main Loop
######################################################

#Choose number of iterations, 30 is suggested
for ix in 1:30

    Random.seed!(100)
    train,test = partitionTrainTest(Df, 0.8)


#Choose columns from the Df for X (inputs) and y (output)
    X=......
    y=......
   
    options = SymbolicRegression.Options(
        binary_operators=(+, *, /, -),
        batching=true,
        npopulations=40,
        maxsize=25,
        maxdepth=10,
        fractionReplaced=0.20f0,
        fractionReplacedHof=0.20f0,
        npop=1000, ns=50,
        ncyclesperiteration=300,
        bin_constraints=[(^)=>(0,10)],
        hofFile="Results/result_$(ix).csv")
    hallOfFame = EquationSearch(X, y, niterations=40, options=options);  
    println("finished run$(ix)")
end




