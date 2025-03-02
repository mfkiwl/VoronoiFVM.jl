# Module for using ReTest during development.
# See  https://juliatesting.github.io/ReTest.jl/stable/#Quick-start
# In an environment with VoronoiFVM, test dependencies and ReTest and Revise
# with `using ReTest; ReTest.load("test/VoronoiFVMReTest.jl")`,
# one can run specific testsets, e.g. `retest("Aqua")`
# `retest(dry=true)` lists all possible testsets.
module VoronoiFVMReTest
using ReTest
include("alltests.jl")
end
