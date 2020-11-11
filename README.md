
# Multi-channel Fitter

This package is designed to fit multiple calorimeters (crystals, LS veto) in multiple channels (single hit, multiple hit) and in multiple energies (low energy, high energy) all simultaneously. In theory this will decrease the parameter space and provide a more accurate description of the background model. 



## Software versions

This package is written in python and uses the PyROOT bindings. Currently this package only works with python-2. I'm planning to make this python-3 compatible but that is not a top priority. This package works with root-v5 and root-v6.



## Data selection

Data selection is done by using the offline good runs list and only selecting run/subruns that are marked as good for analysis. This package is currently using the SET3 data set with the V00-04-15 data production.



