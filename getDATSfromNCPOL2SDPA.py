from __future__ import print_function
from ncpol2sdpa import *
import pickle
# just a basic code to run tri vs bi (or reverse) party NPA optimization
# use inequalities from inequalities.py
# whats up?
         



def __main__():
  X=3
  Y=3
  A=3
  B=3
  P = Probability([A]*X, [B]*Y)
  Objective = 0 # set objective function
# prepare relaxation
  sdpRelaxation = SdpRelaxation(P.get_all_operators(), verbose=0)
  sdpRelaxation.get_relaxation(3, substitutions = P.substitutions)
  sdpRelaxation.set_objective(-Objective)
#  sdpRelaxation.solve(solver="mosek")
# collect data points
#  print(abs(sdpRelaxation.dual))
  sdpRelaxation.write_to_file('chsh.dat-s', 'sdpa')

  print(P.substitutions)
if __name__=="__main__":
   __main__()