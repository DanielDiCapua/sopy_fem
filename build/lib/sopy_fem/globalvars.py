import numpy as np

data = {}
neq = 0
ndof = 1
ncomp = 1
num_unknows = 0
madgln = []
astiff = []
asload = []
u_vec = []
u_known = []
amassmat = []
natFreqVec = []
vibrationModes = []
dynamics_uvec = []
modal_disp = []
react_vec = []
results = {}
Eps_elem_vec = [[], [], []]
Sigma_elem_vec = [[], [], []]
scale_disp = 1.0
dataFileName = ""
