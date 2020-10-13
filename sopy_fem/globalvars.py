import numpy as np

data = {}
neq = 0
ndof = 1
ncomp = 1
num_unknows = 0
madgln = []
astiff = []
asload = []
amassmat = []
natfreq_vec = []
vibrationModes = []
u_vec = []
u_known = []
react_vec = []
results = {}
Eps_elem_vec = [[], [], []]
Sigma_elem_vec = [[], [], []]
scale_disp = 1.0
dataFileName = ""
