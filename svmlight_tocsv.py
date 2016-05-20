import numpy as np
from sklearn.datasets import load_svmlight_file
V,y=load_svmlight_file('U1.lsvm',100)
U,y=load_svmlight_file('M1.lsvm',100)
U=np.array(U.todense())
V=np.array(V.todense())
np.savetxt('U.csv', U, fmt='%f', delimiter=',', newline='\n')
np.savetxt('V.csv', V, fmt='%f', delimiter=',', newline='\n')
