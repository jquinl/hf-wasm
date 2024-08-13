from scipy.linalg import lu,inv,qr,eigh
import numpy as np

f = np.vectorize(lambda x: x if abs(x)>1e-5 else 0.0)

import numpy as np
from typing import Union



def write(arr, n):
    p,l,u = lu(arr)
    lus = l+u - np.eye(arr.shape[0])
    ev,evec = eigh(arr,eigvals_only=False)
    np.savetxt("test_values/csv/{}lu.csv".format(n), f(lus), delimiter=",",fmt='%.6f')
    np.savetxt("test_values/csv/{}inv.csv".format(n), f(inv(arr)), delimiter=",",fmt='%.6f')
    np.savetxt("test_values/csv/{}evec.csv".format(n), f(evec), delimiter=",",fmt='%.6f')
    np.savetxt("test_values/csv/{}eval.csv".format(n), f(ev), delimiter=",",fmt='%.6f')

S1 =np.array([
    [0.9999999908897999, 8.485806523168149e-16, 4.8976976063215463e-11, -8.326316139686609e-11, -8.326316139686609e-11, -8.326316139686609e-11],
    [8.485806523168148e-16, 1.0000000165470366, 0.2367039416499711, 6.311319040196644e-16, 6.311319040196644e-16, 6.311319040196644e-16],
    [4.8976976063215463e-11, 0.23670394164997102, 1.0000000268753377, 0.0, 0.0, 0.0],
    [-8.32631613968661e-11, 6.311319040196645e-16, 0.0, 1.0000000195212955, 0.0, 0.0 ],
    [-8.32631613968661e-11, 6.311319040196645e-16, 0.0, 0.0, 1.0000000195212955, 0.0 ],
    [-8.32631613968661e-11, 6.311319040196645e-16, 0.0, 0.0, 0.0, 1.0000000195212957]],dtype =np.float32)

S2 =np.array([[ 1.000000, 0.000000, 0.000000 ,-0.000000, -0.000000, -0.000000, 0.000002, 0.000255 ,-0.000343 ,-0.000251, -0.000126 ],
    [0.000000, 1.000000, 0.236704, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000 ,0.000000, 0.000000, 0.000000 ],
    [0.000000, 0.236704, 1.000000, -0.000000, -0.000000, -0.000000, 0.000000, 0.000033 ,0.000023, 0.000041, 0.000066 ],
    [-0.000000, 0.000000, -0.000000, 1.000000, 0.000000, 0.000000, -0.000000, -0.000023 ,-0.000006 ,-0.000028, -0.000046 ],
    [-0.000000, 0.000000, -0.000000, 0.000000, 1.000000, 0.000000, -0.000000, -0.000041, -0.000028 ,-0.000041, -0.000083 ],
    [-0.000000, 0.000000, -0.000000, 0.000000, 0.000000, 1.000000, -0.000000, -0.000066, -0.000046 ,-0.000083, -0.000123 ],
    [0.000002, 0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 1.000000, 0.236704, 0.000000, 0.000000, 0.000000 ],
    [0.000255, 0.000000, 0.000033, -0.000023, -0.000041, -0.000066, 0.236704, 1.000000, 0.000001, 0.000000, 0.000000 ],
    [-0.000343, 0.000000, 0.000023, -0.000006, -0.000028, -0.000046, 0.000000, 0.000001, 1.000000, 0.000000, 0.000000 ],
    [-0.000251, 0.000000, 0.000041, -0.000028, -0.000041, -0.000083, 0.000000, 0.000000, 0.000000, 1.000000, -0.000000 ],
    [-0.000126, 0.000000, 0.000066, -0.000046, -0.000083, -0.000123, 0.000000, 0.000000, 0.000000 ,-0.000000, 1.000000 ]],dtype = np.float32)

S3=np.array([[1.000000, 0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 0.000035, 0.036546, -0.048771, -0.035765, -0.017883 ],
[0.000000, 1.000000, 0.236704, 0.000000, 0.000000, 0.000000, 0.000000, 0.002962, 0.001496, 0.002693, 0.004339 ],
[0.000000, 0.236704 ,1.000000, -0.000000, -0.000000, -0.000000, 0.000001, 0.035569, 0.017236, 0.031024, 0.049984 ],
[-0.000000, 0.000000, -0.000000, 1.000000, 0.000000, 0.000000, -0.000001, -0.006672, 0.005238, -0.006534, -0.010527 ],
[-0.000000, 0.000000, -0.000000, 0.000000, 1.000000, 0.000000, -0.000001, -0.012009 ,-0.006534, -0.002893, -0.018948 ],
[-0.000000, 0.000000, -0.000000, 0.000000, 0.000000, 1.000000, -0.000002, -0.019348 ,-0.010527, -0.018948, -0.021659 ],
[0.000035, 0.000000, 0.000001, -0.000001, -0.000001, -0.000002, 1.000000, 0.241137, 0.000000, 0.000000, 0.000000 ],
[0.036546, 0.002962, 0.035569 ,-0.006672, -0.012009, -0.019348, 0.241137, 1.000000, 0.000000, -0.000000 ,0.000000 ],
[-0.048771, 0.001496, 0.017236, 0.005238, -0.006534, -0.010527, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000 ],
[-0.035765, 0.002693, 0.031024, -0.006534, -0.002893, -0.018948, 0.000000, -0.000000, 0.000000, 1.000000 ,0.000000 ],
[-0.017883, 0.004339, 0.049984, -0.010527, -0.018948, -0.021659, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000] ],dtype = np.float32)


S4=np.array([[1.000000, 0.000000 ,0.000000 ,-0.000000 ,-0.000000 ,-0.000000 ,0.000035 ,0.036546, -0.048771, -0.035765 ,-0.017883 ,0.000455, 0.017913, -0.007321, -0.025868, -0.005857 ],
[0.000000, 1.000000, 0.236704 ,0.000000, 0.000000, 0.000000 ,0.000000 ,0.002962 ,0.001496, 0.002693, 0.004339, 0.000000, 0.000000,0.000000 ,0.000000, 0.000000 ],
[0.000000, 0.236704, 1.000000, -0.000000, -0.000000, -0.000000, 0.000001 ,0.035569, 0.017236, 0.031024, 0.049984, 0.000000 ,0.000001 ,0.000001 ,0.000001, 0.000002 ],
[-0.000000, 0.000000 ,-0.000000 ,1.000000, 0.000000, 0.000000, -0.000001, -0.006672, 0.005238 ,-0.006534, -0.010527, -0.000000, -0.000001, -0.000002 ,-0.000001 ,-0.000003 ],
[-0.000000, 0.000000, -0.000000, 0.000000, 1.000000, 0.000000, -0.000001, -0.012009, -0.006534 ,-0.002893 ,-0.018948, -0.000000, -0.000000 ,-0.000001, -0.000000, -0.000001 ],
[-0.000000, 0.000000, -0.000000, 0.000000, 0.000000, 1.000000, -0.000002, -0.019348, -0.010527, -0.018948, -0.021659, -0.000000 ,-0.000001 ,-0.000003, -0.000001 ,-0.000002 ],
[0.000035, 0.000000, 0.000001 ,-0.000001, -0.000001, -0.000002, 1.000000, 0.241137, 0.000000 ,0.000000, 0.000000, 0.000000, 0.006656, 0.013139, -0.002628, 0.002920 ],
[0.036546, 0.002962, 0.035569 ,-0.006672, -0.012009, -0.019348, 0.241137, 1.000000, 0.000000 ,-0.000000 ,0.000000, 0.020798, 0.213670, 0.134812, -0.026962, 0.029958 ],
[-0.048771, 0.001496, 0.017236, 0.005238, -0.006534, -0.010527, 0.000000, 0.000000 ,1.000000 ,0.000000, 0.000000 ,-0.033960 ,-0.313077, -0.176221, 0.054321, -0.060357 ],
[-0.035765, 0.002693, 0.031024, -0.006534, -0.002893 ,-0.018948, 0.000000, -0.000000 ,0.000000 ,1.000000, 0.000000, 0.006792, 0.062615, 0.054321, 0.084520 ,0.012071 ],
[-0.017883, 0.004339, 0.049984, -0.010527, -0.018948,-0.021659, 0.000000 ,0.000000 ,0.000000 ,0.000000 ,1.000000 ,-0.007547, -0.069573, -0.060357, 0.012071 ,0.081972 ],
[0.000455, 0.000000, 0.000000 ,-0.000000 ,-0.000000, -0.000000 ,0.000000, 0.020798, -0.033960 ,0.006792 ,-0.007547, 1.000000, 0.248362 ,0.000000, 0.000000 ,-0.000000 ],
[0.017913, 0.000000, 0.000001 ,-0.000001, -0.000000 ,-0.000001, 0.006656, 0.213670, -0.313077 ,0.062615, -0.069573, 0.248362, 1.000000, -0.000000, 0.000000, 0.000000 ],
[-0.007321, 0.000000, 0.000001, -0.000002, -0.000001, -0.000003, 0.013139, 0.134812, -0.176221 ,0.054321, -0.060357, 0.000000 ,-0.000000, 1.000000, 0.000000 ,0.000000 ],
[-0.025868, 0.000000, 0.000001, -0.000001 ,-0.000000, -0.000001, -0.002628, -0.026962, 0.054321, 0.084520, 0.012071, 0.000000, 0.000000 ,0.000000, 1.000000 ,0.000000 ],
[-0.005857, 0.000000, 0.000002, -0.000003, -0.000001, -0.000002, 0.002920, 0.029958 ,-0.060357, 0.012071, 0.081972 ,-0.000000, 0.000000, 0.000000, 0.000000 ,1.000000 ]],dtype = np.float32)


write(f(S1),1)
write(f(S2),2)
write(f(S3),3)
write(f(S4),4)

