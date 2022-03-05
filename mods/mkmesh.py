import numpy as np

def mkmesh(width, hight, w_size, h_size):
    N_w = int(width/w_size) + 2
    N_h = int(hight/h_size) + 2
    meshed = np.zeros((N_w, N_h))
    return meshed