import scipy
import random
import cupy as np

def generate_random_vector(n, t):
    vector = np.zeros(n, dtype=int)
    vector[:t] = 1
    np.random.shuffle(vector)
    return vector
    
def BitFlipping(h0, h1, s, T, t):
    iter = 0
    n = h0.shape[0]
    u = np.zeros(n, dtype=bool)
    v = np.zeros(n, dtype=bool)
    H = np.hstack([h0, h1])
    flipped_position = np.zeros(2*n, dtype=bool)
    syndrome = np.copy(s)

    lastlasts, lasts = 0, 0
    
    while np.sum(u) != t or np.sum(v) != t and np.sum(syndrome) != 0:
        iter += 1
        sum = np.matmul(syndrome,H) # No modular reduction, values in Z
        
        for i in range(2*n):
            flipped_position[i] = sum[i] >= T
                
        u = np.logical_xor(u, flipped_position[:n])
        v = np.logical_xor(v, flipped_position[n:])
        
        syndrome = np.logical_xor(syndrome, np.matmul(H, flipped_position.T))

        if np.sum(flipped_position) == 0:
            T -= 1
        elif np.sum(syndrome) == lastlasts or np.sum(syndrome) == lasts:
            T = (T+1)*2
            u = np.zeros(n, dtype=bool)
            v = np.zeros(n, dtype=bool)
            syndrome = np.copy(s)

        lastlasts = lasts
        lasts = np.sum(syndrome)
    
    if np.sum((s - np.hstack([h0, h1]) @ np.hstack([u, v]))%2) != 0:
        return None
    else:
        return (u, v)
        
e0 = generate_random_vector(4813, 39)
e1 = generate_random_vector(4813, 39)

h0 = generate_random_vector(4813, 39)
h1 = generate_random_vector(4813, 39)

mh0 = scipy.linalg.circulant(h0)
mh1 = scipy.linalg.circulant(h1)

s = (np.dot(mh0, e0) + np.dot(mh1, e1)) % 2
a = BitFlipping(mh0, mh1, s, 26, 39)