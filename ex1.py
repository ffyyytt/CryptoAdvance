import cupy as np
import time
import random

def print_matrix(matrix):
    print(matrix)

def allocate_matrix(n, m):
    return np.zeros((n, m), dtype=int)

def binary_matrix_inverse(_matrix, n):
    matrix = np.copy(_matrix)
    inverse = np.eye(n, dtype=int)
    # Gaussian
    for i in range(n):
        if matrix[i, i] == 0:
            swap_row = -1
            for j in range(i + 1, n):
                if matrix[j, i] == 1:
                    swap_row = j
                    break
            if swap_row == -1:
                return False, None
            # Rows swap
            matrix[[i, swap_row]] = matrix[[swap_row, i]]
            inverse[[i, swap_row]] = inverse[[swap_row, i]]
        # Remove element below main diagonal
        for j in range(i + 1, n):
            if matrix[j, i] == 1:
                matrix[j, :] = np.bitwise_xor(matrix[j, :], matrix[i, :])
                inverse[j, :] = np.bitwise_xor(inverse[j, :], inverse[i, :])
    # Remove the element above the main diagonal
    for i in range(n - 1, 0, -1):
        for j in range(i - 1, -1, -1):
            if matrix[j, i] == 1:
                matrix[j, :] = np.bitwise_xor(matrix[j, :], matrix[i, :])
                inverse[j, :] = np.bitwise_xor(inverse[j, :], inverse[i, :])
    return True, inverse

def multiply(mat1, mat2, n, k, m, modulo=2):
    res = mat1 @ mat2
    if modulo > 0:
        res %= modulo
    return res

def checkm(m, c, G, n, k, t):
    count = 0
    cp = multiply(m, G, 1, k, n)
    count += np.sum(cp != c)
    return count == t

def generate_random_vector(n, t):
    a = np.zeros(n, dtype=int)
    a[:t] = 1
    np.random.shuffle(a)
    return a

def ISD_encode(m, c, G, n, k, t):
    e = np.expand_dims(generate_random_vector(n, t), axis=0)
    c += multiply(m, G, 1, k, n) + e
    c %= 2

def generate_random_information_set(Gs, G, cs, c, n, k):
    info_set = np.zeros(n, dtype=int)
    indices = random.sample(range(n), k)
    info_set[indices] = 1
    Gs[:,:] = G[:,indices]
    cs[0,:] = c[0,indices]

def ISD_decode(m, c, G, n, k, t, max_iterations=-1):
    cs = allocate_matrix(1, k)
    Gs = allocate_matrix(k, k)
    inv_Gs = allocate_matrix(k, k)
    count = 0
    while True:
        count += 1
        generate_random_information_set(Gs, G, cs, c, n, k)
        inv_exists, inv_Gs = binary_matrix_inverse(Gs, k)
        if inv_exists:
            m[:,:] = multiply(cs, inv_Gs, 1, k, k)
            if checkm(m, c, G, n, k, t):
                return True
        if count == max_iterations:
            break
    return False

def test_ISD(n, k, t):
    G = allocate_matrix(k, n)
    m = allocate_matrix(1, k)
    c = allocate_matrix(1, n)
    mdecode = allocate_matrix(1, k)

    # random values
    for i in range(k):
        G[i, :] = generate_random_vector(n, int(0.6*n))
    m[0, :] = generate_random_vector(k, int(k/2))

    ISD_encode(m, c, G, n, k, t)
    ISD_decode(mdecode, c, G, n, k, t)

def main():
    n = 1000
    k = 500
    t = 10
    start = time.time()
    test_ISD(n, k, t)
    print(f"Time taken: {time.time() - start}s")

if __name__ == "__main__":
    main()