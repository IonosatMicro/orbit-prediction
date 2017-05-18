import numpy as np


def solve(A, b):
    size = A.shape[0]
    B = np.zeros(A.shape)
    C = np.zeros(A.shape)
    y = np.zeros(size)
    X = np.zeros(size)


    for i in range(size):
        B[i, 0] = A[i, 0]
        C[0, i] = A[0, i] / B[0, 0]
    for i in range(size):
        for j in range(size):
            if i >= j > 0:
                B[i, j] = A[i, j] - sum(B[i, k] * C[k, j] for k in range(j))
            if 0 < i <= j:
                C[i, j] = (A[i, j] - sum(B[i, k] * C[k, j] for k in range(i))) / B[i, i]

    y[0] = b[0] / B[0, 0]
    for i in range(1, size):
        y[i] = (b[i] - sum(B[i, k] * y[k] for k in range(i))) / B[i, i]
    X[size - 1] = y[size - 1]
    for i in range(size - 2, -1, -1):
        X[i] = y[i] - sum(C[i, k] * X[k] for k in range(i + 1, size))

    return X


def prediction_coeff(k):

    a = np.zeros((k, k))

    a[0, :] = np.ones(k)
    for i in range(1, k):
        a[i, :] = np.arange(k)
        for j in range(0, k):
            a[i, j] **= i


    b = np.zeros(k)
    for i in range(k):
        b[i] = (- 1)**i / (i + 1)

    return a, b


def corrector_coeff(k):

    a = np.zeros((k, k))

    a[0, :] = np.ones(k)
    for i in range(1, k):
        a[i, :] = np.arange(-1, k - 1)
        for j in range(0, k):
            a[i, j] **= i
    # a[:, 0], a[:, 1] = np.copy(a[:, 1]), np.copy(a[:, 0])

    b = np.zeros(k)
    b[0] = 1
    for i in range(1, k):
        b[i] = (- 1)**i / (i + 1)
    # b[0], b[1] = np.copy(b[1]), np.copy(b[0])

    return a, b


def create_coeff():
    p_coeff = open("Adams_Predictor.txt", 'w')
    c_coeff = open("Adams_Corrector.txt", 'w')
    p_coeff.truncate()
    c_coeff.truncate()
    np.set_printoptions(precision=15, linewidth=2500)
    for order in range(1, 43):
        A, b = prediction_coeff(order)
        S = solve(A, b)
        line = str(order)
        line += str(S)[1:-1] + "\n"
        p_coeff.write(line)
        A, b = corrector_coeff(order)
        S = solve(A, b)
        line = str(order)
        line += str(S)[1:-1] + "\n"
        c_coeff.write(line)
    p_coeff.close()
    c_coeff.close()









