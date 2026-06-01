import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import LinearOperator, eigsh
from numpy.fft import fft, ifft


class FastHankelProduct(LinearOperator):
    """
    Parameters
    ----------
    h: impulse response without the zero lag term, np.ndarray, shape = (n, )
    """

    def __init__(self, h: np.ndarray, dtype='float32'):
        self.h = h
        self.shape = (len(h), len(h))
        self.dtype = np.dtype(dtype)
        self.n = len(self.h)

    def _matvec(self, x):
        """
        Parameters
        x: np.ndarray, shape = (n, )
        """
        chat = np.hstack([self.h[-1], np.zeros(self.n-1), self.h[:-1]])
        Pi = np.eye(self.n)[::-1]    # permutation matrix
        xhat = np.hstack([Pi.dot(x), np.zeros(self.n-1)])
        y = ifft(fft(chat)*fft(xhat))
        return y[:self.n].real


def reduce_fft_truncate(h: np.ndarray, p: int):
    """
    Model reduction using fast hankel matrix vector multiplication
    that can be calculated by FFT

    Parameters
    ----------
    h: array-like, shape = (n, )
    p:  int, the reduced dimension of our model

    Returns
    -------
    lambd: array-like, shape = (p, )
    """
    fft_operator = FastHankelProduct(h)
    evals_large, evecs_large = eigsh(fft_operator, p, which='LM')
    j = np.argsort(np.abs(evals_large))[::-1]
    v, r = la.qr(evecs_large[:, j])
    A = np.dot(v.T[:p, 1:], v[:-1, :p])
    lambd, _ = la.eig(A)
    return lambd


def reduce_svd_truncate(h: np.ndarray, order: int):
    """
    Reduce the system based on the given number of poles
    using svd method
    """
    if h.ndim == 2:
        raise ValueError('h needs to be raveled.')

    n = len(h)
    Gamma = la.hankel(h)
    u, s, v = la.svd(Gamma)

    A = np.dot(v[:order, 1:], v.T[:-1, :order])
    A = np.real(A)

    lambd, _ = la.eig(A)
    return lambd


def toeplog(h: np.ndarray):
    """
    Convert impulse response to the power series coefficients
    of log transfer function.

    Note: in this calculation, h has h_0 which is 1, and a_0
    should be ignored in the reduction step.
    """
    num_ir = len(h)
    d = np.arange(num_ir-1) + 1
    b = d * h[1:]

    c = h[:-1]
    r = np.zeros(num_ir-1)
    r[0] = h[0]
    x = la.solve_toeplitz((c, r), b)
    a_ = x / d

    a = np.insert(a_, 0, np.log(h[0]), axis=0)
    return a


def info_svd_reduce(h: np.ndarray, order: int, rho: float) -> np.ndarray:
    """
    System reduction by minimizing the log transfer function

    Here h is the positive lag impulse response
    """
    if h.ndim == 2:
        raise ValueError('h needs to be raveled.')

    h_ = np.insert(h, 0, 1, axis=0)    # convert to unity zero lag response
    hr = h_ * (rho ** (np.arange(0, len(h_))))
    a = toeplog(h_)
    lambd = reduce_svd_truncate(a[1:], order) / rho
    return lambd


def blkhankel(ir: np.ndarray): 
    """
    Construct block hankel matrix from given 
    impulse response. 
    """
    di, do, l = ir.shape 
    k = int(np.floor((l + 1) / 2))
    H = np.zeros((di*k, do*k))
    for i in range(k): 
        for j in range(k): 
            H[i*di:(i+1)*di, j*do:(j+1)*do] = ir[:, :, i+j]
    return H 


def lschur(A: np.ndarray): 
    # Shur triangularization and convert A to lower 
    T, Z = la.schur(A, output='real')
    L = T[::-1, ::-1]
    Q = Z[:, ::-1]
    return L, Q 


def msvdreduce(ir: np.ndarray, order: int): 
    do, di, num_ir = ir.shape 
    H = blkhankel(ir)

    U, S, Vh = la.svd(H)    # U.dot(S).dot(VT)
    Vh = Vh[:order, :]      # order reduction by truncation 
    A = Vh[:, di:].dot(Vh.T[:-di, :])

    A, Q = lschur(A)
    B = Q.T.dot(Vh[:, :di])
    w = np.diag(A)
    u = np.zeros((di, order, 1))

    for i in range(order): 
        ui = (B[0] / np.linalg.norm(B[0], ord=2)).reshape(-1, 1)
        u[:, i] = ui 
        tmp = np.eye(di) - (1 + np.conj(w[i]))*ui.dot(ui.T)
        B = B[1:, :]
    return w, u[:, :, 0]


