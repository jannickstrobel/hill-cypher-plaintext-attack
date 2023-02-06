import numpy as np
from mpmath.libmp.backend import xrange
from sympy import Matrix


class hill:

    # ABCDEFGHIJKLMNOPQRSTUVWXYZ012
    # THHE
    # UYTK

    def __init__(self):

        self.alphabet = list(input("Input alphabet in order: "))

        input_plain = list(input("Input 4 plaintext characters: "))
        self.plain = np.array([[input_plain[0], input_plain[1]],
                               [input_plain[2], input_plain[3]]])

        input_cipher = list(input("Input 4 ciphertext characters: "))
        self.cipher = np.array([[input_cipher[0], input_cipher[1]],
                                [input_cipher[2], input_cipher[3]]])

        self.plain = np.array([[self.alphabet.index(self.plain[0][0]), self.alphabet.index(self.plain[1][0])],
                               [self.alphabet.index(self.plain[0][1]), self.alphabet.index(self.plain[1][1])]])
        self.cipher = np.array([[self.alphabet.index(self.cipher[0][0]), self.alphabet.index(self.cipher[1][0])],
                                [self.alphabet.index(self.cipher[0][1]), self.alphabet.index(self.cipher[1][1])]])

        self.plain_t = self.plain.transpose((1, 0))
        self.cipher_t = self.cipher.transpose((1, 0))


    def calculate(self):
        matrix = np.concatenate((self.cipher_t, self.plain_t), axis=1)
        matrix = self.modrref(matrix, 29)
        matrix1, matrix2 = np.split(matrix, 2, axis=1)
        A_inv = matrix2.transpose((1, 0))
        #print("A^-1:\n", A_inv)
        A = self.matInvMod(A_inv, len(self.alphabet))
        print("Original key:\n", A)



    #############################

    def modrref(self, M, mod):
        '''
        Computes the row-reduced echelon form of the matrix M modulo mod.
        '''
        r = 0
        while r < M.shape[0]:
            # Ignore non-zero rows.
            try:
                f = self.firstnonzero(M[r])
            except:
                r += 1
                continue

            # Rule 1: Swap with the row above if out of order.
            if r > 0:
                swap = False
                try:
                    g = self.firstnonzero(M[r - 1])
                except:
                    swap = True
                if swap:
                    self.swaprows(M, r, r - 1)
                    continue

            # Rule 2: Normalize each row
            self.normrow(M, r, mod)

            # Rule 3: Subtract it from the others
            self.subrow(M, r)
            r += 1


        return self.mod(M, mod)

    def mod(self, M, mod):
        '''
        returns matrix with elements % mod
        '''
        shape = M.shape
        result = np.zeros(shape, dtype=int)
        for x in range(0, shape[0]):
            for y in range(0, shape[1]):
                result[x, y] = M[x, y] % mod
        return result

    def firstnonzero(self, row):
        '''
        Finds the index of the first non-zero element of the row.
        '''
        for i, r in enumerate(row):
            if r != 0: return i
        else:
            raise Exception('No non-zeros elements found.')

    def swaprows(self, M, i, j):
        '''
        Swaps rows i and j of the matrix M.
        '''
        M[i], M[j] = M[j], M[i]

    def subrow(self, M, i):
        '''
        Subtracts row i from each other row in the matrix M.
        Assumes that the first non-zero element of i is a 1.
        '''
        f = self.firstnonzero(M[i])
        for j in xrange(M.shape[0]):
            if i == j: continue
            M[j] -= M[j, f] * M[i]

    def normrow(self, M, i, mod):
        '''
        Normalizes row i of the matrix M such that the first non-zero element is 1.
        '''
        f = self.firstnonzero(M[i])
        M[i] *= self.modinv(M[i, f], mod)
        M[i] %= mod

    def modinv(self, n, mod):
        '''
        Lazily finds the multiplicative inverse of n modulo mod.
        '''
        for x in xrange(1, mod):
            if (n * x) % mod == 1:
                return x
        else:
            raise ArithmeticError('%i has no multiplicative inverse modulo %i.' % (n, mod))

    def matInvMod(self, vmnp, mod):
        nr = vmnp.shape[0]
        nc = vmnp.shape[1]
        if (nr != nc):
            print
            "Error: Non square matrix! exiting"
            exit()
        vmsym = Matrix(vmnp)
        vmsymInv = vmsym.inv_mod(mod)
        vmnpInv = np.array(vmsymInv)
        print
        "vmnpInv: ", vmnpInv, "\n"
        k = nr
        vmtest = [[1 for i in range(k)] for j in range(k)]  # just a 2-d list
        vmtestInv = vmsym * vmsymInv
        for i in range(k):
            for j in range(k):
                # print i, j, vmtrx2[i,j] % mod
                vmtest[i][j] = vmtestInv[i, j] % mod
        print
        "test vmk*vkinv % mod \n:", vmtest
        return vmnpInv


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    h = hill()
    h.calculate()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
