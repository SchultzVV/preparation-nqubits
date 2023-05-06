import numpy as np
import math
from scipy.linalg import lapack
from src.rdmg import rdm_ginibre

def get_abs_lambda(rho):

    n = np.shape(rho)[0]
    w, u, _ = lapack.zheev(rho, overwrite_a=True)
    psi = 0
    for i in range(len(u)):
        u[:, i] = u[:, i]/np.linalg.norm(u[:, i])
    for i in range(n):
        #rho += w[i] * np.outer(u[:, i], u[:, i].conj()) # se precisar do rho
        psi += math.sqrt(w[i])*np.eye(n)[:, i]
    return psi, u
import sys
from qiskit import QuantumRegister,ClassicalRegister,QuantumCircuit

def get_psi_circuit(psi):
    n_qubits = np.shape(psi)[0]-1
    print(n_qubits)
    qr = QuantumRegister(n_qubits)
    cr = ClassicalRegister(n_qubits)
    print(psi)
    for i in psi:
        print(i)
    #print('qr[0]',[qr[0]])
    if n_qubits == 1:
        lista = [qr[0]]
    else:
        lista = [i for i in qr]#[,qr[1]]
    print('lista',lista)
    #print('qr',[[qr[0], qr[1], qr[2], qr[3], qr[4], qr[5]]])
    #sys.exit()
    #sys.exit()
    qc = QuantumCircuit(qr,cr, name='qc_initialize')
    qc.initialize(psi, lista)
    return qc


if __name__ == "__main__":
    #for i in range(1):
    #    rho = rdm_ginibre(2)
    #    print(get_abs_lambda(rho))
    rho = rdm_ginibre(4)
    psi, u = get_abs_lambda(rho)
    qc = get_psi_circuit(psi)
    print(qc)