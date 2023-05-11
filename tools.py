from qiskit import QuantumCircuit, Aer, execute
import numpy as np
import math
from scipy.linalg import lapack
from src.rdmg import rdm_ginibre
import sys
from qiskit import QuantumRegister,ClassicalRegister,QuantumCircuit

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


def get_psi_circuit(psi, rho):
    n_qubits = int(np.ceil(np.log2(np.shape(rho)[0])))
    d = int(2*n_qubits)
    print(n_qubits)

    qr = QuantumRegister(d)
    cr = ClassicalRegister(n_qubits)
    print(psi)
    for i in psi:
        print(i)
    #print('qr[0]',[qr[0]])
    if n_qubits == 1:
        lista = [qr[0]]
    else:
        lista = []
        for i in range(n_qubits):
            #lista = [i for i in qr]#[:-1]#[,qr[1]]
            lista.append(qr[i])
        lista2 = [qr[0],qr[1]]
    print('lista',lista)
    print(len(qr))
    print(qr)
    print('lista2',lista2)
    #print('qr',[[qr[0], qr[1], qr[2], qr[3], qr[4], qr[5]]])
    #sys.exit()
    #sys.exit()
    qc = QuantumCircuit(qr,cr, name='qc_initialize')

    qc.initialize(psi, lista)
    for i in range(int(n_qubits)):
            qc.cx(i, int(n_qubits)+i)
    return qc

def douglis(n):
    n = 3
    qubits = 2*n
    #qr = QuantumRegister(qubits)
    qc = QuantumCircuit(qubits, qubits, name = 'mistura')

    def qc_CNOTS_mistura():
    #  CNOT application for each pair of qubits
        for i in range(n):
            qc.cx(i, n+i)
       # circuito.cx(n+i, i)
        return qc
    qc_CNOTS_mistura=qc_CNOTS_mistura()
    qc_CNOTS_mistura.draw('mpl')


if __name__ == "__main__":
    #for i in range(1):
    #    rho = rdm_ginibre(2)
    #    print(get_abs_lambda(rho))
    rho = rdm_ginibre(16)
    psi, u = get_abs_lambda(rho)
    qc = get_psi_circuit(psi, rho)
    print(qc)