#!/usr/bin/env python3
# ajustar Kappa para N=4
# sweep_kappa_refined.py
"""
Barrido refinado de kappa 
para N=4 (p=7, spinor_dim=4).
Usa SLQ para estimar d_eff.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import eigh
import math, time, os
import matplotlib.pyplot as plt
import pandas as pd

# ---------------- Parámetros fijos ----------------
p = 7
N = 4
spinor_dim = 4
C = 1.0
phi = (1 + math.sqrt(5)) / 2
alpha = 1.0 / phi

num_probes = 40
m_lanczos = 80
t_grid = np.logspace(-6, 0, 40)   # rango seguro

rng_seed = 123456
kappa_values = np.arange(1.09244, 1.09245, 0.000001)  # refinado

out_dir = "tdp_kappa_refined_N4"
os.makedirs(out_dir, exist_ok=True)

# ---------------- Funciones auxiliares ----------------
def circulant_difference(m):
    row = np.arange(m)
    col_plus = (row + 1) % m
    col_minus = (row - 1) % m
    data = np.concatenate([0.5j*np.ones(m), -0.5j*np.ones(m)])
    rows = np.concatenate([row, row])
    cols = np.concatenate([col_plus, col_minus])
    M = sp.coo_matrix((data, (rows, cols)), shape=(m, m)).tocsr()
    M = 0.5*(M + M.getH())
    return M

def embedding_matrix(m_from, m_to):
    repeats = m_to // m_from
    rows=[]; cols=[]; data=[]
    for i in range(m_from):
        for k in range(repeats):
            j = i + k*m_from
            rows.append(j); cols.append(i); data.append(1.0/math.sqrt(repeats))
    return sp.coo_matrix((data,(rows,cols)), shape=(m_to,m_from)).tocsr()

def build_D_sparse(p=7,N=4,spinor_dim=4,kappa=1.0,C=1.0,alpha=1.0/((1+math.sqrt(5))/2)):
    sizes=[p**n for n in range(1,N+1)]
    if spinor_dim==4:
        gamma = sp.csr_matrix(np.diag([1.0,1.0,-1.0,-1.0]))
    else:
        gamma = sp.eye(spinor_dim, format='csr')
    blocks=[]
    for n,m in enumerate(sizes, start=1):
        s_n = p**(kappa*n)
        diff = circulant_difference(m)
        Dn = sp.kron(diff*s_n, gamma, format='csr')
        blocks.append(Dn)
    Nblocks=len(blocks)
    block_matrix=[[None]*Nblocks for _ in range(Nblocks)]
    for i in range(Nblocks):
        block_matrix[i][i]=blocks[i]
    for i in range(Nblocks-1):
        m_from=sizes[i]; m_to=sizes[i+1]
        E=embedding_matrix(m_from,m_to)
        M_small_to_large=(C/(p**alpha))*sp.kron(E, sp.eye(spinor_dim), format='csr')
        block_matrix[i][i+1]=M_small_to_large.getH()
        block_matrix[i+1][i]=M_small_to_large
    D_big=sp.bmat(block_matrix, format='csr')
    D_big=0.5*(D_big + D_big.getH())
    return D_big

def make_M_linear_op(D):
    n=D.shape[0]; dtype=D.dtype
    def mv(v):
        v=np.asarray(v,dtype=dtype)
        return D.dot(D.dot(v))
    return spla.LinearOperator(shape=(n,n), matvec=mv, dtype=complex if dtype==complex else float)

def lanczos_tridiag(Mop,z,m):
    n=z.size
    q_prev=np.zeros(n,dtype=complex)
    q=z/np.linalg.norm(z)
    alphas=np.zeros(m,dtype=float); betas=np.zeros(max(0,m-1),dtype=float)
    for j in range(m):
        w=Mop.matvec(q)
        alpha=np.real(np.vdot(q,w))
        w=w-alpha*q-(betas[j-1]*q_prev if j>0 else 0.0)
        alphas[j]=alpha
        if j<m-1:
            beta=np.linalg.norm(w)
            if beta<1e-14: betas[j:]=0.0; return alphas,betas
            betas[j]=beta
            q_prev,q=q,w/beta
    return alphas,betas

def slq_trace_expm_of_M(Mop,t_grid,num_probes=32,m_lanczos=80,seed=123):
    rng=np.random.default_rng(seed); n=Mop.shape[0]
    traces=np.zeros(len(t_grid),dtype=float)
    for k in range(num_probes):
        z=rng.choice([1.0,-1.0],size=(n,)).astype(float)
        alphas,betas=lanczos_tridiag(Mop,z.astype(complex),m_lanczos)
        nonzero_beta=np.nonzero(np.abs(betas)>1e-16)[0]
        if nonzero_beta.size>0: last=nonzero_beta[-1]+1; m_eff=min(len(alphas),last+1)
        else: m_eff=1
        T=np.diag(alphas[:m_eff])
        if m_eff>1:
            off=betas[:m_eff-1]
            T+=np.diag(off,1)+np.diag(off,-1)
        evals_T,U=eigh(T)
        e1=np.zeros((m_eff,)); e1[0]=1.0
        w=(U.T@e1)**2
        for i,t in enumerate(t_grid):
            lamb=np.maximum(evals_T,0.0)
            traces[i]+=np.sum(w*np.exp(-t*lamb))
    traces/=num_probes
    return traces

def estimate_d_eff(t_grid,trace_est):
    mask=(t_grid>=1e-5)&(t_grid<=1e-2)
    if mask.sum()<5: mask=np.arange(min(6,len(t_grid)))
    x=np.log(t_grid[mask]); y=np.log(np.maximum(trace_est[mask],1e-300))
    slope=np.polyfit(x,y,1)[0]
    return -2*slope

# ---------------- Main ----------------
def main():
    results=[]
    for kappa in kappa_values:
        print(f"\n=== κ={kappa:.6f} ===")
        D=build_D_sparse(p,N,spinor_dim,kappa,C,alpha)
        Mop=make_M_linear_op(D)
        trace_est=slq_trace_expm_of_M(Mop,t_grid,num_probes=num_probes,m_lanczos=m_lanczos,seed=rng_seed)
        trace_est=np.real(trace_est)
        d_eff=estimate_d_eff(t_grid,trace_est)
        print(f"d_eff = {d_eff:.5f}")
        results.append({"kappa":kappa,"d_eff":d_eff})
    df=pd.DataFrame(results)
    csv_path=os.path.join(out_dir,"d_eff_vs_kappa_refined.csv")
    df.to_csv(csv_path,index=False)
    plt.figure(figsize=(6,4))
    plt.plot(df["kappa"],df["d_eff"],marker='o')
    plt.axhline(4,color='r',ls='--',label="target d_eff=4")
    plt.xlabel("kappa"); plt.ylabel("d_eff")
    plt.title("Barrido refinado d_eff vs kappa (N=4, p=7)")
    plt.legend(); plt.grid(True,ls=':')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir,"d_eff_vs_kappa_refined.png"),dpi=200)
    print("Resultados guardados en",out_dir)

if __name__=="__main__":
    main()
