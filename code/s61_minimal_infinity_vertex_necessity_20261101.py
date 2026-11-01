#!/usr/bin/env python3
"""All M0/M1 infinity points have modular image on the null-correlation curve."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61,X
from s61_nil22_infinity_20261101 import tensor,beta_matrix,sm
ROOT=Path(__file__).resolve().parents[1]
def run():
    m=S61();x,y,z,w=X[1:];g=s.symbols('g0:9');b=s.symbols('b0:16');G=sum(a*c for a,c in zip(g,[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3],strict=True));out=[]
    for label,offset,constant1,constant2 in [('M0',-s.Rational(6,5),-s.Rational(14,5),-s.Rational(32,25)),('M1',-s.Rational(1,5),s.Rational(1,5),-s.Rational(2,25))]:
        Q=tensor([x*x,0,0,-s.diff(G,w),s.diff(G,z),0]) if label=='M0' else tensor([x*x/2,0,0,-x*z-s.diff(G,w),s.diff(G,z),x*x/2])
        M=beta_matrix(m,Q)[70:,:];other=[j for j in range(16) if j not in [0,2]]
        keys=sorted({(j,e) for i in range(64) for j in other for e,c in s.Poly(M[i,j],*g).terms()});index={key:i for i,key in enumerate(keys)};C=s.zeros(len(keys),64)
        for i in range(64):
            for j in other:
                for e,c in s.Poly(M[i,j],*g).terms():C[index[j,e],i]=c
        K=s.Matrix.hstack(*C.nullspace());assert C*K==s.zeros(C.rows,K.cols)
        assert (K.T*M[:,other]).applyfunc(s.expand)==s.zeros(K.cols,len(other))
        N=(K.T*M[:,0]).applyfunc(s.expand);relations=list(dict.fromkeys(f for f in N if f!=0));gb=s.groebner(relations,*g,order='grevlex')
        expected=[g[2],g[3]-offset,g[5],g[6],g[7]];assert list(gb)==expected
        sub={g[2]:0,g[3]:offset,g[5]:0,g[6]:0,g[7]:0,b[0]:1,b[2]:0};F=M*s.Matrix(b)
        row86=s.factor(F[16].subs(sub));row72=s.factor(F[2].subs(sub))
        assert row86==constant1*b[3] and row72.subs(b[3],0)==constant2!=0
        out.append({'label':label,'Q':sm(Q),'beta_matrix':sm(M),'normalization':'b0=u4=1,b2=u12=0','cancelled_beta_columns':other,'constant_row_combinations':sm(K.T),'necessary_polynomials':list(map(str,relations)),'exact_groebner_basis':list(map(str,gb)),'remaining_row86':str(row86),'remaining_row72':str(row72),'row72_after_row86':str(constant2),'contradiction':True})
    inputs=[Path(__file__),ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_infinity_flags.py',ROOT/'code/s61_model.py',m.source/'symplectic_slice_reduced_equations.json']
    return {'id':'S61_MINIMAL_INFINITY_VERTEX_NECESSITY_20261101','all_passed':True,'field':'QQ','g_variables':list(map(str,g)),'beta_coordinate_order':[f'u{i+1}' for i,w in enumerate(m.reduced['parameter_weights']) if w==1],'families':out,'scope':'Set-theoretic necessity tau_beta(image(modular Q))=0 for all M0/M1 solutions of the134 leading equations. No assumptions on beta or discriminant. Does not classify solutions on tau=0 or saturated infinity.','proof_note':'research/notes/S61_MINIMAL_INFINITY_VERTEX_NECESSITY_20261101.md','input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'M0_contradiction':True,'M1_contradiction':True,'output':str(args.output)}))
if __name__=='__main__':main()
