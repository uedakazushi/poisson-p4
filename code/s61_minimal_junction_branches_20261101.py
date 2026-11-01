#!/usr/bin/env python3
"""Actual branches and same-branch length lower bounds at minimal junctions."""
import sys,json,time,argparse,hashlib
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action,minor_witness
ROOT=Path(__file__).resolve().parents[1]
def run():
    m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];src=json.load(open(ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json'));fam=json.load(open(ROOT/'reports/S61_MINIMAL_NILPOTENT_INFINITY_FAMILY_20261101.json'));g=s.symbols('g0:9');k,K,L,U,V=s.symbols('k K L U V');b=s.symbols('b0:6');pars=[k,K,L,U,V]+list(b);bv=s.Matrix([0,0,0,0,k,b[0],U,b[1],b[2],V,b[3],b[4],b[5],2*k+1,K-5*V,5*U-L]);zz,uu,Beta=primitive_beta(m);vec=s.Matrix([Beta[i,j].diff(zz[h]) for (i,j),h in COORDS]);bas=vec.jacobian(uu);rr=list(bas.T.rref()[1]);left=bas[rr,:].inv();actions=s.Matrix.hstack(*[left*infinitesimal_action(vec,H)[rr,:] for H in sp_basis()]);x,y,z,w=X[1:];vv=s.Matrix([x,y,z,w]);out=[]
    for kind,kk,expected in [('M0','1/6',4),('M0','-11/12',4),('M0','-1',4),('M0','-2/3',8),('M1','-11/12',3),('M1','-1',4),('M1','-4/9',5),('M1','-1/3',6)]:
     block=next(a for a in src['families'] if a['label']==kind);qexpr=s.Matrix(list(map(s.sympify,(fam if kind=='M0' else fam['M1_family'])['q_coordinates'])));base_sol={s.sympify(a):s.sympify(v) for a,v in block['baseline_solution'].items()};qmain=qexpr.subs(base_sol,simultaneous=True);sample={k:s.Rational(kk),K:19,L:5,U:2,V:3,**dict(zip(b,[2,3,5,7,11,13]))};q=qmain.subs(sample);beta=bv.subs(sample);p=s.zeros(61,1)
     for i,v in zip(one,beta):p[i]=v
     for i,v in zip(two,q):p[i]=v
     local=calculate(m,True,p);T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);tr=list(T.T.rref()[1]);ti=T[tr,:].inv();cols=local['chart_columns'];fixed=next(i for i in range(45) if i not in cols)
     BD=bv.jacobian(pars).row_join(actions.subs(dict(zip(uu,beta))));assert BD.rank()==16 and BD[:,list(range(1,21))].rank()==16
     shape=[base_sol[gi].subs(sample) for gi in g];G=sum(a*f for a,f in zip(shape,[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3]));Q=s.zeros(4);Q[0,1]=x*x if kind=='M0' else x*x/2;Q[2,3]=0 if kind=='M0' else x*x/2;Q[1,2]=-s.diff(G,w)-(x*z if kind=='M1' else 0);Q[1,3]=s.diff(G,z);Q-=Q.T;qd=[]
     for H in sp_basis():
      dQ=sum(((H*vv)[i]*Q.diff(vv[i]) for i in range(4)),s.zeros(4))-H*Q-Q*H.T;P=X[0]**2*J;P[1:,1:]+=dQ;v=m.normalize(P)['parameters'];v=s.Matrix([v[i] for i in two]);v-=q*v[fixed]/q[fixed];qd.append(v[cols,:])
     DQ=qmain.jacobian(pars).subs(sample);DQ-=q*DQ[fixed,:]/q[fixed];QD=DQ[cols,:].row_join(s.Matrix.hstack(*qd));N=s.Matrix.hstack(*BD.nullspace());C=QD*N;assert C.rank()==1;ix=next(i for i in range(C.cols) if N[0,i]!=0 and C[:,i]!=s.zeros(44,1));tv=C[:,ix].col_join(s.zeros(1,1));vt=ti*tv[tr,:];assert T*vt==tv;branches=[{'kind':'main','tangent':list(map(str,vt)),'lower':3 if kind=='M0' else 2,'dk':str(N[0,ix]),'source_direction':list(map(str,N[:,ix]))}]
     _,F=specialize(m,dict(zip(one,beta)),True);t=s.Symbol('t')
     curves=[('tail',c) for c in block['rational_tails'] if c['k']==kk]+[('free',c) for c in block['baseline_free_lines'] if c['k']==kk]
     for label,c in curves:
      sol={s.sympify(a):s.sympify(v) for a,v in c['coordinate_formulas'].items()};axis=c['axis'];sol[g[axis]]=t;qt=qexpr.subs(sol,simultaneous=True).subs(g[axis],t).subs(sample);center=0 if label=='tail' else base_sol[g[axis]].subs(sample);assert (qt.subs(t,center)-q).applyfunc(s.cancel)==s.zeros(45,1)
      assert all(s.expand(sum(v*s.prod(qt[i] for i in mon) for mon,v in f.items() if len(mon)==2))==0 for f in F)
      dq=qt.diff(t).subs(t,center);dq-=q*dq[fixed]/q[fixed];tv=dq[cols,:].col_join(s.zeros(1,1));vt=ti*tv[tr,:];assert T*vt==tv and vt!=s.zeros(T.cols,1)
      lower=(4 if kind=='M0' else 3) if label=='free' else (3 if kind=='M1' and kk=='-1/3' and axis==2 else 1)
      if label=='free' or lower==3:
       target_path='reports/S61_MINIMAL_BASELINE_THICKENINGS_20261101.json' if label=='free' else 'reports/S61_MINIMAL_CUBIC_TAIL_20261101.json'
       target=json.loads((ROOT/target_path).read_text());known=target['cases'][kind]['point']['local_data'] if label=='free' else target['point']['local_data'];check=p.copy()
       for qi,oldi in enumerate(two):check[oldi]=qt[qi].subs(t,1)
       assert check==s.Matrix(known['point61']).applyfunc(s.Rational)
      branches.append({'kind':label,'axis':axis,'tangent':list(map(str,vt)),'center':str(center),'lower':lower,'q_curve':list(map(str,qt))})
     assert sum(c['lower'] for c in branches)==expected;rank=s.Matrix.hstack(*[s.Matrix(c['tangent']).applyfunc(s.Rational) for c in branches]).rank();assert rank==len(branches)
     result={'kind':kind,'k':kk,'length_lower':expected,'branches':branches,'beta_minor':minor_witness(BD),'restricted_k_beta_minor':minor_witness(BD[:,list(range(1,21))]),'point61':list(map(str,p)),'beta_differential':[[str(a) for a in row] for row in BD.tolist()],'main_q_differential':[[str(a) for a in row] for row in QD.tolist()],'local':local};out.append(result);print(kind,kk,[(a['kind'],a['lower']) for a in branches],'rank',rank,flush=True)

    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json',ROOT/'reports/S61_MINIMAL_NILPOTENT_INFINITY_FAMILY_20261101.json',ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json',ROOT/'reports/S61_MINIMAL_CUBIC_TAIL_20261101.json',ROOT/'reports/S61_MINIMAL_BASELINE_THICKENINGS_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_MINIMAL_JUNCTION_BRANCHES_20261101',all_passed=True,field='QQ',cases=out,scope='Actual fixed-beta branches and generic length lower bounds at eight specified minimal junctions; no upper bound or saturation conclusion in this certificate alone.',proof_note='research/notes/S61_MINIMAL_JUNCTION_BRANCHES_20261101.md',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'cases':len(r['cases'])}))
if __name__=='__main__':main()
