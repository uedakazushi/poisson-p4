#!/usr/bin/env python3
"""Actual C2 Loewy filtration and polynomial eleven-unit cancellation.

Only the accepted finite matrices are used. No actual comparison T2 value.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,product_columns,eye,srec
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {i:F(c) for i,c in v}
def padd(a,b,c=F(1)):return add(a,b,c)
def pmul(a,b):
    o={}
    for m,c in a.items():
        for n,d in b.items():
            p=tuple(x+y for x,y in zip(m,n))
            if p[-1]<2:o=add(o,{p:c*d})
    return o
def prec(p):return [[list(m),str(c)] for m,c in sorted(p.items())]
def vadd(a,b,c=F(1)):
    o={i:dict(p) for i,p in a.items()}
    for i,p in b.items():
        q=add(o.get(i,{}),p,c)
        if q:o[i]=q
        else:o.pop(i,None)
    return o
def vscale_poly(v,p):return {i:q for i,a in v.items() if (q:=pmul(a,p))}
def vrec(v):return [[i,prec(p)] for i,p in sorted(v.items())]
def rank_columns(cols):
    piv={};receipts=[]
    for j,orig in enumerate(cols):
        v=dict(orig);steps=[]
        while v:
            p=min(v)
            if p not in piv:
                c=v[p];piv[p]=scale(v,1/c)
                receipts.append({'source_column':j,'pivot':p,'subtract':[[i,str(d)] for i,d in steps],'divide_by':str(c),'normalized':srec(piv[p])});break
            c=v[p];steps.append((p,c));v=add(v,piv[p],-c)
    return [piv[i] for i in sorted(piv)],receipts
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={}
    def read(p):
        hashes[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
    def write(name,obj):
        p=out/name;p.write_text(json.dumps(obj,separators=(',',':'))+'\n');assert sum(x.stat().st_size for x in out.iterdir())<32*1024*1024;return sha(p)
    r=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/report.json');assert r['all_passed'] and r['actual_C2_dimension']==14
    d=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')
    raw=read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json')
    actions=[[parse(c) for c in A] for A in d['C2_coordinate_actions9_sparse_columns']];L=actions[8];Ar=actions[7]
    assert all(product_columns(A,B)==product_columns(B,A) for A in actions for B in actions)
    assert all(not v for v in product_columns(L,L))
    # Exact filtration: each supplied basis spans the full previous image.
    current=eye(14);dims=[14];filtration=[]
    for power in range(1,6):
        cols=[mv(A,v) for A in actions for v in current]
        nxt,receipts=rank_columns(cols)
        filtration.append({'power':power,'input_columns':[srec(v) for v in cols],
          'column_labels':[[i,j] for i in range(9) for j in range(len(current))],
          'rank':len(nxt),'fixed_rank_receipts':receipts,'basis':[srec(v) for v in nxt]})
        dims.append(len(nxt));current=nxt
        if not nxt:break
    assert dims[-1]==0;hf=[a-b for a,b in zip(dims,dims[1:])];assert sum(hf)==14 and hf[0]==3
    write('actual-C2-Loewy-filtration.json',{'coordinate_order':raw['coordinate_order9'],'dimensions_m_power_C2':dims,'actual_Loewy_HF':hf,'Loewy_length':len(hf),'stages':filtration})
    # Polynomial variables are only(r,ell), ell²=0 in this eleven-unit proof.
    one={(0,0):F(1)};R={(1,0):F(1)};E={(0,1):F(1)}
    raw_indices=d['C2_basis_raw45_indices'];basis33=raw['actual33_basis'];b6=raw['canonical_branch_quotients6'][0]['basis_plane_monomials']
    phys=[]
    for q in raw_indices:
        if q<6:s,a,l=0,*b6[q]
        elif q<12:s,a,l=1,*b6[q-6]
        else:
            t,a,l=basis33[q-12];assert t==0,'all current basis lifts must be pure r/ell monomials';s=2
        phys.append((s,a,l))
    eta=[raw_indices.index(i) for i in [0,6,12]];assert eta==[0,5,10]
    nonphys=[j for j in range(14) if j not in eta];assert len(nonphys)==11
    m0=[{s:{(a,l):F(1)}} for s,a,l in phys]
    Avec=[];Kvec=[];anchors=[]
    for j in nonphys:
        s,a,l=phys[j];ej={eta[s]:F(1)};Ara=dict(ej)
        for _ in range(a):Ara=mv(Ar,Ara)
        assert (mv(L,Ara) if l else Ara)=={j:F(1)},'literal physical basis action binding'
        a2={}
        if l:
            for i,c in Ara.items():a2[i]={(0,0):-c}
        for t in range(a):
            v=dict(ej)
            for _ in range(t):v=mv(Ar,v)
            p={(a-1-t,l):F(-1)}
            for i,c in v.items():a2=vadd(a2,{112+i:scale(p,c)})
        k=vadd({j:one},{eta[s]:{(a,l):F(1)}},-1)
        # d1 has periodic block ellI-L and r block rI-Ar.
        result={}
        for q,p in a2.items():
            if q<14:
                result=vadd(result,{q:pmul(E,p)})
                result=vadd(result,{i:scale(p,-c) for i,c in L[q].items()})
            else:
                assert 112<=q<126;i=q-112
                result=vadd(result,{i:pmul(R,p)})
                result=vadd(result,{a:scale(p,-c) for a,c in Ar[i].items()})
        assert result==k,'literal polynomial d1*a=k'
        constants={i:p.get((0,0),F(0)) for i,p in a2.items() if p.get((0,0),0)}
        assert len(constants)==1 and next(iter(constants.values()))==-1
        anchor=next(iter(constants));anchors.append(anchor);Avec.append(a2);Kvec.append(k)
    assert len(set(anchors))==11
    # Q=selected11 source rows of a. It is -I plus a nilpotent polynomial.
    Q=[[Avec[j].get(i,{}) for j in range(11)] for i in anchors]
    def mm(A,B):
        C=[[{} for _ in B[0]] for _ in A]
        for i,row in enumerate(A):
            for k,a in enumerate(row):
                if a:
                    for j,b in enumerate(B[k]):
                        if b:C[i][j]=add(C[i][j],pmul(a,b))
        return C
    I=[[one if i==j else {} for j in range(11)] for i in range(11)]
    N=[[add(Q[i][j],I[i][j]) for j in range(11)] for i in range(11)]
    term=I;total=[[dict(p) for p in row] for row in I];nil=None
    for power in range(1,12):
        term=mm(term,N)
        if not any(p for row in term for p in row):nil=power;break
        total=[[add(a,b) for a,b in zip(x,y)] for x,y in zip(total,term)]
    assert nil is not None
    inv=[[scale(p,-1) for p in row] for row in total]
    assert mm(Q,inv)==I and mm(inv,Q)==I
    # Target basis [K11,eta3], inverse=nonphysical coordinate rows followed by m0.
    T=Kvec+[{i:one} for i in eta]
    def target_inverse(v):
        outv={i:v[j] for i,j in enumerate(nonphys) if j in v}
        for j,p in v.items():outv=vadd(outv,{11+s:pmul(p,c) for s,c in m0[j].items()})
        return outv
    assert all(target_inverse(v)=={j:one} for j,v in enumerate(T))
    complement=[i for i in range(126) if i not in anchors];assert len(complement)==115
    # The literal physical first differential after cancellation, in variables
    # z0..6,r,ell. G1[:,c]=m0*d1[:,c], c in complement.
    mons9=[tuple(int(i==j) for i in range(9)) for j in range(9)];zero9=(0,)*9
    def emb(p):return {(0,)*7+(a,l):c for (a,l),c in p.items()}
    G1=[]
    for c in complement:
        if c<14:varidx=8;j=c;M=L
        else:varidx=(c-14)//14;j=(c-14)%14;M=actions[varidx]
        result={}
        for s,p in m0[j].items():result[s]=pmul(emb(p),{mons9[varidx]:F(1)})
        for i,a in M[j].items():
            for s,p in m0[i].items():
                q=add(result.get(s,{}),emb(p),-a)
                if q:result[s]=q
                else:result.pop(s,None)
        assert all(not p.get(zero9,0) for p in result.values()),'physical minimal3 target has no unit relation'
        G1.append(result)
    proof={'base':'B=(k[ell]/ell²)[[zA0,...,zA6,rA]], actual9 actions',
      'periodic_free_A_resolution':'... --(ellI+L)--> A14 --(ellI-L)--> A14 -> C2 ->0',
      'augmentation':'u+ell*w maps to u+Lw in actual k-basis14',
      'exactness':'kernelaugmentation=(ellI-L)A14; ker(ellI-L)=im(ellI+L) and conversely, directly from L²=0',
      'total_resolution':'tensor with completed Koszul8(x_iI-A_i), Ai commuting with L andeachother',
      'total_sign':'on periodic degreep/exterior degreeq, d=d_periodic+(-1)^p*d_Koszul',
      'exactness_Koszul':'completed translation by the commuting topologically nilpotent Ai reduces to regular formal variables on C2[[x]]',
      'first_free_ranks':[14,126,518,1302],'after11_unit_cancellation':[3,115,518,1302],
      'physical_target3_basis_positions':eta,'physical_lift_m0_columns':[vrec(v) for v in m0],
      'm0_eta_identity':True,'target_K11_columns':[vrec(v) for v in Kvec],
      'source_a11_columns':[vrec(v) for v in Avec],'all_d1_a_equals_K':True,
      'source_anchors11':anchors,'selected_minor_constant':-1,'selected_minor_unit_proof':'Q=-I+N, N polynomial nilpotent',
      'selected_N_nilpotence_index':nil,'selected_Q11':[[prec(p) for p in row] for row in Q],
      'selected_Q_inverse11':[[prec(p) for p in row] for row in inv],
      'source_complement115_indices':complement,
      'source_basis_clear_formula':'a11 plus (E_C-a11*projection_nonphysical*d1*E_C)',
      'target_basis':'[K11,eta3], inverse is nonphysical coordinate rows followed by m0',
      'reduced_G1_3_by115_poly_columns':[vrec(v) for v in G1],
      'G2_typed':'take115complement coordinates of transformed d2 after source basis change; first11coordinates zero by d1*d2=0',
      'G3_typed':'unchanged periodic/Koszul d3 into F2rank518',
      'actual_V_X_source_gate':'same original fullR0 comparison V3x56,X56x768 restrictedell²; not the ell0 T5 Schur',
      'still_missing':'actual source m,n,h and768cyclelifts through reduced resolution; no coefficient lift follows merely from a finite equality',
      'Vext_optional':'[eta V,K11]:B67->B14, kernel=actualkerV, cokernel=C2; targetsplit proves kernel equality',
      'coarse_T2_recipe_if_maps_supplied':{'generators':518+56,'relations':1302+115+768},
      'unsafe_finite_homology_not_used':True,'T2_or_H0_values_computed':False,
      'Hadamard_order_limit':'C2 action of a regular scalar depends on Taylor degree<Loewy_length; each coordinate division may need one higher coefficient. Formal allorder source lifts must be retained, not inferred from this finite truncation.'}
    write('actual-C2-resolution-and11-unit-source-cancellation.json',proof)
    slots=sum(len(v) for A in actions for v in A)+sum(sum(len(p) for p in v.values()) for v in Avec+Kvec+G1)
    assert slots<150000
    report={'id':STEM,'all_passed':True,'actual_C2_Loewy_dimensions':dims,'actual_C2_Loewy_HF':hf,
      'actual_C2_Loewy_length':len(hf),'exact_source_unit_pairs_cancelled':11,
      'reduced_free_resolution_first_ranks':[3,115,518,1302],'actual_T2_or_H0_computed':False,
      'new_physical_source_values':False,'input_sha256':hashes,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,'stored_named_QQ':slots}
    write('report.json',report)
    print(json.dumps({k:report[k] for k in ['all_passed','actual_C2_Loewy_HF','actual_C2_Loewy_dimensions','exact_source_unit_pairs_cancelled','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
