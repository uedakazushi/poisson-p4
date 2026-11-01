#!/usr/bin/env python3
"""Exact initial ideal and square-zero homogenizer at a binary/mixed junction.

The formal-ring inference is at the specified fixed-beta point only.
This runner uses exact QQ/ZZ algebra and one explicitly labelled modular
non-membership minor; no ignored scratch data are mathematical inputs.
"""
import argparse,hashlib,json,os,subprocess,tempfile,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]


def build():
    m=S61();main=json.loads((ROOT/'reports/S61_BINARY_MIXED_CONNECTION_20261101.json').read_text());b={int(k[1:])-1:s.Rational(v) for k,v in main['point']['sample'].items()};two,F=specialize(m,b,True);fixed=two.index(29);cols=[i for i in range(45) if i!=fixed];J=s.zeros(135,45);terms=[]
    for row,f in enumerate(F):
     ht=[]
     for mon,c in f.items():
      if len(mon)==2:
       a,bb=mon
       if a==fixed and bb!=fixed:J[row,cols.index(bb)]+=c
       if bb==fixed and a!=fixed:J[row,cols.index(a)]+=c
      elif len(mon)==1 and mon[0]==fixed:J[row,44]+=c
      if fixed not in mon:
       a,bb=[cols.index(i) for i in mon]+[44]*(2-len(mon));ht.append((a,bb,c))
     terms.append(ht)
    T=s.Matrix.hstack(*J.nullspace());C=main['point']['jacobian_pivot_columns'];rows=main['point']['jacobian_pivot_rows'];U=s.zeros(45,33);inv=J.extract(rows,C).inv()
    for i,col in enumerate(C):U[col,:]=inv[i,:]
    factors=[s.ilcm(*[s.denom(v) for v in list(J[i,:])+[c for _,_,c in terms[i]]]) for i in range(135)]
    Jn=s.diag(*factors)*J;Un=U*s.diag(*[s.Rational(1,factors[i]) for i in rows]);den=s.ilcm(*[s.denom(v) for v in Un]);Un=(den*Un).applyfunc(s.cancel);td=s.ilcm(*[s.denom(v) for v in T]);Tn=td*T
    assert all(v.q==1 for v in Jn) and all(v.q==1 for v in Un) and all(v.q==1 for v in Tn)
    fmt=lambda M:'matrix{'+','.join('{'+','.join(str(c) for c in row)+'}' for row in M.tolist())+'}'
    text='R=ZZ[v_0..v_11]; T='+fmt(Tn)+'; J='+fmt(Jn)+'; U='+fmt(Un)+'; jetDen='+str(den)+';rows={'+','.join(map(str,rows))+'};\n'
    expr=[];diag=[]
    for row,ht in enumerate(terms):
     expr.append('{'+('+'.join('('+str(c*factors[row])+')*(V_('+str(i)+',0)*W_('+str(j)+',0)+V_('+str(j)+',0)*W_('+str(i)+',0))' for i,j,c in ht) or '0_R')+'}')
     diag.append('{'+('+'.join('('+str(c*factors[row])+')*V_('+str(i)+',0)*V_('+str(j)+',0)' for i,j,c in ht) or '0_R')+'}')
    text+='jetB=(V,W)->matrix{'+','.join(expr)+'};jetQ=V->matrix{'+','.join(diag)+'};\n'
    text+='L1=T*transpose vars R; Q2=jetQ(L1); W2=-U*Q2^rows; A2=jetDen*Q2+J*W2;print("integer twojet");stdio<<flush; Q3=jetB(L1,W2); W3=-U*Q3^rows; A3=jetDen*Q3+J*W3;print("integer threejet");stdio<<flush; Q4=jetB(L1,W3)+jetQ(W2); W4=-U*Q4^rows; A4=jetDen*Q4+J*W4;print("integer fourjet");stdio<<flush;\n'
    text+='S=QQ[v_0..v_11];tjet=(1/jetDen)*sub(W2_(44,0),S)+(1/jetDen^2)*sub(W3_(44,0),S);A2=(1/jetDen)*sub(A2,S);A3=(1/jetDen^2)*sub(A3,S);A4=(1/jetDen^3)*sub(A4,S);R=S;use R;\n'
    assert T.shape==(45,12) and J*T==s.zeros(135,12) and T[-1,:]==s.zeros(1,12)
    assert td % 32003 != 0 and den % 32003 != 0
    free=[19,20,22,23,24,27,32,33,34,36,37,43];assert T[free,:]==s.eye(12)
    # An actual closed double line, not just a tangent-space observation.
    delta=s.zeros(45,1);nil=s.zeros(45,1)
    for i,v in {19:2645,23:4500,24:1587,26:-2484,27:4140,31:972,32:-3240,33:2700}.items():delta[i]=v
    for i,v in {0:1485,19:-1108,24:990,27:378,34:405}.items():nil[i]=v
    assert J*delta==s.zeros(135,1) and J*nil==s.zeros(135,1)
    assert s.Matrix.hstack(delta,nil).rank()==2
    for ht in terms:
     assert sum(c*delta[i]*delta[j] for i,j,c in ht)==0
     assert sum(c*(delta[i]*nil[j]+nil[i]*delta[j]) for i,j,c in ht)==0
    for degree in (2,3,4):
        monomials=['*'.join('v_'+str(i) for i in mon) for mon in combinations_with_replacement(range(12),degree)]
        text+='m'+str(degree)+'=matrix{{'+','.join(monomials)+'}};\n'
    text+=SCHUR_PROGRAM+'\nJ=J4;\n'+PRIMARY_PROGRAM
    return text,dict(sample=main['point']['sample'],tangent_basis=[[str(v) for v in row] for row in T.tolist()],
        chart_free_rows=free,integer_tangent_scale=str(td),implicit_denominator=str(den),
        actual_double_line_tangent=list(map(str,delta)),actual_double_line_nil_direction=list(map(str,nil)))

SCHUR_PROGRAM = """Q=sub(transpose last coefficients(transpose A2,Monomials=>m2),QQ);C=sub(transpose last coefficients(transpose A3,Monomials=>m3),QQ);
Lpol=A2+A3;Hpol=A4;scan(12,i->(Lpol=Lpol||(v_i*A2);Hpol=Hpol||(v_i*A3)));Lo=sub(transpose last coefficients(transpose Lpol,Monomials=>(m2|m3)),QQ);Hi=sub(transpose last coefficients(transpose Hpol,Monomials=>m4),QQ);assert(rank Q==24);assert(rank Lo==226);tc=sub(transpose last coefficients(matrix{{tjet}},Monomials=>(m2|m3)),QQ);assert(rank(Q||tc_{0..77})==24);assert(rank sub(Lo||tc,ZZ/32003)==227);print("t nonzero in degree three checked");print("low matrices ready");stdio<<flush;
C0=transpose(syz transpose Q)*C;rp={84,86,19,117,60,133,17,114,20,22,23,102,61,40,91,85,89,21,92,67,99,43,25,129,219,221,154,252,195,268,63,249,68,101,115,624,626,116,120,759,122,132,894,1029,66,134,121,155,157,77,158,175,59,15,559,98,113,57,289,424,12,694,829,1099,657,1234,1504,78,178,330,88,465,196,600,673,112,654,560,735,870,16,1140,1275,14,1545,1680,354,95,489,39,96,1164,1299,97,55,100,94,33,108,1569,220,109,356,1704,111,384,227,107,34,118,387,522,761,927,37,1197,234,119,403,152,563,896,519,126,924,789,127,1194,538,130,943,1599,41,1213,287,422,557,45,692,827,1734,26,1329,58,1602,962,290,425,642,695,580,202,625,830,562,631,1097,601,607,632,583,965,639,1332,1348,1100,1235,1505,1640,156,561,48,1232,292,1367,697,832,1502,1102,1237,1031,293,1005,736,871,47,760,428,38,698,1637,833,1103,1238,1006,895,52,331,44,1370,226,696,766,831,767,28,30,291,54,629,427,964,160,310,426,445,715,850,966,1120,1255,1101,1236,1371,313,448};cp={6,9,12,13,15,16,17,19,20,27,30,36,44,47,51,54,57,60,61,62,65,69,72,74,84,87,90,91,93,94,95,97,98,105,108,114,117,122,125,129,132,135,136,137,138,139,140,143,147,150,151,152,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,178,179,180,182,183,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,211,212,213,214,215,216,217,218,221,222,224,225,226,227,228,229,231,235,238,241,242,243,244,245,246,248,249,250,251,252,253,255,256,257,258,259,260,261,263,264,267,268,270,271,272,273,280,288,291,295,298,301,302,303,304,305,306,309,313,316,318,322,323,324,327,331,334,337,338,339,340,341,342,345,349,352,353,354,359,362,365,366,367,368,369,370,373,377,380,381,382,386,387,388,389,390,391,394,395,396,397,398,399,400,401,402,403,404,405,406,409,413,416,417,418,423,426,427,428,432,433,434,436,437};rq={1369,1508,1643,1525,1660,718,853,1123,1410,1141,742,877,1434,355,490,1030,237,264,669};aa=(Lo^rp)_cp;schur=((Lo^rq)_cp)*(inverse aa);assert((Lo^rq)-schur*(Lo^rp)==0);H0=Hi^rq-schur*(Hi^rp);print("19 exact quartics ready");stdio<<flush;I=ideal(sub(Q,R)*transpose m2)+ideal(sub(C0,R)*transpose m3);J4=I+ideal(sub(H0,R)*transpose m4);G=gens gb J4;assert(dim(R/J4)==3 and degree(R/J4)==65);print("exact degree",degree(R/J4));print("exact Hilbert",hilbertSeries(R/J4));
"""

PRIMARY_PROGRAM = """B=saturate(J,ideal(v_1));D=saturate(J,ideal(v_7));
assert(degree(R/B)==65 and degree(R/D)==2 and intersect(B,D)==J);print("exact decomposition checked");stdio<<flush;
del={2645,0,0,4500,1587,4140,-3240,2700,0,0,0,0};nn={-1108,0,0,0,990,378,0,0,405,0,0,0};
E=ideal(v_8^2)+ideal apply(select(toList(0..11),i->i!=7 and i!=8),i->v_i-(del#i)*v_7/2700-(nn#i)*v_8/405);
assert(D==E);print("exact double line checked");stdio<<flush;
C=ideal(v_4,v_6,v_7,v_8,v_9,v_11,v_2-(117/110)*v_1,v_5+(3122/99)*v_1,v_10-v_1);
assert((gens B)%C==0);powers=apply(flatten entries gens C,g->first select(toList(1..8),d->g^d%B==0));print("C powers",powers);stdio<<flush;
S=QQ[w_2,w_4,w_5,w_6,w_7,w_8,w_9,w_10,w_11,w_0,w_3,MonomialOrder=>{GRevLex=>9,GRevLex=>2}];ff=map(S,R,{w_0,1_S,w_2,w_3,w_4,w_5,w_6,w_7,w_8,w_9,w_10,w_11}); B1=ff B; GB=gens gb B1;
cc={};scan(flatten entries GB,g->(ee:=take((exponents leadTerm g)#0,9);mm:=product(9,j->S_j^(ee#j));cf:=(sum(select(terms g,tt->take((exponents tt)#0,9)==ee)))//mm;if first degree cf>0 then cc=append(cc,cf/(leadCoefficient cf))));
cc=unique cc;print("QQ leading coefficients",cc);stdio<<flush;den=product cc;assert(saturate(B1,ideal den)==B1);normalInitial=ideal apply(flatten entries GB,g->product(9,j->S_j^((exponents leadTerm g)#0#j)));assert(dim(S/normalInitial)==2 and degree(S/normalInitial)==65);print("QQ C-primary certificate PASS");print("CERTIFICATE PASS");stdio<<flush;
"""

def verify():
    start=time.time();program,data=build()
    binary=Path('/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary')
    assert binary.exists(), 'Macaulay2 runtime required; see the repository setup notes.'
    env=os.environ.copy()
    env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
    with tempfile.TemporaryDirectory(prefix='c3-formal-',dir=ROOT/'work') as scratch:
        source=Path(scratch)/'certificate.m2';source.write_text(program)
        run=subprocess.run([str(binary),'--script',str(source)],cwd=ROOT,env=env,
            capture_output=True,text=True,timeout=180)
    output=run.stdout+run.stderr
    assert run.returncode==0 and 'CERTIFICATE PASS' in output and 'error:' not in output, output[-4000:]
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
        ROOT/'reports/S61_BINARY_MIXED_CONNECTION_20261101.json',
        ROOT/'reports/S61_BINARY_UNIFORM_QUADRATIC_20261101.json',
        ROOT/'reports/S61_BINARY_UNIFORM_CUBIC_20261101.json',
        ROOT/'research/notes/S61_BINARY_MIXED_FORMAL_20261101.md']
    return dict(id='S61_BINARY_MIXED_FORMAL_20261101',all_passed=True,field='QQ, with stated GF(32003) rank lower bound',
        point=data,initial_ideal_dimension=3,initial_ideal_degree=65,
        initial_hilbert_numerator_over_one_minus_z_cubed=[1,9,21,35,-2,1],
        initial_primary_component_dimensions=[3,1],initial_primary_multiplicities=[65,2],
        exact_lower_Macaulay_rank=226,homogenizer_augmented_rank_mod_32003=227,
        actual_double_line_verified=True,completed_homogenizer_nonzero=True,
        completed_homogenizer_square_zero=True,associated_graded_ideal_identified=True,
        source_program_sha256=hashlib.sha256(program.encode()).hexdigest(),macaulay2_program=program,
        macaulay2_output=output,proof_note='research/notes/S61_BINARY_MIXED_FORMAL_20261101.md',
        limitations=['Specified fixed-beta completed local ring only.',
            'No relative family or all-general-beta removal is inferred.',
            'The full completed ring is not asserted to be free rank65 over the binary base.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify()
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps(dict(all_passed=True,initial_degree=65,primary_lengths=[65,2],
        homogenizer_nonzero_square_zero=True,seconds=report['seconds'])))
