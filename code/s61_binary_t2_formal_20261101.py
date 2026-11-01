#!/usr/bin/env python3
"""A binary/T2 junction with one cubic initial socle and square-zero homogenizer.

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
    m=S61();main=json.loads((ROOT/'reports/S61_BINARY_T2_CONNECTION_20261101.json').read_text());b={int(k[1:])-1:s.Rational(v) for k,v in main['point']['sample'].items()};two,F=specialize(m,b,True);fixed=two.index(29);cols=[i for i in range(45) if i!=fixed];J=s.zeros(135,45);terms=[]
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
    for i,v in {19:529,23:675,27:46,32:-2,33:30}.items():delta[i]=v
    for i,v in {0:5,19:-s.Rational(66,5),24:-s.Rational(16,5),27:s.Rational(16,5),34:1}.items():nil[i]=v
    assert J*delta==s.zeros(135,1) and J*nil==s.zeros(135,1)
    assert s.Matrix.hstack(delta,nil).rank()==2
    for ht in terms:
     assert sum(c*delta[i]*delta[j] for i,j,c in ht)==0
     assert sum(c*(delta[i]*nil[j]+nil[i]*delta[j]) for i,j,c in ht)==0
    for degree in (2,3,4):
        monomials=['*'.join('v_'+str(i) for i in mon) for mon in combinations_with_replacement(range(12),degree)]
        text+='m'+str(degree)+'=matrix{{'+','.join(monomials)+'}};\n'
    text+=SCHUR_PROGRAM+'\n'+PRIMARY_PROGRAM
    return text,dict(sample=main['point']['sample'],tangent_basis=[[str(v) for v in row] for row in T.tolist()],
        chart_free_rows=free,integer_tangent_scale=str(td),implicit_denominator=str(den),
        actual_double_line_tangent=list(map(str,delta)),actual_double_line_nil_direction=list(map(str,nil)))

SCHUR_PROGRAM = 'Q=sub(transpose last coefficients(transpose A2,Monomials=>m2),QQ);C=sub(transpose last coefficients(transpose A3,Monomials=>m3),QQ);\nLpol=A2+A3;Hpol=A4;scan(12,i->(Lpol=Lpol||(v_i*A2);Hpol=Hpol||(v_i*A3)));Lo=sub(transpose last coefficients(transpose Lpol,Monomials=>(m2|m3)),QQ);Hi=sub(transpose last coefficients(transpose Hpol,Monomials=>m4),QQ);assert(rank Q==20);assert(rank Lo==214);tc=sub(transpose last coefficients(matrix{{tjet}},Monomials=>(m2|m3)),QQ);assert(rank(Q||tc_{0..77})==20);assert(rank sub(Lo||tc,ZZ/32003)==215);C0=transpose(syz transpose Q)*C;\nrp={84,19,103,60,117,39,133,114,20,102,23,134,40,85,96,68,99,43,21,129,219,154,238,195,252,61,268,249,86,115,116,624,120,122,237,132,91,269,89,155,158,174,63,175,203,59,15,559,113,289,55,424,694,58,829,1099,1234,643,1504,1639,42,178,330,67,465,354,600,657,579,673,654,560,735,870,77,1140,1275,78,1545,1680,489,28,894,1164,88,1569,373,38,95,14,94,98,101,508,220,109,913,1183,759,112,1588,234,106,48,118,384,387,1704,522,119,927,778,121,1197,642,563,1299,1318,126,519,789,128,924,403,131,538,1194,35,943,130,964,1369,127,290,425,1213,1599,695,1329,830,1734,47,1100,1235,1602,580,674,625,1505,264,1640,608,639,583,156,669,1332,1348,561,12,293,428,698,66,1005,52,404,833,1103,1238,54,1508,41,1029,309,231,372,777,124,444,46,714,45,849,1119,912,1047,1182,56,22,984,11,17,965,878,51,1370,291,310,1254,445,715,426,850,1120,1255,696,313,448,831,1101,1236,718,338};cp={9,12,13,15,16,17,18,19,20,27,30,44,47,54,60,62,65,69,72,74,87,90,91,93,94,95,96,97,98,105,108,117,122,125,129,132,135,136,138,140,143,147,150,151,152,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,178,179,180,181,182,183,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,221,222,224,225,226,227,228,229,231,235,238,241,242,243,244,245,246,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,267,268,270,271,272,273,288,291,295,298,301,303,304,306,309,313,316,318,322,323,324,327,331,334,337,338,339,340,341,342,345,349,352,353,354,359,362,365,366,367,368,370,373,377,380,381,382,386,387,388,389,390,391,394,396,398,400,401,402,403,405,406,409,413,416,417,418,423,426,427,428,432,433,434,436,437};rq={1643,1524,1659,1525,1660,853,1123,1410};\n\naa=(Lo^rp)_cp;schur=((Lo^rq)_cp)*(inverse aa);assert((Lo^rq)-schur*(Lo^rp)==0);H0=Hi^rq-schur*(Hi^rp);print("eight exact quartics");stdio<<flush;\nI=ideal(sub(Q,R)*transpose m2)+ideal(sub(C0,R)*transpose m3);J=I+ideal(sub(H0,R)*transpose m4);assert(dim(R/J)==3 and degree(R/J)==65);print("exact Hilbert",hilbertSeries(R/J));stdio<<flush;\n'

PRIMARY_PROGRAM = 'B=saturate(J,ideal(v_1));C=ideal(v_4,v_6,v_7,v_8,v_9,v_11,v_2+(4/3)*v_1,v_5-(76/3)*v_1,v_10-v_1);\ndel={529,0,0,675,0,46,-2,30,0,0,0,0};nn={-66/5,0,0,0,-16/5,16/5,0,0,1,0,0,0};\nD=ideal(v_8^2)+ideal apply(select(toList(0..11),i->i!=7 and i!=8),i->R_i-del#i*v_7/30-nn#i*v_8);\nassert(isSubset(J,D));assert(dim(R/B)==3 and degree(R/B)==65);E=intersect(B,D)/J;assert(dim E==0 and degree E==1 and hilbertFunction(3,E)==1);print("extra cubic socle dimension1");stdio<<flush;\nassert(isSubset(B,C));powers=apply(flatten entries gens C,g->first select(toList(1..8),d->g^d%B==0));print("C powers",powers);stdio<<flush;\nS=QQ[w_2,w_4,w_5,w_6,w_7,w_8,w_9,w_10,w_11,w_0,w_3,MonomialOrder=>{GRevLex=>9,GRevLex=>2}];ff=map(S,R,{w_0,1_S,w_2,w_3,w_4,w_5,w_6,w_7,w_8,w_9,w_10,w_11});B1=ff B;GB=gens gb B1;\ncc={};scan(flatten entries GB,g->(ee:=take((exponents leadTerm g)#0,9);mm:=product(9,j->S_j^(ee#j));cf:=(sum(select(terms g,tt->take((exponents tt)#0,9)==ee)))//mm;if first degree cf>0 then cc=append(cc,cf/(leadCoefficient cf))));cc=unique cc;print("QQ leading coefficients",cc);den=product cc;assert(saturate(B1,ideal den)==B1);normalInitial=ideal apply(flatten entries GB,g->product(9,j->S_j^((exponents leadTerm g)#0#j)));assert(dim(S/normalInitial)==2 and degree(S/normalInitial)==65);print("QQ C-primary PASS");print("CERTIFICATE PASS");stdio<<flush;\n'

def verify():
    start=time.time();program,data=build()
    binary=Path('/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary')
    env=os.environ.copy();env['DYLD_LIBRARY_PATH']=':'.join([str(ROOT/'work/m2-runtime/lib'),'/opt/homebrew/lib','/opt/homebrew/opt/givaro/lib'])
    with tempfile.TemporaryDirectory(prefix='binary-t2-formal-',dir=ROOT/'work') as temp:
        source=Path(temp)/'certificate.m2';source.write_text(program)
        run=subprocess.run([str(binary),'--script',str(source)],cwd=ROOT,env=env,capture_output=True,text=True,timeout=180)
    output=run.stdout+run.stderr
    assert run.returncode==0 and 'CERTIFICATE PASS' in output and 'error:' not in output,output[-4000:]
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
        ROOT/'reports/S61_BINARY_T2_CONNECTION_20261101.json',
        ROOT/'reports/S61_BINARY_UNIFORM_QUADRATIC_20261101.json',
        ROOT/'reports/S61_BINARY_UNIFORM_CUBIC_20261101.json',
        ROOT/'research/notes/S61_BINARY_T2_FORMAL_20261101.md']
    return dict(id='S61_BINARY_T2_FORMAL_20261101',all_passed=True,field='QQ with an explicitly labelled GF(32003) rank lower bound',
        point=data,initial_ideal_dimension=3,initial_ideal_degree=65,
        initial_hilbert_numerator_over_one_minus_z_cubed=[1,9,25,31,-4,4,-1],
        binary_primary_degree=65,actual_double_line_degree=2,
        difference_between_initial_ideal_and_primary_intersection_length=1,
        difference_module_concentrated_degree=3,
        exact_lower_Macaulay_rank=214,homogenizer_augmented_rank_mod_32003=215,
        actual_binary_kernel_torsion_length_upper_bound=1,
        actual_binary_kernel_torsion_filtration_degree_if_nonzero=3,
        actual_double_line_verified=True,associated_graded_ideal_identified=True,
        completed_homogenizer_nonzero=True,completed_homogenizer_square_zero=True,
        source_program_sha256=hashlib.sha256(program.encode()).hexdigest(),macaulay2_program=program,macaulay2_output=output,
        proof_note='research/notes/S61_BINARY_T2_FORMAL_20261101.md',
        limitations=['Specified fixed-beta completed local ring only; relative removal needs a separate argument.',
            'The initial cubic socle does not by itself imply nonzero torsion in the actual local ring.',
            'No free rank65 claim for the full ring or its double-line kernel.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify()
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps(dict(all_passed=True,initial_degree=65,cubic_socle_length=1,seconds=report['seconds'])))
