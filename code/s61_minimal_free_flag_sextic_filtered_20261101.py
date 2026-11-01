"""Full filtered lower tails used to certify the true sixth initial space."""
import json
import sympy as s
from s61_minimal_free_flag_filtered_data_20261101 import ROOT

def multiply(ring, left, right):
    out = [[ring.zero]*len(right[0]) for _ in left]
    sparse = [[(j, value) for j, value in enumerate(row) if value] for row in right]
    for i, row in enumerate(left):
        for k, coefficient in enumerate(row):
            if coefficient:
                for j, value in sparse[k]:
                    out[i][j] += coefficient*value
    return out

def sub(left, right):
    return [[x-y for x, y in zip(row, correction)] for row, correction in zip(left, right)]

def scale(matrix, coefficient):
    return [[coefficient*value for value in row] for row in matrix]

def times_monomial(ring, row, degree, monomial, total, monomials, indices, normal_forms, project=False, columns=None):
    out = [ring.zero]*(1237 if project else len(monomials[total]))
    columns = range(len(row)) if columns is None else columns
    for j, value in zip(columns, row):
        if value:
            column = indices[total][tuple(sorted(monomials[degree][j]+monomial))]
            if project:
                for k, coefficient in normal_forms[column]:
                    out[k] += coefficient*value
            else:
                out[column] += value
    return out

def phase_a(seed, backend, progress):
    R,a,ar,d,mm,ix,nf,C2,M,inverse,selected,N,A3,A4,A5p,free5,G35,A6 = [seed[k] for k in
        ['R','a','ar','d','mm','ix','nf','C2','M','inverse','selected','N','A3','A4','A5p','free5','G35','A6']]
    rp = lambda matrix: [[R.from_expr(x) for x in row] for row in matrix.tolist()]
    mul = lambda left,right: multiply(R,left,right)
    times = lambda row,degree,mon,total,project=False,columns=None: times_monomial(R,row,degree,mon,total,mm,ix,nf,project,columns)
    log = lambda *args: progress(' '.join(map(str,args)))
    G23=mul(rp(inverse),[A3[i]for i in selected]); G24=mul(rp(inverse),[A4[i]for i in selected])
    G25p=mul(rp(inverse),[A5p[i]for i in selected])
    G33=sub(A3,mul(N,[A3[i]for i in selected])); G34=sub(A4,mul(N,[A4[i]for i in selected]))
    G36=sub(A6,mul(N,[A6[i]for i in selected]))
    P3=s.zeros(132,286);P34=[];P35=[];P36=[]
    for i in range(12):
     for j in range(11):
      for k,mon in enumerate(mm[2]):P3[11*i+j,ix[3][tuple(sorted(mon+(j,)))]]+=C2[i,k]
      P34.append(times(G23[i],3,(j,),4));P35.append(times(G24[i],4,(j,),5));P36.append(times(G25p[i],5,(j,),6,True,free5))
    pr=list(P3.T.to_DM().rref()[1]);p3=list(P3.to_DM().rref()[1]);inv3=P3[pr,p3].inv(method='DM');B3=(inv3.to_DM()*P3[pr,:].to_DM()).to_Matrix();B34=mul(rp(inv3),[P34[i]for i in pr]);B35=mul(rp(inv3),[P35[i]for i in pr]);B36=mul(rp(inv3),[P36[i]for i in pr]);GP=[[row[j]for j in p3]for row in G33];G3red=sub(G33,mul(GP,rp(B3)));C3=s.Matrix([[x.exquo(d*d).as_expr()for x in row]for row in G3red]);cr=list(C3.T.to_DM().rref()[1]);rf,pc=C3[cr,:].to_DM().rref();rf=rf.to_Matrix();assert len(pc)==11 and not any(x.free_symbols for x in rf);invC3=C3[cr,list(pc)].inv(method='DM');log('P3 setup')
    C34raw=sub(scale(G34,d),mul(GP,B34));C34=mul(rp(invC3),[C34raw[i]for i in cr]);G44=sub(C34raw,mul(rp(C3[:,list(pc)]),C34));G44+=scale(sub(P34,mul(rp(P3[:,p3]),B34)),d*d);log('G44')
    C35raw=backend.linear(R,G35,GP,B35,d);del G35;C35=mul(rp(invC3),[C35raw[i]for i in cr]);G45=backend.linear(R,C35raw,rp(C3[:,list(pc)]),C35);del C35raw;G45+=scale(sub(P35,mul(rp(P3[:,p3]),B35)),d*d);log('G45',sum(len(x)for row in G45 for x in row))
    C36raw=sub(scale(G36,d),mul(GP,B36));C36=mul(rp(invC3),[C36raw[i]for i in cr]);G46=sub(C36raw,mul(rp(C3[:,list(pc)]),C36));G46+=scale(sub(P36,mul(rp(P3[:,p3]),B36)),d*d);log('G46',sum(len(x)for row in G46 for x in row))
    values=locals()
    return {k:values[k] for k in ['R','a','ar','d','mm','ix','nf','C2','rf','G23','G24','C34','C35','G44','G45','G46']}

def phase_b(previous, backend, progress):
    R,a,ar,d,mm,ix,nf,C2,rf,G23,G24,C34,C35,G44,G45,G46 = [previous[k] for k in
        ['R','a','ar','d','mm','ix','nf','C2','rf','G23','G24','C34','C35','G44','G45','G46']]
    rp = lambda matrix: [[R.from_expr(x) for x in row] for row in matrix.tolist()]
    mul = lambda left,right: multiply(R,left,right)
    linear = backend.linear
    times = lambda row,degree,mon,total,project=False: times_monomial(R,row,degree,mon,total,mm,ix,nf,project)
    log = lambda *args: progress(' '.join(map(str,args)))
    P4=s.zeros(913,1001);P45=[];P46=[]
    for i in range(12):
     for j,mon in enumerate(mm[2]):
      for k,right in enumerate(mm[2]):P4[66*i+j,ix[4][tuple(sorted(mon+right))]]+=C2[i,k]
      P45.append([d*d*x for x in times(G23[i],3,mon,5)]);P46.append([d*d*x for x in times(G24[i],4,mon,6,True)])
    for i in range(11):
     for j in range(11):
      for k,mon in enumerate(mm[3]):P4[792+11*i+j,ix[4][tuple(sorted(mon+(j,)))]]+=rf[i,k]
      P45.append(times(C34[i],4,(j,),5));P46.append(times(C35[i],5,(j,),6,True))
    log('P4 products')
    r4=list(P4.T.to_DM().rref()[1]);p4=list(P4.to_DM().rref()[1]);assert len(r4)==len(p4)==595;inv4=P4[r4,p4].inv(method='DM');B4=(inv4.to_DM()*P4[r4,:].to_DM()).to_Matrix();B45=mul(rp(inv4),[P45[i]for i in r4]);B46=mul(rp(inv4),[P46[i]for i in r4]);SY5=sub(P45,mul(rp(P4[:,p4]),B45));SY6=sub(P46,mul(rp(P4[:,p4]),B46));log('P4 normalize')
    g4p=[[row[j]for j in p4]for row in G44];G4red=sub(G44,mul(g4p,rp(B4)));free4=[j for j in range(1001)if j not in p4];Q4=s.Matrix([[row[j].as_expr()for j in free4]for row in G4red]);rpt4=json.load(open(ROOT/'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUARTIC_20261101.json'));sr4=rpt4['new_quartic_selected_rows'];sp4=rpt4['new_quartic_pivots'];RR4=s.Matrix(rpt4['new_quartic_rref']).applyfunc(s.sympify);assert(Q4.to_DM()-Q4[:,sp4].to_DM()*RR4.to_DM()).is_zero_matrix;CC=Q4[sr4,sp4];invCC=CC.inv(method='DM').applyfunc(s.cancel);N4=(Q4[:,sp4].to_DM()*invCC.to_DM()).to_Matrix().applyfunc(s.cancel);den=s.lcm([s.denom(x)for x in N4]);p=int(s.degree(den,a[1]));assert s.Poly(den,*a).monoms()==[(0,p,0,0,0)];NP=rp(N4.applyfunc(lambda x:s.cancel(x*a[1]**p)));log('Q4 canonical',p)
    N5=linear(R,G45,g4p,B45,d**3);log('N5',sum(len(x)for row in N5 for x in row));N6=linear(R,G46,g4p,B46,d**3);log('N6')
    E5=linear(R,N5,NP,[N5[i]for i in sr4],d**p)+scale(SY5,d**(4+p));log('E5',sum(len(x)for row in E5 for x in row));E6=linear(R,N6,NP,[N6[i]for i in sr4],d**p)+scale(SY6,d**(4+p));log('E6',sum(len(x)for row in E6 for x in row))
    q=int(s.degree(s.lcm([s.denom(x)for x in invCC]),a[1]));IC=rp(invCC.applyfunc(lambda x:s.cancel(x*a[1]**q)));C45=linear(R,[[R.zero]*3003 for _ in range(15)],[[-c for c in row]for row in IC],[N5[i]for i in sr4]);log('C45',q)
    values=locals()
    return {k:values[k] for k in ['R','a','ar','d','mm','ix','nf','C2','rf','RR4','free4','G23','C34','C45','p','q','E5','E6']}
