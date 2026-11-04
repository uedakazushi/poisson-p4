#!/usr/bin/env python3
"""Exact two-coefficient Bezout and global one-cubic-pivot finite recipe."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/S61_NONSQUARE_CUBIC_PROJECTED15_HEAVY_PREFLIGHT_20261104_v63';pr=read(base/'report.json');p=read(base/'certificate.json');assert pr['frozen']and pr['all7_heavy_coefficient_differences_extracted_exactly']
 assert sha(base/'certificate.json')==pr['certificate_sha256']
 for path,h in pr['input_sha256'].items():assert sha(ROOT/path)==h,path
 R,a,e,d=ring('a,e,d',s.QQ);field=R.to_field();dec=lambda ts:R.from_dict({tuple(m):s.QQ(c)for m,c in ts});enc=lambda q:[[list(m),str(c)]for m,c in q.items()]
 C=100*a+18*e**2-9*e-2;F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;J=5*a*d-6*e-4;delta=-64000*C**2
 vals={int(j):dec(v)for j,v in p['unscaled_G61_minus_G20_heavy_coefficients'].items()};z,w=vals[6],vals[3]
 g1=(3*e-2)*d+20;g2=(6*e+1)*d+40;h=C*g1-B;ell=6*e+1;M=600*a+126*e**2-63*e-14;K=M-ell*C
 assert 10000000*z==-a*F*B*g1*g2 and 10000000*w==-F*B*h
 assert 2*ell*h+K*g2==240*C
 assert K*C*g1*g2+(2*ell*B-K*g2)*h==240*B*C
 cZ=K*C;cW=a*(2*ell*B-K*g2);eta=R.from_dict({(0,0,0):s.QQ(192,125)})*a*F*B**2*C**3
 assert cZ*(delta*z)+cW*(delta*w)==eta
 assert vals[0]-3*vals[6]==vals[1]+vals[6]==vals[2]==vals[4]==vals[7]==R.zero
 # The actual two-dimensional target coefficient change and its inverse.
 zz,ww=field(delta*z),field(delta*w);E=field(eta)
 T=[[field(cZ)/E,field(cW)/E],[-ww,zz]];Ti=[[zz,-field(cW)/E],[ww,field(cZ)/E]]
 for left,right in[(T,Ti),(Ti,T)]:assert all(sum((left[i][k]*right[k][j]for k in range(2)),field.zero)==int(i==j)for i in range(2)for j in range(2))
 assert T[0][0]*zz+T[0][1]*ww==1 and T[1][0]*zz+T[1][1]*ww==0
 # Constant heavy-coordinate transport, independent of pure84 columns.
 H=[[0]*7 for _ in range(7)];Hi=[[0]*7 for _ in range(7)]
 H[0][0]=1;H[0][5]=-3;H[1][1]=H[1][5]=H[2][2]=H[3][4]=H[4][6]=H[5][5]=H[6][3]=1
 Hi[0][0]=1;Hi[0][5]=3;Hi[1][1]=1;Hi[1][5]=-1;Hi[2][2]=Hi[3][6]=Hi[4][3]=Hi[5][5]=Hi[6][4]=1
 for left,right in[(H,Hi),(Hi,H)]:assert all(sum(left[i][k]*right[k][j]for k in range(7))==int(i==j)for i in range(7)for j in range(7))
 deps=[ROOT/'reports'/s/'report.json'for s in['S61_NONSQUARE_PRODUCT41_FULL_OPEN_20261103_v63','S61_NONSQUARE_RAW24_PRODUCT41_BRIDGE_20261103_v63_REVIEW','S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63_REVIEW','S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63_REVIEW','S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63_REVIEW','S61_NONSQUARE_ELL_COLON_RANK10_20261104_v63_REVIEW','S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63_REVIEW']];assert all(read(q)['all_passed']for q in deps)
 mpath=ROOT/'reports/S61_NONSQUARE_ELL_CUBIC_RANK10_20261104_v63/certificate.json';minor=read(mpath)['parametric_nonzero34minor'];assert s.Rational(minor['oldpoint_value'])!=0
 encode=lambda v:[enc(v.numer),enc(v.denom)]
 proof={'coefficient_domain':'FULLU=QQ[a,e,d,1/(adFBCJ)]','actual_source_normalized_rows':[43,49],'actual135_rows':[20,61],'remaining57_source_rows':[42,47],
 'g1':enc(g1),'g2':enc(g2),'h':enc(h),'Lscalar':enc(ell),'Kscalar':enc(K),
 'exact_heavy_factors':{'kZ=-aFBg1g2/10000000':enc(z),'kW=-FBh/10000000':enc(w)},
 'Bezout_cZ':enc(cZ),'Bezout_cW':enc(cW),'known_unit_eta':enc(eta),
 'literal_pair_identity':'cZ*delta*kZ+cW*delta*kW=eta=192*a*F*B^2*C^3/125',
 'pair_target_coefficient_change2x2':[[encode(v)for v in row]for row in T],'pair_target_inverse2x2':[[encode(v)for v in row]for row in Ti],
 'heavy7_physical_order':[0,1,2,3,4,6,7],'heavy7_oldX_order':['x0','x1','x2','x3','x4','Z','W'],'constant_heavy_coefficient_change7x7':H,'constant_heavy_inverse7x7':Hi,
 'coefficient_matrix_orientation':'heavy coefficientCOLUMN new=T2*old; rowmatrix usesRIGHT T2transpose. Pure84targetcolumns unchanged.',
 'global_cubic_unit_pivot_recipe':{'source_change':'Replace remaining57 row47 by Cdelta=delta*(row47-row42); other56source rows unchanged. This is invertible withdetdelta.',
 'target_change':'Apply the constant7x7heavychange to the7H-linear columns, thenT2 toHZ/HW coefficientcolumns. Keepall84purecolumns unchanged.',
 'pivot':'Cdelta has heavycoordinates (0,0,0,0,0,1,0), plusitsactualpure84coefficients. Useunit1pivot atnewHnu coordinate.',
 'same_cubic_cokernel_shape':[56,90],'no_expanded56x90_matrix':True,'all41_quadproduct_source_rows_have_zero_newHnu_and_are_unchanged':True},
 'projected15_filtration':{'all_FULLU_exact_sequence':'0->T29->B50->V21->0, usingactualproduct41andactualZWproduct19; noCI/U78needed',
 'Delta34_definition':minor,'Delta34_source_circuit':'reports/S61_NONSQUARE_ELL_COLON_CUBIC11_20261104_v63/circuit.json',
 'on_FULLU_Delta34':'projectedgenuine16rank15withprimitivekernel G61-G20; etaunitmakesitsfullimage nonzero, soall16genuineindependentmodquad41 andACTUALfullcubicrank57',
 'outside_Delta34':'OPEN; full57mayholdtherebutnotprovedbythispacket','Delta34_not_known_global_unit':True},
 'rank57_on_extra_Delta34_open_proved_primary_review_pending':True,'global_actual_cubic56_to90_samecokernel_recipe_proved_primary_review_pending':True,
 'no_allFULLU57_quartic239_genericbeta_or_newkernelcount_claim':True,'no_point_GB_rankjob_fieldsolver_quarticexpansion_or_source26':True}
 (out/'certificate.json').write_text(json.dumps(proof,indent=2)+'\n')
 note=ROOT/'research/notes'/(ID+'.md');paths=[Path(__file__),note,base/'report.json',base/'certificate.json',mpath,*deps]
 result={'id':ID,'all_exact_QQ_heavy_factor_and_Bezout_identities_passed':True,'global_pair_unit_and_target2x2_inverse_verified':True,'actual_global_cubic_unitpivot_samecokernel_shape':[56,90],'actual_full_cubic_rank57_on_FULLU_Delta34':True,
 'Delta34_added_open_not_global_unit':True,'all_FULLU57_or_quartic239_not_proved':True,'knownkernel38_missing42_unchanged':True,'no_newpoint_rank_GB_solver_or_quarticexpansion':True,
 'input_sha256':{str(q.relative_to(ROOT)):sha(q)for q in paths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
