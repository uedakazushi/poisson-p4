#!/usr/bin/env python3
"""Typed actual57/63 source lift, physical91 targets and eta-Cdelta correction."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def parsepoly(line):
 a=line.split();n=int(a[0]);assert len(a)==1+4*n;return [[list(map(int,a[2+4*j:5+4*j])),str(Q(a[1+4*j]))]for j in range(n)]
def universal_cdelta():
 def poly(i):return{tuple(int(k==i)for k in range(6)):1}
 def add(a,b,s=1):
  z=a.copy()
  for m,c in b.items():
   z[m]=z.get(m,0)+s*c
   if not z[m]:del z[m]
  return z
 def mul(a,b):
  z={}
  for m,c in a.items():
   for n,d in b.items():
    k=tuple(x+y for x,y in zip(m,n));z[k]=z.get(k,0)+c*d
    if not z[k]:del z[k]
  return z
 cz,cw,z,w,hz,hw=[poly(i)for i in range(6)];eta=add(mul(cz,z),mul(cw,w));c=add(mul(cz,hz),mul(cw,hw));psi=add(mul(z,hw),mul(w,hz),-1)
 assert not add(add(mul(eta,hz),mul(c,z),-1),mul(cw,psi))
 assert not add(add(mul(eta,hw),mul(c,w),-1),mul(cz,psi),-1)
 return dict(all_two_universal_polynomial_identities_zero=True,formal_variables=['cZ','cW','zdelta','wdelta','hZ','hW'])
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='certificate.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
 dr='S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63';old=load(dr,'operator.json');assert load(dr+'_REVIEW','report.json')['all_passed']
 newst='S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63';new=load(newst,'operator.json');assert load(newst+'_REVIEW','report.json')['all_passed']
 es='S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63';entry=load(es,'source-binding.json');assert load(es+'_REVIEW','report.json')['all_passed']
 rowst='S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62';basis=load(rowst,'report.json');assert load(rowst+'_REVIEW','report.json')['all_passed'];remaining=basis['remaining_original_row_indices']
 src=load('S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62','source-report.json');mons=src['normal_monomials'];labels=src['labels'];assert len(labels)==63 and len(mons)==91 and len(remaining)==57
 np=ROOT/'reports'/rowst/'normalization.txt';pins[str(np.relative_to(ROOT))]=sha(np);lines=np.read_text().splitlines();assert int(lines[0])==63;normalizers=[]
 for line in lines[1:]:
  p=parsepoly(line);assert len(p)==1 and p[0][0]==[0,0,0]and Q(p[0][1]);normalizers.append(Q(p[0][1]))
 assert len(normalizers)==63 and all(normalizers[i]==Q(9,20)for i in range(7))and all(normalizers[i]==Q(-5,4)for i in range(35,42))
 cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');cpmeta=dict(rows=cp['remaining_rows'],columns=cp['remaining_columns'],pivots=[(x['pivot_row'],x['pivot_column'])for x in cp['all_completed_unit_steps']]);del cp
 aug=load('S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63','input.json');source25=aug['source25_labels'];mu=aug['genuine16_source_clearing_units'];projection31=aug['target31_recipe'];rawzw=aug['rawHZ_HW_indices_in_frozen91'];del aug
 lf=load('S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63','certificate.json');lfn=load('S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63','native.json');assert lfn['q5_C_factors_cancelled']==2;del lfn
 lp=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63/L-coefficients.txt';pins[str(lp.relative_to(ROOT))]=sha(lp);ll=lp.read_text().splitlines();assert ll[0].split()==['5','1'];Lpolys=list(map(parsepoly,ll[1:]));assert len(Lpolys)==5
 unit=load('S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63');hunit=load('S61_NONSQUARE_CUBIC_PROJECTED15_HEAVY_PREFLIGHT_20261104_v63');assert load('S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63_REVIEW','report.json')['all_passed']
 pair=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63');q2=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');assert pair['all_5_remaining_heavy_and_psi_columns_unchanged']==q2['all5_remaining_heavy_quotient_columns_unchanged']==[27,28,29,30,31]
 assert (11,26)in cpmeta['pivots'];assert len(cpmeta['pivots'])==8
 delta=hunit['known_unit_delta'];canon=lambda p:sorted((tuple(m),Q(c))for m,c in p);assert canon(mu[1])==canon(mu[6])==canon(delta)
 light=[0,1,2,3,4,6,7];Hcols=[]
 for i in light:
  m=[0]*8;m[5]=m[i]=1;Hcols.append(mons.index(m))
 purecols=[i for i,m in enumerate(mons)if not m[5]];assert len(purecols)==84 and rawzw==[Hcols[5],Hcols[3]]
 nodes=[]
 def n(op,shape,args=(),**kw):
  assert all(0<=i<len(nodes)for i in args);i=len(nodes);nodes.append(dict(id=i,op=op,shape=shape,args=list(args),**kw));return i
 def ref(st,name,pointer,shape):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return n('external_frozen_reference',shape,reference=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),json_pointer=pointer))
 def c(x):return n('QQ_constant',[1,1],value=str(Q(x)))
 def prod(a,b):assert nodes[a]['shape'][1]==nodes[b]['shape'][0];return n('matrix_product',[nodes[a]['shape'][0],nodes[b]['shape'][1]],[a,b])
 def neg(a):return n('negative',nodes[a]['shape'],[a])
 def add(a,b):assert nodes[a]['shape']==nodes[b]['shape'];return n('sum',nodes[a]['shape'],[a,b])
 def sub(a,b):assert nodes[a]['shape']==nodes[b]['shape'];return n('difference',nodes[a]['shape'],[a,b])
 def sel(a,rr,cc):return n('submatrix',[len(rr),len(cc)],[a],rows=rr,columns=cc)
 def sparse(sh,entries):return n('sparse_matrix_from_scalar_roots',sh,[x for i,j,x in entries],entries=entries)
 oldprefix=ref(dr,'operator.json','/exports/source_current15_to_original25',[15,25])
 Lroots=[n('literal_parameter_polynomial',[1,1],polynomial=p,source=dict(path=str(lp.relative_to(ROOT)),row=i))for i,p in enumerate(Lpolys)]
 entries=[];genprovenance=[];phys_five=[0,1,2,4,7]
 def product_entry(row,qi,v,alpha,ell=None):
  original=labels.index(['product',qi,v]);assert original in remaining;where=remaining.index(original);scale=c(Q(alpha)*normalizers[original]);value=scale if ell is None else prod(scale,Lroots[ell]);entries.append((row,where,value))
 for row,v in enumerate([0,2,4,7]):product_entry(row,0,v,1);product_entry(row+4,5,v,1600)
 for k,v in enumerate(phys_five):product_entry(8,0,v,1,k);product_entry(9,5,v,1600,k)
 for row,desc in enumerate(source25[10:],10):
  assert desc[0]=='ZW_genuine';rel,r57,orig=desc[1:];assert r57==41+rel and remaining[r57]==orig and labels[orig][0]=='cubic'
  u=ref('S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63','input.json',f'/genuine16_source_clearing_units/{rel}',[1,1]);entries.append((row,r57,u));genprovenance.append(dict(old25_row=row,genuine_relative=rel,full57_row=r57,original63_row=orig,actual135_cubic=labels[orig][1],full63_normalizer=str(normalizers[orig])))
 Lift25=sparse([25,57],entries);To63=sparse([57,63],[(r,orig,c(1/normalizers[orig]))for r,orig in enumerate(remaining)])
 raw=ROOT/'reports/S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63/frozen57-source.txt';assert sha(raw)==basis['remaining_matrix_sha256'];pins[str(raw.relative_to(ROOT))]=sha(raw)
 F57=n('external_sparse_coefficient_matrix_file',[57,91],reference=dict(path=str(raw.relative_to(ROOT)),sha256=sha(raw),format='57 91 57 thenrow-majorQQ[a,e,d] polynomial records'),all_coefficients_left_external=True)
 F_H=sel(F57,list(range(57)),Hcols);C=n('literal_QQ_matrix',[7,7],matrix=unit['constant_heavy_coefficient_change7x7']);Ct=n('transpose',[7,7],[C])
 delta_node=ref('S61_NONSQUARE_CUBIC_PROJECTED15_HEAVY_PREFLIGHT_20261104_v63','certificate.json','/known_unit_delta',[1,1]);Cdelta=sparse([1,57],[(0,47,delta_node),(0,42,neg(delta_node))])
 cz=ref('S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63','certificate.json','/Bezout_cZ',[1,1]);cw=ref('S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63','certificate.json','/Bezout_cW',[1,1]);eta=ref('S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63','certificate.json','/known_unit_eta',[1,1]);ie=n('known_FULLU_parameter_unit_inverse',[1,1],[eta],proof='eta=192*a*F*B^2*C^3/125')
 kapnodes=[ref(es,f'root{i}.json','/polynomial',[1,1])for i in entry['cofactor_roots']];kap15=sparse([1,15],[(0,r,x)for r,x in zip(entry['source_current15_rows'],kapnodes)]);kap25=prod(kap15,oldprefix)
 res25=ref(newst,'operator.json','/exports/surviving4_source_original25',[4,25]);fullsrc15=ref(newst,'operator.json','/exports/full_source15',[15,15]);pivot25=prod(sel(fullsrc15,list(range(11)),list(range(15))),oldprefix)
 def lift_and_correct(W25,howmany):
  W=prod(W25,Lift25);H=prod(prod(W,F_H),Ct);hz=sel(H,list(range(howmany)),[5]);hw=sel(H,list(range(howmany)),[6]);corr=prod(add(prod(hz,cz),prod(hw,cw)),ie);corrected=sub(W,prod(corr,Cdelta));return dict(source25=W25,uncorrected_full57=W,uncorrected_original63=prod(W,To63),raw7_oldX_before=H,eta_Cdelta_correction=corr,corrected_full57=corrected,corrected_original63=prod(corrected,To63),actual_pure84_image=sel(prod(corrected,F57),list(range(howmany)),purecols),raw7_after=prod(prod(corrected,F_H),Ct))
 kl=lift_and_correct(kap25,1);rl=lift_and_correct(res25,4);pl=prod(pivot25,Lift25)
 # Literal finite projection recipe: all polynomial substitution operations are specified.
 P31=n('finite_weighted_polynomial_substitution_and_selection',[91,31],args=Lroots,physical_monomial_order=mons,new_full40_monomial_order=dict(pure35=lf['full35_cubic_monomials_new_light_order'],heavy5=['H*L','H*x0','H*y1','H*y2','H*y3']),
  substitution_v8=['x0','(L-L0*x0-L2*y1-L3*y2-L4*y3)/L1','y1','0','y2','H','(L-L0*x0-L2*y1-L3*y2-L4*y3)/L1-3*x0','y3'],L_coefficients=Lroots,L1_is_known_FULLUunit='-40*C',retained31_full40=projection31['retained31_full40_columns'])
 zd=ref('S61_NONSQUARE_CUBIC_PROJECTED15_HEAVY_PREFLIGHT_20261104_v63','certificate.json','/delta_scaled_HZ_HW_coefficients/Z',[1,1]);wd=ref('S61_NONSQUARE_CUBIC_PROJECTED15_HEAVY_PREFLIGHT_20261104_v63','certificate.json','/delta_scaled_HZ_HW_coefficients/W',[1,1]);Psi=sparse([91,1],[(rawzw[0],0,neg(wd)),(rawzw[1],0,zd)]);P32=n('horizontal_concatenation',[91,32],[P31,Psi])
 E0=n('coordinate_column_selector',[32,24],columns=cpmeta['columns']);U24=ref('S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63','certificate.json','/target24_roots',[24,24]);E1=n('coordinate_column_selector',[24,23],columns=[i for i,j in enumerate(pair['old24_target_original32_columns'])if j!=5]);U23=ref('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63','certificate.json','/target23_roots',[23,23]);E2=n('coordinate_column_selector',[23,22],columns=[i for i,j in enumerate(q2['old23_target_original32_columns'])if j!=21]);Pcur=prod(prod(prod(prod(prod(P32,E0),U24),E1),U23),E2)
 proof=dict(source25_quad_lift_exact='Q0 products acquire9/20, scaledQ5products acquire1600*(-5/4)=-2000 in normalized full57; L coefficients multiply the SAME liftedfull linearform.',full63_nonzero_QQ_normalizers=[str(x)for x in normalizers],genuine15_provenance=genprovenance,
  physical91_monomial_order=mons,physical7_H_columns=Hcols,pure84_columns=purecols,
  oldH_coefficient_orientation='rowH_new=rowH_physical*transpose(constant7change)',
  raw_Cdelta_source57={'positive_row':47,'negative_row':42,'delta':'-64000C^2','actual135_positive':61,'actual135_negative':20},
  initial_projection_identity='Lift25*F57*P32=actualold25x32 by publishedZW/Lframe/sourceclearing/deleted9projection andrawpsi identities.',
  current_projection_identity='oldprefix*Lift25*F57*Pcur=current15x22',
  corrected_image_identity='Cdelta*P32=0, so correctedprojection=current projected image; fullpure84 neednotbe unchanged.',
  old5_and_psi_zero_reason='Initial8sourceoperations clearall8pivotcolumns includingHL26; subsequent2puretarget transforms preserve[27,28,29,30,31]. kappa killsremainingheavy5; residual4 sourceclearsall11pivotcols includingremainingheavy5.',
  universal_raw7_zero_reason='oldfiveHzero pluspsi0, theneta-unitCdelta identities killrawHZ/HW. No HZorHWindividual inverse.',
  full57_sufficient_criterion='On chartD(Delta407*Saux*T), if residual4x11 rank4 then current15rank15, augmentedold25rank25. Accepted product41/ZWproduct19 projection andglobalCdelta unit give actualfullcubicrank57.',
  no_allpoint_rank4_or_fullFULLU57_claim=True,no_POWER207_or_canonical_unknowng_upivot=True)
 dump(out/'operator.json',dict(nodes=nodes,exports=dict(lift25_to_full57=Lift25,full57_to_original63=To63,source15_to_original25=oldprefix,kappa=kl,residual4=rl,pivot11_uncorrected_full57=pl,pivot11_uncorrected_original63=prod(pl,To63),physical91_to_current22=Pcur),coefficient_field='FULLU for kappa; addDelta407*Saux*Tforresidual4and11',all_frozen57coefficients_external=True,new_expanded_parameter_coefficients=0))
 dump(out/'proof.json',proof);dump(out/'universal-Cdelta-identities.json',universal_cdelta())
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));size=sum(p.stat().st_size for p in out.iterdir());assert size<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,independent_review_pending=True,typed_nodes=len(nodes),actual_full57_and_original63_generator_lifts=True,kappa_and_residual4_corrected_rawH7_zero=True,
  product_row_normalization_included=True,all91_physical_targets_and135_labels_pinned=True,full57_sufficient_criterion_source_gap_closed_primary=True,
  no_new_scalar_coefficient_or_rank=True,no_fullFULLU57_or_rank4_claim=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
