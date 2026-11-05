#!/usr/bin/env python3
"""Independent fixed-source audit; deliberately does not pair the five edge values."""
import ast, hashlib, itertools, json, math, resource, signal, subprocess, sys, time
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NBAR14_KERNEL_SOURCE_20261105_v63_REVIEW'
GATE='AMBIENT_B18_NBAR14_KERNEL_SOURCE_PREFLIGHT_20261105_v63'
PINS={}
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):
    p=Path(p); p=p if p.is_absolute() else ROOT/p
    PINS[str(p.relative_to(ROOT))]=digest(p)
    return json.loads(p.read_text())
def prefix(p,key,n):
    p=ROOT/p; PINS[str(p.relative_to(ROOT))]=digest(p); s=p.read_text()
    token=json.dumps(key); assert s.count(token)==1
    i=s.index(':',s.index(token))+1
    while s[i].isspace(): i+=1
    assert s[i]=='['; i+=1; first=i; dec=json.JSONDecoder(); out=[]
    for _ in range(n):
        while s[i].isspace() or s[i]==',': i+=1
        v,i=dec.raw_decode(s,i); out.append(v)
    return out,{'key':key,'prefix':n,'bytes':[first,i],'full_array_not_decoded':True}
def pairs(x): return {int(i):F(v) for i,v in x if F(v)}
def enc(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,dict): return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return list(map(enc,x))
    return x
def slots(x):
    if isinstance(x,(int,F)): return 1
    if isinstance(x,str):
        try: F(x); return 1
        except ValueError: return 0
    if isinstance(x,dict): return sum(slots(v) for v in x.values())
    if isinstance(x,(list,tuple)): return sum(map(slots,x))
    return 0
def bound(terms):
    # Positive majorant, never the signed functional pairing.
    den=math.lcm(*(v.denominator for v in terms))
    num=sum(abs(v.numerator)*(den//v.denominator) for v in terms)
    return den,num
def sparse(v): return [[i,x] for i,x in enumerate(v) if x]

def main():
    begin=time.monotonic(); signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/STEM; out.mkdir(exist_ok=False)
    gp=f'reports/{GATE}/report.json'; report=load(gp)
    assert report['all_passed'] and not report['actual_edge_scalar_or_kernel_matrix_values_computed']
    for p,h in report['input_sha256'].items(): assert digest(ROOT/p)==h
    source=load(f'reports/{GATE}/five-cache-functionals-and-genuine-fourteen-kernel-source.json')
    assert PINS[f'reports/{GATE}/five-cache-functionals-and-genuine-fourteen-kernel-source.json']==report['evidence_sha256']['five-cache-functionals-and-genuine-fourteen-kernel-source.json']
    paths=report['input_sha256']
    def locate(tail):
        hits=[p for p in paths if p.endswith(tail)]; assert len(hits)==1; return hits[0]
    ordinary={}
    for key,tail in [('m0','phase1-actual-m0-and-section.json'),('first','phase1-actual-fullnine-first-V-columns.json'),('Q2','phase4-all-twenty-actual-quadratic-values.json'),('G3','actual-Cprime9-actions-and-genuine-G3-target.json'),('dual','actual-curve-source-and-one-residual-cubic-functional.json')]: ordinary[key]=load(locate(tail))
    for p in paths:
        if p.endswith('REVIEW/report.json'):
            peer=load(p); assert peer['all_passed']
    curvepath=[p for p in paths if 'CPRIME_CURVE_GA_CUBIC_VALUES_' in p and p.endswith('/report.json') and 'REVIEW' not in p]
    assert len(curvepath)==1 and load(curvepath[0])['actual_residual_scalar']=='0'
    mn,mp=prefix(locate('actual-minimal-Tnil5-849-and-d3-linear-operator.json'),'operator_nodes',9)
    cn,cp=prefix(locate('actual-formal-connecting16-and-ten-source-obstructions.json'),'matrix_operator_nodes',11)
    dif,dp=prefix(locate('actual-reduced-Cbar-three-sector-Koszul.json'),'all8_differentials_sparse_columns',2)
    f=ordinary['m0']; m=[[F(c) for c in row] for row in f['m0_24_by44']]
    E=f['selected_source_columns24']; inv=[[F(c) for c in row] for row in f['selected_minor_inverse']]
    n=[[F(0)]*24 for _ in range(44)]
    for j,i in enumerate(E): n[i]=inv[j]
    assert all(sum(m[i][k]*n[k][j] for k in range(44))==int(i==j) for i in range(24) for j in range(24))
    assert mn[5]==['QQ_matrix_rows',[[str(int(i==j)) for j in E] for i in range(44)]]
    assert mn[6]==['matrix_product',[0,5]] and mn[7][:2]==['origin_unit_matrix_inverse',6] and mn[8]==['matrix_product',[5,7]]
    assert cn[0]==['actual_external_matrix','minimal','matrices.m']
    assert cn[1]==['actual_external_matrix','minimal','matrices.direct_unit_section_n']
    assert cn[7]==['matrix_product',[1,2]] and cn[8]==['horizontal_blocks',[6,7]] and cn[10]==['matrix_product',[8,9]]
    heads=[42,43,50,68,70]
    assert cn[9]==['QQ_matrix_rows',[[str(int(i==j)) for j in heads] for i in range(128)]]
    assert all(row[43]==0 for row in m)
    c0=[F(i==42)-sum(n[i][j]*m[j][42] for j in range(24)) for i in range(44)]
    assert pairs(source['head42_origin_cycle44_sparse'])==dict(sparse(c0))
    assert source['source_m0_24_by44']==enc(m) and source['source_n0_44_by24']==enc(n)
    comb=list(itertools.combinations(range(8),2)); d2r=[]
    for k,j in enumerate([6,24,26]):
        assert comb[j]==([0,4,5][k],7)
        v=[F(0)]*24
        for row,mon,c in dif[1][j]:
            assert sum(mon)>0
            if mon==[0]*7+[1]: v[row]+=F(c)
        assert v==[F(-int(i==[0,4,5][k])) for i in range(24)]
        d2r.append(v)
    assert all(sum(mon)>0 for col in dif[1] for row,mon,c in col)
    # The actual16 basis identifies exactly the four retained first classes.
    con=load(locate('actual-formal-connecting16-and-ten-source-obstructions.json'))
    basis=load(con['source_namespace_bindings']['basis'])['explicit_sixteen_filtered_basis_lifts']
    assert [(v['target'],v['monomial']) for v in basis[5:10]]==[(0,[5]),(0,[7]),(1,[7]),(3,[7]),(4,[7])]
    assert all(v['degree']>=2 for v in basis[10:]) and con['genuine_kernel_dimension']==14
    linear=[pairs(c) for c in ordinary['first']['all44_actual_linear_columns27_sparse']]
    V=[[[col.get(t*9+v,F(0)) for col in linear] for t in range(3)] for v in range(9)]
    Q2=[[F(c) for c in row] for row in ordinary['Q2']['all44_actual_V2_pureplane_nine_values']]
    actions=[[[F(x) for x in col] for col in mat] for mat in ordinary['G3']['actual_R0_actions9_by9_columns']]
    assert ordinary['G3']['actual_genuine_G3_action_scalars9']==['0']*7+['-11/16','0']
    W=[[-actions[v][t][7] for t in range(3)] for v in range(9)]
    phi=pairs(ordinary['dual']['residual_actual_cubic_dual495_sparse'])
    mons=list(itertools.combinations_with_replacement(range(9),3))
    phirow=lambda mon:[phi.get(t*165+mons.index(mon),F(0)) for t in range(3)]
    b=phirow((7,7,7)); a=phirow((7,7,8)); c=phirow((7,8,8))
    assert a==[F(-1,75),F(0),F(-14,55)] and b==[F(92,20625),F(0),F(3756,75625)] and c==[F(0),F(0),F(1)]
    assert all(sum(a[t]*V[8][t][j]+b[t]*V[7][t][j] for t in range(3))==0 for j in range(44))
    zscale=[F(2,25)*W[z][0]+W[z][2]/75 for z in range(7)]
    assert all(W[z][1]==0 for z in range(7))
    assert all(sum(W[z][t]*V[8][t][j] for t in range(3))==zscale[z]*V[7][0][j] for z in range(7) for j in range(44))
    rho=[sum(b[t]*Q2[j][3*t]+a[t]*Q2[j][3*t+1]+c[t]*Q2[j][3*t+2]-W[7][t]*V[8][t][j] for t in range(3))+sum(zscale[z]*V[z][0][j] for z in range(7)) for j in range(44)]
    assert pairs(source['remaining_pure_r_first_cycle_adjoint44_sparse'])==dict(sparse(rho))
    assert pairs(source['lambda0_source_row44_sparse'])==dict(sparse(V[8][1]))
    assert [i for i,v in enumerate(rho) if v]==[0,4,21,22,24,28,32,36,40,41]
    assert rho[22]==F(11807,48400) and rho[20]==0
    assert all(Q2[43][i]==0 for i in range(9))
    alpha=F(11,12)*F(ordinary['dual']['gA_changed_r_ell_coefficients'][0]); assert alpha==F(1375,288)
    s1=[[-alpha*n[i][23] for i in range(44)]]+[[-n[i][j] for i in range(44)] for j in [0,4,5]]
    assert [pairs(v) for v in source['all4_other_head_pure_r_first_cycles44_sparse']]==[dict(sparse(v)) for v in s1]
    # Majorants use ABSOLUTE products only: no five signed dots here.
    Amaj=bound([V[8][1][i]*c0[i] for i in range(44)])
    Bmaj=[bound([rho[i]*v[i] for i in range(44)]) for v in s1]
    ratio_bound=max(max((16*nb*Amaj[0]).bit_length(),(11*db*Amaj[1]).bit_length()) for db,nb in Bmaj)
    assert ratio_bound<=82
    custody=load('reports/AMBIENT_B18_NBAR14_KERNEL_PREPARED_WORKER_20261105_v63/report.json')
    driver=ROOT/custody['prepared_code']; PINS[custody['prepared_code']]=digest(driver)
    assert digest(driver)==custody['prepared_code_sha256']=='32620bc12a49d40cf81e6400204c8a455396fa29272c0430cd2f8556e115b452'
    tree=ast.parse(driver.read_text()); imports=[n for n in ast.walk(tree) if isinstance(n,(ast.Import,ast.ImportFrom))]
    assert all(not isinstance(n,ast.ImportFrom) or n.module in ['fractions','pathlib'] for n in imports)
    assert not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id in ['eval','exec','__import__'] for n in ast.walk(tree))
    helpcall=subprocess.run([sys.executable,str(driver),'--help'],capture_output=True,text=True,timeout=5)
    assert helpcall.returncode==0 and '--source-peer-report' in helpcall.stdout
    assert not (ROOT/'reports/AMBIENT_B18_NBAR14_KERNEL_CACHED_VALUES_20261105_v63').exists()
    owned=slots([source,ordinary,m,n,c0,basis,V,Q2,actions,W,rho,s1,con,mn,cn,dif])
    assert owned<60000 and slots(source)*4+5000<30000
    data={'independent_direct_section_and_connecting_heads':True,'all576_point_section_identities':True,'actual_three_wedge_derivatives':d2r,'four_head_s1r':list(map(sparse,s1)),
          'independent_sigma2_rows':[a,b],'all44_sigma2_identities':True,'independent_sigma1z_coefficients':zscale,'all308_sigma1z_identities':True,'independent_full_unnormalized_rho44':sparse(rho),
          'native_curve43_is_not_raw_f63':True,'curve_order3_kills_V2_all_directions':'accepted curve source and curve scalar peer; actual pureplane zeros independently checked',
          'source_bounds_without_pairing':{'a_positive_majorant':Amaj,'four_b_positive_majorants':Bmaj,'kernel_coefficient_ratio_bit_upper':ratio_bound},
          'exact_basis_first5':basis[5:10],'source_only_no_five_edge_pairing':True,'formal_beta_never_zero_filled':True,'free_I14_and_edge_kernel_universal_proof':'retain rows1,2,3,4,5,7,8,9,10..15; colj=e_j-(b_j/xi)e6 for top1..4; edge first row ae0 and second beta e0+sum bje_j+xi e6. Every col has e0 coefficient0. a*xi nonzero by surjectivity/cyclicity; xi=-11a/16.',
          'selective_prefixes':[mp,cp,dp],'source_and_prepared_help_checked':True,'independent_owned_numeric_upper':owned,'future_static_numeric_bound':30000,'hard60000':True,'raw_source_custody_bytes':sum((ROOT/p).stat().st_size for p in report['input_sha256'])+len(json.dumps(source)),
          'accepted_scope':'genuine N/ellN barcycle kernel only, no full R0 nil lifts, generators or products; beta remains formal'}
    ep=out/'independent-five-functional-source-and-symbolic-kernel-proof.json'; ep.write_text(json.dumps(enc(data),separators=(',',':'))+'\n')
    result={'id':STEM,'all_passed':True,'status':'SOURCE_PASS_exact_five_cached_functionals_and_actual14_kernel_recipe','input_sha256':PINS,'code_sha256':digest(Path(__file__)),'note_sha256':digest(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:digest(ep)},'actual_five_edge_values_paired':False,'prepared_code_sha256':digest(driver),'seconds':time.monotonic()-begin,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n'); signal.alarm(0); print(json.dumps(result))
if __name__=='__main__': main()
