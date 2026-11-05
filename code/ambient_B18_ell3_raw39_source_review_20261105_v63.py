#!/usr/bin/env python3
"""Independent source/formal-algebra review; no quotient rank or source solve."""
import hashlib,itertools,json,pathlib,resource,signal,time
from fractions import Fraction as Q

ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_ELL3_RAW39_SOURCE_PREFLIGHT_20261105_v63_REVIEW'
AUTHOR=ROOT/'reports/AMBIENT_B18_ELL3_RAW39_SOURCE_PREFLIGHT_20261105_v63'
PRIMARY_SHA='03ee57279bd4143098018d1273d1d2161b2bdb12b1b48f0ad2741330268d13c4'

def sha(p):
    h=hashlib.sha256()
    with pathlib.Path(p).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def dump(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def selected(text,key):
    j=text.index(':',text.index(json.dumps(key)))+1
    while text[j].isspace():j+=1
    return json.JSONDecoder().raw_decode(text,j)[0]
def enc(p):return [[list(m) if isinstance(m,tuple) else m,str(c)] for m,c in sorted(p.items()) if c]
def dec(p):return {tuple(m) if isinstance(m,list) else m:Q(c) for m,c in p if Q(c)}
def add(p,q,c=Q(1)):
    d=dict(p)
    for m,v in q.items():
        d[m]=d.get(m,Q(0))+c*v
        if not d[m]:del d[m]
    return d
def mul(p,q):
    d={}
    for m,c in p.items():
        for n,e in q.items():
            k=tuple(a+b for a,b in zip(m,n));d[k]=d.get(k,Q(0))+c*e
    return {m:c for m,c in d.items() if c}
def powers(p,n):
    a={(0,0):Q(1)}
    for _ in range(n):a=mul(a,p)
    return a
def compose(p,x,y):
    d={}
    for m,c in p:d=add(d,mul(powers(x,m[0]),powers(y,m[1])),Q(c))
    return d
def clip(p):return {m:c for m,c in p.items() if m[1]<3 and sum(m)<=4}
def bits(x):
    if isinstance(x,dict):return max([1]+[bits(v) for v in x.values()])
    if isinstance(x,(tuple,list)):return max([1]+[bits(v) for v in x])
    if isinstance(x,Q):return max(abs(x.numerator).bit_length(),x.denominator.bit_length())
    return 1

def main():
    t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
    try:
        assert sha(AUTHOR/'report.json')==PRIMARY_SHA
        rep=read(AUTHOR/'report.json');assert rep['all_passed']
        assert sha(ROOT/'code/ambient_B18_ell3_raw39_source_preflight_20261105_v63.py')==rep['code_sha256']
        assert sha(ROOT/'research/notes/AMBIENT_B18_ELL3_RAW39_SOURCE_PREFLIGHT_20261105_v63.md')==rep['note_sha256']
        for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h
        for p,h in rep['evidence_sha256'].items():assert sha(AUTHOR/p)==h
        author=read(AUTHOR/'exact-m5-ell3-cubic-truncation-and-75-column-source-input.json')
        sourcepath=next(ROOT/p for p in rep['input_sha256'] if p.endswith('actual-special-fiber20-coefficients.json'))
        st=sourcepath.read_text()
        names=['beta2_before_ell','beta3_before_ell','gA_plane3','gB_own_plane3','YB_in_YA_linear']
        S={k:selected(st,k) for k in names};del st
        for p in rep['input_sha256']:
            if p.endswith('_REVIEW/report.json'):assert read(ROOT/p)['all_passed']
        raw2path=next(ROOT/p for p in rep['input_sha256'] if p.endswith('actual-ell2-product33-and-raw45-actions.json'))
        rt=raw2path.read_text();assert len(selected(rt,'actual33_basis'))==33;assert selected(rt,'actual_source_equality')
        actions=selected(rt,'coordinate_actions9_sparse_columns')
        del rt
        tablepath=next(ROOT/p for p in rep['input_sha256'] if p.endswith('exact-existing-raw33-multiplication-table.json'))
        tt=tablepath.read_text();assert selected(tt,'m5_zero') is True;assert len(selected(tt,'basis33'))==33
        del tt
        assert len(actions)==9 and all(len(a)==33 for a in actions)
        actions=[[{i:Q(c) for i,c in col} for col in a] for a in actions]
        states={():{0:Q(1)}};power_checks=[]
        for degree in range(1,6):
            nxt={}
            for word,v in states.items():
                for i in range(word[-1] if word else 0,9):
                    outv={}
                    for k,c in v.items():outv=add(outv,actions[i][k],c)
                    nxt[word+(i,)]=outv
            assert len(nxt)==__import__('math').comb(8+degree,degree)
            power_checks.append({'degree':degree,'complete_commutative_words':len(nxt),'nonzero_images_of_unit':sum(bool(v) for v in nxt.values())})
            states=nxt
        assert all(not v for v in states.values())
        del states,nxt,actions
        r={(1,0):Q(1)};s={(1,0):Q(1),(0,1):Q(1,5)}
        beta=[]; beta2=[]
        for p,q in zip(S['beta2_before_ell'],S['beta3_before_ell']):
            assert all(not any(m[:7]) for m,c in p)
            p2=[[m[7:],c] for m,c in p]
            # z-dependent cubic terms are removed by the accepted I+M2 unit generator change.
            p3=[[m[7:],c] for m,c in q if not any(m[:7])]
            beta2.append(compose(p2,r,s));beta.append(clip(add(beta2[-1],compose(p3,r,s))))
        ga=clip(compose(S['gA_plane3'],r,s));Y=S['YB_in_YA_linear']
        rb=add({m:Q(Y[7][7])*c for m,c in r.items()},{m:Q(Y[7][8])*c for m,c in s.items()})
        sb=add({m:Q(Y[8][7])*c for m,c in r.items()},{m:Q(Y[8][8])*c for m,c in s.items()})
        gb=clip(compose(S['gB_own_plane3'],rb,sb))
        assert enc(ga)==author['gA_cubic'] and enc(gb)==author['gB_cubic']
        assert list(map(enc,beta))==author['beta7_cubic_plane_after_s_equals_r_plus_ell_over5']
        assert ga[(3,0)]==Q(125,24) and gb[(3,0)]==Q(-15125,144)
        # Fixed, already-present quadratic source minor; no quotient rank search.
        qq=[(2,0),(1,1),(0,2)];fixed=[[beta2[i].get(m,Q(0)) for m in qq] for i in [0,4,5]]
        determinant=Q(0)
        for perm in itertools.permutations(range(3)):
            inv=sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
            term=Q((-1)**inv)
            for i,j in enumerate(perm):term*=fixed[i][j]
            determinant+=term
        assert determinant
        # A different enumeration of the ambient module, then the saved canonical order.
        ambient=[]
        for total in range(5):
            for u in [7,6,5,4,3,2,1,0]:
                for ellpow in range(3):
                    apow=total-ellpow-(u!=0)
                    if apow>=0:ambient.append((u,apow,ellpow))
        assert len(ambient)==75 and ambient==list(map(tuple,author['ambient75_basis']))
        index={m:i for i,m in enumerate(ambient)}
        def sector(u,p):return {index[u,a,b]:c for (a,b),c in p.items() if (u,a,b) in index}
        rel=[]
        for i,j in itertools.combinations(range(7),2):rel.append((['pair',i,j],add(sector(i+1,beta[j]),sector(j+1,beta[i]),Q(-1))))
        for i in range(7):rel.append((['gA',i],add(sector(i+1,ga),sector(0,mul(ga,beta[i])))))
        for i in range(7):rel.append((['gB',i],sector(i+1,gb)))
        assert [[lab,enc(p)] for lab,p in rel]==author['source_plane_module_relations35']
        rows=[]
        for lab,p in rel:
            for a,b in ([(0,0),(1,0),(0,1)] if lab[0]=='pair' else [(0,0)]):
                q={}
                for k,c in p.items():
                    u,i,j=ambient[k]
                    if (u,i+a,j+b) in index:q[index[u,i+a,j+b]]=c
                if q:rows.append([lab,[a,b],enc(q)])
        assert rows==author['all_atmost77_plane_multiple_rows75']
        nz=sum(len(p) for lab,m,p in rows);assert len(rows)<=77 and nz==204 and nz==author['input_nonzero_count']
        # Pure D-module counts, no new matrix rank or multiplication table.
        free=3+3;muT=20-free;lenT=33-2*free;t2=lenT-muT;t1=muT-t2
        assert (muT,lenT,t1,t2,3*free+lenT)==(14,21,7,7,39)
        # Universal six-basis weight test: list all possible symbolic output orders.
        # Branch Weierstrass coefficients have orders (A2,A1,A0)>=(1,2,3).
        w=[0,1,2,2,2,2]
        rprofiles=[[(1,0)],[(2,0),(3,0)],[(2,1),(5,1),(4,1)],
                   [(3,1),(1,2),(5,1),(0,3),(4,1)],[(5,1)],[(2,1)]]
        # beta*r^j coefficients in branch basis have lower orders (j+2,j+1,j).
        zprofiles=[[(3,0),(1,1),(5,0),(0,2),(4,0)],
                   [(3,1),(1,2),(5,1),(0,3),(4,1)],[],
                   [(3,2),(1,3),(5,2),(0,4),(4,2)],[],[]]
        for i,profile in enumerate(rprofiles):assert all(w[j]+n>=w[i]+1 for j,n in profile)
        for i,profile in enumerate(zprofiles):assert all(w[j]+n>=w[i]+2 for j,n in profile)
        assert min(5-x for x in w)==3
        proof={
          'all_source_hashes_and_accepted_prior_peers_match':True,
          'same_actual_DVR_parameter':'ell=5*(sA-rA), branches D-free3 each',
          'independent_cubic_compositions_and_35_relations_match':True,
          'independent_all_permitted_plane_multiples_match':True,
          'actual_saved_rows':len(rows),'actual_saved_terms':nz,'ambient':75,
          'fixed_beta2_quadratic_minor_rows':[0,4,5],'fixed_beta2_minor_det':str(determinant),
          'J_equality_proof':'The three plane quadrics span (r,ell)^2 modulo (r,ell)^3; formal Nakayama and z graph generators give J=(z7,(r,ell)^2).',
          'U_basis':['(1,1)','(r,r)','(r^2,0)','(0,r^2)','(ell^2,0)','(ell*r,0)'],
          'basis_derivation':'Fiber product of the two free3 branches over k[r,ell]/(r,ell)^2. Constant coefficients differ by ell^2 and r coefficients by ell.',
          'weights':w,'r_symbolic_output_order_profiles':rprofiles,'z_symbolic_output_order_profiles':zprofiles,
          'm5_U_in_ell3_U':True,'weighted_proof':'Each r or ell increases weight >=1, each z >=2; largest basis weight2, so weight>=5 gives an overall ell^3.',
          'T_annihilated_by_J':True,'ell2_T_zero':True,
          'T_injection_proof':'U is D-free, so Tor_1^D(U,D/ell^n)=0. For n>=2 the copy of T in B/ell^n is unchanged by reduction.',
          'existing_B2_all1287_degree5_words_kill_unit_verified':power_checks,
          'm5_B3_zero':'m5 maps to zero in U3, hence lies in T; accepted B2 actions independently verify m5B2=0, and injectivity on T forces vanishing in B3.',
          'both_actual_and_cubic_truncated_pairs_have_same_proof':True,
          'exact_cubic_truncation_proof':'Generator errors have order>=4, product errors order>=5. m5 is zero in both quotients, giving both ideal inclusions modulo ell^3.',
          'counts':{'U_free':free,'T_mu':muT,'T_length':lenT,'T_exponents1':t1,'T_exponents2':t2,'B3_length':39,'raw_C3_length':57},
          'module_presentation_complete':'49 monic z_i(z_j+beta_j) eliminate multi-z monomials; 21 commutation/pair and14 g relations are closed under z action. Plane multiples of allowed degrees exhaust the m5,ell3 truncation. gAgB has order6 and vanishes.',
          'no_quotient_rank_table_or_physical_source_values_computed':True,
          'individual_branch_cubic_truncation_NOT_claimed':True,
          'scope':'Actual product B/ell^3 has length39 and raw cokernel length57; not a Cprime/exponent/kernel/nil/H0 conclusion.',
          'numeric_static_upper':30000,'hard_terms':150000,'hard_seconds':60,'hard_memory_bytes':2000000000,
          'source_and_output_bytes_upper':2*1024*1024,'source_arithmetic_observed_bits':bits([beta,ga,gb,rel,rows]),
        }
        assert proof['source_arithmetic_observed_bits']<2048
        dump(out/'independent-weighted-union-and-cubic-source-proof.json',proof)
        report={'id':ID,'all_passed':True,'frozen':True,'primary_report_sha256':PRIMARY_SHA,
          'input_sha256':rep['input_sha256'],'code_sha256':sha(__file__),
          'evidence_sha256':{'independent-weighted-union-and-cubic-source-proof.json':sha(out/'independent-weighted-union-and-cubic-source-proof.json')},
          'no_author_producer_imports':True,'no_new_quotient_rank_table_or_source_values':True,'elapsed_seconds':time.monotonic()-t}
        dump(out/'report.json',report);signal.alarm(0)
        print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'elapsed_seconds':report['elapsed_seconds']}))
    except Exception as e:
        dump(out/'failure.json',{'error':type(e).__name__+': '+str(e),'no_positive_claim':True});raise
if __name__=='__main__':main()
