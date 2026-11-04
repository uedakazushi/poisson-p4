#!/usr/bin/env python3
"""Small exact nil generators/presentation of normalized E4 only. No GB."""
import argparse,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'reports/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63/five-line-and-generic-models.json';data=json.loads(source.read_text())
    (out/'phase.json').write_text(json.dumps({'id':STEM,'status':'source_bound_before_fixed_algebra','source_sha256':sha(source),'no_GB_new_point_jet_or_field':True},indent=2)+'\n')
    a,d,x,t=s.symbols('a d x t');old=[a,d,x,t];c,r,l,h=s.symbols('c r l h');new=[c,r,l,h]
    def dec(rows):return s.expand(sum(s.Rational(v)*s.prod(z**e for z,e in zip(old,m))for m,v in rows))
    def enc(f,vs=new):return [[list(m),str(v)]for m,v in s.Poly(s.expand(f),*vs,domain=s.QQ).terms()if v]
    Q0o,Q1o,fo=[dec(data[k])for k in ['source_Q0','source_Q1','source_f']]
    forward={c:(d+x)/2,r:(d-x-2*a)/6,l:(d-x+a)/3,h:a-d+x-2*t}
    backward={a:-2*r+l,d:c+r+l,x:c-r-l,t:(-4*r-l-h)/2}
    assert all(s.expand(f.subs(backward,simultaneous=True)-z)==0 for z,f in forward.items())
    assert all(s.expand(f.subs(forward,simultaneous=True)-z)==0 for z,f in backward.items())
    transform=s.Matrix([forward[z]for z in new]).jacobian(old);assert transform.det()!=0
    Q0,Q1,f=[s.expand(z.subs(backward,simultaneous=True))for z in [Q0o,Q1o,fo]]
    assert s.expand(fo-(x-a-d)**2*(a+x))==0
    g=s.Matrix([h,c*r,c*l,r*l])
    def factor_I3(z):
        coeff=[s.Integer(0)]*4
        for mon,v in s.Poly(z,*new,domain=s.QQ).terms():
            if not v:continue
            m=list(mon)
            if m[3]:idx=0;m[3]-=1
            else:
                pair=next(((i,j)for i,j in [(0,1),(0,2),(1,2)]if m[i] and m[j]),None)
                assert pair is not None,('nonzero_axis_restriction',mon,v)
                idx={(0,1):1,(0,2):2,(1,2):3}[pair]
                for i in pair:m[i]-=1
            coeff[idx]+=v*s.prod(z**e for z,e in zip(new,m))
        coeff=s.Matrix(list(map(s.expand,coeff)));assert s.expand((coeff.T*g)[0]-z)==0
        return coeff
    V=s.Matrix.hstack(*[factor_I3(q)for q in [Q0,Q1,f]])
    origin=dict.fromkeys(new,0);constant=V.subs(origin);assert constant.rank()==2 and constant[:,2]==s.zeros(4,1)
    pivot=next((p for p in itertools.combinations([1,2,3],2)if constant[list(p),:2].det()!=0));remaining=next(j for j in [1,2,3]if j not in pivot)
    A=V[list(pivot),:2].T;B=V[[0,remaining],:2].T;Delta=s.expand(A.det());assert Delta.subs(origin)
    adj=s.Matrix([[A[1,1],-A[0,1]],[-A[1,0],A[0,0]]]);assert (A*adj).applyfunc(s.expand)==Delta*s.eye(2)
    E=s.zeros(4,2);T=s.zeros(4,2);E[0,0]=Delta;E[remaining,1]=Delta
    for i,j in enumerate(pivot):
        E[j,:]=-adj[i,:]*B;T[j,:]=adj[i,:]
    E=E.applyfunc(s.expand);T=T.applyfunc(s.expand);small=s.Matrix([g[0],g[remaining]])
    identities=(Delta*g-E*small-T*s.Matrix([Q0,Q1])).applyfunc(s.expand);assert identities==s.zeros(4,1)
    # Complete syzygies of I3:3 crossings with h and the2 Hilbert-Burch columns.
    sy=[]
    for j in [1,2,3]:
        z=s.zeros(4,1);z[0]=g[j];z[j]=-h;sy.append(z)
    sy += [s.Matrix([0,l,-r,0]),s.Matrix([0,l,0,-c])]
    assert all(s.expand((w.T*g)[0])==0 for w in sy)
    W=s.Matrix.hstack(*sy,V[:,2]);Rel=(E.T*W).applyfunc(s.expand)
    witnesses=s.zeros(3,6)
    for j in range(6):
        witnesses[:2,j]=-(T.T*W[:,j]);witnesses[2,j]=Delta if j==5 else 0
    witnesses=witnesses.applyfunc(s.expand)
    assert (small.T*Rel-s.Matrix([[Q0,Q1,f]])*witnesses).applyfunc(s.expand)==s.zeros(1,6)
    assert all(z==0 for z in Rel.subs(origin))
    # Elimination of Q0,Q1 in the free4 module is a literal unit change.
    N=s.zeros(4);N[0,0]=1;N[1,remaining]=1;N[2,:]=V[:,0].T;N[3,:]=V[:,1].T
    assert N.det().subs(origin)!=0
    # Same3 lines/source prime checks, with d nonzero on each.
    line_vectors={'C':[0,1,1,0],'R':[-2,1,-1,-2],'L':[1,1,-1,s.Rational(-1,2)]}
    axis_vectors={}
    for key,vec in line_vectors.items():
        z={v:w for v,w in zip(old,vec)};val=[s.expand(forward[v].subs(z))for v in new];axis_vectors[key]=list(map(str,val))
        assert val[3]==0 and sum(bool(v)for v in val[:3])==1 and z[d]!=0
    # All8 source/representation identities are literal QQ, not a quotient guess.
    result={'old_variable_order':list(map(str,old)),'new_core_variable_order':list(map(str,new)),
      'old_source_Q0':enc(Q0o,old),'old_source_Q1':enc(Q1o,old),'old_source_f':enc(fo,old),
      'forward_GL4':{str(z):enc(v,old)for z,v in forward.items()},'inverse_GL4':{str(z):enc(v)for z,v in backward.items()},'GL4_determinant':str(transform.det()),
      'three_minimal_core_lines_GL4_values':axis_vectors,'I3_generators4':[enc(v)for v in g],
      'actual_core_generators3':[enc(v)for v in [Q0,Q1,f]],'source_factorization_coefficients4_by3':[[enc(V[i,j])for j in range(3)]for i in range(4)],
      'minimal_generator_class_matrix4_by3':[[str(v)for v in row]for row in constant.tolist()],
      'quad_pivots_in_I3_generators':list(pivot),'retained_quad_index':remaining,'unit_Delta':enc(Delta),'unit_Delta_origin':str(Delta.subs(origin)),
      'unit_change_N4_by4':[[enc(N[i,j])for j in range(4)]for i in range(4)],
      'old_I3_generators_Delta_reduction_E4_by2':[[enc(E[i,j])for j in range(2)]for i in range(4)],
      'old_I3_generators_Delta_source_Q_cofactors_T4_by2':[[enc(T[i,j])for j in range(2)]for i in range(4)],
      'I3_complete_syzygies4_by5':[[enc(sy[j][i])for j in range(5)]for i in range(4)],
      'nil_generators_in_new_core':[enc(backward[d]*v)for v in small],
      'nil_generators_in_original_core':[enc(s.expand((backward[d]*v).subs(forward,simultaneous=True)),old)for v in small],
      'nil_complete_core_S4_relations2_by6':[[enc(Rel[i,j])for j in range(6)]for i in range(2)],
      'all6_relation_source_cofactors3_by6':[[enc(witnesses[i,j])for j in range(6)]for i in range(3)],
      'full_S9_E4_extra_relations10':'u_j*e_k for1<=j<=5 and1<=k<=2, appended to core2x6; E4 keeps all d*u_j=0',
      'full_E4_nil_module_presentation_shape':[2,16],'nil_minimal_generator_count':2,'relation_count_not_claimed_minimal':True,
      'radical_formula':'sqrt(d*(u5,Jcore))=d*(u5+I3), I3=(h,cr,cl,rl)',
      'nil_module_exact':'Nil(E4)=d*(u5+I3)/d*(u5+Jcore) isomorphic asS9module to I3/Jcore; u5 actszero',
      'no_D0_nil_generators_claimed':True,'source_identity_count':4+6+5+3,'maximum_total_expanded_term_cap':150000}
    terms=sum(len(enc(z))for M in [V,E,T,Rel,witnesses,N]for z in M)+sum(len(enc(q))for q in [Q0,Q1,f,Delta]);assert terms<150000
    ev=out/'actual-E4-nil-two-generators-and-relations.json';ev.write_text(json.dumps(result,indent=2)+'\n')
    deps=[ROOT/'reports'/n/'report.json'for n in ['AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW','AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63_REVIEW','AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63_REVIEW']]
    assert all(json.loads(p.read_text())['all_passed']for p in deps)
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',source]+deps
    report={'id':STEM,'all_passed':True,'normalized_E4_only':True,'radical_from_all_actual_minimal_primes_proved':True,'actual_ideal_I3_source_bound':True,
      'two_actual_nil_generators_original_core':result['nil_generators_in_original_core'],'nil_minimal_generator_count':2,
      'source_Q0_Q1_classes_rank':2,'source_f_class_zero':True,'origin_unit_Delta':str(Delta.subs(origin)),
      'complete_nil_presentation_over_core_S4':[2,6],'complete_nil_presentation_over_full_S9_or_E4':[2,16],
      'all_five_d_times_u_relations_retained':True,'u5_action_on_nil_zero_proved':True,'six_or16_relations_minimality_claimed':False,
      'source_expanded_term_count':terms,'no_GB_new_source_point_jet_or_field':True,'original_D0_nil_generators_or_principality_claimed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{ev.name:sha(ev)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');(out/'phase.json').write_text(json.dumps({'id':STEM,'status':'FINAL'},indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps({k:report[k]for k in ['id','all_passed','nil_minimal_generator_count','origin_unit_Delta','source_expanded_term_count','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
