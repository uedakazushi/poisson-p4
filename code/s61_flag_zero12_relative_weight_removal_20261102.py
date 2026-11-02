"""Relative nilpotence from a finite kernel over a nonpositively weighted base."""
import argparse,hashlib,json,os,subprocess,tempfile,time
from pathlib import Path
import sympy as s
from s61_model import S61,BASIS,X
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_FLAG_ZERO12_RELATIVE_WEIGHT_REMOVAL_20261102'
NOTE='research/notes/'+STEM+'.md'

def verify():
    start=time.monotonic();m=S61()
    p0='reports/S61_FLAG_ZERO12_ZERO_NUISANCE_FORMAL_20261102.json'
    pd='reports/S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101.json'
    r=json.loads((ROOT/p0).read_text());d=json.loads((ROOT/pd).read_text())
    assert r['all_passed'] and d['all_passed']
    assert r['all135_full_polynomial_substitutions_checked']
    assert r['implicit_coordinate_degree_bound']==6 and r['residual_weight_bound']==10
    point=s.Matrix(r['point61']).applyfunc(s.Rational);u=s.symbols('u1:62');sub=dict(zip(u,point))
    K=s.Matrix(d['kernel']).applyfunc(s.sympify).subs(sub)
    C=s.Matrix(d['actual_C_tangent']).applyfunc(s.sympify).subs(sub)
    N=K[:,d['squarezero_normal_columns']];B=C.row_join(N);cut=[19,26]
    H=B*s.Matrix.hstack(*B[cut,:].nullspace());Cc=C*s.Matrix.hstack(*C[cut,:].nullspace())
    assert H.rank()==5 and Cc.rank()==1
    qw=[0 if i==25 else r['raw_coordinate_weights'][i-(i>25)]for i in range(45)]
    W=s.diag(*qw);assert H.row_join(W*H).rank()==5 and Cc.row_join(W*Cc).rank()==1
    def characters(M):return {w:M[[i for i,a in enumerate(qw)if a==w],:].rank()for w in sorted(set(qw))if M[[i for i,a in enumerate(qw)if a==w],:].rank()}
    hc=characters(H);cc=characters(Cc);assert hc=={1:3,2:2} and cc=={1:1}
    # Exact full polynomial algebra, not a finite jet initial ideal.
    def expr(a):
        return '+'.join('('+c+')'+''.join('*x_'+str(i)+(('^'+str(n))if n!=1 else '')for i,n in enumerate(e)if n)for e,c in a)or'0'
    with tempfile.TemporaryDirectory(prefix='zero12-relative-hilbert-')as tmp:
        script=Path(tmp)/'verify.m2'
        program='R=QQ[x_0..x_6,Degrees=>{1,1,2,2,1,2,2}];I=ideal('+','.join(expr(f)for f in r['polynomial_residuals'])+');H=hilbertSeries(R/I,Reduce=>true);A=ring numerator H;z=A_0;assert(numerator H==1+2*z+7*z^2+12*z^3+18*z^4+9*z^5-34*z^6-10*z^7);assert(denominator H==1-z);D=1+2*z+2*z^2;E=5*z^2+17*z^3+35*z^4+44*z^5+10*z^6;assert(numerator H-D==(1-z)*E);print("HILBERT EXACT PASS");'
        script.write_text(program)
        env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'}
        rr=subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary','--script',str(script)],capture_output=True,text=True,timeout=35,env=env)
        assert rr.returncode==0 and 'HILBERT EXACT PASS'in rr.stdout,(rr.stdout,rr.stderr)
    # Raw parameter weights, including the nonpositive weights of every H14 base parameter.
    torus=[0,1,-1,0,0];weights=[]
    for j,index in enumerate(m.certificate['raw_free_indices']):
        pair,mon=BASIS[index];ep=s.Poly(mon,*X).monoms()[0]
        raw=sum(torus[i]for i in pair)-sum(a*b for a,b in zip(torus,ep))
        weights.append(m.reduced['parameter_weights'][j]-raw)
    one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1]
    h14=[i for i in one if i not in [3,11]]
    baseweights={str(i+1):-weights[i]for i in h14};assert max(baseweights.values())==0
    assert {int(i)for i,w in baseweights.items()if w==0}=={5,13,14,49}
    z,b,beta=primitive_beta(m);v=s.Matrix([beta[i,j].diff(z[k])for(i,j),k in COORDS])
    basis=v.jacobian(b);rows=list(basis.T.to_DM().rref()[1]);left=basis[rows,:].inv()
    action=s.Matrix.hstack(*[left*infinitesimal_action(v,a)[rows,:]for a in sp_basis()])
    normal=action[[one.index(3),one.index(11)],:]
    atzero=normal.subs(dict(zip(b,[point[i]for i in one])))
    assert atzero[:,[5,6]]==s.diag(-1,1)
    assert atzero[:,[4,9]]==s.zeros(2) and atzero[:,7]==s.zeros(2,1)
    # Remaining C direction is A6; the canonical u12 coordinate changes nontrivially.
    assert atzero[1,6]==1
    assert r['orbit_minor']=='2'
    nuisance=[24,33,34,35,36,47,50]
    assert [weights[i]for i in nuisance]==[3,3,1,2,2,2,1]
    sample=point.copy()
    for k,i in enumerate(nuisance):sample[i]=k+1
    Dact=action.subs(dict(zip(b,[sample[i]for i in one])))
    units=s.eye(16)[:,[one.index(i)for i in nuisance]];source=units.row_join(Dact)
    cols=list(source.to_DM().rref()[1]);assert len(cols)==16
    det=source[:,cols].det(method='domain-ge');assert det!=0
    paths=[str(Path(__file__).relative_to(ROOT)),NOTE,p0,pd,
        'code/s61_model.py','code/primitive_jet_elliptic_curve_20261101.py','code/weight1_structure.py',
        'research/notes/S61_FLAG_ZERO12_ZERO_NUISANCE_FORMAL_20261102.md',
        'research/notes/S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101.md',
        'research/notes/S61_FLAG_UNIFORM_OPEN_20261101.md',
        'research/notes/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.md',
        'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',
        'reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json',
        'reports/S61_FLAG_UNIFORM_OPEN_20261101.json']
    paths+=list(r['input_sha256'])+list(d['input_sha256'])
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in dict.fromkeys(paths)}
    for old in [r,d]:
        for p,h in old['input_sha256'].items():assert hashes[p]==h,p
    return dict(id=STEM,all_passed=True,field='QQ',slice_hilbert_numerator=[1,2,7,12,18,9,-34,-10],slice_hilbert_denominator='1-z',
        D_cut_tangent_characters=hc,C_cut_tangent_characters=cc,D_cut_hilbert_numerator=[1,2,2],
        actual_special_kernel_HF={'2':5,'3':17,'4':35,'5':44,'6':10},actual_special_kernel_length=111,
        relative_base='completion of u4=u12=0 at u5=u13=1, other beta=0',relative_base_weights=baseweights,
        t_weight=5,kernel_generator_weight_bounds=[2,6],strict_weight_inequality='2+5>6',
        H14_normal_action=[[str(a)for a in row]for row in normal.tolist()],compensation_columns=[5,6],compensation_minor='-1',
        actual_C_ordinary_coordinate_derivative='1 (up to inverse-action sign)',two_orbit_cut_minor='2',
        nuisance_weights=[3,3,1,2,2,2,1],H7_source_sample_beta=[str(sample[i])for i in one],H7_source_minor_columns=cols,H7_source_minor=str(det),
        relative_t_kernel_zero_by_finite_weight_lemma=True,relative_full_H14_t_squared_zero_at_basepoint=True,
        all_H7_relative_section_points_t_squared_zero_by_contraction=True,
        all_general_beta_zero12_flag_points_ambient_t_squared_zero=True,
        all_flags_or_affine_exhaustion_claimed=False,full_relative_polynomial_membership_expanded=False,
        proof_note=NOTE,input_sha256=hashes,seconds=time.monotonic()-start)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=verify()
    with args.output.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds'],'source_minor':r['H7_source_minor']}))
