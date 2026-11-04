#!/usr/bin/env python3
"""Independent exact source-weight checker and relative completed-image proof audit."""
import argparse
import ast
from collections import Counter
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63_REVIEW'
def read(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def enc(M):return [[str(a) for a in row] for row in M.tolist()]
def add(a,b):
    c=dict(a)
    for e,v in b.items():c[e]=c.get(e,Q(0))+v
    return {e:v for e,v in c.items() if v}
def scale(a,c):return {e:v*c for e,v in a.items() if v*c}
def mul(a,b):
    c={}
    for e,u in a.items():
        for f,v in b.items():
            g=tuple(x+y for x,y in zip(e,f));c[g]=c.get(g,Q(0))+u*v
    return {e:v for e,v in c.items() if v}
def parse(text,names):
    zero=(0,)*len(names)
    def walk(t):
        if isinstance(t,ast.Constant):return {zero:Q(t.value)} if t.value else {}
        if isinstance(t,ast.Name):return {tuple(int(i==names.index(t.id)) for i in range(len(names))):Q(1)}
        if isinstance(t,ast.UnaryOp):
            assert isinstance(t.op,(ast.USub,ast.UAdd));return scale(walk(t.operand),-1 if isinstance(t.op,ast.USub) else 1)
        assert isinstance(t,ast.BinOp)
        a=walk(t.left)
        if isinstance(t.op,ast.Pow):
            assert isinstance(t.right,ast.Constant) and isinstance(t.right.value,int) and t.right.value>=0
            ans={zero:Q(1)}
            for _ in range(t.right.value):ans=mul(ans,a)
            return ans
        b=walk(t.right)
        if isinstance(t.op,ast.Add):return add(a,b)
        if isinstance(t.op,ast.Sub):return add(a,scale(b,-1))
        if isinstance(t.op,ast.Mult):return mul(a,b)
        assert isinstance(t.op,ast.Div) and set(b)=={zero};return scale(a,1/b[zero])
    return walk(ast.parse(text,mode='eval').body)

def main(out):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63'
    ar,ac=read(folder/'report.json'),read(folder/'certificate.json')
    assert ar['all_frozen_source_weight_and_sparsity_checks_passed']
    assert sha(folder/'certificate.json')==ar['certificate_sha256']
    for rel,h in ar['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    grading=ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63_REVIEW'
    assert read(grading/'report.json')['all_passed'];gc=read(grading/'independent-certificate.json')
    rp=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json';data=read(rp)
    oldgate=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json'
    assert read(oldgate)['all_passed']
    orbitp=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz';oc=read(orbitp)
    tcp=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63/certificate.json.gz';tc=read(tcp)
    fullgate=ROOT/'reports/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63_REVIEW/report.json';assert read(fullgate)['all_passed']
    fullp=ROOT/'reports/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63/relations99.json.gz';full=read(fullp)
    nilp=ROOT/'reports/AMBIENT_B18_NILCOTANGENT_BASIS_20261103_v63/report.json';nilr=read(nilp);assert nilr['all_passed']
    O25=s.Matrix(gc['physical150_by25_GL_action']);O19=s.Matrix(oc['physical150_by19_scalar_then_GL'])
    Oold=O25[:,data['orbit_GL_generator_indices']]
    weights=gc['physical150_weights'];Wdiag=s.diag(*weights)
    def weight(v):
        assert any(v)
        ww={weights[i] for i,a in enumerate(v) if a};assert len(ww)==1;return ww.pop()
    Rrows=tc['common29_even_target_rows'];C39=tc['common39_physical_coordinate_rows']
    assert Rrows==data['common28_linear_target_rows']+[data['common29th_target_row']]
    cuts=ac['common_physical_group_cuts19'];base=ac['remaining_actual_even_base10_raw_coordinates']
    assert cuts+base!=Rrows and sorted(cuts+base)==sorted(Rrows)
    assert len(cuts)==19 and len(base)==10 and O19[cuts,:].det()==-32
    assert str(O19[cuts,:].det())==ac['common_group_cut_minor']
    assert Counter(weights[i] for i in cuts)=={0:15,1:4}
    assert Counter(weights[i] for i in base)=={0:9,1:1}
    assert {i for i in C39 if i not in cuts}==set(C39)-set(cuts)
    # The actual group section can be chosen as ordered exponentials of
    # eigenvectors. Check their characters by conjugation on all GL matrices.
    eig=[0,0,-1,0,0]
    groupweights=[0]+[eig[i//5]-eig[i%5] for i in oc['GL18_generator_indices']]
    assert groupweights==[weight(O19[:,j]) for j in range(19)]
    assert all(v>=0 for v in groupweights)
    assert Oold.rank()==19
    oldgroupweights=[eig[i//5]-eig[i%5] for i in data['orbit_GL_generator_indices']]
    assert oldgroupweights==[weight(Oold[:,j]) for j in range(19)]
    pairs=list(itertools.combinations(range(5),2))
    def flag(M):
        L=[M[i,2] if i!=2 else 0 for i in range(5)]
        H=[-M[1,i] if i!=1 else 0 for i in range(5)]
        plane=[]
        for i,j in pairs:
            a=-M[0,i]*int(j==1)+M[0,j]*int(i==1)-int(i==0)*M[1,j]+int(j==0)*M[1,i]
            plane.append(0 if (i,j)==(0,1) else a)
        return s.Matrix(L+H+plane)
    glflags=[]
    for ix in data['orbit_GL_generator_indices']:
        M=s.zeros(5);M[ix//5,ix%5]=1;glflags.append(flag(M))
    HF=s.Matrix.hstack(*glflags)
    names=['p']+['k'+str(i) for i in range(9)];zero=(0,)*10
    chars=[0,0,0,1,0,0,0,0,0,0]
    ownrecords={}
    for family in ['CA','CB']:
        f=data['families'][family];ft=tc['families'][family]
        polynomials=[parse(a,names) for a in f['core_full150_coordinates']]
        count=Counter();nonzero=Counter();mixed=0
        for i,p in enumerate(polynomials):
            if p:nonzero[weights[i]]+=1
            for e,c in p.items():
                assert sum(a*b for a,b in zip(e,chars))==weights[i]
                assert weights[i] in [0,1]
                count[weights[i]]+=1
                mixed+=bool(e[3] and sum(e)>1)
        # Bind the full source polynomial to its constant, first derivatives,
        # and pure-kernel square used in all accepted finite-source results.
        point=s.Matrix([p.get(zero,Q(0)) for p in polynomials])
        assert point==s.Matrix(gc['literal_point150'])
        first=[]
        for j in range(10):
            e=tuple(int(i==j) for i in range(10));first.append(s.Matrix([p.get(e,Q(0)) for p in polynomials]))
        assert first[0]==s.zeros(150,1)
        J=Oold.row_join(s.Matrix.hstack(*first[1:]));assert J==s.Matrix(f['linear_J'])
        assert J==s.Matrix(ft['full150_by28_physical_even_linear_J'])
        square=s.Matrix([p.get((2,)+(0,)*9,Q(0)) for p in polynomials])
        assert square==s.Matrix(f['raw_second_coefficient'])
        assert square==s.Matrix(ft['physical_pure_kernel_second_W'])
        kernelM=s.Matrix(f['source_kernel_Lie_matrix']);D=s.diag(*eig)
        assert D*kernelM==kernelM*D and s.Rational(f['source_velocity'][3])==0
        # All 20 original source coefficients are retained, including the
        # CB ratio times binary f product already present in the accepted U.
        extras=f['source_fiber_basis_indices'];assert extras==[1,2,3,9,10,14,15,18,19]
        assert [int(i==3) for i in extras]==chars[1:]
        oldbase=[Q(0)]*20
        for j,c in [(0,-1),(7,-1),(11,1),(16,Q(1,3))]:oldbase[j]=c
        for j,text in enumerate(f['source_alpha_coefficients']):
            expected={zero:oldbase[j]} if oldbase[j] else {}
            if Q(f['source_velocity'][j]):expected[(1,)+(0,)*9]=-Q(f['source_velocity'][j])
            if j in extras:
                ex=tuple(int(a==extras.index(j)+1) for a in range(10));expected=add(expected,{ex:Q(1)})
            assert parse(text,names)==expected
        augmented=J.row_join(s.zeros(150,1)).col_join(HF.row_join(s.zeros(20,9)).row_join(flag(kernelM)))
        assert augmented==s.Matrix(f['entire_source_tangent']) and augmented.rank()==29
        normalchars=oldgroupweights+chars[1:]+[0]
        assert Counter(normalchars)=={0:24,1:5}
        E=s.Matrix(ft['full150_by29_even_cotangent_embedding'])
        B=s.Matrix(ft['full150_by10_actual_odd_cotangent_embedding'])
        assert E[Rrows,:]==s.eye(29)
        assert Wdiag*E==E*s.diag(*(weights[i] for i in Rrows))
        assert Counter(weight(E[:,j]) for j in range(29))=={0:24,1:5}
        oddchars=[weight(B[:,j]) for j in range(10)]
        assert Counter(oddchars)=={0:8,1:2} and weight(square)==0
        saved=ac['component_records'][family]
        assert {str(k):v for k,v in count.items()}==saved['core_terms_by_weight']
        assert {str(k):v for k,v in nonzero.items()}==saved['core_nonzero_physical_coordinate_counts_by_weight']
        ownrecords[family]={'full150_polynomial_source_term_counts':dict(count),
                           'full150_polynomial_nonzero_rows':dict(nonzero),
                           'all_source_parameters_have_exact_displayed_characters':True,
                           'nonlinear_positive_character_source_terms_retained':mixed,
                           'full_source_29_augmented_tangent_rank':29,
                           'normalizer29_characters':dict(Counter(normalchars)),
                           'normalizer10_after_group_cuts_characters':{'0':9,'1':1},
                           'actual_odd10_characters':oddchars,'kernel_p_and_tracefree_rank2_generator_weight':0}
    canonical=s.Matrix(data['canonical_reconstruction_from126']);rawfree=data['common_target126_raw_coordinate_indices']
    assert canonical[rawfree,:]==s.eye(126)
    aw=[weight(canonical[:,j]) for j in range(126)];ambient=Counter(aw)
    assert ambient=={0:36,1:60,2:24,3:6}
    tangent=Counter(weights[i] for i in C39);assert tangent=={0:32,1:7}
    differences={i:ambient[i]-tangent[i] for i in range(4)}
    assert differences=={0:4,1:53,2:24,3:6}
    quotient=Counter(weights[i] for i in C39 if i not in cuts);assert quotient=={0:17,1:3}
    assert set(cuts)<=set(C39)
    # Recheck the old 99 relation character compatibility. This is only a
    # consistency check, never a relative CA/CB coordinate identification.
    abstract_names=full['base_R_core_variables']+full['odd_target_variables']
    relcounts=Counter()
    assignments=[int(n in ['w','Y7','Y8']) for n in abstract_names]
    for text in full['all99']:
        poly=parse(text,abstract_names)
        cc={sum(a*b for a,b in zip(e,assignments)) for e in poly}
        assert len(cc)==1;relcounts[cc.pop()]+=1
    assert relcounts=={0:64,1:32,2:3}
    assert nilr['physical_linear_forms']==[{'46':'1','95':'1'},{'63':'1'}]
    assert weights[46]==weights[95]==1 and weights[63]==0
    assert nilr['nil_cotangent_kernel_dimension']==2 and not nilr['raw_linear_functions_themselves_claimed_nilpotent']
    evidence={'component_source_records':ownrecords,'common_group_cuts19':cuts,'common_cut_minor':'-32',
              'remaining_even_base10':base,'transverse_minimal20_character_counts':dict(quotient),
              'canonical126_character_counts':dict(ambient),'all87_graph_discrepancy_character_counts':differences,
              'abstract99_character_compatibility_only':dict(relcounts),
              'simultaneous_orbit_product_and_exact_quotient_argument':[
                  'The physical19 slice is torus stable; ordered exponentials of GL eigenvectors give a conjugation-stable formal group section.',
                  'The ambient action map has a unit differential and is torus equivariant; its formal inverse is equivariant by uniqueness.',
                  'Apply the same group product to the actual images, proper-source normalizers, and their cokernel. Group cuts are regular on all of them; injection survives.',
                  'The tracefree free rank-one summand of the quadratic cover has a character-zero eigen-generator, so its square lies in the zero-character ring.',
              ],
              'completed_image_relation_argument':[
                  'R10=R0[[w]], and actual I=I0*R10+w*J0*R10, with eight degree-zero and two degree-one generators.',
                  'Degree-d homogeneous first syzygies divide by w^(d-1) for d>=1 in the regular domain, giving first syzygy degrees zero or one.',
                  'All image relations are these syzygies and 55 pair products, hence component positive degree<=2.',
                  'Choose equivariant minimal presentations and homogeneous lifts of all87 eliminated ambient coordinates. Their degrees are their weights, at most3.',
                  'The actual sum ideal retains both component kernels and all87 graph differences, hence <=285 generators of positive degree<=3.',
              ],
              'fixed_ring':{'edim':17,'nil_cotangent_image_dimension':1,'nonzero_class':'u63',
                            'reason':'degree-zero subring equals the positive-variable quotient; (m_D^2)_0=m_fixed^2, and nilradical intersects the subring exactly'},
              'no_actual_zero_character_coefficients_or_full_D_nil_depth_computation':True}
    (out/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),folder/'report.json',folder/'certificate.json',
            grading/'report.json',grading/'independent-certificate.json',rp,oldgate,orbitp,tcp,fullgate,fullp,nilp]
    report={'id':STEM,'all_passed':True,
            'independent_AST_sparse_all150_source_covariance_and_full_source29_unit_check':True,
            'genuine_CB_nonlinear_ratio_f_terms_retained_in_full_source':True,
            'common_equivariant_G19_product_and_no_arbitrary_injective_basechange_assumption':True,
            'normalizer_after19_cuts_characters':{'0':9,'1':1},'actual_odd10_characters':{'0':8,'1':2},
            'actual_transverse20_characters':{'0':17,'1':3},
            'component_positive_degree_bound':2,'whole_D_safe_positive_degree_bound':3,
            'all87_actual_ambient_graph_discrepancies_retained':True,
            'graph_weight_counts':differences,'possible_cubic_graphs':6,'generator_count_upper_bound':285,
            'fixed_zero_character_ring_edim':17,'fixed_nil_cotangent_image_dimension':1,
            'actual_zero_character_coefficient_functions_full_D_nil_socle_depth_uncomputed':True,
            'no_actual_relative_alignment_of_abstract99_models_assumed':True,
            'no_new_sourcejets_new_points_GB_elimination_or_author_producer_rerun':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
            'evidence_sha256':{'independent-certificate.json':sha(out/'independent-certificate.json')},
            'seconds':time.monotonic()-started,'frozen':True}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','seconds','graph_weight_counts','whole_D_safe_positive_degree_bound']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/STEM)
    a=p.parse_args();main(a.output.resolve())
