#!/usr/bin/env python3
"""Independent reverse source AD and octic-base dominance at one old point."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

from q_firstjet_spin4_mixed_discriminant_review_20261103_v63 import mat,mm,tp,det,inv,add

ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63_REVIEW'


def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def zero(n,m):return [[Q(0)for j in range(m)]for i in range(n)]
def identity(n):return [[Q(i==j)for j in range(n)]for i in range(n)]
def pivots(A):
    A=[row[:]for row in A];r=0;cols=[]
    for j in range(len(A[0])):
        i=next((i for i in range(r,len(A))if A[i][j]),None)
        if i is None:continue
        A[i],A[r]=A[r],A[i];p=A[r][j];A[r]=[x/p for x in A[r]]
        for i in range(r+1,len(A)):
            c=A[i][j];A[i]=[x-c*y for x,y in zip(A[i],A[r])]
        cols.append(j);r+=1
        if r==len(A):break
    return cols


def main(output):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    output.mkdir(parents=True,exist_ok=False)
    author=ROOT/'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63/report.json'
    certpath=author.parent/'certificate.json.gz';ar=read(author);cert=read(certpath);assert ar['all_passed']
    mappath=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz'
    mapreview=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63_REVIEW.json'
    pairpath=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/certificate.json.gz'
    pairreview=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW/report.json'
    pointpath=ROOT/'reports/Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63/certificate.json.gz'
    pointreview=ROOT/'reports/Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63_REVIEW/report.json'
    spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    spinreview=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63_REVIEW/report.json'
    normalreview=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_NORMAL_FORM_20261103_v63_REVIEW/report.json'
    primitivepath=ROOT/'reports/PRIMITIVE_JET_ELLIPTIC_CURVE_20261101.json'
    origpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    payload,pair,oldpoint,spin,primitive,orig=map(read,[mappath,pairpath,pointpath,spinpath,primitivepath,origpath])
    checked={}
    for path in [author,mapreview,pairreview,pointreview,spinreview,normalreview]:
        data=read(path);assert data['all_passed']
        for name,want in data['input_sha256'].items():assert sha(ROOT/name)==want,name;checked[name]=want
    assert pair['same_point_index']==oldpoint['sample_index']==cert['same_existing_point_index']==0
    point=oldpoint['accepted_source_point']['point'];assert point==cert['same_source_point']
    names=payload['input_names'];nodes=payload['nodes'];assert len(names)==18 and len(nodes)==6079
    values=[];local=[]
    for i,(op,argument)in enumerate(nodes):
        if op=='q':value=Q(argument);partials=[]
        elif op=='v':value=Q(point[names[argument]]);partials=[]
        elif op=='+':
            assert all(0<=j<i for j in argument);value=sum((values[j]for j in argument),Q(0));partials=[(j,Q(1))for j in argument]
        elif op=='*':
            assert all(0<=j<i for j in argument)
            prefix=[Q(1)]
            for j in argument:prefix.append(prefix[-1]*values[j])
            suffix=Q(1);partials=[]
            for k in reversed(range(len(argument))):
                partials.append((argument[k],prefix[k]*suffix));suffix*=values[argument[k]]
            value=prefix[-1]
        elif op=='inv':
            assert 0<=argument<i and values[argument];value=1/values[argument];partials=[(argument,-value*value)]
        else:raise AssertionError(op)
        values.append(value);local.append(partials)
    outputs=payload['output_u1_through_u61'];U=[values[i]for i in outputs]
    assert list(map(str,U))==cert['all61_output_values']
    sourceJac=[]
    #61 independent reverse sweeps, instead of the author's18-source forward AD.
    for target in outputs:
        cotangent={target:Q(1)};row=[Q(0)]*18
        for i in range(len(nodes)-1,-1,-1):
            weight=cotangent.pop(i,Q(0))
            if not weight:continue
            if nodes[i][0]=='v':row[nodes[i][1]]+=weight
            for child,derivative in local[i]:
                value=weight*derivative
                if value:cotangent[child]=cotangent.get(child,Q(0))+value
        sourceJac.append(row)
    assert sourceJac==mat(cert['full61x18_source_Jacobian'])
    vp=values[payload['normalization']['value_pivot_node']];pf=values[payload['normalization']['value_pfaffian_node']]
    assert vp==Q(cert['value_pivot'])==Q(2051,12)and pf==Q(cert['value_pfaffian'])==Q(117727,2592)
    # Full primitive conversion from the independent saved S61 coefficient map.
    assert primitive['all_passed']
    for name,want in primitive['source_sha256'].items():assert sha(ROOT/name)==want,name;checked[name]=want
    us=list(map(s.Symbol,primitive['primitive_jet_variables']));zs=list(map(s.Symbol,primitive['affine_variables']))
    one=[int(str(u)[1:])-1 for u in us];assert one==cert['primitive61_weight1_indices']
    betaU=[]
    for i in range(4):
        for j in range(i+1,4):
            expression=s.sympify(primitive['primitive_beta'][i][j])
            for z in zs:
                coefficient=s.expand(expression).coeff(z)
                row=[Q(str(coefficient.coeff(u)))for u in us]
                assert s.expand(coefficient-sum(s.Rational(q.numerator,q.denominator)*u for q,u in zip(row,us)))==0
                betaU.append(row)
    B=mat(orig['primitive24_by16_basis']);br=pivots(tp(B));assert len(br)==16
    conversion=mm(inv([B[i]for i in br]),[betaU[i]for i in br])
    assert mm(B,conversion)==betaU and conversion==mat(cert['primitive_S61_to_beta16_linear_map'])
    beta=[sum(c*U[i]for c,i in zip(row,one))for row in conversion]
    assert list(map(str,beta))==cert['beta16']==oldpoint['accepted_source_point']['beta16']
    betaJac=mm(conversion,[sourceJac[i]for i in one]);assert betaJac==mat(cert['full_beta16x18_source_Jacobian'])
    C,D=mat(pair['actual_C4']),mat(pair['actual_D4']);At=list(map(mat,pair['pair_coordinate_basis4']))
    pairJac=mat(pair['jacobian20x16']);quadratic=spin['quadratic_monomials']
    ccoef=mat(spin['spin_adjoint_C_coefficients16x136']);dcoef=mat(spin['vector_adjoint_D_coefficients16x136'])
    def matrix_derivatives(coefficients):
        out=[zero(4,4)for k in range(16)]
        for entry,row in enumerate(coefficients):
            assert sum(c*beta[i]*beta[j]for c,(i,j)in zip(row,quadratic))==[x for r in(C if coefficients is ccoef else D)for x in r][entry]
            for c,(i,j)in zip(row,quadratic):
                out[i][entry//4][entry%4]+=c*beta[j]
                out[j][entry//4][entry%4]+=c*beta[i]
        return out
    dC,dD=matrix_derivatives(ccoef),matrix_derivatives(dcoef)
    words=['CCCC','CCCD','CCDD','CDCD','CDDD'];assert words==cert['trace8_basis_words']
    def derivative(word,deltaC,deltaD):
        value=identity(4);jet=zero(4,4)
        for letter in word:
            M=C if letter=='C'else D;delta=deltaC if letter=='C'else deltaD
            jet=add(mm(jet,M),mm(value,delta));value=mm(value,M)
        return sum(jet[i][i]for i in range(4))
    pairGradient=[[derivative(w,E,zero(4,4))for E in At]+[derivative(w,zero(4,4),E)for E in At]for w in words]
    assert pairGradient==mat(cert['trace8_gradient5x20_pair'])
    betaGradient=[[derivative(w,dc,dd)for dc,dd in zip(dC,dD)]for w in words]
    assert betaGradient==mat(cert['trace8_gradient5x16_beta'])==mm(pairGradient,pairJac)
    composed=mm(betaGradient,betaJac);assert composed==mat(cert['trace8_source_derivative5x18'])
    groupcols=[i for i,name in enumerate(names)if name.startswith('g')]
    assert len(groupcols)==10 and all(composed[i][j]==0 for i in range(5)for j in groupcols)
    chosen=cert['selected_source_columns'];assert [names[i]for i in chosen]==cert['selected_source_names']==['n0','n1','n2','scale','base1']
    minor=[[row[j]for j in chosen]for row in composed];determinant=det(minor)
    assert minor==mat(cert['exact5_minor'])and determinant==Q(cert['exact5_minor_determinant'])and determinant
    assert abs(determinant.numerator).bit_length()==461 and determinant.denominator.bit_length()==329
    assert len(pivots(composed))==5
    change=mat([[0,0,'-1/2',0,-1],['-1/2',-1,1,0,2],['1/2',0,-2,'-1/2',-2],['-1/4',0,'1/2',0,0],['1/32',0,'1/16',0,'1/8']])
    assert change==mat(cert['trace8_to_pencil5_linear_change'])and det(change)==Q(-1,128)==Q(cert['linear_change_determinant'])
    # Full octic coefficient identity in the already independently accepted
    # coefficient tables; no evaluation of degree24 polynomials is needed.
    traces=[{m:Q(c)for m,c in p}for p in pair['trace_polynomials'][:5]]
    pencil=[{m:Q(c)for m,c in p}for p in spin['canonical_adjoint_pencil']['determinant_coefficients']]
    keys=set().union(*traces,*pencil)
    for k in range(5):
        assert all(sum(change[k][j]*traces[j].get(m,Q(0))for j in range(5))==256*pencil[k].get(m,Q(0))for m in keys)
    evidence={
        'same_point_index':0,'input_names':names,'source_point':point,
        'AD_method':'61 reverse sweeps with all6079 local derivatives; multiplication handles zero factors',
        'all61_values':[str(q)for q in U],
        'independent_full61x18_source_Jacobian':[[str(q)for q in row]for row in sourceJac],
        'independent_primitive_conversion16x16':[[str(q)for q in row]for row in conversion],
        'independent_beta16x18_source_Jacobian':[[str(q)for q in row]for row in betaJac],
        'independent_trace8_direct_beta_gradient5x16':[[str(q)for q in row]for row in betaGradient],
        'independent_trace8_source_derivative5x18':[[str(q)for q in row]for row in composed],
        'selected_source_columns':chosen,'exact_minor':str(determinant),
        'all10_Sp4_source_columns_zero':True,'all_universal_octic_change_coefficients_checked':True,
        'constant_change_determinant':'-1/128',
        'dominance_and_full_characteristic_base_proof_in_review_note':True,
    }
    evpath=output/'independent-certificate.json.gz';evpath.write_bytes(gzip.compress((json.dumps(evidence,indent=2)+'\n').encode(),mtime=0))
    note=ROOT/'research/notes'/(ID+'.md')
    files=[Path(__file__),note,author,certpath,mappath,mapreview,pairpath,pairreview,pointpath,pointreview,spinpath,spinreview,normalreview,primitivepath,origpath,
           ROOT/'code/q_firstjet_spin4_mixed_discriminant_review_20261103_v63.py',
           ROOT/'research/notes/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63.md',
           ROOT/'code/q_new_family_octic_base_dominance_20261103_v63.py']
    report={
        'id':ID,'all_passed':True,
        'independent6079_node_reverse_AD_all61x18_entries_match':True,
        'full_primitive_beta_binding_and_direct_quadratic_covariant_trace_derivatives_match':True,
        'same_existing_point_only_no_new_beta_points':True,
        'actual_Inv8_source_rank_exact5_by_QQ_minor':True,
        'all10_Sp4_source_columns_exactly_zero':True,
        'all_universal_Inv8_to_pencil_coefficients_and_det_minus1_over128_verified':True,
        'NEW25_dominates_A5_pencil_and_full_irreducible_characteristic_base':True,
        'all_proper_spectral_base_only_divisors_excluded_as_NEW25_equation':True,
        'generically_finite_multisection_scope_accepted':True,
        'multisection_degree_genus_rationality_equation_or_global_finiteness_claim':False,
        'no_author_code_no_degree24_or_higher_polynomial_probe':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},
        'checked_dependency_sha256':checked,
        'evidence_sha256':{evpath.name:sha(evpath)},
        'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'seconds':time.monotonic()-started,
    }
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0);print(json.dumps({'all_passed':True,'source_octic_rank':5,'seconds':report['seconds']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())
