#!/usr/bin/env python3
"""Independent whole-coefficient review of the three fixed Theta descendants."""
import argparse,json,os,resource,signal,subprocess,time
from pathlib import Path
from fractions import Fraction as F
from s61_nonsquare_p_branch_chain_review_20261104_v63 import ROOT,sha,read,dec,add,mul,sc,shift
STEM='S61_NONSQUARE_THETA_FIXED_DEFECTS_20261104_v63_REVIEW'
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    stems=['RAW_H_GENERIC_RANK','THETA_RATIONAL_COORDINATES','RAW_H_THETA_UNIT_PROFILE','RAW_H_THETA_ONE_DEFECT','THETA_DELTA_LAMBDA_PROFILES','THETA_RAW12_DEFECT','THETA_REMAINING_DEFECTS']
    dirs=[ROOT/'reports'/f'S61_NONSQUARE_{s}_20261104_v63' for s in stems];certs=[];inputs=[]
    for folder in dirs:
        report=read(folder/'report.json');certs.append(read(folder/'certificate.json'));assert sha(folder/'certificate.json')==report['certificate_sha256']
        inputs +=[folder/'report.json',folder/'certificate.json']
        for p,h in report['input_sha256'].items():assert sha(ROOT/p)==h,p
    g,c,t,o,d,r,w=certs
    prior=ROOT/'reports/S61_NONSQUARE_RAW_H_P_BRANCH_CHAIN_20261104_v63_REVIEW/report.json';pr=read(prior);assert pr['all_passed']
    for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h,p
    inputs.append(prior)
    B=[[dec(x) for x in row] for row in g['actual_block7x3']];Delta=dec(g['nonzero_selected_minor_Delta']);piv=g['selected_rows7'];rows=t['full_source_rows'];assert piv==[0,5,2] and rows==[7,9,10,11,12,13,14]
    lamb=[[dec(x) for x in row] for row in g['cleared_upper_row_coordinates7xr']]
    source_checks=0
    for i in [1,4,6]:
        for j in range(3):
            rhs={}
            for k,b in enumerate(piv):rhs=add(rhs,mul(lamb[i][k],B[b][j]))
            assert mul(Delta,B[i][j])==rhs;source_checks+=1
    images={x['name']:dec(x['polynomial']) for x in d['actual_images']}
    S={x['name']:dec(x['polynomial']) for x in r['actual_transformed_S']}
    new={x['name']:dec(x['polynomial']) for x in w['new_actual_source_images']};S.update({k:v for k,v in new.items() if k.startswith('S_')})
    substitutions=[('Delta',Delta,images['Delta'])]
    substitutions +=[(f'lambda_{i}_{j}',lamb[i][j],images[f'lambda_{i}_{j}']) for i in [1,4,6] for j in range(3)]
    substitutions +=[(f'S_raw{raw}',dec(t['fixed_7x4_matrix'][rows.index(raw)][3]),S[f'S_raw{raw}']) for raw in [7,9,10,12,13,14]]
    substitutions +=[('old_V',B[0][0],new['old_V'])]
    # Clear only powers of d before generic library composition, unlike the
    # producer's per-term cached power/summation algorithm.
    one={(0,0):F(1)};D={(1,0):F(1)}
    E=shift(dec(c['target_coordinates']['e']),(1,0));A=shift(dec(c['target_coordinates']['a']),(2,0))
    text=[]
    def emit(p):
        text.append(str(len(p)))
        text.extend(str(v)+' '+' '.join(map(str,m)) for m,v in sorted(p.items()))
    for p in [D,E,A]:emit(p)
    text.append(str(len(substitutions)));clearings={}
    for name,p,expected in substitutions:
        n=max([2*a+e-d for d,e,a in p]+[0]);clearings[name]=n
        q={(d+n-e-2*a,e,a):v for (d,e,a),v in p.items()};want=shift(expected,(n,0));assert all(min(m)>=0 for m in q) and all(min(m)>=0 for m in want)
        text.append(name);emit(q);emit(want)
    ip=out/'independent-compose-input.txt';ip.write_text('\n'.join(text)+'\n')
    cpp=Path(__file__).with_suffix('.cpp');exe=ROOT/'work/Theta-fixed-defects-review-compose'
    prefix=Path('/workspace/cloud-setup/prefix/usr');lib=prefix/'lib/x86_64-linux-gnu';env=os.environ.copy();env['LD_LIBRARY_PATH']=str(lib)+':'+env.get('LD_LIBRARY_PATH','')
    cmd=['g++','-O2','-std=c++17',str(cpp),'-I'+str(prefix/'include'),'-L'+str(lib),'-lflint','-lmpfr','-lgmp','-o',str(exe)]
    subprocess.run(cmd,check=True,timeout=20,env=env);np=out/'independent-compose-checks.json';subprocess.run([str(exe),str(ip),str(np)],check=True,timeout=30,env=env)
    assert read(np)['all_passed'];delta=images['Delta'];assert len(delta)==56 and new['old_V']
    assert dec(r['transformed_Delta'])==delta and dec(r['transformed_lambda0'])==images['lambda_4_0'] and dec(r['transformed_lambda2'])==delta
    assert not images['lambda_4_1']
    p12=[mul(delta,S['S_raw12']),mul(images['lambda_4_0'],S['S_raw7']),mul(delta,S['S_raw10'])]
    assert p12==[dec(p) for p in r['three_product_polynomials']]
    assert not add(add(p12[0],p12[1],F(-1)),p12[2],F(-1)) and not dec(r['defect'])
    defects={};stats=[]
    for i,raw in [(1,9),(6,14)]:
        prods=[mul(delta,S[f'S_raw{raw}'])]+[mul(images[f'lambda_{i}_{j}'],S[f'S_raw{rows[piv[j]]}']) for j in range(3)]
        defect=prods[0]
        for p in prods[1:]:defect=add(defect,p,F(-1))
        rec=next(x for x in w['two_fixed_defects'] if x['raw_row']==raw)
        assert prods==[dec(p) for p in rec['transformed_product_polynomials']] and defect==dec(rec['defect']) and (not defect)==rec['defect_zero']
        defects[raw]=defect;stats.append({'raw_row':raw,'product_terms':list(map(len,prods)),'defect_terms':len(defect)})
    assert not defects[9] and len(defects[14])==269 and not w['all_four_nonpivot_defects_zero'] and w['sharp_consequence_if_allzero'] is None
    assert not o['defect_remainder'] and o['defect_identically_zero']
    assert dec(w['old_V_actual_image'])==new['old_V'] and len(new['old_V'])==8
    # Accepted original polynomial V*S11-U*S7 and Delta*U-lambda*V
    # identities plus the domain and nonzero V imply the fourth coordinate
    # defect is zero without declaring V a global unit.
    evidence={'exact_source_B_Cramer_identities_rebuilt':source_checks,'generic_composition_all_coefficients_checked':len(substitutions),
      'composition_context':'source(d,e,a),target(d,z),DEGREVLEX; generic compose after known-d clearing',
      'known_d_clearing_exponents':clearings,'Delta_nonzero_terms':len(delta),'old_V_nonzero_terms':len(new['old_V']),
      'raw12_products_whole_QQ_checked':True,'remaining_defects':stats,
      'raw11_zero_follows_by_domain_cancellation_only':True,
      'generic_Theta_augmented_rank':4,'generic_Theta_heavy_rank':5,'generic_Theta_original_raw_H_rank':7,
      'pointwise_rank5_heavy_on_Brank3_and_D14_nonzero':True,
      'all_Theta_rank_drop_or_Delta_D14_unit_or_full_cubic_conclusion':False}
    ep=out/'independent-fixed-defects.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    inputs +=[Path(__file__),cpp,ROOT/'research/notes'/f'{STEM}.md']
    report={'id':STEM,'all_passed':True,'independent_all17_source_compositions_and_all_products_QQ':True,
      'all_three_descendants_reviewed':True,'raw14_nonzero_not_all_Theta_rank_drop':True,
      'generic_heavy5_rawH7_accepted_only_on_actual_open':True,'no_point_prime_ranksearch_GB_or_unknown_global_inverse':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'evidence_sha256':{p.name:sha(p) for p in [ip,np,ep]},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['all_passed','seconds','max_rss_kib']}))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)
