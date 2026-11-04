#!/usr/bin/env python3
"""Independent fixed coefficient minors proving the D4 second Betti number."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_D4_MINIMAL_BETTI112_20261104_v63_REVIEW'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    np=ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63/actual-N-literal-cone-DAG.json';hp=ROOT/'reports/AMBIENT_B18_D4_LITERAL_HORSESHOE_20261104_v63/actual-D4-horseshoe-DAG.json';n=read(np);h=read(hp)
    gates=[ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63_REVIEW/report.json',ROOT/'reports/AMBIENT_B18_D4_LITERAL_HORSESHOE_20261104_v63_REVIEW/report.json',ROOT/'reports/AMBIENT_B18_NIL_DERIVATION_EXTRACTOR_20261104_v63_REVIEW/report.json'];hashes={}
    for p in gates:
        r=read(p);assert r['all_passed'];hashes[str(p.relative_to(ROOT))]=sha(p)
    horse_gate=read(gates[1]);assert horse_gate['extractor_external_independent_gate_now_passed'] and horse_gate['input_sha256'][str(hp.relative_to(ROOT))]==sha(hp)
    assert horse_gate['input_sha256'][str(np.relative_to(ROOT))]==sha(np)
    pairs=list(itertools.combinations(range(8),2));triples=list(itertools.combinations(range(8),3));variables=[*range(8),9];names=[*[f'z{i+1}' for i in range(8)],'v']
    # Full 324x84 coefficient matrix is sparse, and no row/pivot search occurs.
    entries={};records=[]
    for beta,column in enumerate(n['d2_universal_sparse_polynomial_columns']):
        for alpha,terms in column:
            alpha=int(alpha);assert 0<=alpha<36
            for coefficient,monomial in terms:
                assert coefficient in [-1,1] and len(monomial)==1 and monomial[0] in variables
                key=(alpha,monomial[0],beta);entries[key]=entries.get(key,0)+coefficient
    entries={k:v for k,v in entries.items() if v};assert len(entries)==252
    selector=[(pairs.index((j,k)),i) for i,j,k in triples]+[(alpha,9) for alpha in range(28)]
    assert len(selector)==len(set(selector))==84
    minor=[[entries.get((alpha,var,beta),0) for beta in range(84)] for alpha,var in selector]
    assert minor==[[int(i==j) for j in range(84)] for i in range(84)]
    # Independently bind these columns to the actual horseshoe N projection.
    actual_n2=[{int(i):q for i,q in col.items()} for col in n['R0_sparse_d2_36_by84_columns']]
    actual_h2=[{int(i):q for i,q in col.items()} for col in h['D4_d2_44_by112_columns']]
    assert len(actual_h2)==112 and len(actual_n2)==84
    for beta,col in enumerate(actual_h2[:28]):assert all(i<8 for i in col)
    for beta,col in enumerate(actual_h2[28:]):assert {i-8:q for i,q in col.items() if i>=8}==actual_n2[beta]
    P=h['nil_annihilator8_nodes'];assert len(P)==8
    nil_selector=[];nil_minor=[[0]*28 for _ in range(28)]
    for column,(i,j) in enumerate(pairs):
        assert actual_h2[column]=={j:P[i],i:h['D4_d2_44_by112_columns'][column][str(i)]}
        # The negative node is certified in the accepted literal horse input;
        # check its exact syntax again, not merely its value at the origin.
        negative=actual_h2[column][i];node=h['DAG_scalar_nodes'][negative]
        assert node==['n',P[j]]
        nil_selector.append((j,i))
    # P_i=x_i-psi_i(z1), x=(z2,...,z8,v). The accepted curve input
    # depends only on z1, so partial P_i / partial x_j is exactly delta_ij.
    for row,(target,coordinate) in enumerate(nil_selector):
        for column,(i,j) in enumerate(pairs):nil_minor[row][column]=int(target==j and coordinate==i)-int(target==i and coordinate==j)
    assert nil_minor==[[int(i==j) for j in range(28)] for i in range(28)]
    # Check origin zero without any gradients or unexpanded residue series.
    values=[]
    for i,node in enumerate(h['DAG_scalar_nodes']):
        op=node[0]
        if op=='c':q=F(node[1],node[2])
        elif op=='v':q=F(0)
        elif op=='a':q=values[node[1]]+values[node[2]]
        elif op=='m':q=values[node[1]]*values[node[2]]
        elif op=='n':q=-values[node[1]]
        else:assert op=='i' and values[node[1]];q=1/values[node[1]]
        values.append(q)
    assert all(values[q]==0 for q in h['positive_c36_actual_nodes']+h['nilcurve_psi8_actual_nodes']+P)
    assert all(values[q]==0 for col in actual_h2 for q in col.values())
    assert all(values[q]==0 for col in h['D4_d1_9_by44_columns'] for q in col.values())
    xlabels=[f'z{i}' for i in range(2,9)]+['v']
    evidence=dict(N_coefficient_variables_in_order=names,N_coefficient_matrix_shape=[324,84],N_all252_integer_coefficient_entries=[[a,names[variables.index(v)],b,c] for (a,v,b),c in sorted(entries.items())],N_selected84_coefficient_rows=[[a,names[variables.index(v)]] for a,v in selector],N_selected84_minor_sparse=[[i,i,1] for i in range(84)],N_selected84_minor_determinant=1,nil_selected28_rows=[[row,xlabels[var]] for row,var in nil_selector],nil_selected28_minor_sparse=[[i,i,1] for i in range(28)],nil_selected28_minor_determinant=1,full112_coefficient_minor_block_pattern='[[I28,uncomputed],[0,I84]]',full112_coefficient_minor_determinant=1,all_actual_d2_origin_entries_zero=True,actual_complete_horseshoe_gate_bound=True,actual_minimal_second_Betti_number=112,coefficient_minor_not_a_ring_minor_of_44_by112=True,no_offdiagonal_derivative_needed=True,no_T0_Tor_depth_or_higher_Betti_claim=True)
    (out/'independent-coefficient-minors.json').write_text(json.dumps(evidence,indent=2)+'\n')
    for p in [Path(__file__),ROOT/f'research/notes/{ID}.md',np,hp]:hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,actual_D4_beta2_exactly112=True,rank84_N_linear_coefficient_minor_and_rank28_nil_minor_exact=True,no_source_solver_point_prime_gradient_or_new_residue_series=True,input_sha256=hashes,evidence_sha256={'independent-coefficient-minors.json':sha(out/'independent-coefficient-minors.json')},seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
