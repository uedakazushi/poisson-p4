#!/usr/bin/env python3
"""Read-only arithmetic-method inventory; no coefficient DAG evaluation or quartic retry."""
import argparse,collections,hashlib,importlib.util,inspect,json,subprocess,time
from pathlib import Path
import sympy
from sympy.polys.domains import QQ
from sympy.polys.fields import FracElement
from sympy.polys.rings import PolyElement
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_H1_NATIVE_ARITHMETIC_PREFLIGHT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 ip=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json';x=json.loads(ip.read_text());a=x['scalar_circuit']['atoms_num_den'];nodes=x['scalar_circuit']['nodes']
 den={json.dumps(d,separators=(',',':')):d for n,d in a}
 funcs={'FracElement.__mul__':FracElement.__mul__,'FracElement.__add__':FracElement.__add__,'FracElement.new':FracElement.new,'PolyElement.cancel':PolyElement.cancel}
 snippets={name:inspect.getsource(f)for name,f in funcs.items()}
 # This reads configured tool metadata only; no install, network, polynomial solve or benchmark.
 flint=subprocess.check_output(['bash','-lc','source /workspace/cloud-setup/activate.sh\npkg-config --modversion flint'],text=True).strip()
 record={'sympy_version':sympy.__version__,'QQ_ground_dtype':str(QQ.dtype),'python_flint_module_available':importlib.util.find_spec('flint')is not None,'existing_activated_native_FLINT_version':flint,
  'scalar_atoms':len(a),'scalar_nodes':len(nodes),'node_operations':dict(collections.Counter(n[0]for n in nodes)),
  'unique_atomic_denominators':len(den),'nonmonomial_atomic_denominators':sum(len(v)>1 for v in den.values()),'max_atomic_denominator_terms':max(map(len,den.values())),
  'max_atomic_denominator_totaldegree':max(sum(m)for v in den.values()for m,c in v),'max_atomic_denominator_rational_bits':max(max(abs(int(z)).bit_length()for z in c.split('/'))for v in den.values()for m,c in v),
  'static_callchain':'FracElement arithmetic -> new(numer,denom) -> PolyElement.cancel -> QQ denominator clearing and multivariate polynomial cofactors over ZZ',
  'runtime_hotspot_profile_NOT_measured':True,'static_method_sources':snippets,
  'next_finite_gate_NOT_executed':'Exact factor all228 distinct atomic denominators by a,d,F,B,C,J and record literal reconstruction; if all knownunits, use native polynomial numerator plus knownunit exponent-vector instead of generic multivariate GCD at each node.',
  'no_DAG_coefficient_evaluation_quartic_retry_GB_point_prime_or_install':True}
 (out/'certificate.json').write_text(json.dumps(record,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),ip,ROOT/'reports/S61_NONSQUARE_H1_KNOWN_UNIT_CONTINUATION_20261104_v63/report.json',ROOT/'reports/S61_NONSQUARE_H1_CACHE_UNITS_20261104_v63/report.json']
 result={'id':ID,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'all_static_counts_and_library_method_sources_recorded':True,'no_runtime_profile_or_new_coefficient_claim':True,'seconds':time.monotonic()-start,'frozen':True}
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
