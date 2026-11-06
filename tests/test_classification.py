"""Control-plane unit tests only; no claim about Poisson mathematics."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('classification_control', ROOT/'scripts/classification.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class ClassificationTests(unittest.TestCase):
    def setUp(self):
        self.plan=json.loads((ROOT/'research/classification/PLAN.json').read_text())
        self.sources=json.loads((ROOT/self.plan['source_registry']).read_text())
        self.temp=tempfile.TemporaryDirectory(); self.root=Path(self.temp.name)
        (self.root/'evidence.txt').write_text('synthetic evidence, not a mathematical proof\n')
    def tearDown(self): self.temp.cleanup()
    def valid(self): return m.validate(self.plan,self.sources,self.root,check_inputs=False)
    def test_packaged_plan_valid(self): self.assertEqual(len(self.valid()),10)
    def test_first_task(self): self.assertEqual(m.available(self.valid())[0]['id'],'C0')
    def test_suspended_not_selected(self): self.assertNotIn('S',[t['id'] for t in m.available(self.valid())])
    def test_conditional_not_selected(self): self.assertNotIn('B',[t['id'] for t in m.available(self.valid())])
    def test_cycle(self):
        self.plan['tasks'][0]['depends_on']=['G']
        with self.assertRaises(m.PlanError): self.valid()
    def test_unknown_dependency(self):
        self.plan['tasks'][0]['depends_on']=['ZZ']
        with self.assertRaises(m.PlanError): self.valid()
    def test_duplicate_id(self):
        self.plan['tasks'].append(copy.deepcopy(self.plan['tasks'][0]))
        with self.assertRaises(m.PlanError): self.valid()
    def test_false_global_completion(self):
        self.plan['mathematical_status']['global_exhaustion']=True
        with self.assertRaises(m.PlanError): self.valid()
    def test_false_generic_completion(self):
        self.plan['mathematical_status']['generic_five_exhaustion']=True
        with self.assertRaises(m.PlanError): self.valid()
    def test_done_requires_evidence(self):
        self.plan['tasks'][0]['status']='done'
        with self.assertRaises(m.PlanError): self.valid()
    def test_done_with_references_unlocks_g(self):
        t=self.plan['tasks'][0];t.update(status='done',proof='evidence.txt',review='evidence.txt')
        self.assertEqual(m.available(self.valid())[0]['id'],'G')
    def test_invalid_bound(self):
        self.plan['mathematical_status']['global_upper_bound']=57
        with self.assertRaises(m.PlanError): self.valid()
    def test_activation_without_reason(self):
        self.plan['tasks'][-1]['status']='ready'
        with self.assertRaises(m.PlanError): self.valid()
    def test_paths(self):
        for p in ('../x','/tmp/x','x/../y','x\\y','./x',''):
            with self.subTest(p=p),self.assertRaises(m.PlanError): m.relative_path(p)
    def job(self):
        return dict(task_id='C0',target_statement='audit',implication_to_classification='cover all rank4',
          coefficient_ring='QQ',hypotheses='none',inverted_elements=[],complements=[],
          success_certificate='a proof map',independent_check='read the actual maps',
          negative_result_next_step='isolate a missing bridge',
          exact_inputs_and_hashes=[dict(path='evidence.txt',sha256=m.digest(self.root/'evidence.txt'))],
          resource_limit=dict(timeout_seconds=30,memory_mb=256))
    def test_job_valid(self): m.validate_job(self.job(),self.valid(),self.root)
    def test_job_bad_hash(self):
        j=self.job();j['exact_inputs_and_hashes'][0]['sha256']='0'*64
        with self.assertRaises(m.PlanError): m.validate_job(j,self.valid(),self.root)
    def test_job_localization_requires_complement(self):
        j=self.job();j['inverted_elements']=['f']
        with self.assertRaises(m.PlanError): m.validate_job(j,self.valid(),self.root)
    def test_job_suspended(self):
        j=self.job();j['task_id']='S'
        with self.assertRaises(m.PlanError): m.validate_job(j,self.valid(),self.root)
    def test_job_resource_limit(self):
        j=self.job();j['resource_limit']['timeout_seconds']=0
        with self.assertRaises(m.PlanError): m.validate_job(j,self.valid(),self.root)
    def test_replan_gate(self):
        self.plan['tasks'][0]['nonclosing_attempts']=2
        with self.assertRaises(m.PlanError): m.validate_job(self.job(),self.valid(),self.root)
        j=self.job();j['route_change_reason']='Switch from repeated minors to a support argument.'
        m.validate_job(j,self.valid(),self.root)
    def test_strict_missing_inputs(self):
        with self.assertRaises(m.PlanError): m.validate(self.plan,self.sources,self.root,check_inputs=True)
    def test_git_blob_hash(self):
        import hashlib
        data=(self.root/'evidence.txt').read_bytes()
        self.assertEqual(m.digest(self.root/'evidence.txt',git_blob=True),hashlib.sha1(f'blob {len(data)}\0'.encode()+data).hexdigest())

if __name__=='__main__': unittest.main()
