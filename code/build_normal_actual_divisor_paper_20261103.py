#!/usr/bin/env python3
"""Build the standalone Japanese normal-actual-divisor paper and inspect its PDF.

Requires existing LuaLaTeX/Japanese fonts.  PyMuPDF, when available, adds all-page
geometry/widget/text checks.  This does not run any archived algebra certificate.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def build(folder,destination,report):
    start=time.monotonic();source=ROOT/'papers/normal_actual_divisor_recognition_ja.tex'
    tex=source.read_text();labels=re.findall(r'\\label\{([^}]+)\}',tex)
    refs=re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',tex)
    assert not set(refs)-set(labels)
    assert all(n==1 for n in Counter(labels).values())
    assert not re.search(r'\\(?:input|include)\s*\{',tex)
    citations=set()
    for value in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]+)\}',tex):citations.update(value.split(','))
    assert not citations-set(re.findall(r'\\bibitem\{([^}]+)\}',tex))
    folder=folder.resolve();folder.mkdir(parents=True,exist_ok=False)
    latex=shutil.which('lualatex');assert latex
    for index in range(1,4):
        with (folder/f'pass-{index}.stdout').open('w') as stream:
            subprocess.run([latex,'-interaction=nonstopmode','-halt-on-error',f'-output-directory={folder}',str(source)],
                           cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=180)
    stem=source.stem;log=(folder/f'{stem}.log').read_text()
    critical=[line for line in log.splitlines() if any(value in line for value in
              ['LaTeX Warning:','Package hyperref Warning:','Overfull','Missing character:',
               'Undefined control sequence','multiply defined'])]
    assert not critical,critical
    pdf=folder/f'{stem}.pdf';assert pdf.exists()
    qa={}
    try:
        import fitz
        document=fitz.open(pdf);outside=[];question_pages=[];widgets=0;blank=[]
        for index,page in enumerate(document):
            widgets+=len(list(page.widgets() or []));text=page.get_text()
            if not text.strip():blank.append(index+1)
            if '??' in text:question_pages.append(index+1)
            rect=page.rect
            for block in page.get_text('dict')['blocks']:
                if 'lines' not in block:continue
                for line in block['lines']:
                    for span in line['spans']:
                        x0,y0,x1,y1=span['bbox']
                        if x0 < -1 or y0 < -1 or x1>rect.width+1 or y1>rect.height+1:
                            outside.append(dict(page=index+1,bbox=span['bbox'],text=span['text']))
        assert not outside and not question_pages and not widgets and not blank
        qa=dict(engine='PyMuPDF',pages=len(document),all_pages_text_inside_page=True,
                unresolved_question_mark_pages=question_pages,blank_pages=blank,widgets=widgets,
                bookmarks=len(document.get_toc()),page_size_points=list(document[0].rect))
    except ImportError:
        qa=dict(engine='unavailable',limitation='PyMuPDF all-page geometry inspection not performed')
    destination=destination.resolve();destination.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(pdf,destination)
    proof_reports=[ROOT/'reports'/name for name in [
        'NORMAL_QUINTIC_RECOGNITION_AUDIT_20261103_v62.json',
        'NORMAL_ACTUAL_INACTIVE_LINE_SCOPE_AUDIT_20261103_v62.json',
        'NORMAL_ACTUAL_RNC4_COVER_AUDIT_20261103_v62.json',
        'NORMAL_ACTUAL_E6_E8_GEOMETRY_AUDIT_20261103_v62.json',
        'NORMAL_ACTUAL_EIGHT_SUPPORT_AUDIT_20261103_v62.json']]
    for p in proof_reports:assert json.loads(p.read_text())['all_passed']
    result=dict(id='NORMAL_ACTUAL_DIVISOR_PAPER_20261103',all_passed=True,
                source=str(source.relative_to(ROOT)),pdf=str(destination.relative_to(ROOT)),
                source_sha256=sha(source),pdf_sha256=sha(destination),build_script_sha256=sha(Path(__file__)),
                latex_engine=subprocess.check_output([latex,'--version'],text=True).splitlines()[0],passes=3,
                no_recursive_tex_inputs=True,labels=len(labels),citations=len(citations),
                undefined_references=[],duplicate_labels=[],critical_log_warnings=critical,pdf_qa=qa,
                proof_record_sha256={str(p.relative_to(ROOT)):sha(p) for p in proof_reports},
                scope='Standalone theorem/proof paper; normal ACTUAL divisor closed-point classification and geometric-fibre-normal relative recognition.',
                unchanged_old_algebra_certificates_rerun=False,
                history_independent_model_definitions=['all150quadraticcoefficientordering andRREFconstraints',
                                                      'explicitfreecolumnindices for211triple,22,31',
                                                      'RN36+28candidate coordinatebasisrules',
                                                      'generalfirstjetH,B reconstruction andalllocalizedexceptions'],
                limits=['Historical novelty and journal submission are not claimed.',
                        'No full Poisson component exhaustion or classification of nonnormal/nonreduced actual divisors.',
                        'The relative normal condition concerns geometric fibres; total schemes over nilpotent bases are not asserted normal.'],
                seconds=time.monotonic()-start)
    report.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-directory',type=Path,required=True)
    p.add_argument('--pdf-destination',type=Path,default=ROOT/'papers/normal_actual_divisor_recognition_ja.pdf')
    p.add_argument('--report',type=Path,required=True);a=p.parse_args()
    r=build(a.output_directory,a.pdf_destination,a.report);print('PASS',r['pdf_qa'],r['seconds'])
