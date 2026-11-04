#!/usr/bin/env python3
"""One typed/pinned10x17+seven pure-kernel charts; NO polynomial expansion."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PURE_KERNEL9_CHARTS_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parsepoly(s):
    z=s.split();n=int(z[0]);assert len(z)==1+4*n
    return [[list(map(int,z[2+4*k:5+4*k])),str(Q(z[1+4*k]))]for k in range(n)]
class SLP:
    def __init__(self):self.nodes=[];self.cache={};self.polys={};self.explicit_terms=0;self.const(0);self.const(1)
    def put(self,key,value=None):
        if key not in self.cache:
            self.cache[key]=len(self.nodes);self.nodes.append(list(key)if value is None else value)
        return self.cache[key]
    def const(self,c):return self.put(('QQ',str(Q(c))))
    def neg(self,a):
        if not a:return 0
        n=self.nodes[a]
        if n[0]=='QQ':return self.const(-Q(n[1]))
        if n[0]=='neg':return n[1]
        return self.put(('neg',a))
    def add(self,a,b):
        if not a:return b
        if not b:return a
        if self.neg(a)==b or self.neg(b)==a:return 0
        if self.nodes[a][0]==self.nodes[b][0]=='QQ':return self.const(Q(self.nodes[a][1])+Q(self.nodes[b][1]))
        return self.put(('add',min(a,b),max(a,b)))
    def sub(self,a,b):return 0 if a==b else self.add(a,self.neg(b))
    def mul(self,a,b):
        if not a or not b:return 0
        if a==1:return b
        if b==1:return a
        aa,bb=self.nodes[a],self.nodes[b]
        if aa[0]==bb[0]=='QQ':return self.const(Q(aa[1])*Q(bb[1]))
        if aa==['inv',b]or bb==['inv',a]:return 1
        if aa[0]=='neg':return self.neg(self.mul(aa[1],b))
        if bb[0]=='neg':return self.neg(self.mul(a,bb[1]))
        return self.put(('mul',min(a,b),max(a,b)))
    def inv(self,a):
        assert a
        n=self.nodes[a]
        if n[0]=='QQ':return self.const(1/Q(n[1]))
        if n[0]=='inv':return n[1]
        return self.put(('inv',a))
    def poly(self,p,ref):
        if not p:return 0
        p=tuple(sorted((tuple(m),str(Q(c)))for m,c in p if Q(c)))
        if len(p)==1 and p[0][0]==(0,0,0):return self.const(p[0][1])
        key=hashlib.sha256(json.dumps(p,separators=(',',':')).encode()).hexdigest()
        if key not in self.polys:
            self.explicit_terms+=len(p);assert self.explicit_terms<=150000
            self.polys[key]=p
        return self.put(('parameter_polynomial',key),['parameter_polynomial',dict(reference=ref,canonical_SHA256=key,terms=len(p),type='QQ[a,e,d]')])
    def classify(self,p,known):
        if not p:return 'ZERO_INVALID'
        p={tuple(m):Q(v)for m,v in p}
        if len(p)==1 and next(iter(p))==(0,0,0):return 'nonzero_QQ_unit'
        if len(p)==1 and next(iter(p))[1]==0:return 'known_FULLU_a_d_monomial_unit'
        mi=min(m[0]for m in p);mk=min(m[2]for m in p);r={(m[0]-mi,m[1],m[2]-mk):v for m,v in p.items()}
        def projective(q):
            z=sorted(q.items());c=z[0][1];return tuple((m,str(v/c))for m,v in z)
        key=projective(r)
        return 'known_FULLU_literal_'+known[key] if key in known else 'UNKNOWN_chart_factor_NOT_proved_nonunit'
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    lb=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63';lf=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63'
    ab=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63';ub=ROOT/'reports/S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63';hb=ROOT/'reports/S61_NONSQUARE_CUBIC_PROJECTED15_HEAVY_PREFLIGHT_20261104_v63'
    a=read(ab/'input.json');unit=read(ub/'certificate.json');heavy=read(hb/'certificate.json');lc=read(lb/'certificate.json')
    for st in ['S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63_REVIEW','S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63_REVIEW']:
        assert read(ROOT/'reports'/st/'report.json')['all_passed']
    o=SLP();known={}
    def know(p,name):
        d=sorted((tuple(m),Q(v))for m,v in p);first=d[0][1];known[tuple((m,str(v/first))for m,v in d)]=name
    know([[[1,0,0],'500'],[[0,2,0],'90'],[[0,1,0],'3'],[[0,0,0],'-24']],'F')
    know([[[1,0,0],'100'],[[0,2,0],'18'],[[0,1,0],'-9'],[[0,0,0],'-2']],'C')
    know([[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,1,0],'60'],[[0,0,1],'-2'],[[0,0,0],'-40']],'B')
    know([[[1,0,1],'5'],[[0,1,0],'-6'],[[0,0,0],'-4']],'J')
    for i,p in enumerate(a['genuine16_source_clearing_units']):know(p,'accepted_source_Cpower')
    know(unit['known_unit_eta'],'accepted_eta')
    ref=lambda path,**r:dict(path=str(path.relative_to(ROOT)),**r)
    mu=[o.poly(p,ref(ab/'input.json',key='genuine16_source_clearing_units',index=i))for i,p in enumerate(a['genuine16_source_clearing_units'])]
    columns=[[0,0,0]for _ in range(26)]
    for i in range(16):
        psi=o.poly(a['psi16_coefficients'][i],ref(ab/'input.json',key='psi16_coefficients',index=i))
        zw=[o.mul(mu[i],o.poly(a['raw16_HZ_HW_coefficients'][i][j],ref(ab/'input.json',key='raw16_HZ_HW_coefficients',index=[i,j])))for j in range(2)]
        columns[10+i]=[psi,*zw]
    assert columns[11][0]==columns[16][0]
    assert [i for i in range(16)if columns[10+i][0]==0]==[0,3,7]
    rawlines=(lf/'genuine16x31.txt').read_text().splitlines();assert rawlines[1+31*6:1+31*7]==rawlines[1+31:1+31*2]
    phi=[[int(i==j)for j in range(26)]for i in range(26)]
    prefix=[];inverse=[];localizers=[];counts={'ELIMINATE':0,'NORMALIZE':0};paramterms=0;phase='initial_normalizations';prefixbytes=None
    def loc(p,node,record,phase):
        nonconstant=not(len(p)==1 and p[0][0]==[0,0,0]);
        if nonconstant:localizers.append(dict(node=node,classification=o.classify(p,known),introduced_phase=phase,**record))
    with(lb/'operations.txt').open()as f:
        while True:
            offset=f.tell();h=f.readline().split()
            if not h:break
            if h[0]=='ELIMINATE'and int(h[2])==14 and int(h[3])==25:prefixbytes=offset;break
            if h[0]=='NORMALIZE':
                t=int(h[1]);pos=f.tell();p=parsepoly(f.readline());c=Q(f.readline());assert c
                g=o.poly(p,ref(lb/'operations.txt',byte_offset=pos,format='one_polynomial_record'));gc=o.mul(g,o.const(c));q=o.inv(gc)
                loc(p,g,dict(kind='NORMALIZE_g',row=t,operation=len(prefix),byte_offset=pos),phase)
                columns[t]=[o.mul(q,x)for x in columns[t]];phi[t]=[o.mul(q,x)for x in phi[t]]
                prefix.append(dict(op='NORMALIZE',row=t,g=g,QQ_scalar=str(c),inverse_scale=q,phase=phase));inverse.append(dict(op='multiply_row',row=t,scale=gc))
                counts['NORMALIZE']+=1;paramterms+=len(p)
            else:
                assert h[0]=='ELIMINATE';t,p,c=map(int,h[1:]);phase='quad_echelon'if t<10 else 'genuine_reduce_quad10'if p<10 else 'heavy5'
                posu=f.tell();pu=parsepoly(f.readline());posv=f.tell();pv=parsepoly(f.readline());assert pu
                u=o.poly(pu,ref(lb/'operations.txt',byte_offset=posu,format='one_polynomial_record'));v=o.poly(pv,ref(lb/'operations.txt',byte_offset=posv,format='one_polynomial_record'))
                loc(pu,u,dict(kind='ELIMINATE_u',row=t,pivot_row=p,pivot_column=c,operation=len(prefix),byte_offset=posu),phase)
                columns[t]=[o.sub(o.mul(u,x),o.mul(v,y))for x,y in zip(columns[t],columns[p])]
                phi[t]=[o.sub(o.mul(u,x),o.mul(v,y))for x,y in zip(phi[t],phi[p])]
                prefix.append(dict(op='ELIMINATE',row=t,pivot_row=p,pivot_column=c,u=u,v=v,phase=phase));inverse.append(dict(op='inverse_eliminate',row=t,pivot_row=p,coefficient_own=o.inv(u),coefficient_pivot=o.mul(o.inv(u),v)))
                counts['ELIMINATE']+=1;paramterms+=len(pu)+len(pv)
            assert len(o.nodes)<150000
    assert counts=={'ELIMINATE':206,'NORMALIZE':231}and paramterms==30417 and len(localizers)==136 and prefixbytes==733511
    assert columns[16][0]==0 and all(phi[i][j]==0 for i in range(10)for j in range(10,26))
    remaining=lc['remaining11_genuine_relative_rows'];keep=[i for i in range(11)if i!=4]
    assert remaining[4]==6;relative=[remaining[i]for i in keep];rids=[10+i for i in relative]
    R=[];pureterms=0;rowprofiles=[]
    with(lb/'completed-R11x16.txt').open()as f:
        assert f.readline().split()==['11','16']
        for i in range(11):
            row=[];nz=nt=0
            for j in range(16):
                p=parsepoly(f.readline());nt+=len(p);nz+=bool(p);pureterms+=len(p)
                n=o.poly(p,ref(lb/'completed-R11x16.txt',row11=i,column16=j,format='matrix_polynomial_entry'));row.append(n)
            rowprofiles.append(dict(row11=i,genuine_relative=remaining[i],nonzero=nz,terms=nt))
            if i==4:assert row==[0]*16
            else:R.append(row)
    assert pureterms==95376 and len(R)==10
    stage=read(lb/'current-stage.json');pivots=stage['pivots'][:15];assert set(r for r,c in pivots[:10])==set(range(10))and set(r-10 for r,c in pivots[10:])==set([1,2,8,10,11])
    pivotset=set(map(tuple,pivots));pivatoms=[]
    with(lb/'current-matrix.txt').open()as f:
        assert f.readline().split()==['26','31']
        for i in range(26):
            for j in range(31):
                line=f.readline()
                if(i,j)in pivotset:
                    p=parsepoly(line);assert p;n=o.poly(p,ref(lb/'current-matrix.txt',row26=i,column31=j,format='matrix_polynomial_entry'))
                    pivatoms.append(dict(row=i,column=j,node=n,classification=o.classify(p,known),phase='quad10'if i<10 else'heavy5'))
    assert len(pivatoms)==15
    purecols=[j for j in range(31)if j not in [c for r,c in pivots]];assert len(purecols)==16 and max(purecols)<26
    psi=[columns[r][0]for r in rids];zz=[columns[r][1]for r in rids];ww=[columns[r][2]for r in rids]
    Z=[i for i,x in enumerate(psi)if not x];P=[i for i,x in enumerate(psi)if x];assert Z==[0,1,4]and len(P)==7
    cz=o.poly(unit['Bezout_cZ'],ref(ub/'certificate.json',key='Bezout_cZ'));cw=o.poly(unit['Bezout_cW'],ref(ub/'certificate.json',key='Bezout_cW'))
    eta=o.poly(unit['known_unit_eta'],ref(ub/'certificate.json',key='known_unit_eta'));ie=o.inv(eta)
    zd=o.poly(heavy['delta_scaled_HZ_HW_coefficients']['Z'],ref(hb/'certificate.json',key=['delta_scaled_HZ_HW_coefficients','Z']))
    wd=o.poly(heavy['delta_scaled_HZ_HW_coefficients']['W'],ref(hb/'certificate.json',key=['delta_scaled_HZ_HW_coefficients','W']))
    charts=[]
    for q in P:
        vectors=[];images=[];hzproof=[]
        order=[('pure',i)for i in Z]+[('pair',i)for i in P if i!=q]
        for kind,j in order:
            if kind=='pure':
                v=phi[rids[j]][:];z=zz[j];w=ww[j];image=R[j][:];psizero=psi[j];basis=[int(k==j)for k in range(10)]
            else:
                v=[o.sub(o.mul(psi[q],x),o.mul(psi[j],y))for x,y in zip(phi[rids[j]],phi[rids[q]])]
                z=o.sub(o.mul(psi[q],zz[j]),o.mul(psi[j],zz[q]));w=o.sub(o.mul(psi[q],ww[j]),o.mul(psi[j],ww[q]))
                image=[o.sub(o.mul(psi[q],x),o.mul(psi[j],y))for x,y in zip(R[j],R[q])]
                psizero=o.sub(o.mul(psi[q],psi[j]),o.mul(psi[j],psi[q]));basis=[0]*10;basis[j]=psi[q];basis[q]=o.neg(psi[j])
            assert psizero==0
            correction=o.mul(o.add(o.mul(cz,z),o.mul(cw,w)),ie)
            v[16]=o.sub(v[16],correction);v[11]=o.add(v[11],correction)
            vectors.append(dict(kind=kind,remaining10_source_index=j,source_basis_coefficients10=basis,
                initial26_coefficients=v,original57_genuine_coefficients=[o.mul(v[10+i],mu[i])for i in range(16)],
                before_Cdelta_HZ=z,before_Cdelta_HW=w,Cdelta_correction=correction))
            images.append(image);hzproof.append(dict(psi_zero_root=psizero,
                exact_raw7_Hzero_reason='oldfiveHzero byacceptedprefix; subtract(cZ*HZ+cW*HW)/eta*(G6-G1). Universal2x2identity gives rawHZHWzero.'))
        assert len(vectors)==9
        selected=[j for kind,j in order]
        for i,(kind,j)in enumerate(order):
            expected=1 if kind=='pure'else psi[q]
            assert [vectors[i]['source_basis_coefficients10'][k]for k in selected]==[expected if k==i else 0 for k in range(9)]
        charts.append(dict(pivot_remaining10_index=q,pivot_genuine_relative=relative[q],localizer=psi[q],inverse_localizer=o.inv(psi[q]),
            source9=vectors,pure9x16_projected_roots=images,kernel_Hzero_proofs=hzproof,
            source_basis_minor='psi_pivot^6 (three1diagonal entries andsixpsi_pivotentries afterdeletingpivotcolumn)',
            exact_kernel_basis_on_this_chart=True,chart_nonempty_NOT_proved_by_AST=True))
    assert o.explicit_terms<=150000 and len(o.nodes)<150000
    classification={}
    for x in localizers+pivatoms:classification[x['classification']]=classification.get(x['classification'],0)+1
    fieldscope=dict(field='QQ(a,e,d), and specified localizationof FULLU only',
        all136_old_nonconstant_inverse_occurrences_retained=True,all15_LITERAL_pivot_localizers_retained=True,
        no_unknown_factor_declared_global_unit=True,classification_is_conservative='UNKNOWN means not recognized by literalnamed-unit tests; it does NOT prove the factor is a true bad divisor.',
        chart='Omega=FULLU intersect D(product136u/g atoms * product15pivot atoms); then intersect H7 and the selected psi chart.',
        complement='Omega complement may meet H7+Bgood. It is NOT claimed contained in knownBbad13 orHdropcurves.',
        eta_known_unit='192*a*F*B^2*C^3/125',eta_inverse_root=ie,
        source_primitive_Cdelta_line='Initialnormalized G6-G1 has rawheavy(0five,zdelta,wdelta), pi=psi=0. cZ*zdelta+cW*wdelta=eta, acceptedwholeFULLU identity.',
        full_cubic_exact_criterion='Withproduct41 andrawH7, fullsource57 iff actualpurekernel9->pure43 is injective. Our9x16 is itsprojection, so injectivity is SUFFICIENT ONLY.',
        exact_projected_rank='rank(aug25x32)=15+rank(10x17)=16+rank(kernel9x16) onOmega,H7,psi chart; thusrank(old17x24)=8+rank(kernel9x16).',
        full57_scope='kernel9x16rank9 =>aug25rank25=>fullcubic57. FailuredoesNOTprovefullcubicrankloss orquarticfailure.',
        canonical_rank10='Old10minor firstgavegenericLOWER10; exactLframezeroG6row givesUPPER10 forcanonical11 ONLY afterQ10+H5Schur. No upperboundforfixedpureprojectionofallG16.',
        no_new_rank_minor_point_prime_GB_or_expanded_parameter_products=True)
    proof=dict(source_operator='Orderedprefix normalizations and ELIM row_t=u row_t-v row_p; full26x26Phi roots saved. inverseops inreverseorder saved.',
        universal_inverse_templates=[{'op':'NORMALIZE','forward':'row/(g*c)','inverse':'row*(g*c)'},{'op':'ELIMINATE','forward':'u*x-v*y','inverse':'(new_x+v*y)/u'}],
        appended_columns='Exactlytheacceptedpsi16 andsame-scaledrawHZHW; sameparameter-linear rowops preserve psi=-wdelta*HZ+zdelta*HW.',
        duplicate_zero='InputG6=G1 inpiandpsi; sourceunitnormalizationsmatch; literalELIM16,11,26 hasu=v=297/21250, so appendedpsi16rootzero.',
        remaining_three_zero_psis='G0,G3,G7 initiallyhaveALL7heavyentrieszero; QsourcehaszeroH, hencecolumnsremainzero.',
        Cdelta_universal_identities=['eta*hZ-(cZ*hZ+cW*hW)*zdelta=-cW*psi','eta*hW-(cZ*hZ+cW*hW)*wdelta=cZ*psi'],
        rank9_basis='Inremaining10sourceframe,thethreepurevectorsand6pairvectorshaveinvertibleminorpsi_pivot^6. Phi invertible,Q10sourcequotientandCdelta liftgiveactualrawH7kernel9.',
        pure_projection='Cdelta pi=0, so correctiondoesnotchange thesaved16pureimagecoefficients.',
        old_source_all_entries_binding_reused='ELL_LFRAME_BORDERS_REVIEW bindsallrowoperations, all15pivots andcompleteR11. Weonlyappend3columns, no oldmatrixregeneration.')
    # Tiny independent-variable identities, never an expansion in(a,e,d).
    def va(i):return {tuple(int(j==i)for j in range(6)):1}
    def fa(a,b,s=1):
        z=dict(a)
        for m,c in b.items():
            z[m]=z.get(m,0)+s*c
            if not z[m]:del z[m]
        return z
    def fm(a,b):
        z={}
        for m,c in a.items():
            for n,d in b.items():
                k=tuple(x+y for x,y in zip(m,n));z[k]=z.get(k,0)+c*d
        return {m:c for m,c in z.items()if c}
    z,w,cZ,cW,hZ,hW=[va(i)for i in range(6)]
    ee=fa(fm(cZ,z),fm(cW,w));pp=fa(fm(z,hW),fm(w,hZ),-1);rr=fa(fm(cZ,hZ),fm(cW,hW))
    assert not fa(fa(fm(ee,hZ),fm(rr,z),-1),fm(cW,pp))
    assert not fa(fa(fm(ee,hW),fm(rr,w),-1),fm(cZ,pp),-1)
    assert not fa(fa(fa(fm(z,hZ),fm(w,hW),-1),fm(w,hW)),fm(z,hZ),-1)
    proof['both_Cdelta_Hzero_templates_all_Z_coefficients_checked']=True
    proof['all7_source_basis_9minor_diagonal_coefficients_checked']=True
    cost=dict(prefix_operations=counts,prefix_bytes=prefixbytes,prefix_parameter_terms_counting_repeats=paramterms,
        materialized_pure_R10_terms=pureterms,existing_psi16_terms=sum(map(len,a['psi16_coefficients'])),
        raw_HZHW32_terms=sum(len(p)for row in a['raw16_HZ_HW_coefficients']for p in row),
        active_unique_explicit_parameter_terms=o.explicit_terms,scalar_DAG_nodes=len(o.nodes),
        no_parameter_polynomial_products_expanded=True,no_POWER207_atoms=True)
    data=dict(scalar_nodes=o.nodes,scalar_type='rationalfunctions ONLY in a,e,d; inverseparameters explicitlylocalized',
        source26_forward_Phi=phi,source_prefix_operations=prefix,source_inverse_operations_reversed=list(reversed(inverse)),
        source10_genuine_relative=relative,source10_initial26_rows=rids,pure16_original31_columns=purecols,
        actual10x17_roots=[R[i]+[psi[i]]for i in range(10)],psi10_roots=psi,HZ10_roots=zz,HW10_roots=ww,
        pure_three_positions=Z,seven_possible_psi_positions=P,raw_Cdelta_data=dict(cZ=cz,cW=cw,eta=eta,zdelta=zd,wdelta=wd),
        seven_kernel9_charts=charts,field_scope=fieldscope)
    forward_roots=[x for ch in charts for row in ch['pure9x16_projected_roots']for x in row]+[x for ch in charts for row in ch['source9']for x in row['initial26_coefficients']]
    ancestors=set();todo=forward_roots[:]
    while todo:
        n=todo.pop()
        if n in ancestors:continue
        ancestors.add(n);z=o.nodes[n]
        if z[0]in ['add','mul']:todo+=z[1:]
        elif z[0]in ['neg','inv']:todo.append(z[1])
    inverse_nodes={n for n in ancestors if o.nodes[n][0]=='inv'}
    normalizer_inverse_nodes={p['inverse_scale']for p in prefix if p['op']=='NORMALIZE'}
    assert inverse_nodes<=normalizer_inverse_nodes|{ie}
    normalizer_localizers=[x for x in localizers if x['kind']=='NORMALIZE_g']
    assert len(normalizer_localizers)==38
    forward_scope=dict(forward_inverse_nodes=sorted(inverse_nodes),all_only_normalizer_or_eta=True,
        normalizer_nonconstant_occurrences=38,known_eta_inverse=ie,
        forward_only_sufficient_domain='FULLU intersect H7 intersect D(product38NORMALIZE_g), plus selectedpsi chart. Keeporiginal136+15forrankequivalence.',
        proof='ForwardPhi identitiesandzeroheavycomponents remainvalid whereverNORMALIZEdenominatorsareunits, evenifELIMuorpivotvanishes. Sourcekernelcofactorsandprojectedimagesarethen defined. Ifnineprojectedimagesareindependent, nineactualrawHkernelvectorsareindependent; rawH7+product41 givesfull57. No basis/rankequivalenceassertedoutsideoldOmega beforethisranktest.',
        no_unknown_complement_removed=True)
    fieldscope['larger_forward_only_sufficient_chart']=forward_scope
    proof['forward_only_denominator_ancestor_gate']=forward_scope
    for name,x in [('circuit.json',data),('localizer-inventory.json',dict(all136_inverse_occurrences=localizers,all15_pivots=pivatoms,classification_counts=classification)),('symbolic-proof.json',proof),('cost-and-source-profile.json',dict(cost=cost,pure_row_profiles=rowprofiles))]:
        payload=json.dumps(x,separators=(',',':'))+'\n';assert len(payload.encode())<30*1024*1024;(out/name).write_text(payload)
    paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',lb/'report.json',lb/'certificate.json',lb/'operations.txt',lb/'completed-R11x16.txt',lb/'current-stage.json',lb/'current-matrix.txt',lf/'report.json',lf/'genuine16x31.txt',ab/'report.json',ab/'input.json',ub/'report.json',ub/'certificate.json',hb/'report.json',hb/'certificate.json',ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63_REVIEW/report.json',ROOT/'reports/S61_NONSQUARE_CUBIC_PURE_HEAVY_KERNEL_SCOPE_20261104_v63/report.json',ROOT/'reports/S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63_REVIEW/report.json']
    report=dict(id=ID,all_passed=True,typed_source10x17_and_seven_kernel9x16_constructed=True,
        all10_psi_roots_and_Cdelta_actual7Hzero_lifts_saved=True,nonconstant_inverse_occurrences_retained=136,pivot_localizers_retained=15,
        active_explicit_parameter_terms=o.explicit_terms,scalar_nodes=len(o.nodes),
        no_global_chart_coverage_rank_minor_GB_or_parameter_expansion_claim=True,independent_review_pending=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:main(a.output)
    except Exception as e:
        if a.output.exists():(a.output/'failure.json').write_text(json.dumps(dict(all_passed=False,error=str(e),no_auto_retry=True),indent=2)+'\n')
        raise
