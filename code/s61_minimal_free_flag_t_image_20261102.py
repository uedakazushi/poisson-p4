import argparse,hashlib,json,time,tarfile,numpy as np
from pathlib import Path
from itertools import combinations_with_replacement as cwr
from s61_minimal_free_flag_transverse_sources_20261102 import echelon
ROOT=Path(__file__).resolve().parents[1]

def image_layer(Rfull,Tfull,degree):
    p=32003
    mons=[m for d in range(degree)for m in cwr(range(3),d)]
    idx={m:i for i,m in enumerate(mons)}
    n=141*len(mons)
    E,piv=echelon(Rfull[:,:n])
    for i in range(len(piv)-1,-1,-1):
        ids=np.flatnonzero(E[:i,piv[i]])
        if len(ids):E[ids]=(E[ids]-E[ids,piv[i],None]*E[i])%p
    free=[j for j in range(n)if j not in piv]
    W=np.zeros((n,n),dtype=np.int64)
    for a,m in enumerate(mons):
        for b,mon in enumerate(mons):
            if len(m)+len(mon)<degree:
                c=idx[tuple(sorted(m+mon))]
                W[a*141:(a+1)*141,c*141:(c+1)*141]=Tfull[:,b*141:(b+1)*141]
    W=(W[:,free]-W[:,piv]@E[:,free])%p
    dimension=len(echelon(W)[1])
    times_u=len(echelon(W[141:])[1])
    return dict(r=degree,ambient_length=n-len(piv),image_length=dimension,base_times_image_length=times_u,minimal_generators=dimension-times_u)

def run():
    p=32003;st=time.time()
    with tarfile.open(ROOT/'research/wip/S61_MINIMAL_FREE_FLAG_THIRD_BASE_QUOTIENT_20261102.tar.gz')as tf:
     lines=tf.extractfile('closure.txt').read().decode().splitlines();actionlines=tf.extractfile('actions.txt').read().decode().splitlines()
    R=np.zeros((294,2820),dtype=np.int64)
    for i,line in enumerate(lines[1:295]):
     z=list(map(int,line.split()));R[i,z[1::2]]=z[2::2]
    E,piv=echelon(R)
    for i in range(293,-1,-1):
     ids=np.flatnonzero(E[:i,piv[i]])
     if len(ids):E[ids]=(E[ids]-E[ids,piv[i],None]*E[i])%p
    free=[j for j in range(2820)if j not in piv];T=np.zeros((141,2820),dtype=np.int64)
    for j,line in enumerate(lines[298:439]):
     z=list(map(int,line.split()));T[j,z[1::2]]=z[2::2]
    mons=[m for d in range(4)for m in cwr(range(3),d)];idx={m:i for i,m in enumerate(mons)};W=np.zeros((2820,2820),dtype=np.int64)
    for a,m in enumerate(mons):
     for b,n in enumerate(mons):
      if len(m)+len(n)<4:
       c=idx[tuple(sorted(m+n))];W[a*141:(a+1)*141,c*141:(c+1)*141]=T[:,b*141:(b+1)*141]
    W=(W[:,free]-W[:,piv]@E[:,free])%p
    powers=[]
    for d in range(5):powers.append(len(echelon(W[sum(len(m)<d for m in mons)*141:])[1])if d<4 else 0)
    print('powers',powers,flush=True)
    U,upiv=echelon(W[141:]);C=W[:141].copy()
    for row,j in zip(U,upiv):
     ids=np.flatnonzero(C[:,j]);C[ids]=(C[ids]-C[ids,j,None]*row)%p
    _,gens=echelon(C.T);ng=len(gens);assert ng==12
    D=W[[a*141+j for a in range(20)for j in gens]].copy();A=D.copy();tr=np.eye(240,dtype=np.int64);pc=[];h=0
    for j in range(2526):
     ids=np.flatnonzero(A[h:,j])
     if not len(ids):continue
     k=h+int(ids[0]);A[[h,k]]=A[[k,h]];tr[[h,k]]=tr[[k,h]];iv=pow(int(A[h,j]),-1,p);A[h]=A[h]*iv%p;tr[h]=tr[h]*iv%p;ids=np.flatnonzero(A[h+1:,j])+h+1
     if len(ids):
      co=A[ids,j,None].copy();A[ids]=(A[ids]-co*A[h])%p;tr[ids]=(tr[ids]-co*tr[h])%p
     pc.append(j);h+=1
    assert h==76 and not np.any(A[h:]);K=tr[h:];assert not np.any(K@D%p);assert not np.any(K[:,:12]);assert len(echelon(K)[1])==164
    UK=np.zeros((492,240),dtype=np.int64)
    for j in range(3):
     for a,m in enumerate(mons):
      if len(m)<3:
       b=idx[tuple(sorted(m+(j,)))];UK[j*164:(j+1)*164,b*12:(b+1)*12]=K[:,a*12:(a+1)*12]
    ukrank=len(echelon(UK)[1]);assert ukrank==126
    print('image relations',164,'uK',ukrank,'minimal',164-ukrank,flush=True)
    def solve(v):
     v=v.copy();out=np.zeros(240,dtype=np.int64)
     for row,(j,trow)in enumerate(zip(pc,tr[:76])):
      c=int(v[j])
      if c:v=(v-c*A[row])%p;out=(out+c*trow)%p
     assert not np.any(v)
     return out
    acts=[]
    for var in range(8):
     mat=[]
     for gen in gens:
      z=list(map(int,actionlines[1+var*141+gen].split()));target=np.zeros(2526,dtype=np.int64)
      for i,c in zip(z[1::2],z[2::2]):target=(target+c*W[i])%p
      coeff=solve(target);assert np.array_equal(coeff@D%p,target);mat.append(coeff)
     acts.append(mat)
    print('all96 actions',time.time()-st,flush=True)
    # All induced actions preserve the finite-image presentation, via its independently computed decoder.
    expanded=[]
    for var,mat in enumerate(acts):
     X=np.zeros((240,240),dtype=np.int64)
     for a,m in enumerate(mons):
      for b,n in enumerate(mons):
       if len(m)+len(n)<4:
        c=idx[tuple(sorted(m+n))];X[a*12:(a+1)*12,c*12:(c+1)*12]=np.array(mat)[:,b*12:(b+1)*12]
     assert not np.any((K@X%p)@D%p)
     expanded.append(X)
    for i in range(8):
     for j in range(i+1,8):
      assert not np.any(((expanded[i][:12]@expanded[j]-expanded[j][:12]@expanded[i])%p)@D%p)
    print('stability and commutators PASS',time.time()-st,flush=True)
    sparse=lambda r:[[int(i),int(c)]for i,c in enumerate(r)if c]
    result=dict(prime=p,scope='image(t:S/a^4 -> S/a^4)',length=76,own_base_ideal_power_lengths=powers,own_Hilbert_function=[powers[i]-powers[i+1]for i in range(4)],minimal_generators=12,generator_B_indices=gens,relation_dimension=164,base_times_relation_rank=ukrank,minimal_relations=164-ukrank,relations=[sparse(r)for r in K],normal_actions=[[sparse(r)for r in a]for a in acts],seconds=time.time()-st)
    result["ambient_image_layers"]=[image_layer(R,T,d)for d in range(1,5)]
    assert [x["image_length"]for x in result["ambient_image_layers"]]==[4,18,40,76]
    assert [x["minimal_generators"]for x in result["ambient_image_layers"]]==[4,7,8,12]
    result.update(id="S61_MINIMAL_FREE_FLAG_T_IMAGE_20261102",all_passed=True,actual_tS_minimal_generator_lower_bound=12,relation_stability_checks=8*164,commutator_generator_checks=28*12,limitations=["Finite image tS/(tS intersect a^4 S), not tS/a^4(tS).", "No upper bound12 on actual tS generators, no full-S t-nilpotence, and no QQ or relative lift."])
    paths=[Path(__file__),ROOT/"reports/S61_MINIMAL_FREE_FLAG_THIRD_BASE_QUOTIENT_20261102.json",ROOT/"research/wip/S61_MINIMAL_FREE_FLAG_THIRD_BASE_QUOTIENT_20261102.tar.gz",ROOT/"code/s61_minimal_free_flag_transverse_sources_20261102.py",ROOT/"research/notes/S61_MINIMAL_FREE_FLAG_T_IMAGE_20261102.md"]
    result["input_sha256"]={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    result["seconds"]=time.time()-st
    return result

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    result=run()
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",result["seconds"])
