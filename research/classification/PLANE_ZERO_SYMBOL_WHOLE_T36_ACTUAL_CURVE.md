# 零symbol全39源の実T36曲線による全点被覆

2026-11-08、PRIMARY E、producer plane_collision_review。
**全零symbol laneは既知閉成分 T36 に入る。新成分増分は零。**
既受理の全39-source一実軌道定理を使い、一つの prescribed actual curve を true completed T36 から構成する。
旧 local normal六変数や固定Cartan pencilを再実行したものではない。

## 1. 真のT36全式と実座標変換

入力は `reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T36.json` の full completed source。
その八つの非零括弧は、座標 (x0,x1,x2,x3,x4)=(a,b,c,u,v) として

```
Pab=t0*u², Pac=t1*v²,
Pbc=p0*b*c+10*t2*t3*u*v/p0,
Pbu=3*p0*b*u/5+t2*v², Pbv=p0*b*v/5,
Pcu=-p0*c*u/5, Pcv=-3*p0*c*v/5+t3*u²,
Puv=-p0*u*v/5.
```

他二括弧は零。この表示には Pbc の必要 correction を含む。
全sourceがPoissonであることは既 exact completed-source certificate を再使用する。
period open は 3p0²/25 非零である。

ε≠0 で

```
p0=ε, t0=t1=3/ε, t2=t3=1
```

とし、sourceからtargetへの実GL5を

```
A=a+(u+v)/ε,
B=b+5(v-u)/ε,
C=c+5(v-u)/ε,
u=u, v=v
```

とする。六shear係数は (1,1,-5,5,-5,5)/ε。detは一、逆は三shiftの符号反転。
追加tensor scaleは一、Euler補正は零である。
source period/smoothing parametersとこのmapの全ては同じ一つの ε∈Gm 上で定義される。

全 raw pushforward を M P(M^{-1}X) M^t と計算すると、ε^{-1} の全係数が消え、
正確に

```
P_ε=ε L+R
```

となる。L の非零括弧は

```
L_BC=BC, L_Bu=3Bu/5, L_Bv=Bv/5,
L_Cu=-Cu/5, L_Cv=-3Cv/5, L_uv=-uv/5.
```

R は次の全十括弧である。

```
R_AB=-B(3u+v)/5, R_AC=C(u+3v)/5,
R_BC=B(2u-4v)+C(4u-2v),
R_Au=uv/5, R_Av=-uv/5,
R_Bu=3u²-2uv+v², R_Bv=2uv-v²,
R_Cu=-u²+2uv, R_Cv=u²-2uv+3v², R_uv=0.
```

これはprojective classだけの一致ではなく、全150係数のraw恒等式である。
従って P_ε は ε=0 まで正則な実多項式Poisson曲線である。
ε≠0でPoissonなのは元completedsourceと実GL変換から従い、元Jacobiは多項式なのでε=0にも延びる。
新しいJacobi operatorを構成する必要はない。

## 2. endpointの全actual quinticと全zero-symbol仮定

E∧R∧R/2 の signed Euler-Pfaffian は、座標順(A,B,C,u,v)で

```
F_R=-3A(u²-uv+v²)²+(3/5)B*u⁴-(3/5)C*v⁴.
```

全126係数の恒等式である。特に F_R は恒等零でなく、Rのprojective rankは四。
C[u,v][A,B,C]でbase variablesに関する一次式で、B/C係数は 3u⁴/5、-3v⁴/5。
そのgcdは一なのでprimitiveであり、C(u,v)上の一次式の既約性とGauss lemmaからC上絶対既約である。
よって被約性も従う。generic chosen parameter の議論ではない。

R はbase-base行でnormal degree一、base-normal行でnormal degree二、normal-normal行は零。
したがってplane u=v=0 の restrictionは零、normal-weight-zero symbolは全零、P2は零である。
実Dのplane沿い位数は正確に四。
既受理 `PLANE_ZERO_SYMBOL_WHOLE_SINGLE_ORBIT.md` の全39-source theoremの全仮定を満たす。
anchor共通零、abelian kernel、split extension、repeated extension cubic、CE cocycle の全zero補集合は
その全源定理と別作者reviewの結論をそのまま使い、今回のendpointだけから一般化しない。

## 3. 一endpointから全39源への正当な橋

ε≠0では P_ε は真のcompletedT36 sourceの実PGL像なので、既知の閉PGL不変成分 T36 に入る。
P_ε は非零projective tensorの全A1曲線であり、P_0=Rも非零。
従って閉性により R∈T36。

旧全源定理は、指定された全零symbol39-sourceの任意rank四／被約既約quintic点が
一つのactual PGL軌道へ入ることを証明している。
今回 R がそのdomainの点であるため、その同じ全軌道が閉PGL不変 T36 に入る。
従って **全零symbol39-sourceの全点がT36**。
P2の任意cohomology係数、全frame、全base変換は旧定理の真の有限shear／block mapで戻る。
単にbundle orbitが一つ、又はdivisorが一致するというだけの推論ではない。

orbit23／接空間29という旧ambient診断からsmoothness又は唯一親を推論しない。
ここで証明するのは全閉成分所属であり、T36以外の親への所属・交差成分・局所scheme構造は主張しない。

## 4. finite exact readbackと実装失敗の区別

rootが凍結契約 `jobs/plane_zero_symbol_t36_actual_curve.json` を審査し、ONE prescribed tiny readbackを承認した。
最終再生は **0.189590秒／RSS49012KiB、全PASS**。
全150 Laurent push、全raw逆、detと左右逆、全126 endpointPf、全curve Pf共変性、
primitive gcd、全zero-symbol supportを確認した。報告には実行source全文とSHAを埋め込んだ。
sourceの元350Jacobiは既証明を再使用し、新operator/profile/point/GB/jetは零。

初回は全150正逆／全126／gcdまで成功後、最後のnormal-weight loopで
SymPy Poly(0).terms()がliteral zeroを一項返す実装guardにより停止した。
初回sourceと契約exactbytes／errorを別IMPLEMENTATION_FAILURE reportに不変保存した。
変更は零coefficientをskipする一行だけで、同じ入力／同じ曲線／同じresource予算を再生した。
失敗runを成功と呼ばず、最終reportを全PASSの根拠とする。

会話中の別diagnosticに「四つのdouble-zero全てが22型」と述べた点も訂正する。
元canonical Fのplane内四点は、一つが22型3z²w²、他三点は
F4(pζ)=(ζ²/2)(z-ζ²w)^4、ζ³=1という4型である。
本T36曲線及び全source所属証明はその誤った全22主張に依存しない。
保存済み有理Casimir／四点noteは全local同型型を主張していないので、そのbytesは不変に保つ。
