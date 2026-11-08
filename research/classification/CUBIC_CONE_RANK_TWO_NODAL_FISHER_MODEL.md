# Nodal quinticの明示Pfaffianモデルと多項式Fisher認識

2026-11-08、PRIMARY E。前の `NODAL_QUINTIC_T39_BRIDGE` のsmoothability・五quadricsを使う。
Casimir-kernel operatorは実行せず、既受理の多項式Fisher構成による、より直接な陽性認識を選ぶ。
元の二次余核零という結果は変更せず、local高次計算も行わない。

## 1. 全Weierstrass族のliteral alternating model

前noteの変数(z0,z1,z2,z3,z4)、parameter eを用い、次のlinear 5×5交代行列Phiを取る。
非対角上三角は順に

```
Phi01=−z3−e*z0, Phi02=−z3, Phi03=−z4, Phi04=z2,
Phi12=z2, Phi13=z3, Phi14=−z1,
Phi23=−z1, Phi24=z0, Phi34=0.
```

各iを除いた4×4Pfaffianに(−1)^iを掛けた五成分p_iは、前noteの五g_iの
正確に−一倍である。特にPhi*p=0。
e=0でnodal elliptic quintic、generic eでsmooth elliptic normal quinticを与える。
smooth locusの非空性は明示cubicのdiscriminant e(4+27e)≠0で確認出来る。
五quadric basis又はnormalization-flatの無証明仮定は用いない。

## 2. 既受理Fisher式の正確な使用

`research/notes/KNOWN_FIVE_SMOOTH_LOCUS_20261101.md` §2とその独立reviewは、
Fisherのquadratic Omegaとsmooth elliptic normal quinticのq(5,1)を同定している。
独立source `code/known_five_smooth_locus_review_20261101.py` の32–50行と同じ規約で、
a<b、残三添字c<d<f、epsilon=sign(a,b,c,d,f)に対して

```
Omega_ab=epsilon * [
  (d_c p)^T (d_d Phi) (d_f p)
 +(d_f p)^T (d_c Phi) (d_d p)
 +(d_d p)^T (d_f Phi) (d_c p) ].
```

これはPhiの全五十係数の多項式であり、discriminant、Jacobian、curve parameterの
非定数を一切除さない。smoothモデル上で既受理q5,1=T39構成である。
上のexplicit smooth Weierstrass族だけにこのsmooth解釈を使えば十分である。
e=0にも同じpolynomial式が定義され、Omega(Phi_0)が非零ならそのclassはclosed T39に入る。
singular curveにsmooth-only theoremを直接適用したのではなく、regularな族の閉性である。

## 3. 残す有限義務と全sourceへの戻し

一つの有限jobで全eを保ってPhiの五Pfaffian／Omega、全ten Jacobi／canonical divergence、
全126Euler-PfaffianとJacobian secant式を直接再構成する。
前noteのconstant GL5のdetと正逆、e=0での全150係数のcanonical代表への非零scalar一致も調べる。
新M1又は全150Casimir operatorは不要で、実行しない。

全150一致が通れば、canonical代表はactualにT39へ入る。
既凍結 `CUBIC_CONE_RANK_TWO_CRITICAL_ORBITS.md` の実shear／Cartan全逆により
semisimple ab≠0全源が同じ軌道へ戻る。ab=0とnil全六fiberはその軌道閉包内なので
全λ3/5 sourceがT39となる。λ2/5は既T38、他の全familiesの既約quintic空性も保持する。
この陽性包含では新component数を増やさず、generic rank・因子条件外の別sourceを消さない。

もしliteral一致が失敗すれば結果をそのまま保存し、同時に別のkernel operatorを走らせない。
判定の前提と独立reviewは新finite packetで固定し、旧shared runnerと全旧source bytesは不変とする。
