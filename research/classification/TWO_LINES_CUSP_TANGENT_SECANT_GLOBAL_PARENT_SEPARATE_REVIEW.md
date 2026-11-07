# cusp tangent/secant 全源の既知三親除外：独立全係数・proper flag監査

2026-11-07、PRIMARY E、`/root/log_manuscript`。
producer `TWO_LINES_CUSP_TANGENT_SECANT_GLOBAL_PARENT_EXCLUSION.md`
（SHA256 `9028e7636bf43714890eabf809a1b3b69ec3654ef8f776daa7601117a9986905`）を別作者として読解し、
全parameter算術を独立Fractionで再構成して **PASS** とする。
新成分、他親不存在、I=0のC39/C40外性は結論しない。

## 1. 不変入力と一回の独立算術

入力入口は `inputs/two_lines_cusp_tangent_secant_global_parent_manifest_20261107.json`
（SHA256 `68f5fe5517a1b7ad292461e9a9588fb979e9946c3d93bfee7becc54c909ae3b9`）。
原全源は `CLASSIFICATION_TWO_LINES_PLANE_CUBIC_20261107.json` の
`cusp_cusp_tangent_cusp_secant` であり、全六磁気係数を保持した。
旧producer/source/report/jobは改変せず、旧二次法線診断は依存に含めない。

新しい `classification_two_lines_global_parent_review.py` は、過去の独立
`classification_minimal_quartic_review.py` / `classification_quadratic_cofactor_review.py`
のsparse Laurent/Fraction演算だけを使う。producer、S61 engine、CASをimportしない。
有限job `jobs/two_lines_global_parent_separate_review.json` をgate PASS後に一回実行し、
**0.718秒**で次を確認した。実測RSS・実行source全文・入力hashは新reportに保存する。

* 全括弧、全六磁気上のstrict十式と元projective Jacobi全350係数。
* 全六磁気translation、det一、正逆全150係数。同じ全六源からのI非零fiber尺度とtensor尺度の正逆。
* 全λの同一rank四点、Pf恒等値−9、Darboux・perspective・全PGLと逆のpolynomial identity、det恒等値−1/9。
* 全61値からの元150係数の復元、Euler quotientの正逆両方向、元134式、復元tensorの全350式。
* 両全45 graph残式をQQ[λ]上で比較。T39は11本、T40は13本非零。
  非零定数1/4、−2/3をそれぞれ含むので、各残式idealは厳密に(1)。
* T21の全五parameter raw source、全Euler gauge、両intrinsic field条件。
* 原全六源のactual D、全modular fieldとlog-Hamiltonian、条件違反の恒等式。

λの一点標本・generic inverse・追加anchor/frameは使わない。
QQ[λ]環のparserではparameter逆元を許しておらず、全61値・全変換・両45残式に非定数分母は無い。
I=r5だけを反転する別Laurent環での磁気尺度と混同しない。

## 2. 磁気全源とC39/C40外性の正しい範囲

原括弧はh-Casimirで

```
{x,z}=−hx,       {x,w}=3hx−3x²,
{y,z}=3(h−x)²−3hy, {y,w}=−9xy,
{z,w}=3xz+3(r0h²+r1hx+r2hy+r3x²+r4xy+r5y²).
```

translation
`ell=(r3+r2)x/2+r4y/4`、
`kappa=(−3r1−3r3/2+9r2/2)x−r2y`、
`z'=z+ell,w'=w+kappa`
の逆は同じbase上で差を取る式であり、全A6上で正則である。
全150の直接変換は磁気部分を `3(Jh²+Iy²)`、I=r5、J=r0+r2へ送る。
I非零でfiberをIで割りtensorをI倍すると、全anchorは元に戻り、磁気だけが `3(λh²+y²)`、λ=J/Iとなる。

既受理の `S61_GLOBAL_COMPONENT_EVALUATION_FIBER.md` と独立reviewの全scopeを再読した。
connected Sp4を持つ評価associated bundleから全J-fiberは既約、定数det36の全perspective商はA4積、
元55消去は全域polynomial三角消去なので、各25成分の全normalized交叉は既全A16 graphに等しい。
したがって今回の全rank四評価と全正逆変換は、このglobal必要十分判定の仮定を満たす。
両残式idealが(1)であることは、全有限λ、従ってI非零の全六源がC39/C40双方の外であることを証明する。
**陰性を閉境界I=0へ延長する論理は無い。** I=0は元sourceのまま残す。

## 3. 二単純Hのintrinsic条件とproper flag

rank四tensorに相異なる単純Poisson超平面H,Wがあるとする。
Hを無限遠に取るquadratic affine代表は、そのh-Casimir homogeneous liftとして一意に定まる。
H-Casimir gauge、quadratic代表の発散Z_H、正則linear field
`X=Π#dlog(W/H)` は独立な二linear formsの各frame chartで正則に計算される。
Wによる除算は各該当括弧の正確な多項式除算で、W非零への局所化ではない。
affine volumeの定数倍とhomogeneous Euler差は、以下のprojective条件を変えない。

```
(a) Z_H+X はW上でprojectively零、
(b) Z_H+4X はH上でprojectively零。
```

各々は、fieldとEulerのwedgeを指定超平面へ制限した全係数の消滅である。
したがって独立H,Wのflag open上で閉条件となる。
volume変換 `vol_W=(W/H)^(-5) vol_H` とdivergenceの符号から
`Z_W=Z_H+5X`、`X_new=−X`。
よってswapは(a),(b)を交換し、両条件の同時成立は順序に依存しない。

ここでpair flagのopen自体をproperと呼ばない。
全射影係数積 `P(H)×P(W)×P(C3)` の乗法はproperであり、rank四openではtensorから実Dへの写像も正則である。
一般sourceのtensorとfactor flagのgraphをこのproper積内で閉じれば、全rank四boundary点はflagの極限を持つ。
target Dが依然二つの相異なる単純Hと既約C3なら、UFDにより極限の二linear factorsは正確にこの二Hである。
対角H=W、cubicへのlinear因子の紛れ込みはtargetの単純性・既約性に反する。
したがってtarget近傍は独立flag openに戻り、そこでの閉条件(a),(b)が両方延長される。
非proper openで単に「閉」と言うだけの議論でも、H=Wの境界へ本条件を無理に適用する議論でもない。

## 4. 完成T21全五parameterのliteral計算

`MPS_S61_RAW_CATALOGUE_20261103_v62/T21.json` の全p0,p1,t0,t1,t2を戻す。
座標A,B,C,H,Wに対し
`L=(0,−p1B/3,p1C/3,0,p0W)` を用いた `P_H=P−E∧L` を全括弧で再構成した。
H行は零、残りは

```
f=(p1ABC+t0C³−t1B³+t2A³)/3,
P_H=Jac(f)+p0 W ∂W∧E_ABC.
```

period又はedgeを反転せず、発散とlog-Hamiltonianは

```
Z_H=p0 E_ABC−3p0 W∂W,  X=−p0 E_ABC,
Z_H+X=−3p0 W∂W,
Z_H+4X=−3p0(E−H∂H).
```

これが(a),(b)の全parameter恒等式である。
この完成sourceのPGL飽和閉包がT21全成分という依存は、既受理
`MPS_S61_CATALOGUE_20261103_v62.md` §§1–2の全period/edge、像次元、MPS局所成分同定定理に置く。
今回「八括弧がJacobiを満たす」「同じ次元」というだけで新たに成分IDを同定したのではない。
従って前節のproper flagから、指定二単純H+既約三次を持つT21全boundaryも(a),(b)を満たす。

## 5. 原全六源の違反はI=0でも消えない

全A6 sourceについて直接

```
D=9hx((h−x)³−h²y),
Z_h=−4h∂z+(3h−12x)∂w,
X=h∂z+(−3h+3x)∂w,
Z_h+X=−3h∂z−9x∂w,  Z_h+4X=−9h∂w.
```

磁気項は二linear fieldsに寄与しない。cubic `(h−x)³−h²y` はyについてprimitiveな一次式で、
二係数のgcd一から絶対既約。h,xは相異なる単純因子である。
x零上の第一fieldとEulerの(h,z)成分は恒等的に **−3h²** なので、(a)は全sourceで失敗する。
swapは両条件を交換するだけで、片方だけ成立しても救済されない。
したがって **I=0を含む原全六源はT21全成分の外**。

新reportは `CLASSIFICATION_TWO_LINES_GLOBAL_PARENT_SEPARATE_REVIEW_20261107.json`。
全算術と全pure proofの含意を合わせた受理範囲は、I非零でのC39/C40外性、全A6でのT21外性だけである。
T36/T37/T38、Ga9 cuspidal24又は他の親、未知成分性、全E/G、成分数へ外挿しない。
