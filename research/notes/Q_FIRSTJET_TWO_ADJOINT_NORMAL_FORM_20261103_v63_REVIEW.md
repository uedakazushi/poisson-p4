# 二 adjoint の七座標切片：独立監査

2026-11-03。対象 `Q_FIRSTJET_TWO_ADJOINT_NORMAL_FORM_20261103_v63`。
四関係の scheme-theoretic な graph 消去と、同じ既存一点の chart 非空性を受理する。
この七変数・一方程式は周囲の6次元切片であり、NEW25 の未知5次元 divisor の方程式ではない。

## 1. 四式と全77係数を別算術で再構成した

author module と GB は使っていない。全一般 sp4 行列 C と対角 D を、Fraction 係数の
Laurent 多項式辞書で作り直した。全 symplectic 条件と trace の全係数を検査し、
四式 f1,f2,f3,f4、H=(f3+f4)/12、Newton による determinant 置換恒等式を比較した。
4角 determinant は subset 展開であり、author の symbolic determinant 消去を再実行していない。

Delta=x²−y²、M=gj−hi、L=(a+x)²−y²として、b=e=1 の gauge 上で

```
d0=−(ax+x²+y²)/y,
c0=M−L,     f0=M−x²L/y²
```

を独立に作った。z1=d−d0,z2=c−c0,z3=f−f0 に関する三つの forward 恒等式を全係数で確認した。
f2=2y*z1、残り二つの係数行列が [[1,1],[x²,y²]] である。
したがって y と Delta を単元とすれば、保存された reverse の二つの線形 cofactor によって
z1,z2,z3 を元の三式から解ける。Delta=0 へこの消去を延長しない。

det C−Delta² の graph への代入を d,c,f の順に行う際の三つの difference quotient は、
有限幾何和 (t^n−v^n)/(t−v)=sum t^(n−1−r)v^r で独立に計算した。
全 telescoping cofactor と全恒等式が一致した。
最後に y²を掛けると非負指数の整数77項となり、通常全次数8、七座標の全保存係数と一致する。
graph 上で f1=f2=0、f3=−4F/y²、f4=4F/y² も全係数で確認した。

これで、b=e=1 後の QQ[x,y,a,c,d,f,g,h,i,j,1/(xyDelta)] において

`(f1,f2,f3,f4)=(d−d0,c−c0,f−f0,F)`

という両方向の ideal 等号を受理する。標本での一致や radical の等号だけではない。
F の斉次性を仮定しておらず、77項の一部を「不要因子」として捨てることもしていない。

## 2. gauge が有限に取れる範囲

regular semisimple D の固有値を (x,−x,y,−y) と順序付け、symplectic eigenbasis を取る。
必要な開条件は xyDelta!=0 であり、これは四固有値が相異なる条件である。
残る centralizer は二次元 torus diag(t,t^-1,u,u^-1)。共役により b は t²b、e は u²e へ移る。
b,e が単元なら t²=b^-1、u²=e^-1 の有限 etale 拡大で b=e=1 とできる。
この平方根の選択を有理的な gauge と呼ばない。

残る四つの torus 符号の同時反転は中心であり、effective な相対符号は
全 cross-block 座標 (g,h,i,j) を同時に反転する。その全行列係数も検査した。
regular D の順序の違いは有限 Weyl 群8通りである。
bc*ef!=0 なら全ての固有値順序で b,e が非零に取れ、一般 orbit に対する正規形は高々16個となる。
generic stabilizer による重なりはあり得るので、必ず16個とも述べない。

## 3. chart 非空性は同じ既存点で確認した

新しい beta 点は一つも使わない。前の独立二 adjoint 監査で受理した同じ sample[0] の C,D を使う。
S=tr D²/2、P=det D、d=S²−4P を別の Fraction matrix 演算で再計算した。
P,d はともに非零で、全16係数の D^4−S D²+P I=0 も確認した。
したがって、この実点の D は regular semisimple である。

QQ[eta]/(eta²−d) を有理数対で独立に実装し、
X=(S+eta)/2,Y=(S−eta)/2、PX=(D²−Y I)/eta を再構成した。
PX²=PX、PX+conjugate(PX)=I、tr PX=2、および X*Y=P を検査した。
この計算は一般には二次 etale algebra の計算で足り、体であるとの仮定を要しない。
今回の具体的 d は負なので、実点では真の二次体でもある。

二つの block product をそれぞれ独立に

`bc=tr(PX C PX C)/2−tr(D C PX)^2/(4X)`

とその Y-projector 版から求めた。後者は前者の共役と一致し、積は正確に QQ の非零値である。
完全な有理数を独立 certificate に保存した。295/282 bit の非零は QQ 上の完全値で判断し、
補助的な有限体還元から昇格していない。

この式は対角基底の第一blockに対して tr(Cblock²)/2=a²+bc、tr(D C PX)=2xa と計算すれば得られる。
従って bc と ef がともに非零で、四つの b,c,e,f は全て非零である。
有限平方根 gauge が適用でき、七座標切片は実際の像成分と交わる。
同じ点が既に受理済み NEW25 source の点なので、NEW25 との交わりも非空である。
固有ベクトルの数値近似や別の beta 点の探索は不要である。

## 4. 次元と適用範囲

regular D を対角化する8次元の群方向と torus の2方向により、一般 pair orbit は10次元である。
torus Lie algebra の元を diag(t,−t,u,−u) とすれば、Cとのcommutatorの01、23成分は2tb、2ue。
b,e!=0 から、pair stabilizer Lie algebra は零となる。
この開で一般 orbit の分類を七座標 hypersurface の選ばれた成分の有限正規形へ移せる。
先の四gradient rank4 と、今回の unit graph 消去は、その像成分に対応する切片が6次元であることと整合する。

NEW25 の first-jet image は既知の15次元既約 hypersurface で、この開にも会うため、
一般 NEW25 orbit の次元も10、その rational quotient は5次元である。
したがって、その未知条件は周囲の6次元切片の選ばれた成分上でさらに codimension1となる。
有限 cover への引戻しが既約であるとは仮定せず、未知 divisor の方程式は構成していない。

F が全像を尽くすこと、全 hypersurface の既約性・正規性、quotient の有理性、
全退化 orbit の正規形は受理範囲外である。とくに F 自体を NEW25 方程式とは呼ばない。

入口 `code/q_firstjet_two_adjoint_normal_form_review_20261103_v63.py`。
60秒・2GB cap、全係数と同じ一点だけの検査であり、凍結primary・原稿・台帳を変更しない。
