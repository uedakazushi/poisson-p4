# trace/traceless零・K非零の全二源：別作者の手証明監査

2026-11-07。reviewer `exceptional_coverage`、producer `generic_exhaustion`。
入力 `PLANE_RANK_ZERO_ZERO_TRACE_BINARY_SYMBOL.md`
SHA256 `36d4f8e73dac3b4037e5d169abf3563e4d49a688e23ef0151459c59555ac2b88`。
**forcing・全逆変換・projective compatibility・全source rank≤2を独立に確認して PASS。**
root の全350/30/9 Fraction 再構成は別の計算証明書であり、この手監査と区別する。

## 全sourceと forcing を直接確認した

normal GL2と非零tensor scalarで、全非零binary quadratic Kはzw又はz²へ移る。
逆行列・逆scalarで元に戻り、全30/9 sourceを保つ。rank零planeでの可換conormalを仮定し、
trace/tracelessの両coHiggsが零であることが正確な入口である。
Gでの全重複度への適用は既受理の被約実D conormal補題に依存する。

Π0=K partial_z wedge partial_w、全Π1は元ノートの L,M,B,C,D,E、全Π2は九係数とする。
normal二次場Vについて

```
L_V(K partial_z wedge partial_w)
 =(V(K)-K*div_normal(V))*partial_z wedge partial_w
```

を使い、zL,wMの Hamiltonian 項も加えると、全 base²normal² forcingは

```
K=zw: z³ Eb wedge D+w³ Eb wedge C-z²w L-zw² M,
K=z²: z³(-Eb wedge E-L)+z²w(Eb wedge B-M)+2zw² Eb wedge C.
```

normal monomialが独立なので、前者は C=D=L=M=0、B,E全六自由。
後者は C=0,L=−Eb wedge E,M=Eb wedge B、B,D,E全九自由。
Eb wedge constant=0はconstant自体が零を意味する。
残る base³normal の二成分は Eb wedge L,Eb wedge Mであり、これらの解では零。
四vectorのindexにはbase²normal²とbase³normalしか無いので、未処理の別成分はない。
特に z²の九次元kernelを六へ減らす仮定は使っていない。

## 有限shearとcanonical Euler representative

実座標変更 a'=a+zY+wZ は可逆、逆はminusである。
Π0のpushforwardのweight一は K(Y wedge partial_w−Z wedge partial_z)、
weight二は全九係数内の K Y wedge Zとなる。
元のΠ1の座標置換・微分による追加もweight二だけで、Π2を狭めない。

K=zwでは B=−Z,E=Yなので全六kernelを消せる。
K=z²では raw shear の −z²Z wedge partial_z を Euler gaugeへ直す。
正確なcanonical式との差は z Z wedge Euler_totalであり、
D=Y,E=Z,L=−Eb wedge Z,B=M=0となる。
従ってD,Eの六係数だけを消し、残る全sourceは

```
Pi1=w Eb wedge B+zw B wedge partial_z
```

である。Bの三係数を捨てていない。

## weight二injectivity、compatible representativeとrank

任意のΠ2のbase²normal²成分は全weight二作用で −2KΠ2となる。
これはbinary coefficientのEuler公式 z*f_z+w*f_w=2f から直接従う。
K!=0と多項式環の整域性により全D2がinjectiveである。
K=zwではΠ1=0なのでΠ2=0が従う。

K=z²では

```
R=B wedge W,  W=2zw partial_z+w² partial_w,
R-Pi1=w B wedge Euler_total.
```

よって真に同じprojective bivectorである。
全coefficientsは二次、gauge vector wBは一次なので、Euler商でのPoisson括弧とも整合する。
Bはconstant-base、Wはnormal変数だけに依存するから [B,W]=0、従って[R,R]=0。
また W(z²)=4z²w、div_normal W=4wより L_W Π0=0、[Π0,R]=0。
従ってΠ1自身のprojective Jacobiは零であり、weight二はD2Π2=0へ戻る。
injectivityによりΠ2=0となる。

残る代表Π0+Rは span(B,partial_z,partial_w) の高々三次元空間に入る。
skew rankは偶数なのでambient rank≤2、Euler商へ降りてもrankは増えない。
全ての有限shear、normalGL2、非零scalarを戻しても同じrankである。
以上で全二源の全higher係数を含む主張を受理する。

## 射程

T=V=K=0にはD2injectivityを適用しない。
その全39-sourceは別のLie algebroid証明で扱う。
非零trace scalar源、非零traceless dependent globalV源も本命題の対象外である。
新CAS、局所環、higher jet、特別点scanは行っていない。
