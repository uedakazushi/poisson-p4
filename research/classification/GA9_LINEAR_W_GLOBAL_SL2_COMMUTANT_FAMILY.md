# Ga9 linear-w：global SL2と全九次元commutantによる有限族の準備

2026-11-08、PRIMARY E。**純数学と有限契約の準備だけ。NO GO、未実行。**
Poisson族の存在、T37帰属、他親の非所属又は新成分性は主張しない。

## 1. endpointとglobal SL2

Ga9 canonical四倍tensorをP=P2+P3に分ける。P2は

```
Pxy=4h²−y², Pxz=−hx+4w²+2yz,
Pxw=−4wy, Pyz=hy, Pzw=−4hw,
```

で、P3は

```
P3=16hy∂x∧∂w+4y²∂z∧∂w.
```

N=y∂x+h∂z、T=x∂y+z∂h、H=[T,N]と置く。
H=x∂x+z∂z−y∂y−h∂hで、wはこのSL2のtrivial座標。
L_N P2=L_N P3=0、L_H P2=−2P2、L_H P3=−3P3。
従ってglobal二次係数表現において、P2が生成するmoduleはSym²、P3が生成するmoduleはSym³である。
特にL_T³P2=0、L_T⁴P3=0。Euler quotientでも同じで、canonical Euler projectionはSL2と可換である。

PがPoissonであり、三つのSchouten項のH-weightが−4,−5,−6と異なることから、P2/P3はそれぞれprojective Poissonかつ互いにcompatible。
またSym³全体はV_c(h,x,y,z)∧∂wの形で、V_cはwに依らないので、どの元もself Jacobi零である。

一方、Sym²全体はPoissonではない。B=L_T P2のw-Casimir Euler liftを取り、w=1で直接計算すると、strict Jacobiの(h,x,z)成分は

```
−4z(hx−yz+5).
```

Sym²のself-obstructionは非零のSL2-trivial成分を持つ。この小moduleだけでsemisimple qへ動かす方法は停止し、同じ部分族の係数を増やして追わない。

## 2. 全projective SL2-invariant二次二ベクトルのbasis

二列行列を `M=[[x,z],[y,h]]` と見る。左SL2と可換な右GL2の四fieldを

```
D1=x∂x+y∂y, D2=z∂z+h∂h,
R+=x∂z+y∂h, R−=z∂x+h∂y
```

とする。次の九つを**全て**保持する。

1. 四field D1,D2,R+,R−の全六wedge。
2. `w²∂x∧∂y`。
3. `w²∂z∧∂h`。
4. `w²(∂x∧∂h+∂z∧∂y)`。

それぞれがSL2-invariantであり、w-rowは零。この九つはprojective classで独立で、全不変部分空間を尽くす。理由を記す。

V=W⊗U⊕Cw、dim W=dim U=2とすると

```
Λ²V = (Sym²W⊗Λ²U) ⊕ (Λ²W⊗Sym²U) ⊕ (W⊗U),
Sym²V* = (Sym²W*⊗Sym²U*) ⊕ (Λ²W*⊗Λ²U*)
          ⊕ (W*⊗U*)w ⊕ Cw².
```

raw invariant次元は3+6+4=13。
Euler map A↦E∧Aのinvariant imageは、End(V)^SL2の次元5からscalar kernel一を引いた四次元。
SL2の有限次元表現は完全可約なので、quotientのinvariant次元は13−4=9。

右GL2四fieldはdet M≠0上でframeなので、その全六wedgeは定数係数の組として独立。残三つはw²係数だけを持ち、他六つから独立である。
w-Casimir raw liftでEuler classが零なら、w A_i−x_i A_w=0からAはscalar Euler場で、lift自体が零。
従って表示九元のclassは独立、上の次元と一致して全不変classを尽くす。

これはw-rowを零にした都合のよいsubsourceへの制限ではない。
以下のK=Σk_i K_iは、全projective invariant classの任意の元である。

## 3. split eigenrootsと真のcubic periodを保つ一つの族

binary variablesを(s,t)と書く。SL2 moduleへの同定ρ2,ρ3は

```
ρ2(s²)=P2, ρ2(st)=L_TP2/2, ρ2(t²)=L_T²P2/2,
ρ3(s³)=P3, ρ3(s²t)=L_TP3/3,
ρ3(st²)=L_T²P3/6, ρ3(t³)=L_T³P3/6
```

と固定する。generic rootsをs±εtとし

```
qε=(s+εt)(s−εt),
cε=Σ_(j=0)^3 r_j (s+εt)^(3−j)(s−εt)^j,
Σ r_j=1,
Πε=ρ2(qε)+ρ3(cε)+εK.                                  (1)
```

四r_jは同じmoving eigenframeの**実cubic係数**で、無関係なTaylor係数を独立に足していない。
ε≠0では四productがSym³のbasisであり、全三自由periodを保つ。ε=0ではc0=s³なのでΠ0=P2+P3が正確なGa9 endpointとなる。
反転する必要のあるparameterは無い。(1)は最初から元全150係数空間のpolynomial族である。

qεのdiscriminantは定数倍ε²であるため、そのmodified self-obstructionをinvariant self-bracketで補う最初の候補はεKである。
**このεKという次数制限は一つの明示ansatzの定義であり、全T37 arcへの必要次数定理ではない。**
失敗時にε²K等を足す自動延長はしない。

endpointはSL2-weight二／三の和であるから、この族をGa9を固定するCstar-equivariant graphと呼ばない。
εの有限次数は(1)そのものから従い、targetのtorus stabilizer又は形式IFTから仮定したものではない。

## 4. 一回の有限な全係数課題

未知定数は九k_iと四r_j（Σr_j=1）。(1)はε次数三以下のglobal homogeneous quadratic bivectorなので、全projective Jacobi350式はε次数六以下である。
全係数を保持した有限QQ代数方程式を一度だけ作ることが出来る。normal slice、M1、新ambient rank、追加point又は高次jetは必要ない。

分類へ十分な成功証明書は以下である。

* 一つのexactなQQ又は明示有限代数拡大上のconstant solution。
* 全150係数のpolynomial正則性、Π0=Pのliteral identity、全350 Jacobiのε係数0–6が全零。
* 全五Pfaffian／全126 signed quintic係数と、generic actualD非零。
* generic Πが既受理のknown closed familyへ入る**独立の全源bridge**。例えば、true completed T37のactual period／可逆frame恒等式、又はgeneric Dの絶対既約性と適用範囲を確認した既full pointwise theoremが使える。

最後のbridgeが無ければ、Poisson族が見つかってもGa9 E帰属を宣言しない。
known有限closed和へgeneric像が入れば、その閉性と正確なendpointからGa9全軌道閉包が入る。

成功解が無い場合、又は代数方程式が予算内で未決なら、記録するのは**(1)という全九不変補正ansatzの不適合／未決**だけ。
T37全非所属、Ga9成分性又は全Eは結論しない。別点、追加normal、ε次数増加のseriesへ移らずglobal幾何へ戻る。

## 5. 現在の実行状態

本稿はbasis、全parameterとdegree bound、成功／失敗の意味を固定した準備であり、まだJacobi係数のtool算術やalgebraic solveを実行していない。
rootのwhole source読解と明示GOの後にのみ、180秒／3GiBの一回契約をgateして実行する。
既Ga9六normal二次支持、固定N unfolding、単独unipotent又は1PS陰性のbytesは変更しない。
