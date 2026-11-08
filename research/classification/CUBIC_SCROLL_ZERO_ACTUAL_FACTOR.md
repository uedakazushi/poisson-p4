# Smooth cubic scroll rank零の全 source：actual D は必ず可約又は零

2026-11-08、PRIMARY `E_reduced_nonnormal`。
全 reduced source の純還元と、一回の全350／126係数評価による producer exact PASS。
**smooth cubic scroll 上で零となる指定全 rank四 source に、被約既約 actual quintic は存在しない。**
別作者全元係数 review は独立 receipt で管理し、それ以前に受理済み status や新成分数を増やさない。
可約 D を持つ rank四 Poisson 構造の不存在は主張しない。

## 1. 全 source と四つの閉源

既独立受理の abelian-conormal source は、元三十 coefficient から十二係数 r0,…,r11 への
完全30×12 matrixで与えられる。元350 Jacobi rowspace の rank二十四 basis を用いる。
`CUBIC_SCROLL_ZERO_FOUR_WHOLE_SOURCES.md` は全元二十四式を手で還元し、
全 reduced point support が L6 と三つの Segre cone F/U/S の和に等しいと証明した。
η、v、a−2f、f、u、s の全反転零補集合を保持する。

L6 の全六 parameter と、三coneの rank≤1 matrix の全外積 map を使う。
外積は rank零 matrix を含め全 cone を覆い、chosen nonzero chart や generic factorization のみではない。
四つの r-map は同純証明の§5と、実行 report の各 case に元順序のまま保存した。

旧 generic dimension≤6 や leading-ideal certificate を E帰属としたのではない。
今回の全 source split と full actual factor 恒等式が新しい含意である。

## 2. 元実 tensor からの全係数評価

reuse runner `classification_cubic_scroll_abelian_conormal.py` へ isolated mode
`--smooth-four-sources` を追加した。
直前の current code03e053…と全三契約 exact bytes を preappend snapshotへ保存し、
既 compute/main/dimension 三 FunctionDef の AST は全不変。
旧 source report、旧計算 sourcebytes、current契約の binding を黙って上書きしていない。

active finite gate PASS 後、一回の batch、**1.76058秒／RSS53,216KiB**。
新 native、GB、rank、minor、matrix profile、parameter point 又は jet 計算は無い。

各全 map で、元24 generator と元350保存行の全 parameter係数が恒等零。
さらに原30×12 sourceから ten actual bracket を独立に組み、
全元 projective Jacobi E∧[Π,Π] を全物理変数と全 parameter 上で恒等零と確認した。
sourceは necessary linearization の chosen sliceでなく、全Poisson actual polynomial family である。

各元 tensor の signed Euler Pfaffian を定義から作り、
**全126 actual D coefficients**が元十二 source の保存二次係数式を mapした結果と完全一致した。
以下の全 polynomial factor product も元 actualD と展開一致した。
正式 report は `CLASSIFICATION_CUBIC_SCROLL_ZERO_FOUR_SOURCES_20261108.json`。
全 ten brackets、126係数、因子、parameter map、実行時source全bytes/SHAを保存した。

## 3. 全 actual factor の短い式

物理座標は x0,…,x4、scroll ideal は

`q0=x0*x3−x1*x2, q1=x0*x4−x1*x3, q2=x2*x4−x3²`。
以下の F は full signed Euler Pfaffian Σ_i(−1)^i x_i Pf(Πのi補集合)である。

### L6

純証明の namespace `(P,A,S,Q,V,R)` をここでは `(l0,…,l5)` とする。

```
H = l1*x3−l2*x2+l4*x4,
K = l2*q0²−l0*q0*q2−l1*q0*q1−l3*q1*q2−l4*q1²−l5*q2².
```

全源で **F_L=H*K**。両因子の物理次数は一と四。

### 三つの Segre cone

全五外積parameterを `(A,B,C,D,E)` とする。

```
L_F = (4A²D−2ACE)*x0 + (−2ACD+C²E)*x1 +4ABE*x2
      +(−4ABD−2BCE)*x3 +2BCD*x4,

L_U = 2ACE*x0−C²E*x1−4A(AD+BE)*x2
      +2C(2AD+BE)*x3−C²D*x4,

L_S = AE*x1+2BE*x3+(−AC+2BD)*x4.
```

全 polynomial ring 内で

```
F_F =  (1/2)*q2*(E*q0−D*q1)*L_F,
F_U = −(1/2)*q2*(E*q0−D*q2)*L_U,
F_S =  (A/2)*q2*(C*q2+D*q1+E*q0)*L_S.
```

全三式の物理因子次数は二・二・一。S の A は parameter-only factor である。
report の factor list の展開を上の q0,q1,q2 表示へ書き直したもので、
scalar又はEulerのnormal symbolだけの因子ではない。

## 4. 全 parameter specialization、全 inverseへの含意

これらは source parameterを一つも反転しない polynomial identity。
positive homogeneous physical factor は、任意の parameter specialization後、
零になるか、同じ positive physical degreeを保つ。
零因子又は parameter-only factor の零は F=0であり、rank四・quintic対象にはならない。
非零Fなら、L6では少なくとも一次×四次、三coneでは二次×二次×一次の proper factorを持つ。
一因子が別因子と一致した場合も、被約既約quinticにはならない。

外積 rank零、各 parameter零、coneの特異点を generic chartへ外挿していない。
四源の全像を覆い、全反転補集合を純証明で戻したため、
**全十二 reduced Poisson source に被約既約 actual quintic は無い**。
元scrollの全 actual PGL normalizer inverse も factor degree／可約性を保つ。
従って指定 smooth cubic-scroll rank零・被約既約 quintic branch は全空である。

twisted cubic cone、scroll上のrank二制限、他の導手曲面又は可約Dの全E帰属は別である。
新成分数は零。本結論の採用は独立 whole-source／全係数 review後、統合担当が行う。
