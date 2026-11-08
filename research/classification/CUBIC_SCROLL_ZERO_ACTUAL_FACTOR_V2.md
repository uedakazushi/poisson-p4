# Smooth cubic-scroll rank零の全空性：修正 L6＋全Segre五の actual factor

2026-11-08、PRIMARY `E_reduced_nonnormal`。修正全支持と有限全係数 evaluation の producer PASS。
**指定 smooth cubic-scroll rank零・被約既約actualquintic source は全空**という新証明を得た。
受理は別作者の正しい whole-cover と全係数 review 後に行う。
可約actualDを持つrank四source、twisted cubic cone又は他のE branchesの不存在は主張しない。

## 1. 旧誤りを保持した正しい入力

旧四源 whole-cover は g11 coefficient の誤りで撤回した。
`CUBIC_SCROLL_ZERO_FOUR_COVER_ERRATUM.md` と別作者 rejection receipt がその正確な射程を記録する。
旧 frozen report の wholeempty boolean を今回の証明として用いない。
四つの個別 source の fullJacobi／Pfaffian factor arithmetic は有効であり、
今回再利用するのは別作者も確認した **individual L6** 一つだけである。

正しい全支持は `CUBIC_SCROLL_ZERO_TWO_WHOLE_SOURCES_V2.md`。
元二十四式と正しい八式の定数両方向、全parameter零補集合、
row1非零／零による missing二minorの消滅から、
全 reduced support は L6 と全3×3 rank≤1 Segre cone五に等しい。
fu非零を保持し、旧counterfamilyも新外積へ入る。

## 2. 一回の全Segre五 coefficient evaluation

reuse runner に isolated `--smooth-two-sources-v2` を追加した。
直前のCBC257 fullbytesと全四契約を preappend snapshotへ保存し、旧四 FunctionDef AST は完全保持。
現sourceは0b7f7997…で固定。旧 report／旧 sourcebytes／契約は上書きしない。
active gate PASS後、一回の新全Segre五評価で **0.88410秒／RSS53,312KiB**。
L6は元 frozen individual case を明示再利用し、同じsourceを再評価していない。

全六外積parameter `(A,B,C,D,E,F)` に対し、元30×12 source matrixから全十 actual bracketを作った。
元二十四generator、保存全350行、原tensorから作る全projective350Jacobiは恒等零。
元signedEulerPfaffianの全126 coefficientは保存された全十二源の二次式のmapと完全一致した。
さらに以下の fullpolynomial factor product を元actualDと展開一致確認した。
finite field、point sample、minor、M1、GB、native又はjet計算は無い。

正式 producer report は `CLASSIFICATION_CUBIC_SCROLL_ZERO_TWO_SOURCES_V2_20261108.json`。
全 r-map、全 ten bracket、全126 coefficient、factor、実行時source全bytes/SHAを保存する。

## 3. 全actualfactorと全specialization

物理座標は x0,…,x4、

`q0=x0x3−x1x2, q1=x0x4−x1x3, q2=x2x4−x3²`。

### L6

全六parameter l0,…,l5 に対し、既正しい全源式は

```
H_L = l1*x3−l2*x2+l4*x4,
K_L = l2*q0²−l0*q0*q2−l1*q0*q1−l3*q1*q2−l4*q1²−l5*q2²,
F_L = H_L*K_L.
```

物理次数は一と四。元tensor全係数確認と別作者個別factor確認済み。

### 全3×3 Segre五

全outer product `(A,B,C)^t(D,E,F)` に対する元座標mapは

```
r=(CD+2AF, 2AE, AD, −CE−BF, 0, AE,
   −CF, BF, AF, BE, 2AE+BD, 2AD).
```

次の全一次形式を置く。

```
H_S = (4A²E−2ABD)*x0 + (B²D−2ABE)*x1
      +4A(AF+CD)*x2
      +(−4ABF−4ACE−2BCD)*x3
      +(B²F+2BCE)*x4.
```

全 source 上で、符号も含めて

`F_S=(1/2)*q2*(D*q0−E*q1−F*q2)*H_S`。

物理次数は二・二・一であり、これはnormal symbol又はoldF/U/S subfamilyへの制限だけの式ではない。
全六外積parameterの元actualPfaffian恒等式である。

どの式にも source parameter の反転は無い。
任意のspecialization後、positive homogeneous physical factorは零になるか同じpositive degreeを保つ。
因子零ならactualD零、非零ならproper factorを持ち、同一因子の重複も被約既約quinticを与えない。
外積は全rank≤1 matrixを覆い、rank零も含む。
したがって L6 と全Segre五の全specializationに対象既約quinticは無い。
正しい whole support と元sourceの全PGL逆により、全 smooth cubic-scroll rank零 branch が空となる。

新成分数は零。旧falsecoverを履歴から消さず、修正V2とその独立receiptを優先して読む。
