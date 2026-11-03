# actual42 Euler delta28の独立検証

2026-11-03。`S61_AFFINE_ACTUAL42_EULER_DELTA28_20261103_v63` を監査する。
producerは再実行せず、新しい算術方法で次数28と正の非零証明を確認する。
全42変数の四次・五次多項式展開、membership探索、GBは行わない。

## exact Euler簡約と実代数

既証明の全QQ cross恒等式 `z_k*n_j=P_kj(n)` から、正の次数dの斉次成分fについて

    z_k*f = sum_j P_kj*partial_j(f)/d

が元の商代数内で成立する。Kには定数項がなく、Pもnormal原点で零なので、各段階で
正次数を保持する。homogeneous成分へ分けて外側のoperator順を固定すれば、有限のQQ
多項式recipeとなる。plane cofactor恒等式に代入したS1,S2,Hはそれぞれ次数4,4,5以下であり、
既承認のglobal quinticの二つの逆写像にそのまま使用できる。旧凍結Redの値を改変する必要はない。
新しいEuler版は同じactual Aと同じnormal mapを表す。R42上の自由性は仮定しない。

## 次数28の上界と最高部分

全134 remainderを読み、normal重み1、S1/S2重み4、H重み5を直接適用する。
h0係数の重みは8以下、h1–h3は5以下、h4は4以下でHを含まない。
一般betaでも、元式の二次性・代入のh次数4以下・monic quintic divisionにより同じ型となる。
companion bottom4 determinantを全24 permutation項の和として独立に展開し、最高重みが28、
最高部分が正確に

    H² L3² L4² − H L3³ Q − L4 Q³

であることを確認する。自由formal45での最高次数7だけを代入して実次数を推測しない。

## 根付き三階微分とQQへの正の移行

QQ cross inverseの全15876係数、Pの全119070係数を保存されたexact出力から独立に読み、
分母が32003のunitであることとPの全有限体係数への一致を確認する。
このexact出力の両側inverse恒等式と原134式へのbindingは既承認の整数監査を依存証拠とし、
その全入力hashを照合する。三つのplane cofactor恒等式自体もQQで再確認し、全分母を確認する。

今回の根付き評価はproducerのHessian/Jacobian公式を使わない。
`u²=v²=w²=0` の八次元係数環で、二次vector fieldの三つの一次flowを順番に合成する。
`f(Flow_c(w,Flow_b(v,Flow_a(u,e0))))` のuvw係数が `D_a D_b D_c f(e0)` である。
全六組の二階値・全十組の三階値をこれで求め、Euler分母6と24を適用する。

    S1_4=15338, S2_4=11727, H5=31338,
    Q8=30427, L3_5=13746, L4_4=2943,
    delta28=31492 != 0 (mod32003).

全演算は既にQQで定義されたp-integralな有限circuitの像であり、この**非零**はQQ上の
実最高係数の非零を証明する。a=0の係数が非零なのでgeneric moment parameterでも非零。
一般betaへの移行はp-unitの有限inverse/cofactor recipe上の非零特殊化による。
有限体での零・unit・空性をQQへ移す議論ではない。巨大なQQ値の展開保存も主張しない。

## scalar criterionと限界

既承認のnormal idealは `J0=m intersect I_B=ker(R42->A)`、
生成元数は3+2583+210+42=2838、全次数4以下である。
Euler版でもnormal mapは同じ。delta_Eは原点で零であり、既承認のsharp moment議論は
このmonic quintic presentationにも適用され、B(a)上でregularとなる。
従って `B=0 iff delta_E belongs J0` は正しい。
逆向きにはmembershipがdelta_EによるBのannihilationを意味し、regular性からB=0となる。

targetのnormal次数は正確に28である。従って生成元次数4以下の多項式cofactor表示があるなら、
少なくとも一つのcofactorのnormal次数が24以上となる。この下界を承認する。
cofactorの上界、次数24以上の非存在、actual membership、残余の空性、generic五点尽くしは
いずれも証明していない。

独立コードは `code/s61_affine_actual42_euler_delta28_review_20261103_v63.py`。
新しい出力先を指定して既設Cloud環境で実行する。全入力hashと独立flow証拠は
`reports/S61_AFFINE_ACTUAL42_EULER_DELTA28_20261103_v63_REVIEW/` に保存する。
