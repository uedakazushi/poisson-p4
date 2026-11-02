# 自由 M0 旗端点：横断原方程式39本と相対 source reduction

2026-11-02。指定された従来の有限体点、係数体 k=F_32003 に限定する。
通常の正規化五変数 slice の一般点についての結論ではない。
原135式の形式的陰関数消去後の環を S、法線変数を x_0,...,x_7、
底変数を u=(v_8,v_9,v_10)、a=(u)、極大 ideal を m とする。
既証明の S/a は長141、最高非零次数7である。

## 1. 全 r の有限 cutoff の改善

既存 `S61_MINIMAL_FREE_FLAG_R_ADIC_CUTOFF_20261102` は、初期下側 ideal J10
について (k[v]/J10)/(u) の Hilbert 関数
(1,8,24,38,38,23,8,3,1)、長144、最高次数8を確認した。
従って gr_m(S) は k[u] 上次数8以下で生成される。
特に m^9=a m^8+m^10 である。有限 S 加群
N=m^9/(a m^8) は N=mN を満たすから、Nakayama により
m^9=a m^8。帰納的に

    m^(r+7)=a^(r-1)m^8,  r>=1.

別の既証明である actual 横断商の最高次数7から m^8⊂a を使うと、
全 r>=1 について m^(r+7)⊂a^r を得る。
J10 の長144だけから m^8⊂a を推論していないことが重要である。
従って S/a^4 は通常次数10以下で厳密に計算できる。
この改善は同じ指定有限体点での主張である。

## 2. 原135式から39式を選ぶ

既存の actual A3 商の八つの乗法行列を底0へ特殊化し、A0=S/a の
Koszul 複体を別に作る。d1 の階数140、d2 の階数949なので

    dim H1(x;A0)=8*141-140-949=39.

正則形式環 k[[x_0,...,x_7]] において、これは横断 ideal I0 の
最小生成元数 dim I0/(x)I0 である。
原135式の八次 jet は元の符号付き135式と固定接枠から fresh に再生成する。
各 Fi を、各単項式の最初の変数で割って Fi=sum x_j r_ij と書き、
(r_ij mod I0) を Koszul cycle として計算する。
135個の class の階数も39で、保存された39行が基底をなす。
八次以後を捨てても class が変わらないのは、m_x^8⊂I0 より
m_x^9⊂m_x I0 だからである。Nakayama によりこの39個の**実際の形式式**が
I0 を生成する。

相対環で残り96式を捨ててよいとは主張しない。それらの底次数正の
残差は求める torsion relations の一部である。

## 3. Source DAG と89個の monic normal rewrites

選択39式の八次切断に、既存の local standard-basis 演算を適用する。
167個の非零行が得られ、最小先頭単項式は89個、その staircase は
既存 A0 の141基底と完全一致する。各非零行を追加する際、元式または
S-pair、実際の各 reduction、最後の非零規格化係数を全て DAG に記録する。
別の native verifier が 1,861,027 reductions を再実行し、全167行の
全八次係数を比較する。任意の pair criterion の完全性を用いて、元式への
所属を推論する必要はない。

各 DAG を八次切断せず実際の相対形式式 Fi(u,x) に適用することにより、
89個の形式式 h_b が定義される。その source coefficient matrix A(x)
は89×39で、h=A(x)F_selected である。
独立 replay は A(0) も計算し、その39×39 minor が9252 mod32003と
非零であることを確認する。従ってこの minor は k[[u,x]] で単位であり、

    (h_1,...,h_89)=(F_selected1,...,F_selected39)

が**相対形式環でも**成立する。h_b の底0八次先頭は89個の monic normal
border である。高次項の切捨てを相対 ideal の等号へ取り違えていない。
残る96原式を加えると元の相対 ideal になる。

これにより、次の base-degree3 計算を、出所不明の abstract deformation
ではなく元135式に戻せる固定 normal rewrite によって進められる。
本証明書はその有限入力を確定する。全相対 presentation の計算は未完了である。

## 再現と範囲

    .venv/bin/python code/s61_minimal_free_flag_transverse_sources_20261102.py \
      --output reports/S61_MINIMAL_FREE_FLAG_TRANSVERSE_SOURCES_20261102.json

二つの native helper を一時ディレクトリにコンパイルし、原八次 jet、
Koszul class、39式 standard-basis source、独立 DAG replay、単位 minor を再生成する。
生成された大きい trace は保存済みデータを入力せず再構成され、報告には digest を記録する。
入力 monomial 順は既存 producer の total-degree/CWR 順に固定されている。

係数体0への lift、全 S の t 冪零、相対 ambient 飽和除去、全体 exhaustion は
いずれも主張しない。既存 A3 商での t²=0 とも範囲を区別する。
