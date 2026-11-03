# 実axis polarizationと普遍G6 first-w-colonの独立監査

2026-11-03。`S61_AFFINE_AXIS_POLARIZATION_SPLIT` と `S61_AFFINE_G6_PARTIAL_COLON` の二packetを承認する。
これは実ideal操作による二段階のw-torsion処理であり、全w飽和や未知アフィンの除去を完了したという主張ではない。

## 全134 quadratic chartと実idempotent

元S61の全134式をweight-two q次数2へhomogenizeし、保存axisの非零q座標を固定する45変数chartを独立に再構築した。
各式を対称二次行列 `x^T A_i x+b_i^T x+c_i` に分け、midpoint rowを
`A_i x+A_i p+b_i` として作った。作者のderivative midpoint実装は呼ばない。
全134でF_i(p)=0と `F_i(x)=M_i(x)(x-p)` を全QQ係数で確認し、保存45×45行列と比較した。
元134の選択45行だけからのdeterminantは厳密に

    d=-124102513209883356037120000000000000000/59049 != 0.

e_p=det(M_selected(x))/dとするとadjugate恒等式からe_p(x_j-p_j)が元ideal Iに入る。
e_p(p)=1なので評価差の因数分解からe_p²-e_pもIに入り、p_w=0からe_p*wもIに入る。
I+(1-e_p)は正確に点ideal(x-p)である。従って

    I=(I+(1-e_p)) intersect (I+(e_p)),  sum=O

という実scheme分解となる。D(w)ではe_p=0なので補idealは元Iと同じであり、未知の点・正次元成分・nilpotentsを保つ。
determinantの巨大展開を計算したと主張しない。

## 選択finite etale factorへのtrace

先行する独立RAW_PURITY_SPLIT reviewは、元134のglobal45minor Delta、Eulerのradial関係、
`(I+(w)):Delta^infinity` がgeneric K上で選択するclopenなfinite etale因子を検査済みである。
Delta!=0ならminorで削除したq座標も非零なので、この選択因子全体を同じchartで扱える。
必要なbeta-only normを逆にして、普遍点pにおけるd(p)をfinite etale代数Eの単元にする。

分離閉包上で異なる二点p,qのe_p e_qは、ある非零座標差p_j-q_jによって零になる。
従って各e_pは互いに直交し、その和 `Tr_(E/K)(e_p)` は選択因子全体のidempotentになる。
wとの積も零である。chart外へ零でglueでき、先行するannihilator補sheafと同じ分解を与える。
これはgeneric Kまたは十分縮めた共通beta開上の有限recipeである。
**Eの次数・basis・trace係数は未計算**であり、全特殊fiberでの一様分解とは主張しない。

## 全61変数の実w-colon恒等式

保存された125項のq一次cofactorから、元70本weight4式を用いてG6を独立に再構築した。
全61変数QQ上で保存9257項と一致し、実q次数はG6が2、cofactorが1である。
両辺を別々にhomogenizeし、全係数で

    w G6^h = sum_(i<70) K_i(q) F_i^h

を確認した。従って実idealとして `I134 subset I135 subset I134:w` であり、両者のw飽和は同じである。
modpやbeta標本を一般恒等式へ外挿していない。

idempotentは任意のquotientへdescendする。従ってG6を加えてからその像で選択axis因子を除いても、
D(w)の全schemeとw-saturated closureは保たれる。保存minimal axisではg7=0で、既承認の全135 source恒等式によりG6最高式も零である。
ただしDeltaが他の孤立因子も選ぶ場合、それら全てでG6が零だとは仮定しない。
G6非零の選択因子が先に消えても、残りのidempotentによるcompositionは正しい。

## 同じquintic表示への正確な五列追加

元9257項から同じT01 fiberへ独立に代入・中心化し、保存790項と照合した。
全45qを実frameとformal46変数(h,n42,S1,S2,H)へ代入し、独立のflat polynomial divisionで
保存quotientとremainderに全係数一致することを確認した。
さらに元cofactorK_iの同じ代入を使い、`Rem_G(sum K_i*r_i)=r_G6` を直接確認した。
従ってこのremainderはformal45環でも旧134のmonic quotient idealの帰結である。
h^j倍の五つの正確なremainder(j=0..4)を追加証拠として保存した。
5×675表示は同じ実代数を与え、新しいPoisson条件やrank上界を加えていない。

全w-torsion、colon GB、trace basis、残余のunit、一般五点尽くしは未解決である。
review code: `code/s61_affine_w_torsion_pair_review_20261103_v63.py`。
receipt/追加五列: `reports/S61_AFFINE_W_TORSION_PAIR_20261103_v63_REVIEW/`。
