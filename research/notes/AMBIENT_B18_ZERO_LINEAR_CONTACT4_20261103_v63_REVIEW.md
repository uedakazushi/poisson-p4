# B18：N の零一次関係と R 曲線の接触 h=4

2026-11-03。独立検証。generic B18 の既受理な同時 G19 分解と、実際の被約 N、R を扱う。
本稿は **h=4、3≤k≤4** を証明する。k の3と4の選択、実際の三次ジェット、
全 D 完成環の同定は行わない。

## 入力と独立算術

`AMBIENT_B18_N_R_SECOND_CONTACT_20261103_v63` の source は、全次数で
固定 N19 空間の条件 C2=0 を保つ実際の N の形式曲線である。その正規化 source
の28次元 étale 性は既受理の `ODD_FIRSTJETS_N_PROJECTION` に含まれる。
`N_FULL_COMPLETION` の全環と全81関係、`NILRADICAL_LINEAR_VISIBILITY` の
DVR 被約和、`ACTUAL_INTERSECTION_ORBIT_PRODUCT` の同時 slice を使用する。

今回の検証コードは producer を再実行せず、通常の有理数と物理150係数の添字演算だけを使う。
保存された実際の R の多項式曲線から R0,R1,R2 を別に抽出する。
固定 source の速度を F、外側の GL 速度を M、補償された kernel 速度を K とすると、
添字から独立に作った `L_M(P)=MP+PM^T−(Mx)·∇P` により全150係数で

    N20=½ L_M²(P0)+L_M(F),  H=L_K(F),  W=−½ L_K²(P0)

を検査する。これは新しい三次ジェットではない。
J27 の保存逆行列も掛け戻し、Δ=R2−N20、H、W の剰余の全150成分を検査した。
結果は

    Δ = A W + B H  modulo J27,
    A=4800/23, B=300/23, B²−A=−20400/529≠0.

W,H の剰余の行27,61の小行列式は −138 である。従って H は J27 と W の和に入らない。
Δ が実際の N の37次元接空間に入ることも保存した全150行で確認する。
この最後の事実により、任意の N 関係の R への制限は少なくとも三次となる。

## 正規化の trace と二次係数の対応

既受理の **実際の** N 完成環は

    R=C[[g19,ν,s1,...,s6,w,v]], I=(ν,s1,...,s6,wv,v²),
    A_N=R⊕I t ⊂ B=R[t]/(t²−ν).

このモデルでは元の判別式 ν を保って座標を持ち上げている。
正規化 B は正則28次元で、原点で dν=0 である。
source 微分 J27 は even R の `dν=0` 超平面を張り、odd 出力 `Y_i=b_i t`
の微分は全て零である。従って source の一次 kernel k 上で dt(k)≠0。
非零定数で t と ν=t² を同時に正規化し dt(k)=1 としてよい。
これは ν を単元倍する座標変更に過ぎず、I と全像 A_N を変えない。
古い R の座標を新しい別成分の R と同一視する議論は使わない。

baseline の N 曲線の trace 一次係数を d とし、even 係数 b_i の一次係数を l_i とする。
R 曲線の一次方向は baseline と等しいのでこれらは共通である。
ν の一次係数は零、core `wv,v²` の一次係数も零である。
baseline では

    ν2=d²,   (Y_i)2=l_i d.

source の kernel を c 倍加えると trace の一次係数は d+c となる。
even source 一次速度は変わらないため、pure kernel 曲率 W の
`(ν2,Y2)` は `(1,0)`、混合曲率 H は `(2d,l)` となる。
source の二次加速度 J27 はこれらの出力には寄与しない。
物理的な Δ の分解から実際の R 曲線では

    ν2=d²+A+2dB,   (Y_i)2=l_i(d+B).

even の接空間は J27 と W で尽くされる。従って H がその和に入らないことは
odd 成分 l≠0 を意味する。ν と二つの core の l は零なので、ある s_i の l_i が非零である。
この記述は trace の初期座標に依存しない。dt(k)=q のままなら後の係数全体が
非零因子 q² 倍されるだけである。

その i に対する実際の関係

    Y_i²−ν b_i²

は ambient 一次項を持たず、R 上の四次係数は

    l_i²[(d+B)²−(d²+A+2dB)] = l_i²(B²−A) ≠0.

従って零一次関係の最小制限次数 h は高々4である。
単なる「正規化に二次持ち上がらない」という事実からこの結論を出していない。

## 全関係について h≥4

最小 ambient は R と九つの Y の37変数である。全 kernel は、
36本の Taylor pair syzygy と45本の `Y_iY_j−ν b_i b_j` で生成される。
独立コードは81本全てを lcm から再構築し、保存済みの全 kernel と係数ごとに比較する。

R 曲線上では、七つの regular b は次数≥1、二つの core b は次数≥2である。
全 Y は次数≥2、core Y は二次係数が零なので次数≥3である。
七 regular の pair syzygy の三次係数は
`l_j l_i(d+B)−l_i l_j(d+B)=0`。
regular/core の14本は次数≥4、短い core/core の一本
`vY_wv−wY_v²` も次数≥4である。全45 product も次数≥4である。
これで全81を覆う。検証コードは一般記号の有限多項式演算でも、この三次以下の消滅と
上記四次係数の d の相殺を確認する。未知の三次係数を零と推測する計算ではなく、
それらは各係数に必ず一次以上の因子が掛かるため四次以降にしか寄与しない。

大きい実 ambient に戻すと、独立な一次項を持つ graph 方程式 G_a が加わる。
全 I_N は `(G_a,81関係)` で生成される。
`F∈I_N∩m_ambient²` の表示で G_a の係数の定数項は全て零でなければならない。
一方、上の曲率により **全** G_a の R 制限は次数≥3である。
従って graph の寄与 `m_ambient*(G_a)` も次数≥4である。
これは graph 関係同士の一次項の相殺や非線形係数を含めた議論である。
よって全零一次関係が次数≥4となり、**h=4** が従う。

## k と被約接空間に関する残る二択

同時 G19 slice 上で被約和は

    F=A_N,slice ×_{C[[T]]/(T^k)} C[[T]],   edim(A_N,slice)=18.

actual D slice の edim は20である。曲率から k≥3、上の関係から k≤4である。
R 曲線は orbit に横断する一次速度を持つため T はこの DVR の正則 parameter である。

F→A_N,slice の kernel は κ=(0,T^k) で生成される主 ideal である。
κ が m_F² に入ることと、zero-linear N 関係が R 上でちょうど次数 k を持つことは同値である。
前向きには κ を ambient m² の元で持ち上げる。
逆向きにはその関係の R 制限 `T^k*unit` の単元を F→C[[T]] の全射で持ち上げて除く。
従って κ∈m_F² ⇔ h=k。

以上から、まだ決まっていない二択は

* k=3：被約 slice の edim=19、actual slice の線形可視な nil 次元=1。
* k=4：被約 slice の edim=18、actual slice の線形可視な nil 次元=2。

である。cone の被約 edim はそれぞれ38、37である。
この稿はどちらが実際に起きるかを決定せず、具体的な物理座標の四次式も展開していない。
主張は形式的な generic B18 open orbit に限り、全特殊点や完全な D の環型には拡張しない。

検証コード：`code/ambient_B18_zero_linear_contact4_review_20261103_v63.py`。
flat receipt：`reports/AMBIENT_B18_ZERO_LINEAR_CONTACT4_20261103_v63_REVIEW.json`。
