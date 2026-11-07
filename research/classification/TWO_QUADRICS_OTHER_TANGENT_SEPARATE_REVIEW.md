# 四重接触conic対・別の接線：別作者の三次支持と相違の確認

2026-11-07。`E_reduced_nonnormal`。
対象はcontact4_other_tangentの全六磁気源である。
sourceの完全性は `TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md` の全Saito還元と
その既存独立確認を入力とする。局所唯一性や正確な非被約厚さを結論しない。

独立コードは `code/classification_plane_conic_cubic_separate_review.py`、証明書は
`CLASSIFICATION_TWO_QUADRICS_PLANE_OTHER_TANGENT_SEPARATE_REVIEW_20261107.json`。
producer又はその作者のcubic checker、CASをimportしていない。
既存の別作者fullmagnetic再構成と、以前の独立Fraction/S61係数loop算術だけを再利用した。
実行3.121秒、最大RSS18,816 KiB、メモリ上限1 GiB。

## 1. 全実sourceと全周囲

conicは `q1=hx-y²,q2=hx+x²-y²`。両方とも滑らかで相異なる。
全六磁気parameterを残した全raw Jacobiと、全parameter上の実Pfaffian恒等式を再構成した。
六磁気方向は表示された全sourceを実際に微分して作り、保存値と一致した。

全150二次係数から全350projective Jacobi行を独立生成した。
非零101minor、実source方向48、全核49=F48+R、全D・核積の零を確認した。
従って元の全空間のD階数は正確に101である。
一法線Rだけを持つことは、小さいsliceだけの接空間計算から仮定していない。

## 2. 真正sliceにおける固定Tと全三次

101minorのstandard座標列が張るUとRが作る102次元空間Sは、実F48の線形補空間。
これは元の全150係数空間における真正な48本の線形cutである。
選択した101行の陰関数定理に対応する二次補正Tを、Uの中でFraction消去により
独立に一意に解き、保存されたTの全150成分と一致した。

`Q=Jac_projective(R)` とすると、選択行だけでなく全350行で

    D*T+Q=0

を確認した。従って元の全Jacobiはこの補正で二次まで一致する。
次に `P+tR+t²T` の全350三次係数Cを、疎多項式の実展開
`Jac(R+T)-Jac(R)-Jac(T)` により再構成し、保存された全成分と一致した。
実350行上の関数ϕについて

    ϕD=0 （全150列）, ϕQ=0,
    ϕC=-533355840000 != 0

を正確に確認した。
さらにFの独立48基底とRからなる**全49次元核**について `ϕ[R,K]=0` を検証した。
したがって二次補正を任意の全核方向だけ変えても、この三次pairingは変わらない。
固定されたS内の補正を使う証明にはこの追加不変性は不要だが、補正曖昧さも直接監査できた。

三次以降のUの補正はDの像を加えるだけで、非零三次pairingを消せない。
従って真のslice残式には非零t³初項があり、その被約支持は原点、次元は零。
Krullの不等式で全raw周囲の局所次元は48以下となり、実sourceの下界48と一致する。
全六磁気源のPGL・垂直・scalar飽和像閉包は既約でwitnessを含むので、
その像閉包自体が全空間成分となる。射影次元は `48-24-1=23`。
元の全六磁気源の全点がこの成分へ入る。

## 3. 基準90および受理済みbitangent成分との相違

凍結した基準90projectionと、既に別作者受理済みのcommon-bitangent成分を合わせ、
比較対象は論理的に91成分である。91の台帳更新やGit保存が済んだと仮定しない。
基準90との比較はbitangent別作者監査と同じで、rank、次元、本質変数数三以下の
閉条件によるSNC不在、残る一般因子型の相違で全てを区別できる。
特にCPLUS_EXTは24次元だから、23次元の本成分とは次元だけで異なる。

二つの新成分は、ともに一般因子型が被約1+2+2で本質変数数は正確に三。
したがって、退化式の偏微分を消す二次元の空間から共通vertex P¹と平面商が
一意に復元される。唯一の一次因子Hと、二つのunorderedな既約conicも固有に決まる。
平面conic対の交点divisorの重複度分割はPGL不変である。

本源のq1の正規化を

    (h,x,y)=(s²,t²,st)

とすれば `q2|q1=t⁴`。すなわち交点は一点で重複度4。
bitangent源のconicを前の別作者noteの座標で書き、そのq1を

    (h,x,y)=(s²,t²,t(s-t)/4)

と正規化すると `q2|q1=t³(s-t)/4`、相異なる二点で重複度3と1となる。
両restrictionをFraction疎多項式で実際に確認した。

滑らかで相異なるconic対の開集合上では、restriction binary quarticが四乗になる条件は
Veronese像という閉条件である。同値に、正規化上の一点でcontact order≥4となる
proper incidenceの像として定義できる。因子を交換しても条件は変わらず、
上の一意な平面商からintrinsicに復元される。
本源の全一般像はこのcontact4条件に入り、bitangentの全一般像は外にある。
二つの像閉包が同じ成分なら、それぞれの稠密constructible像が含む非空開が交わり、
不可能な二つの交点分割を同時に持つことになる。従って成分は相異なる。

これで本成分も追加受理でき、両新成分を合わせた下界92、一般rank四86が得られる。
他の三平面配置、非平面1+2+2、全分類の完了は主張しない。
