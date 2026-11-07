# 本質三変数既約五次の全八23次元成分：全源と最弱支持

2026-11-07、PRIMARY E。全八curve/sourceのproducer算術と、残七canonical源の全周囲算術はPASSした。
最初の(4,5)源は既採用の104番目である。本稿の残七源の採録は独立reviewと統合担当の受理を区別する。
全P4分類の完了や一般五点の尽くしは主張しない。

## 対象の全支持

実Dが既約被約・本質変数三のrank四Poisson構造の全支持は、
`ESSENTIAL_THREE_QUINTIC_FINITE_CURVES.md` の全curve分類と
`ESSENTIAL_THREE_QUINTIC_EIGHT_SOURCES.md` の全Saito算術から、八つのactual PGL5軌道で尽くされる。
全磁気inverseによって各固定Fの全frame・六磁気項を一軌道へ移せる。
各軌道は23次元であり、本稿はその全閉包がそれぞれ全Poisson空間の成分になることを証明する。

| curve型 | 全D階数 | 実作用M階数 | true normal | 最弱支持 |
|---|---:|---:|---|---|
| (4,5) |101|48|−1一本|既済の非零二次、再計算しない|
| A12 |102|48|なし|smooth23|
| E8+A4 |101|48|−1一本|非零二次|
| E6+A6 |101|48|−1一本|非零二次|
| E6+A4+A2 |101|48|−1一本|非零二次|
| A8+A4 |102|48|なし|smooth23|
| 3A4 |102|48|なし|smooth23|
| A6+3A2 |102|48|なし|smooth23|

各行はその全源が一actual軌道であるという全source theoremから選ぶ唯一のcanonical anchorである。
恣意的な点・素数の標本を増やした検査ではない。
三次以上の障害や正確な非被約局所環は計算していない。

## 全150/350と真正切片

各新FのSaito核のstandard frame、磁気R=0を取り、全tensorの分母を整数scaleで払った。
全350×150 projective Jacobi微分D、Euler25列、GL5二十五列、scalar一列のMを保持する。
全DM=0、全階数minor、全kernel/actual actionのweight blockを保存した。

C*の座標重みは全源で(0,0,0,1,1)、tensor character+1であり、
g_t pushforward/tは正当なprojective stabilizer作用である。
全Dがrelative weightsを保つことを全非零成分で確認した。
実Mの独立48方向に非零な48coordinate cutを選び、その上の全normalを求め、

    C150=actual48 ⊕ true normals ⊕ independent implicit coordinates

の全直和determinantを非零と確認した。
これはEuler quotientや選択部分空間から全支持を推測したnormalではない。

最初の(4,5)点の全Dと二次は旧不変reportをそのまま参照する。
残る七点だけを新たに計算し、各全150/350行列は公開report内にある。

## smooth四源

A12、A8+A4、3A4、A6+3A2ではrankD=102なので、全raw接次元は48である。
実orbit/vertical/scalar源が48次元の下界を与えるため、局所次元も48である。
全Poisson schemeはこの点でregular、射影化とEulerを戻すとsmooth23である。
actual23軌道閉包はその唯一の全空間成分と一致する。

## 一法線の三源

E8+A4、E6+A6、E6+A4+A2ではrankD=101で、真正48cut上のnormalは一本である。
選択101Jacobi式のimplicit blockと全exact inverseを保存した。
normal Rの全二次Jacobi Qから一意のimplicit補正T2を求め、全350残差DT2+Qを検査した。

各caseにλD=0を**全150列で**満たすfunctionalがあり、λQは次の通りである。

| curve型 | λQ |
|---|---:|
| E8+A4 |−47546093/13602816|
| E6+A6 |−913230352584/7808230609235|
| E6+A4+A2 |−3431736/3246473|

全て非零である。従って真の一変数IFT残式は非零二次初項を持ち、
被約切片の局所支持は原点だけである。
48linear cutのKrull不等式から全raw局所次元は48以下、実源48と一致する。
よってそれぞれの全23軌道閉包は全空間成分である。
PGL5が連結で全成分を保つこととactualorbit23から、点を通る被約成分の唯一性も従う。
schemeのsmoothnessをこの三源には主張しない。

## prior104全成分との相違

比較入力はimmutable commit
`6ab28a327b64543ba0a7bf8cd302bea9def517e6`、COMPONENTS blob
`f7ec76f1b2bb1850f2e3314e36f673f218558c1c` の全104行である。
最小projection `inputs/essential_three_quintic_prior104_20261107.json` に固定した。

最初の103成分の比較は、独立受理済み
`CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_FIRST_PRIOR103_20261107.json` を再利用する。
その全103行の理由のうち必要仮定は、今回の各全源にも共通する次の二つだけである。

1. 実Dは絶対既約である。
2. 全source閉包は実Dの偏微分spanが三以下という閉determinantal条件に入る。

従って異rank・異次元は自動的に異なり、同次元MPS二十九成分と既存二つの既約五次成分は
保存済みのactual essential-rank五の実memberと区別される。
他の同次元成分は受理済みgeneric因子分割のproper closed multiplication像に入るため、絶対既約Dと異なる。
旧actualSNCや旧二vectorは今回再計算せず、その全証明書・独立receiptをそのまま保持する。

104番目 `R4_PLANE_UNICUSPIDAL_QUINTIC_45_23` は、一つの重複度四cuspを持つcurveのactual23軌道閉包である。
残七curveの最大尖点重複度は三以下なので、PGL3同値ではない。
本質三変数の部分空間は実Dから内在的に一意に復元される。
従ってPGL5でPoisson構造が共役なら、その平面curveもPGL3同値でなければならない。
双方がdense actual23軌道を持つ同次元成分なので、異なるcurve orbitを持つ同一成分は存在しない。

## 八源相互の区別

最大重複度、尖点数、各尖点の型による八つのconfigurationは全て異なる。
全28組について、このintrinsic curve型とdense23軌道の理由をreportへ記録した。
同一のirreducible23閉包に二つの異なる23軌道があれば、どちらもdense openとなって交わり矛盾する。
従って全八軌道閉包は相互に異なる。
閉包のboundary同士が交わらないとは主張しない。

## 証明状態と再現

- 全source: `CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_EIGHT_SOURCES_20261107.json`。
- 全七new ambient/最弱二次/条件付き全比較:
  `CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_EIGHT_AMBIENT_20261107.json`。
- 契約: `jobs/essential_three_quintic_eight_ambient.json`。
- runner: `code/classification_plane_quintic_saito.py` のisolated `--eight-ambient-input` mode。

全周囲batchはgate後約37.506秒、RSS約75MiBでPASSした。
各実行producer bytesはreportに固定し、以前の全functionはmain以外AST不変である。
derived sourceは公開八源reportのJSON pointerからinlineで作り、workファイルを入力にしていない。
全curve分類・全source・全component算術の独立review結果は別receiptへ接続する。
台帳の数と新canonical IDsは統合担当だけが更新する。

この結果は本質三変数・既約被約五次という全subtypeのE帰属を閉じるもので、
本質四・五変数の既約非正規五次、他の未帰属source、全P4分類の尽くしは残る。
