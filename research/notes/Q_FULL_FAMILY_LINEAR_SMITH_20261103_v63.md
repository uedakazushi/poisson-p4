# 全Q_a族の一次Jacobiと軌道：全パラメータのSmith証明書

2026年11月3日。Qaの全canonical126係数のJacobi一次写像と、GL5軌道写像をQQ[a]上で計算した。
選んだ109変数切断が失敗する値も、全126係数で検査している。

結論は次の表である。射影次元はPGL商を取る前のもの。

|a|全Jacobi rank|canonical affine tangent|projective tangent|projective orbit|横断接次元|
|---|---:|---:|---:|---:|---:|
|0|73|53|52|16|36|
|2/9|74|52|51|16|35|
|その他|75|51|50|16|34|

軌道のaffine rankは **全aで17**。非零aで全Jacobi rankが落ちる値は正確に2/9だけである。
親模型の群minorの9a+4因子や、選択した陰関数blockの追加因子を、内在的な軌道・接次元の跳びと扱わない。
ただしa=0ではPoisson自体の一般階数が二となり、全Jacobiも二つ分落ちる実際の現象がある。
群minorのa因子が座標frame由来であることは、この別の現象を消さない。

## 全写像の再構築

divergence零の定数有理基底Cは150×126である。
元の150個の単項式方向それぞれについて、物理的な

    Σ_l (Qa_il*∂l R_jk+R_il*∂l Qa_jk)+cyclic

の全350係数を直接再構築し、Cを掛けて350×126行列J(a)を得た。
既存producerの108列専用helperを126列へ拡張したと仮定していない。
別に全25個のGL5方向から150×25軌道行列O(a)を作った。
J*O=0、QaがPoisson、weighted diagonal方向が5Qaであることも全係数で確認した。

Qaのtorus weightは5であり、各入力の相対重みd=5−WによりJは12blockへ分かれる。
各blockをQQ[a]上のSmith形D=S*Jblock*Tへ移した。
S,Tの全係数を保存し、積の等号、diagonal性、detS,detTが非零有理定数であることを直接確認した。
したがって任意のaへの特殊化でも、同じ可逆な行・列操作で階数が決まる。
標本点のランクから例外の尽くしを推論したものではない。

全75個の非零Smith invariantのうち、非単元は定数倍を除いて

    weight6: 9a−2、weight7: a、weight8: a

の三つだけである。他の72個は全て単元。
全GL5軌道の七つの重みblockでは、17個の非零invariantが **全て単元** となる。
さらに、最初に固定する七つの非正重み係数にも、全aで単元なrank7の群frameがある。

## 意味と限界

この結果は全Q_a族の一次変形と軌道を尽くす。
有限な高次ideal・枝・平坦性は一次行列だけでは決まらない。
例えばT40のvertical成分はa=−8/9にあり、全Jacobi rankがgeneric75のままでも実全環はaについて非平坦である。

単元Smith frameは同じ重みの座標だけを混ぜる。
したがって正重み群補空間を多項式で選び、72単元の陰関数を全aで有限多項式消去できる。
この構成を実装すれば37変数模型、aを反転して二つのa方向も消せば35変数模型を得る。
**その高次構成は本証明書そのものの結論として実行済みとは扱わず、別のproducerで検証する。**

コードは `code/q_full_family_linear_smith_20261103_v63.py`。
全行列と全S,Tは `reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz`、
射程は同directoryのreport.json。
QQ[a]全係数と特殊0,2/9,−4/9,1の照合を約3.7秒で検証し、work出力を同一bytesでpublicへ保存した。
