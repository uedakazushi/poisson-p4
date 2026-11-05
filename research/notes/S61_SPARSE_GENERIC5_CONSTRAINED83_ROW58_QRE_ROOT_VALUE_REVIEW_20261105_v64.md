# 固定 row58/col56 の正確な既知単元指数と第10段候補

受理済み第9段 col56 の実際のrow58分子 N9 は951項である。固定Q=5q0+1、R=q0+1、E=32s0+15q0+35に対し、**N9=Q^4 R^2 E F10/(353894400000 H t0)** が正確に成り立つ。F10は317項のprimitive整数多項式であり、真のLaurent単項式H^(-1)t0^(-1)以外に単項式contentはない。

主計算d1f1940fは10個の整数擬除算の商と剰余を保存した。独立agentであるROOTの検査fff391b8は、除算engineもproducerのcodecもimportせず、独自のstream decoderで317項を読み、Q^4 R^2 Eを直接掛けて**元の951係数すべて**を復元した。さらにF10へq0=−1/5、q0=−1、s0=−(15q0+35)/32を記号代入し、それぞれ52・69・109項の非零多項式を得た。したがって上記の3指数は正確である。別途81c59418は同一agentの異なるQQ直接積で全10個の保存恒等式と終端剰余を照合しており、独立agent検証とは明確に区別する。

入力N9のbinary SHAは5774ebc24c84576c4d7b21f5e245644a4e1a6ad224447d9ab82ec0d1bf2a035b。F10のbinary SHAは526ec2470b0ceb344ac7a837e9e45e3bfbf99e60e0a646831c5a2ae6ad917a68、整数value SHAは31a9d7277afd02c1a52e5e4f50d0b0394eea2bddef5cea0aabca2f52839644a7である。全試行の商・非零剰余と元データを保存した。

共通表現D9=Q^4 R^2 E F9から、現在の実際のpivotは

p10 = F10/(353894400000 H t0 F9)

となる。F9は除算していない。F10の既約性、F9との互いに素、F10=0や既存chartの補集合の被覆は未証明である。第10段の追加chartはD(F10)だけと記述できるが、第10段の他の行列係数はまだ計算していない。全168×15、87追加RHS、66制約行、P290一回と全復元式を保持し、全体逆行列やexhaustionを主張しない。

独立検証reportは reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_QRE_ROOT_VALUE_REVIEW_20261105_v64/report.json、SHA fff391b83d5f1ff72549982310a7d37a04a875c4d4c00936f1d6c639a79bf4f0。ROOTのSOURCE/code gateはa03e3d78、元951項の別QQ検証は788daf79である。旧凍結report/code/noteには変更を加えていない。
