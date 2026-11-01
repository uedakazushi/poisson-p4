# AB 正規化環の次数三生成性の独立点検

2026年11月1日。[主証明](AB_ORIGIN_DEGREE_THREE_20261101.md)を別担当者が点検し、
D3=B3、dim D3=6691、H1(E3)=0 を受け入れる。全次数の生成性・原点の正規性は含まない。

## 1. Filtration と cohomology の独立再構成

proper birational な正規 incidence と、vertical 座標による P-stable filtration の使用は
[二次の独立点検](AB_ORIGIN_DEGREE_TWO_REVIEW_20261101.md)と同じである。
同点検で一次資料本文を確認した
[Fonarev, Theorem 3.2](https://www.mathnet.ru/php/archive.phtml?jrnid=sm&option_lang=eng&paperid=9880&wshow=paper)
の U* への双対 weight 変換を用いる。

別コードは主コードの Gröbner 基底を使わず、core の次数二の9個の minor と、
その8座標倍の有理線形包から次数0〜3の weight ごとの商を計算した。
関係空間の次元は順に0、0、9、56で、core の商次元は1、8、27、64となる。
これを vertical の対称冪と掛け、GL2 分解と Bott を独立に再計算した。
結果は h0(gr E3)=7285、h1(gr E3)=594。
唯一の H1 型は W(9,1) と W(8,0) の各一コピーで、H0 における重複数は一、二となり、主証明と一致した。

## 2. P-cyclicity と評価の重複数制約

P の Levi と unipotent radical は、選んだ lowest vector に対する反対 Borel を含む。
したがって既約表現は lowest vector から P によって生成される。
W(9,1) の場合、y8(y∧w) は Sym8(V)⊗Λ²0(V) の Cartan 成分の lowest vector である。
幾何学的には、旗 (line y, plane〈y,w〉) の P 軌道は、q と横断的な Lagrangian 旗の開集合を含む。
任意の transverse 平面を symmetric shear で作り、その中の任意の線を Levi で動かせる。

この議論を別コードでも確認した。q=(x,z)、p=(y,w) 順で P の7個の行列を構成し、
Sym8 の y8 から165次元、Sym8⊗Λ² の y8(y∧w) から429次元の巡回空間を厳密に得た。
後者は W(9,1) の Weyl 次元429に一致する。

非零 G-部分加群の固定平面での評価は零写像にならない。
もし零なら、G の推移性と部分加群の G-安定性によって全点で零となるからである。
従って G 型の H0(E3) 内の重複数は P-map 空間の次元で抑えられ、
cyclicity によって P-map は lowest vector の像だけで決まる。

## 3. Fiber における必要条件

W(9,1) の weight (−9,−1) の候補は h0g0² のみである。
w→w+uz は y8(y∧w) の二階微分を零にするが、候補の二階微分は288h0³となる。
h0³ は core 上で非零なので、候補の係数は零。従って W(9,1) の H0 重複数は零である。

Sym8 の weight (−8,0) には h0²g1、h0h1g0、h0s0g0 の三候補がある。
独立コードは rank-one trace-zero chart 上でこれらの独立性を確認した。
二つの shear による annihilator 条件も同じ chart 上で階数二となり、
候補係数は (1,−1/2,1/24) の直線に限られる。
ここから得るのは P-map 次元の上界一であり、これらの必要条件だけから
新しい生成元や P-map の存在を主張していない。

## 4. Euler 差と係数下界

associated graded の cohomology が H0,H1 に限られるので、E3 にも H2 以上はない。
H1 に現れ得る型は上の二つだけである。
W(9,1) の Euler 重複数は0、H0 重複数上界は0なので H1 も0。
W(8,0) の Euler 重複数は1、H0 重複数上界は1なので、H0=1、H1=0である。
これにより H1(E3)=0、h0(E3)=7285−594=6691 が従う。
複数段の filtration の接続写像を、根拠なく分裂として扱っていない。

主スクリプトを新しい出力先へ再実行し、6691次 modular minor の値20093を再現した。
この minor は有理階数の下界6691を与え、上の特性零の上界と一致する。
よって実際の一次係数の三次積が D3 全体を張る。
6691次行列全体を別実装で再生成したとは主張しない。

独立証明書は `code/ab_origin_degree_three_review_20261101.py`、
`reports/AB_ORIGIN_DEGREE_THREE_REVIEW_20261101.json`。
core の線形関係、巡回空間、fiber の条件、Bott の点検は主導出を import せずに行った。

```
.venv/bin/python code/ab_origin_degree_three_review_20261101.py --output /tmp/ab-origin-degree-three-review.json
```

次数四以上の全次数生成性はこの点検の範囲外である。
