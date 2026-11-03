# J=5ad−6e−4=0 上の product 階数欠損は有限

2026-11-03 / v62。U={aFBC≠0} のうち d≠0、J=5ad−6e−4=0 を扱う。ここでは六つの light 二次式の light 一次倍42行の階数が41未満となる集合は、QQ 上で次元0以下である。全開上の全点階数41、J≠0 の三 core、四次・ambient 冪零性はこの結果単独では結論しない。

## 二 core の正確な表示

J=0 では d=(6e+4)/(5a) であり、d≠0 は 3e+2≠0 と同値である。定数座標

    (x0,x1,x2,x3,x4,x5,x6)=(v0,v2,v4,v6,v7,Z,W),
    Z=3v0−v1+v6, W=v3

を使う。ℓ の x2 係数は正確に −32(3e+2)/27 となるため、ℓ=0 で x2 を消去できる。残る六 light 変数で、四つの二次式 q0,q1,q3,q5 の一次倍は24×56行列になる。

新 producer は元の六二次式からこの変換・代入・消去をすべて QQ(a,e) で再構成する。分母と除く行 content は

    a, F=500a+90e²+3e−24,
    Bs=300ae+54e³+9e²−24e−4,
    C=100a+18e²−9e−2, 3e+2

の単元と非零有理数だけである。元の B はこの面上で 2Bs/(5a) に等しい。未知の core 因子を除かない。既知単元 pivot を12個、exact な行列演算で除いて12×44残余を得た。

## 四つの QQ minor と gcd

残余の四つの実際の12×12 minor を FLINT の QQ 多項式 Bareiss で計算した。全中間割算の exact 性と符号を検査し、得た四多項式の gcd を正確に計算する。結果は非零有理数倍を除いて

    a^12 C^3 (3e+2)^31

である。全係数と四つの列選択を [report](../../reports/S61_BINARY_NONSQUARE_PRODUCT_J_ZERO_20261103_v62/report.json) に保存した。minor 計算そのものは48.82秒。さらに (a,e)=(1,2),(2,1) の二点で、独立 SymPy の QQ 行列式と四多項式の評価がすべて一致することを検査した。

従って既知単元を反転した QQ[a,e] の UFD では、これらの全 minor を含む高さ1の素 ideal はない。二 core なので rank24 欠損集合は次元0以下となる。幾何学的な曲線成分が存在すればその有限個の共役から QQ の共通因子が生じるため、代数閉包への拡大後も同じ結論である。有限体で gcd=1 だったことから QQ を推測した証明ではない。

## rank41 への線形代数

六二次式は、a を単元として

    (q0,q1,q3,q5, ℓW,ℓZ)

へ可逆に変換できる。[pair cover](S61_BINARY_NONSQUARE_PAIR_COVER_20261103_v62.md) により Z,W は独立で、q0,q1,q3,q5 を Z=W=0 に制限した四行は U の全点で独立である。

四二次式の一次倍が ℓ を法として独立な点で、42本の関係に線形 syzygy

    Σ li qi + uℓW + vℓZ=0

があれば、ℓ を法として各 li は定数倍 ciℓ となる。ℓ で割り Z=W=0 とすると四行の独立性から全 ci=0。残る uW+vZ=0 は、WZ−ZW の定数倍一つだけである。従って syzygy 次元は1、42行の階数は41となる。これにより J=0 の指定開で product-rank41 欠損も次元0以下である。

## 再現

    source /workspace/cloud-setup/activate.sh
    .venv/bin/python code/s61_binary_nonsquare_product_J_zero_20261103_v62.py \
      --output /tmp/product-J-zero-new

今回の集約では、先に同セッションで実行した四 minor の入力が byte 単位で完全一致することを確認して結果を再利用し、座標表示と単元 pivot の再構成・独立 QQ 検算を再実行した。再利用元を report に記録した。上記コマンドは再利用なしで全 minor を計算する。
