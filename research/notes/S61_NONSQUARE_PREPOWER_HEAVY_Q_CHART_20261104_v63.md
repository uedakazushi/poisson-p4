# prePOWER の実 H 小行列による九つの純三次判定チャート

本ノートは既存の係数値を展開せずに、実在する 15×22 入力から型付きの有限演算手順を固定する。新しい小行列式、階数、有限体点は計算しない。主対象は FULLU 内の未知の純三次欠損であり、その消滅をここで主張するものではない。

## 入力と局所化の範囲

二つの大域的 Q 行消去の受理済み入力
`S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json`
を使う。源行は Q9 と十四の genuine 行、標的は十七の pure 座標と五つの H 商座標である。Q9 の九つの小係数 h_i だけでなく、残る八つの pure 尾も全て保持する。

受理済み `GENUINE14_H5_PREFLIGHT` の全七十の恒等式は

    H_actual = diag(rho_i) N diag(1,1,1,1,gamma)

を与える。rho_i と gamma は a,d,F,B,C,J の既知単元である。源 G_i を rho_i で割り、最後の H 標的を gamma で割った後、五つの重い列は実多項式行列 N となる。pure 部分は P_i/rho_i となる。

今回は順序付き行 I=(1,3,6,12,9)、順序付き列 pi=(0,4,1,3,2) だけを選ぶ。この 5×5 行列 A の行列式 Delta は、受理済み `RAW_H_GENERIC_RANK` の原行 (2,4,7,13,10)、原列 (0,5,2,4,3) の既存の 1752 項非零多項式そのものである。新しく行列式を展開していない。Delta は一般係数体の単元であり、FULLU 全域の単元とはしない。

Q9 の九つの実係数 h_i は全体で 373 項である。既存の `LAST_Q9_POWER_BEZOUT` の受理済み恒等式 sum H_i h_i=a^29 により、D(h_i) は FULLU を覆う。ただし高次数の H_i は評価しない。H 側は今回の D(Delta) のみであり、その補集合を除去しない。D(p) 上の三十五チャートと、p=0 にも H7 が存在し得るという警告を証明書に引き継ぐ。

## 10×17 と九つの 9×16

u が I 外の九行を走るとき

    lambda_u = N[u,pi] adj(A),
    K_u = Delta P_u/rho_u - sum_i lambda_ui P_i/rho_i

と置く。源係数も同じ式で保存する。恒等式 A adj(A)=adj(A) A=Delta I により全五つの H 成分が零である。源 G14 内で対応する九行の小行列式は Delta^9 であり、Q9 を加えて H5 核の 10×17 表示を得る。実パラメータでの積は展開せず、25 個の独立形式変数について両随伴恒等式全五十成分を検査する。

Q9 のある pure 列 q の係数を h とすると、D(h) 上で

    h K_u - K_u[q] Q9

を取る。q 成分は全八十一行で文字通り零である。残る十六列を保存したものが九つの 9×16 判定行列である。源 fifteen-row frame の九行小行列式は

    h^9 Delta^9 product_{u outside I}(rho_u^-1)

なので、指定されたチャート上で九行は H5 と Q9 の消去から得られる基底である。

受理済み二つの大域 Q 完成と八つの既知単元ピボットを合わせると、表示されたチャートで

    rank(old17×24) = 8 + rank(new9×16),
    rank(aug25×32) = 16 + rank(new9×16)

となる。ここでは new9×16 の階数を計算していない。

## 実 raw H7 への持ち上げ

源十五行の係数ベクトルを原二十五行へ戻す演算は完全に固定する。消した Q5 の位置へ零を挿入して受理済み source16 を右から掛け、消した Q4 の位置へ零を挿入して source17 を右から掛ける。その後、残る十七行を二十五行へ埋め、八ピボットの消去を逆順に戻す。行正規化・消去の共変係数の向きを明示し、全 15×25 演算根を保存する。

前七ピボットは Q 行であり genuine 源係数と H に影響しない。最後の H*L ピボットは

    Gi_mu - (c_i/p_HL) G1_mu,
    p_HL = 2 F B^2 C^2 / 15625

を与える。この genuine 部分が全十五源ベクトルで保存された有限演算と一致することを検査する。二つの大域 Q 完成も genuine 間を混ぜない。

H5 核は実 raw H7 核そのものとは区別する。受理済み Cdelta=G6_mu-G1_mu は pure 射影と psi が零であり、raw 二成分 (zdelta,wdelta) を持つ。既知単元 eta=192 a F B^2 C^3/125 に対する

    cZ zdelta + cW wdelta = eta

を用い、持ち上げたベクトルから

    ((cZ HZ+cW HW)/eta) Cdelta

を引く。H5 零から psi=-wdelta HZ+zdelta HW=0 が従い、二つの随伴恒等式により raw HZ,HW も零となる。証明書は各八十一行の修正前 HZ/HW、修正係数、全原二十六源係数、原五十七源中の genuine 十六係数を保存する。

したがって、今後いずれかの 9×16 が階数九であると証明されれば、実九次元 H7 核の pure 像の独立性を得て full cubic57 が従う。ただしこの十六座標は全四十三座標の pure 商の射影である。射影の階数降下は full cubic 欠損と同値ではなく、四次欠損にも結論しない。

## 再現・資源・未決事項

新コードは `code/s61_nonsquare_prePOWER_heavy_Q_chart_20261104_v63.py`、証明書は同名報告ディレクトリの `circuit.json` と `symbolic-proof.json`。全入力の SHA256 と既存レビューを `report.json` に固定する。初期の `phase.json` は再起動前から保存されており、実行はその状態から継続する。

実行例:

    .venv/bin/python code/s61_nonsquare_prePOWER_heavy_Q_chart_20261104_v63.py \
      --output reports/S61_NONSQUARE_PREPOWER_HEAVY_Q_CHART_20261104_v63

上限は 60 秒、アドレス空間 2GB、成果物合計 32MiB、明示パラメータ項 150000。事前見積りは 25000 項未満で、pure 行列と二つの源変換は固定済み DAG 参照のまま扱う。旧四つの正規化因子も POWER207 係数も導入しない。実際の追加局所化は Delta と h_i だけである。H チャートの補集合、純三次判定の階数、全 FULLU の full57、全四次239、一般アフィン尽くしはいずれも未決である。
