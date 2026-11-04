# P100 の H15 有限代数上で全四次239を調べるための入力・資源監査

これは source-only preflight である。実係数、normal form、階数、小行列式、点は評価しない。既存素数32003における H15 の因数一覧は別の凍結報告 `S61_NONSQUARE_P100_H15_MODP_FACTOR_INVENTORY_20261104_v63` を参照する。因子は次数4,3,3,2,1,1,1、全て重複度一である。

## 完全な入力

`QUARTIC_319_ROWS` の retained catalogue は六つの二次式と重み二単項式の積174行、十六の実三次式と七 light の積112行、原 genuine quartic33行を全て保持する。既存の型付き `P100_PURE90_TYPED_SPLIT/operator.json` の full_Gamma319x98 を入力とする。未構成の「70+15+4」核や特定の部分集合の完全性は仮定しない。

各物理座標多項式をまず joint 座標へ移し、受理済み NF117 を適用する。残る H2,H1,pure 部分を f2,f1,f0 とすれば

    Gamma=f0-(f2/U2)*chi0-(f1-(f2/U2)*chi1)*Section,
    U2=(50aBC)^2,
    Section=Lactual^t * Ddelta7 * PiOrig,
    delta=-64000 C^2

である。PiOrig は112個の元 cubic-product の pure tail である。Lactual と組み合わせる前に実 delta7 源変換を必ず適用する。chi0 と全 pure tail を零と置かない。

十六三次の pure84 は全1344係数・139746 core 項として既に materialize されている。元の圧縮57×91を逐次読むことで全源を取得でき、巨大な source JSON は保持しない。実三次ラベル・正規化除算は凍結 FrameAdapter と pure-cubic source certificate に固定されている。

quartic 側は33の原多項式根を保持する方法が最も直接的である。Oracle は1373ノード、係数原子3484個で、Qhat_l=(50/B)^2 F_l(T34 diag(Tden) v) が正しい物理 frame である。H2を消した21行のラベルと演算も記録するが、21×pure210の新しい係数根が既に計算されたとはしない。21行は psi77/131,17,22,38,48,55,58,63,67,78,91,96,98,101,109,112,118,119,121,126,128 である。

## 有限体上の後続計算案と現段階の範囲

各既約因子 f に対し F_32003[z]/f の元を次数 d 以下の係数配列で表す。最大 d=4 なので98×98階数基底は38416個の prime-field slot で足りる。全319×98を同時に保持しない。係数処理と階数処理を分け、前者は一行ずつ書き出して源・NFキャッシュを解放し、後者は逐次行消去を行う方式を見積もる。selected98×98元行と消去基底の二重保持も明示して見積もる。

preflight は Oracle の単項式 support のみを伝播する。係数の打消しを仮定せず、各 quartic 単独と全33根共有の場合の last-use live 上界、全 NF122 の scalar/linear closure、source の literal QQ 分母を数える。これは値の計算ではなく、将来の係数計算の厳密な support 上界である。元139746項の三次入力も一係数ずつ走査する。

五つの曲線極 z,z+20,z+60,203z+2700,z+12 と全 Hbar 因子との共通因子を記録する。線形因子は剰余定数で、次数二以上の既約因子は次数で判定できる。他のパラメータ分母の有限代数上での値や逆元はまだ計算しない。QQ係数分母の p-unit 検査だけでパラメータ分母の可逆性を結論しない。

全因子で Gamma の98列が生成されれば、CRT により F_p[z]/Hbar 上で全射となる。モニック有限自由モデル Z_(p)[z]/H と全必要分母の単元証明を用いれば有限 Nakayama により標数零の H15 有限候補全体で full239 が従う。この正の方向だけを後続判定として認める。ある因子での欠損、極、時間超過は標数零の欠損を意味せず、因子を勝手に除去してはならない。

source/support builder:

    .venv/bin/python code/s61_nonsquare_P100_H15_full239_source_preflight_20261104_v63.py \
      --output reports/S61_NONSQUARE_P100_H15_FULL239_SOURCE_PREFLIGHT_20261104_v63

資源上限は60秒、2GB、32MiB。実係数・NF・階数の仕事にはこの測定結果に基づく別の承認が必要であり、ここでは自動的に実行しない。
