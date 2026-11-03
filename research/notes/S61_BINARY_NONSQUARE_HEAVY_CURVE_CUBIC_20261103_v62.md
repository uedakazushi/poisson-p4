# Heavy 三次欠損曲線の一般点では全三次階数57

2026-11-03 / v62。三 core 開 `U={aFBC≠0}`、nuisance 零 section の正確な weighted chart を扱う。既知の15×7 heavy-light 行列を N、元135式から再構成された全三次行列の21既知単元 pivot 後の42×70行列を C とする。今回の結論は

    dim{ x∈U : rank N(x)<7 かつ rank C(x)<36 } ≤ 0  （QQ上）

である。この minor が与える下界は57以上であり、既存の [全Uでの三次階数上界57](S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102.md) と合わせると、N の欠損曲線はどの成分もその一般点で全三次階数が正確に57となる。16行の元 heavy block の欠損集合は N の欠損集合に含まれるので、同じ結論が成り立つ。**四次全射、横断長65、ambient t 冪零性はまだ結論しない。** 全三次階数の不足曲線が他にないとも主張しない。

## 1. 正確な QQ 曲線模型

[既存 heavy producer](../../code/s61_binary_hyperplane_heavy_cubic_codimension_20261102.py) は、元135式から16×7 template の全成分を再構成し、行尺度と列分母が a,F,B,C の単元だけであることを検査している。その template に既知の階数15の有理定数行変換を掛ける。新 producer は行 content のうち既知単元因子と非零有理定数だけを除く。例えば残る `e+1/2` の content は除かずに保持した。従って未知の core 因子を反転していない。

R=QQ[a,e,d,z]、D=aFBC とし、行列

    E = transpose(N) | (zD−1) I7

の余核 M と J=Ann_R M を計算する。M の支持は正確に `rank N<7` の U 上の表示であり、z は D の逆数である。QQ の exact module Gröbner 計算と exact annihilator 計算で、

    dim R/J=1,   degree R/J=16

を得た。J の Gröbner basis 25本を全整数係数で保存した。全25本の grevlex 先頭係数は p=32003 と素であり、係数も p 整である。全計算の M2 入力と出力は [report directory](../../reports/S61_BINARY_NONSQUARE_HEAVY_CURVE_CUBIC_20261103_v62/report.json) にある。

## 2. 特殊 fiber の全曲線を覆う

ここから使う特殊模型は、QQ で得た同じ25多項式の mod p 還元 Jbar である。有限体側で別途計算した annihilator を同じものと仮定しない。

Fp(a)[e,d,z] 上でこの25式の商は長さ4である。Gröbner basis と元25式の間の**実際の変換行列**も計算し、積がその basis に等しいことを検査した。basis の係数だけでなく変換行列の全係数の分母を取り、次数39の h(a) にまとめた。h(1)≠0 である。従って h≠0 上ではこの商は Fp[a,1/h] 上自由4である。basis の係数だけから h を採って、生成 ideal の同一性に必要な分母を落とすことはしていない。

h は相異なる10個の既約因子を持つ。各因子について実際の剰余体を構成し、25式の fiber を exact Gröbner 計算した。次元はすべて0以下で、長さは順に

    5,3,3,3,3,0,4,4,4,5

だった。各既約因子と重複度、結果を report に保存した。これで a が定数の曲線成分はない。従って特殊 fiber の全曲線成分は a を優越し、h≠0 上の有限な模型を通じて a=1 に会う。

## 3. 同じ一つの36 minor

a=1 では二つの二次拡大体

    Fp[s]/(s²+235s+1065),   Fp[s]/(s²−2903s−7388)

を使う。各定義式の判別式が非平方であることを検査した。report の e,d,z 座標に代入すると25式すべてが零になる。e の最小多項式が二つの組で異なることも検査したので、これは相異なる二つの Galois 対、計4点である。自由4の fiber の全幾何点を覆っている。

C の全成分は [Cloud 再構成](S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62.md) の QQ 行列そのものであり、分母はすべて p 単元である。同じ36行・36列を選んだ**一つの minor**は、二つの二次体でそれぞれ

    8483+1549s,   7601+9450s

となり、共に非零である。選択行列と全座標を report に保存した。体演算と Gaussian determinant を独立の短い Python 実装で exact に行う。これは有理点の標本から既約曲線一般点を推測する議論ではない。全4点で非零だから、有限自由4の模型上のこの元の norm は a=1 で非零となる。従ってこの minor は特殊 fiber の**全曲線成分**の一般点で非零である。

## 4. 非 proper な局所化を使わない QQ への移送

V=Z_(p) とする。25本の固定した QQ 多項式を、その実際の全次数で P4_V の同次座標 (a,e,d,z;w) に斉次化する。mod p の最高次数部分が生成する ideal L は、exact 計算で `dim Fp[a,e,d,z]/L=1` である。従ってこの射影模型の w=0 の部分は次元0である。25本の先頭係数の p 単元性は、ここで次数が還元によって落ちないことの追加確認でもある。

選択36×36 minor の各成分を、全成分の QQ 全次数の最大値 D0 に揃えて斉次化する。その行列式 m_h は V 上の次数36D0の実際の同次多項式であり、w=1 では元の minor である。巨大な行列式を展開する必要はない。還元後に余分な w の因子が生じても、w=0 の ann 模型は既に次元0なので曲線は生じない。

25斉次式と m_h が定める P4_V 内の閉部分は proper である。その特殊 fiber のアフィン部分は第2・3節によって次元0以下、無限遠も次元0以下である。proper morphism の fiber 次元の上半連続性により一般 fiber も次元0以下となる。そのアフィン部分に QQ の `J+(minor)` が含まれるので、冒頭の次元上界を得る。局所化した Z_(p)[a,e,d,z] の特殊 fiber だけから QQ の次元を推測する手順は使わない。斉次化模型の平坦性も仮定しない。

## 5. 再現と次の計算

    source /workspace/cloud-setup/activate.sh
    .venv/bin/python code/s61_binary_nonsquare_heavy_curve_cubic_20261103_v62.py \
      --output /tmp/heavy-curve-cubic-new

公開済み環境で実行済みの module/ann 計算を今回の集約時に再利用した。再利用に際しては、新 producer の生成した M2 入力が既存の入力と byte 単位で完全一致することを必須にし、再利用元を report に明記した。特殊 fiber、全分母、全10 fiber、4点、同じ minor は新 producer で再実行した。再利用なしの上記コマンドはすべてを最初から計算する。

次は、これらの曲線の一般点で**四次を含む実際の横断長65上界**を与える必要がある。全三次階数57はその一つの必要条件に過ぎない。全開上の quadratic-product 階数41は Fp では確認済みだが QQ の600秒の計算は時間切れであり、別の未完了課題として保存した。
