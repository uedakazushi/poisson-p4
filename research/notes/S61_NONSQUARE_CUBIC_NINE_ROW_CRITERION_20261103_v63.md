# product41 と heavy7 の後に残る実9×43判定

2026-11-03。元135式から受理済みの63×91全三次行列と、その自由57行基底を用いる。
これは新しい全点 full cubic rank57 定理ではない。新しい pair proof が product41 を
全域で与えても、残る9方向の特殊化階数を別に扱う必要があることを明確にし、
実際の保存行列の小さな有限 circuit を与える。

## 実 block の exact sequence

重み三の91単項式は pure light 84 個と heavy-light 7 個に分かれる。
六二次式×七 light の42行は全て pure 側にある。既存の六つの unit 関係により、
全三次 row module M は41個の product 基底と16個の genuine cubic 基底からなる。
これは U0={aFBC!=0} 全域の自由57基底であり、d や J を逆にしない。
抽象商 M/P は自由16だが、特殊化した周囲91空間への単射性とは異なる。

product rank41 が成立する開では、その像 P は周囲の pure84 束の subbundle となる。
周囲を P で割った束は rank50 である。genuine16 の heavy7 projection が全射の開では、
heavy の核は rank9、pure 商の rank は84−41=43となる。従って実際に残る判定は

`rank9 source -> rank43 pure quotient`

の単射性である。heavy rank7 だけが保証する全三次下界は41+7=48であり、57ではない。

## 元行列を参照する387個の bordered minor

保存された正規化57×91行列 A を使う。最初の41行は product、残り16行は genuine cubic。
pure41列 I を選び product minor DeltaP が非零、genuine7行 E を選び heavy7 minor
DeltaH が非零とする。行 (product41,E)、列 (I,heavy7) の48角行列 D は block 三角形で、
`det D=DeltaP*DeltaH` となる。

残り9行、残りpure43列に対し、D に一行一列を足した49角 bordered determinant を
Sij と定義する。387成分を展開しなくても、これは元の多項式係数に対する有限 circuit である。
det D を単元とする chart で Schur 消去は可逆なので、厳密に

`rank A = 48 + rank S`。

よってこの chart の full57 判定は S の rank9 判定そのものである。
`DeltaP*DeltaH` は追加の **chart 単元** であり、全 core の既知単元とは呼ばない。
product41 かつ heavy7 の全開は、固定41 product 行の全41列 minor と、16 genuine 行の
heavy7全7行 minor による有限 chart 被覆で扱える。

今回の code は、元の全5190成分を保持する圧縮行列を再読し、product の heavy 成分が
全て零であること、84+7 の単項式分割を exact に確認する。既存の QQ 点 (1,2,3) だけで
非空な chart を選ぶ。FLINT の rational matrix 演算で product41、heavy7、48pivot、
残余9、全57を確認し、387個の rational Schur 成分と非零9minorを保存する。
選んだ indices と全 parametric 成分の定義は `circuit.json` に固定する。
これは点での非零を全点に読み替える検査ではなく、有限 circuit の chart が非空という追加確認である。
大きなパラメータ行列、generic function field、有限体、GB は新しく展開・起動しない。

## 誤った合成に対する明示的 module 反例

R=QQ[x,y,z]、周囲 R91=pure84⊕heavy7 とする。product は pure の固定41座標、
genuine の7行は heavy の固定7座標、さらに8行は未使用のpure8座標とし、最後の行を
`x*e49+y*e50` とする。57行は R 上独立で、product は全点rank41、heavy は全点rank7、
抽象商は自由16である。しかし x=y=0 の曲線上で全 rank は56に落ちる。
余次元二の full defect、一般rank57、全点上界57という既存の条件も満たす。
これは Jacobi を満たす実例ではなく、**既存の抽象自由性と二つの rank だけでは9方向が
消えない** という論理上の障害を証明する反例である。

## Heavy 曲線・J0・四次との接続の限界

既存の heavy-curve cubic proper proof は、heavy 欠損 locus Y に限って
`dim(Y∩{full cubic rank<57})<=0` を与える。新 product41 と合わせても、Y の外の
9×43 residual の欠損曲線は除外できない。heavy quartic proper proof も Y 上の全曲線に
沿う有限性であって、これを任意の cubic/quartic defect へ流用することはできない。

新 pair の P=ad 表示は d!=0 を要する。J=0 では既存
`S61_NONSQUARE_PRODUCT_J_ZERO_FINITE` が同じ ad 開で product41 defect の有限性を
与えるが、full57 defect の有限性とは異なる。d=0 はこの二命題の範囲外である。

次の実際の有限目標は、上の chart cover で S の rank9 欠損に曲線がないことを、
元行列からの最大 minor / annihilator / proper model のいずれかで証明することである。
その後も四次は独立であり、既存の **354×239 全重み四行列** またはそれと同じ ideal の
partial135 **428×300** 判定を使う必要がある。full57 だけから次の239列が全射とは言えない。
full57 の pivot chart で下位次数の帰結を四次へ掛けた後に再び Schur を取る、という
有限 recipe はあるが、その四次 pivot 数は未計算であり、小残余のサイズを創作しない。
ambient 長65・nuisance収縮・joint像<=15の適用には、三次だけでなく四次の実上界が必要である。

再現時は既設の `source /workspace/cloud-setup/activate.sh` を先に読み込む。
その後 `python code/s61_nonsquare_cubic_nine_row_20261103_v63.py --output <新しい出力>`
を実行する。公開済み FLINT の pkg-config 設定を再利用し、依存を導入し直さない。

出力: `reports/S61_NONSQUARE_CUBIC_NINE_ROW_CRITERION_20261103_v63/`。
この packet は具体的な9×43 criterion と不正な rank 合成への障害を凍結するもので、
新しい full57/ambient 除去定理を登録するものではない。
