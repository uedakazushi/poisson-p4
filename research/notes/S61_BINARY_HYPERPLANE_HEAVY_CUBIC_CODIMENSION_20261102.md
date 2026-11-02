# K=0 三core・heavy cubic block の階数欠損は余次元二以上

2026-11-02。既証明 `S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102` の
nuisance-zero三core `(a,e,d)`、指定開 `a F B C !=0` を用いる。

\[
F=500a+90e^2+3e-24,\quad
B=50ad+18de^2-9de-2d+60e-40,\quad
C=100a+18e^2-9e-2.
\]

正重み八法線座標は `(1,1,1,1,1,2,1,1)` であり、heavy座標は `v5`。
本稿は重み三の式に含まれる七つの `v5*vj (j!=5)` の係数ブロックを扱う。
このブロックの `rank<7` 集合は指定三core開内で高々一次元である。
全開上のrank7、pure cubicの階数、重み四次、ambient飽和は主張しない。

## 1. 元135式から必要な七列を直接作る

局所元135式を `J y+H(y,y)` とし、凍結したlinear sectionを `T,U` とする。
`T` の各法線列の共通分母を `t_j`、`U` の全成分の共通分母を `delta` とし、
`Tn_j=t_j T_j, Un=delta U` を多項式化する。
producerは各分母が既知 `a,F,B,C` の積と非零有理数だけからなることを検査し、
再構成した元Jについて `J Tn=0` と `J_selected Un=delta I34` も全係数で検査する。
実際の列分母は `(a,a,a,1,a,B/50,a,a)`、`delta=a²FB/25000` である。
これは指定開内の可逆な法線座標の尺度変更である。

二次通常次数の各heavy-light係数を

\[
 Q_{i,j}=[v_5v_j]H_i(T_n v,T_n v),\qquad
 W_2=-U_n Q_{selected},\qquad A_2=\delta Q+J W_2
\]

で作る。具体的な各 `H_i` の単項式 `c*y_r*y_s` の寄与は
`c*(Tn[r,5]*Tn[s,j]+Tn[r,j]*Tn[s,5])`。
重みを保つ陰関数解で、通常三次以上の補正はこの通常次数二係数を変えない。
従ってこれが正確なheavy cubicブロックであり、raw四次展開のキャッシュは不要である。
元135式のうちheavy成分非零は16行で、残りはこのブロックでは零。

## 2. 小さいtemplateへの可逆な比較

探索で得た16×7多項式templateを
`code/data/s61_binary_hyperplane_heavy_template_20261102.json` に保存する。
これは計算結果を無条件に仮定する入力ではない。
producerは元式から直接得た16行と、templateの対応する各行とが七成分すべてで
正確に比例することをQQ(a,e,d)上で検査する。
その比例係数の分子・分母から既知 `a,F,B,C` の因子を完全に払い、残りが非零定数で
あることを確認する。各係数の定数と既知因子指数はreportに保存する。
実際の比は非零有理数と `a^r F B (r=1,2,3)` の積である。
従って指定開上で両行加群は一致する。template各行の七成分の全QQ多項式GCDも
非零定数であると検査し、未知の行contentは残っていない。

このtemplateは以前の大きいweighted入力から共通単元F*Bを除いて得たものだが、
今回の証明はその入力、raw四次、元の行contentの採り方には依存しない。
全比例恒等式と単元条件を直接再検証するためである。

保存した15×16のQQ定数行列でさらに行結合する。その定数行列の階数15を検査し、
新15行を実際にtemplateとの積として構成する。
以下のminorは従って直接heavy行の局所化した行加群に属する行から得られる。

## 3. 三つの最大minorと支持

新15×7行列の次の三組の行を使う（零始まり）。

- `(13,14,12,11,7,9,1)`
- `(7,4,14,3,8,5,0)`
- `(5,13,11,12,2,0,7)`

native FLINTのfraction-free Bareissで三minorの全QQ[a,e,d]係数を計算する。
各中間除算はexact polynomial divisionとして検査し、行・列交換の符号も保持する。
三多項式の正確なGCDは

\[
  a^4 C/100
\]

である。全minorをreportへ保存する。
補助的に `(a,e,d)=(1,2,3),(2,1,1)` の二点で、独立SymPyの有理数7×7行列式と
保存多項式の評価が一致することも確認する。この標本検算を全多項式GCD証明の
代用にはしない。

`QQ[a,e,d,(aFBC)^-1]` はUFDであり、三minorのGCDは同環で単元となる。
従って全最大minor idealを含む高さ一素idealは存在せず、階数欠損集合は
余次元二以上、すなわち三core内で高々一次元である。
GCDが単元でも多変数idealが単位idealとは限らないため、全開の全点でrank7とは言わない。

## 4. 再現と依存

```sh
.venv/bin/python code/s61_binary_hyperplane_heavy_cubic_codimension_20261102.py \
  --output /tmp/heavy-cubic-codimension.json
```

必要な元模型・凍結linear section・小template・native source・本noteのhashを主reportに
保存する。出力は新規作成のみ。82MBの探索rowsファイルも通常三次・四次raw展開も読まない。
既知weighted chartと元S61模型を依存に保持し、新たな全K=0 sourceの飽和除去を
この補題だけから主張しない。
