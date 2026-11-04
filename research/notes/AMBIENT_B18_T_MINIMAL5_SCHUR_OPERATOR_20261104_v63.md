# 実比較加群 T の五生成 Schur 表示

## 範囲と実入力

基底環は実 CA-F9 座標の
`Rbar=QQ[[zA0,...,zA6,rA]]=R0/(ell)`, `ell=5(sA-rA)` である。
独立監査済みの全116生成936関係の表示を出発点とし、凍結済みの実 X0
の27 pivot と、実 X1 から得た84 wedge pivot だけを選ぶ。新しい source
方程式、jet、正則商の再除算、348次元上の scalar 行列は求めない。

出力は、元の全形式函数を係数とする **五生成825関係の完全な表示**である。
係数は二つの単位行列逆を持つ factored typed operator として保存する。
これは冪級数係数を展開したという主張ではなく、825関係が極小という主張
でもない。五生成は全関係の定数項が零なので極小である。

## 第一の27単位消去

元の target は84 wedge と32自由 target に分かれる。実 X0 の27独立列を
凍結済みの定数 column 操作で正規化する。同じ操作により残る741 source
列 `Xprime` は全56成分で定数項零となる。32自由 target の定数 row 操作
は、pivot27成分とその商五成分に分ける。五成分の projection は、凍結済み
の27正規化列を pivot 番号の昇順に消去する操作そのものである。その元の
V56 source row ラベルは `[22,23,29,30,55]` となる。

この frame で、168本の literal d3 と source 列の全行列は

```
             d3    pivot27   source741
 wedge84      D      Wp         Wr
 pivot27      0       A          C
 free5        0       B          Y
```

と書ける。`A(0)` は単位であり、`B,C,Y` は最大 ideal `m` に入る。
特に `B(0)=0` は選んだ定数 quotient projection の定義による。
27対を全次数で消去した結果は

```
W_hat = Wr - Wp A^(-1) C,
Y_hat = Y  - B  A^(-1) C.
```

`B A^(-1) C` は m二乗に属するので `Y_hat_1=Y_1` である。
また `C(0)=0` より `W_hat_0=Wr_0` である。`Wp(0)` や A の一次係数を
知らなくても、この二つの低次数結論は全次数の Schur 恒等式から従う。
元の d3 の自由32成分は全次数で零なので、この段階で D は変わらない。

## 第二の84単位消去

凍結済みの W0 消去の **定数** column 操作を `(W_hat,Y_hat)` の両方に施す。
84本の正規化 pivot と657本の残余に分けて

```
             d3       pivot84      residual657
 wedge84      D          U            Wrest
 free5        0        Ystar          Yrest
```

を得る。`U(0)` は実84 pivot の正規化された可逆行列であり、`Ystar`,
`Yrest`, `Wrest` はすべて m に入る。ここで `Wrest(0)=0` には、選択列を
並べ替えるだけでなく、凍結済みの定数 column 消去を全部施したことが必要
である。

残る完全な関係行列は、全次数で正確に

```
R_min = ( -Ystar U^(-1) D  |  Yrest - Ystar U^(-1) Wrest )
```

という五対825行列である。第一段階27対と第二段階84対はともに形式単位
の Schur 消去であり、元の116対936の cokernel を保つ。実 source の全768
列、全168 d3、旧正則 r/g 係数の参照を残す。Euler operator は元の形式
函数に作用したものを参照し、商環の normal form に落としてから微分しない。

## 既存 X1 だけで確定する低次数

全 source657列の一次関係は `Yrest_1` そのものである。第二の補正積は
`Ystar,Wrest` がともに m に入るため m二乗から始まる。第一段階でも先述
の `B,C` の消滅がある。このため未知の m1、X2、W_p の定数、単位 U の
一次係数は一次関係へ寄与しない。

literal d3 は定数項が零である。その消去後の列は m二乗から始まり、二次
初期形は正確に

```
-(Ystar)_1 (U_0)^(-1) D_1
```

である。これも全次数の行列積の次数だけで従い、新しい係数計算を要しない。
一方、657本の source 関係の二次項がすべて分かったとは主張しない。

一次関係が五対八の40次元空間で rank37なら
`length(T/m²T)=5+3=8` となる。これら一次関係と八変数との積、および上記の
168二次初期形が五対 `Sym²(QQ^8)` の180次元を全部張れば、全次数の実関係
から `m²T=m³T` を得て Nakayama により `m²T=0` となる。この場合に限り
この十分条件は全長8を与える。部分的な二次 span や単なる接空間の値は、
未知の残余二次関係を零とする理由にはならず、全長・元の深さを決めない。

本 packet 自身は rank を再探索せず、A0 の定数操作・source label・低次数
interface の全固定式を再生する。深さ19対20と元の全 nil 三乗の決定は、
全長を確定する別の証明が必要である。

## 再現と資源

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_T_minimal5_Schur_operator_20261104_v63.py
```

既存の外部形式函数 DAG を参照し、追加の matrix operator と4125個の scalar
readback root を保存する。実 source、全 constant column 操作、二つの単位
逆の origin 証明、全 relation order を hash で結ぶ。60秒・2GB・32MiB・
150,000参照の guard 内であり、旧凍結ファイルは変更しない。
