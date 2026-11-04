# 非平方 cubic11：L 座標の多項式入力 preflight

これは六つの bordered determinant の再試行ではない。凍結済み
`S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63` の60秒中断を保持し、
同じ実 source の座標変更と入力サイズだけを確認する。rank・行列式・GB・
追加標数・追加 core 点・四次係数は計算しない。

旧五 light は `(x0,x1,x2,x3,x4)=(v0,v1,v2,v4,v7)`、
`v3=0, v6=x1-3*x0`。実二次 source の四行は `(q0,q1,q3,q5)`。
`C=100a+18e²-9e-2` として

```
L=-(405/4) (∂q1/∂x1)/C,   L_x1=-40C,
q1=L²/8100,   q3=L(eL-400a*x0)/8100.
```

従って `aC!=0` 上で `(L,x0,y1,y2,y3)`, `y=(x2,x3,x4)` は
可逆な light 座標であり、`(q1,q3)=(L²,L*x0)`。
source を直接読み、この三つの恒等式を core 多項式で検証する。
分母は既知単元 `-40C` だけである。元 genuine cubic16 行を全て保持し、
各行へ `(-40C)^3` を掛けてから全係数に共通な C の因子だけを厳密に除く。
q0 は x1 を含まないので無変更、q5 は `(-40C)^2` で分母を払う。
除いた C の個数は行ごとに保存する。未知因子は除かない。

次数3の `(L²,L*x0)` には九つの monomial pivot がある。
残る26 monomial を「L を含まない四変数の三次20本、
L 掛ける y 三変数の二次6本」の順にする。
`q0,q5` の light 五倍10行は

```
       no-L cubic20    L*y-quadratic6
         A(8x20)           B(8x6)
          0(2x20)        Cpair(2x6)
```

となる。後の二行は q0*L,q5*L。Cpair の全12係数は既存の
全域 pair dual source と列順・既知 C 単元を含め一致させる。
この dual をここで再計算しない。genuine16 は26 pure 列に H*light5 を
加えた16x31の多項式入力として保存する。

従来の20x35二次倍行から19個を密に消去する代わりに、九つを単項式で
消し、残る二次10行を8+2の block に分けた正確な入力である。
既知 pair dual と A の structured split を実際の cubic 行へ作用させる
計算は今回実施しない。canonical R11 の係数も未構成であり、次の
determinant 計算の時間をこの入力サイズだけから保証しない。

実行：

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/s61_nonsquare_ell_Lframe_preflight_20261104_v63.py
```

全域 rank57、全 quartic239、generic rank10上界、新しい kernel 方向の数は
変えない。generic cubic11 rank は既存の受理済み証明により10または11。
