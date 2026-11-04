# 独立監査：四つの g 差二次値に必要な source 五次入力

本監査は source/order/input の preflight であり、114未知係数や g2値を求める
計算ではない。既存 own CA/CB57 source4 を一次元へ縮約し、degree4までの
方程式、固定Jacobian逆、signed numerator と分母の source同一性だけを
検証する。実行する環は QQ[r]/r五乗であり、source5値は生成しない。

## 座標と source 57 の同一性

CAの実common曲線は own F9 の零7、最後二座標(r,r)である。CBの同じ曲線
の own F9 座標は保存105 sourceの読戻しにより

```
(YB7,YB8)=(17/12,7/12)*r + (0,0)*r² + O(r³).
```

二つの nodeID が cached CB57 F9 の最後二個と一致することを確認し、単に
傾きの数値が同じというだけでは扱わない。両系の元57未知・元57方程式・
frameを以前の full9 source2、独立監査済み source4 と照合する。

全 saved source4係数をその傾きへ普通の多項式係数として代入し、元scalar
DAGを別の五係数 series arithmeticで再評価する。全57式は最初55式が零、
最後二式が slope*r にdegree4まで一致する。全scalar単位逆の積も確認する。
固定57 Jacと保存逆の積は恒等である。新しい Gaussian elimination はない。

## regular quotient と符号

四つの旧CB quotient は、元CB sourceの bbar/g を full unit inverse14へ
写した関数である。四つの新CA quotient は、元CA sourceの bbar/g を
zero-prefix0 inverseへ写した関数である。元の numerator/denominator rootを
実bindingで構造的に置換し、保存された old/new numerator/denominatorの
expressionと一致することを検証する。

この検査はnode表示名を使わず、variableの実nodeIDを保つ。後のDAGにある
重複negative aliasはcanonical expression signatureで処理し、異なる
variableを名前で同一視しない。従って元 gdiff は正確にoldCB quotientから
newCA quotientを引いた値であり、別の符号やunitを勝手に追加しない。
下流SchurのXp最後成分はその差の負号であることと混同しない。

## 必要次数の理由

regular cofactor hの存在は元の既知I0-membershipにより証明されている。
保存された四つの bbar **全own二変数**のdegree3以下の係数が零であり、
gのdegree3 homogeneous部分は非零である。従って関係 bbar=g*h のdegree3を
取れば h0=0が従う。選択line上だけの b3=0から推測するのではない。

両直線に沿う分母は

```
g_CA=(125/24)*r³+0*r⁴+O(r⁵),
g_CB=(-15125/144)*r³+0*r⁴+O(r⁵).
```

ordinary係数で比較すると

```
bbar5=g3*h2+g4*h1+g5*h0=g3*h2.
```

よって h2=bbar5/g3である。r三乗を形式的に取り除いた後の分母は非零有理数
を定数とする単位であり、零のrを反転していない。g5の値は不要である。

actual curved pullbackの二次係数は、own straight-slope h2と、own h1に
actual YBの二次curvatureを掛けたものの和である。curvatureが零なので両者
は一致する。YB三次以上の値やcommon105のdegree4/5を使う必要もない。
Hessianのfactorial規約は使わない。

## 将来の一回の114係数計算

各familyは同じ57単位sourceを使い、saved u0..4に未知u5*r五乗だけを加える。
元57方程式のdegree5残差に、固定Jac0*u5が加わる。Jac(u)-Jac0による項は
degree6以上なので

```
u5=-Jac0^(-1)*residual5
```

で足りる。両familyで57ずつ、計114個の普通の係数だけが新入力である。
実行が許可された場合でも、補正後は元57式全てをdegree5まで検証する。
元の分母、旗、非線形ratio依存を省いてはならない。

source closureはCA9735、CB9978で逐次処理でき、既存4次評価時の最大live
nodeは373/381である。six-coefficient workerの保守的live/record上界は
18,605 QQ slotであり、2GB・32MiB・150000係数のguard内を狙える。新しい
五次係数の高さや実時間を計算前に保証したとは主張せず、60秒guardと失敗
prefix保存を維持する。本監査自体はnumeric workerを許可・実行しない。

## 再現

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_T_gdiff_second_univariate5_preflight_review_20261104_v63.py
```

新しい source5、114未知の求解、619実値、HF3、T長さ、nil三乗の結論はない。
