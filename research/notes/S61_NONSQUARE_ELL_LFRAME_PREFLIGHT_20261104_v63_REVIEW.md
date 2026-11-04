# ell L-frame 入力 preflight の独立監査

2026-11-04。対象は凍結済み `S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63`。
これは実入力の座標・row normalization・単項式削除の監査であり、進行中の別 job の
structured split、rank、border identity を再計算しない。

## 元 source と座標

元の quad 行は source 上の **q0,q1,q3,q5**、元 genuine cubic は凍結した57×91の行41,…,56。
同じ Tn frame の Z=3v0−v1+v6=W=v3=0 制限を使う。
旧五 light は (x0,x1,x2,x3,x4)=(v0,v1,v2,v4,v7)、H=v5。
公開された元8変数 monomial metadataから、反復整数 convolution で制限写像を独立に再構成する。
author の selection.txt と forward-coordinate helper は使わない。

C=100a+18e²−9e−2、L=−(405/4)(∂q1/∂x1)/C とする。
実 q1 の全係数から derivative identity と

    L_x1=−40C,
    q1=L²/8100,
    q3=L(eL−400a*x0)/8100

を確認する。よって aC≠0 上で新 light (L,x0,y1,y2,y3)、
(y1,y2,y3)=(旧x2,旧x3,旧x4) は可逆であり、(q1,q3) は既知 unit row operation により
(L²,L*x0) になる。q3 は元 source の第三採用行で、q2 と取り違えない。
row操作の分母は非零有理数、a、C だけで、FULLU の未知 factor を反転しない。

## 全 polynomial identity の独立な逆向き検証

レビュー native は別実装で新座標 L を **L=Σ Lj*旧xj と戻す**。
author の分母を払った forward substitution と C division を再実行するのでなく、
以下の逆向き polynomial identities を全係数で照合する。
den=−40C、q5 の C cancellation 数 k5=2、保存された各 genuine 行の cancellation 数 kr∈{1,2} として、

    Qnew0(L(oldx),…)=Qold0(oldx),
    Qnewi(L(oldx),…)=den² Qoldi(oldx)  (i=1,3),
    C^k5 Qnew5(L(oldx),…)=den² Qold5(oldx),
    C^kr Gnew_r(L(oldx),H,…)=den³ Gold_r(oldx,H).

すなわち四 quad の4×15=60係数、genuine16の16×40=640係数をすべて検証する。
元 quad20×35の全700係数も、元四 quad の五 light 倍として独立に再構成する。
逆向き検証では polynomial division は不要である。既知 aC unit 上の literal equality
なので、nilpotent baseを含む basechange でもこの source/frame 一致は保たれる。
新 genuine 16×40は **pure35 + H*light5** の全列を保持し、H部分にも同じ den³倍率を使う。

## 九つの固定 monomial pivot と残りの全入力

L² の五 light 倍と L*x0 の四 non-L light 倍は、cubic monomial 順の最初の九列を
ちょうど張る。これは既に正規化した unit monomial 行による商であり、
残る quad 行から未知 minor を探す操作ではない。
商の pure26列は、no-L の (x0,y1,y2,y3) 三次20本、L*(y1,y2,y3) 二次6本の順。
genuine16×31にはこの全26列と **H*L,H*x0,H*y1,H*y2,H*y3** を残す。

元 q0,q5 の全 light倍から、この商の10×26を独立に再構成する。
先の八行は q0,q5 の四 non-L倍、後二行は q0*L,q5*L。
後二行の最初の20係数が全て零なので、保存された

    [ A(8×20)  B(8×6) ]
    [    0     Cpair(2×6) ]

に一致する。全ての10×26、8×20、8×6、2×6、16×31係数を照合する。
quad10×26は317 core terms、pair2×6は40 terms、genuine16×31は50085 terms・最大 core degree19。
この count は新 source の展開からの独立な数え直しで、階数の主張ではない。

Cpair は既存 pair dual 入力の列反転 j↔5−j により、全十二係数について

    Cpair[0,j]=oldPair[0,5−j],
    C*Cpair[1,j]=(den²/C^k5)*oldPair[1,5−j]

を満たす。oldPair 第二行が **C*q5|L=x0=0** であることを保持する。
レビューでは分母を払った等価式を検証する。dualを再発見したり、新たな rank test をしない。

## 受理範囲

同名 REVIEW report の全チェックが PASS のとき、source-bound L frame、既知 unit row scalings、
九単項式削除、後続計算用の全 block input を受理する。
canonical R11、structured8+2+heavy5 elimination、六 border identity、generic rank10上界、
全点 cubic57、全四次239、追加 kernel count はこの packet の結論に含めない。
保存点 rank10と一般点 rank≥10の既存結果も変更しない。

実行（新しい空出力 directory）：

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_ell_Lframe_review_20261104_v63.py --output work/marked-square-v63/Lframe-review-REPLAY
```
