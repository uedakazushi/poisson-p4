# 既存二次値からの Cprime/m³ の全作用の復元

2026-11-05。これは source の代数的縮約であり、以下の作用行列の値や genuine connecting 値はまだ計算していない。Cprime の真の HF は次数 2 まで (3,3,3) と独立検証済みである。V=V1+V2+… は実 R0³ presentation の 44 列で、V0=0 である。

V1 の固定 24 独立列を B、固定 pivot 行を R、対応小行列の逆を I とする。一次 quotient の三自由行の標準 embedding を E1、標準 quotient を q1 とすれば、

    identity_27 = E1 q1 + B I restriction_R

である。真の二次 quotient は NF3=pi3 であり、その標準三行の embedding を E2 とする。全二次空間では pi3 E2=identity_3 である。

二次 homogeneous 関係の既知空間 K2 は、一次式の一次倍と、実二次 kernel 列 18 の initial の和であり、pi3 の核に等しい。一次式の一次倍は、実関係 V の一次倍に等しい modulo m³ ので、K2 の各元は Cprime/m³ で実際に零になる。

従って V2_j を、pi3(V2_j) と同じ二次代表へ置き換えてよい。今回 pi3 は pure-plane 三方向しか読まないため、保存済み 44 列の pure-plane 二次値だけで十分である。零埋めした source 近似を使うのではなく、捨てる部分が実関係空間 K2 に入ることを用いる。

head 三個、一次 E1 三個、二次 E2 三個を九基底とする。一次単項式行 e の reduction は

    [ 0 ; q1(e) ; −pi3(V2_piv) I restriction_R(e) ]

であり、二次単項式行 f の reduction は [0;0;pi3(f)] である。全 44 元関係 V1_j+V2_j がこの reduction で零になることは、固定 V1 factorization と、完全 kernel 二十列の pi3(V2 K20)=0 によって直接確認できる。head を変数倍した結果を前式で還元し、一次 E1 を変数倍した結果を後式で還元し、二次 E2 の変数倍を零とすれば、全九変数の実作用が得られる。

この九次元加群は head 三個から生成される。全元関係とそれらの一次倍を満たし、既知の長さ 9 と一致するので、Cprime/m³ そのものである。作用の可換性、全 44 元関係、および関係の全一次倍と全三次単項式作用の消滅は固定行列恒等式として検査できる。必要な値はすべて既存 V1、固定 24 小行列逆、44 個の pi3(V2_j)、および固定 pi3 にある。新 physical source・gauge・三次 jet は不要である。

さらに、別の既知 initial 三次上界試験が m³(Cprime/ell³)=0 を証明した場合に限り、ell³Cprime⊂m³Cprime と逆包含が一致する。このとき上の実九次元加群は自然に Cprime/ell³ でもあり、全 torsion 指数が 3 未満という別の DVR 条件と合わせて

    Cprime[ell] ≅ ker(ell:Cprime/ell³→Cprime/ell³) / ell²(Cprime/ell³)

を使える。三次上界未証明の段階では Cprime/m³ と Cprime/ell³ を同一視しない。

これは finite target の追加 source 計算を省く縮約であり、formal numerator V_full sigma を ell で割る edge の source 次数を省略する根拠ではない。真の二次元 target が m² により消えることから、その Cprime への埋込みが filtration strict だとは推論しない。connecting の必要係数は別途、元 formal source 上で感度を抽出する。
