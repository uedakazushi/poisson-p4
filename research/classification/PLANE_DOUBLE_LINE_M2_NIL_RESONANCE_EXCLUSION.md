# p=−2/7 nil normal-double 全 fiber：一つの元 Jacobi coefficient で空

2026-11-08、producer 全源証明。全 p の necessary source 還元とは分ける。
入力は immutable 三十四源 report の `HA3_NIL1_WHOLE_P`。
保存済み generic complete K1 十一を p=−2/7 へ評価する。
この値は generic minor／complete-kernel denominator の wall ではなく、
その nonzero minor と free-coordinate identity を正確に評価して完全性を確認した。
新しい M1 profile は再計算していない。

全 source は k0,…,k10 を持つ P1 十一と、**FULL P2 九** q0,…,q8。
homogeneous K2 三に制限しないので未知の inhomogeneous P2 particular も全て含む。
元 bad normal functional は `m=coefficient(P1_cu,v²)=k5`、
`T=coefficient(P1_av,u²)=k6`。
全 source の raw 150 coefficients は resonance report に保存済み。

## 1. 正規化しない元 coefficient

strict Jacobi を `[Π,Π]`、その projective equation を `E∧[Π,Π]` とする。
reusable SymPy helper は `[Π,Π]` を cyclic expression の二倍で実装する。
元三百五十 row のうち

```
quad=(a,b,c,u), coordinate monomial=a*v³
```

の coefficient（zero-based row 34）は、**正規化しない全二十 source variables** 上で正確に

`J34=−14 k5²/9=−14 m²/9`。

全 other k、全 q、T に依存しない。
single cyclic expression を使う独立実装では coefficient は −7m²/9 であり、同じ零条件になる。
元 row ordering を変えたり、他の未証明 equation で割ったものではない。

別の手確認：全 source の faithful c-gauge で Q/W notation を使えば
この coefficient は `2m(2W−V)`。
保存済み全 kernel の global coefficient identities は p=−2/7 で
`V=0,W=−7m/18` を与え、同じ −14m²/9 に一致する。
原 source を直接代入した coefficient extraction でも確認した。

従って複素 Poisson point では必ず m 零。
これは reduced support に対する結論であり、非被約 scheme 上で m² 零から m 零と等置しない。

## 2. actual divisor の u 因子

元 P0 の全 u-row は u の倍数。
full K1 と FULLP2 を戻した元 quadratic tensor を u 零へ制限すると、その u-row は

`(Π_au,Π_bu,Π_cu,Π_vu)|u=0=(0,0,m v²,0)`。

P2 九のどの係数も normal column を持たない。
m 零の全 Poisson point では full u-row が u 倍となるので、
H={u=0} で tensor は三次元 tangent space 内に入り rank 四を持てない。
従って actual quintic D は u を因子に持つ。
固定 P0 の actual leading coefficient は `(−2/7−4)a³u²=−30a³u²/7` 非零。
したがって D 非零の既約 reduced quintic source は全空である。

これは **T 零／非零、m 零／非零を全て含む p=−2/7 whole fiber** の空性。
先に valid に証明した Tm-unit normalizer を使う必要すらない。
unit normalization 後の 1=−9J34/14 という exact report はそのまま保存するが、
最終の弱い十分証明は元全 source の一係数と actual divisor 因子だけである。
GB、higher jet、追加 profile、選んだ source 点又は generic density をこの固定値で用いない。

## 3. 全 p との接続と射程

`PLANE_DOUBLE_LINE_M2_LAST_NIL_ORE_REDUCTION.md` は全 p の補集合を保持し、
T 零は essential-three、m 零は actual u 因子、p−二も u 因子へ送り、
未処理 Tm-unit branch の必要値を p=−2/7 へ縮約した。
本 note がその残 whole fiber を全空にする。
両 proof の source 完全性／係数 convention／全 wall への延長は独立 review で別に確認する。
本固定値の証明だけで全 rank二 plane、全 E帰属又は全分類の完了は主張しない。
