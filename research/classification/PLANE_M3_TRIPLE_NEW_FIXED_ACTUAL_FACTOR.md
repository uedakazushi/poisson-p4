# 正しい triple 表の新 fixed 源：元 Jacobi ideal からの actual z 因子

2026-11-08、producer 全源証明。正しい V2 の
`E2_LOG_ETA-1_4`、すなわち

`f=a²b,h=a²,p=dlog a,q=2dlog a−dlog b,r=dy/x`

を扱う。旧 plus-curvature coverage に依存しない。
全 first-symbol inverse と fullK/P0 は corrected-V2 report に保存済み。
元 P0 は K 零、c 非依存で

`Q0=3az ∂b∧∂w, W0=(a²,−2ab,az,4aw)`。

## 1. 全 source と二つの不成立な Ore 入口

元 full weight-one＋actual D3-line の全 kernel 十五を p0,…,p14、
**FULL weight-two source 九**を q0,…,q8 とする。
P2 を homogeneous Jacobi kernel 五に制限しない。
従って weight-two の inhomogeneous solution を含む全 actual source は
`P0+P1(Kp)+P2(q)` に入る。
full original three hundred fifty Jacobi equations を J0,…,J349 とする。

全 source を c-Ore とする試みは Qab の c-derivative p8 z を残し、
b-Ore とする試みは Qac の b-derivative p5 z を残した。
両方の covering criterion は false のまま別 report に保存される。
p5 p8 又はその平方が元 Jacobi rowspace にあるという予想も不成立。
以下はそれらの仮定を一切使わず actual divisor を直接扱う。

## 2. 三つの元 Jacobi row 結合

full original support report の exact QQ row reduction の三式を

```
A = p14(p14−p9),
B = p14 p8,
C = q2 − (p14−p9)(p10+4p2−p4)/15 − p3 p8/3
```

と書く。A,B,C はそれぞれ元 J0,…,J349 の **QQ-linear combination** である。
focused factor certificate は三つの全 row lifts を保存し、
元 source の polynomial として各 identity を正確に再確認する。
これらは必要 equation であり、全 Jacobi support の分類を意味しない。

## 3. actual quintic の制限に対する exact ideal identity

全 source に対し、actual divisor は convention

`D=Σ(−1)^i xi Pf(Π omit i)`

で直接計算する。Euler gauge の変更 E∧L はこの polynomial を変えない。
元 full support report の D を z 零へ制限すると、全二十四 source coordinate を残したまま

`D|z=0 = a A w^4 + 2b B w^4 + p14 q2 w^5`。

一方、定義から全 polynomial ring 上で

`p14 q2 = p14 C + (p10+4p2−p4) A/15 + p3 B/3`。

従って

```
D|z=0 = [a w^4+(p10+4p2−p4)w^5/15] A
       +[2b w^4+p3 w^5/3] B
       +p14 w^5 C.
```

これは source coordinate を一つも逆元にしない、元 Jacobi ideal 内の exact identity。
focused certificate は A,B,C の QQ row lifts を代入した
**D|z=0 の元三百五十式への polynomial lift** も保存する。
generic rank、選んだ点、有限 prime、radical membership の推定は使わない。

全 actual Poisson point で D|z=0 は零なので global homogeneous polynomial D は z の倍数。
rank 四なら D 非零であり、degree 五の既約 reduced divisor とは両立しない。
よってこの新 fixed first-symbol の **全既約 quintic source は空**。
通常 Ore でない点を消去したり homogeneous P2 に制限したりした結論ではない。

## 4. 射程

corrected-V2 九十源の direct criteria、正しい表からの旧 valid literal source への実包含、
旧 valid 六源の十八全 affine module、及び本 exact factor を合わせれば、
非被約 R／位数三／generic triple normal cubic の全 first symbolsについて
未知 β優越成分の一般元を排除する G obligation が閉じる。
direct dimension-bound で閉じた source の E帰属はこれとは別であり、
全 generic triple の全点既知成分被覆又は全研究の完了は主張しない。
whole-source と独立 review の受理／global claim 登録は統合担当の別工程である。
