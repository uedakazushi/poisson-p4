# Shepherd-Barron の一次文献による W16/Sp4 の複素有理性

2026-11-04 UTC。新しい文献適用記録。以前の stable level≤2 の凍結記録は変更しない。
以下は既知文献の正確な対象への適用であり、新規性・優先権の主張ではない。

## 一次資料と確認箇所

N. I. Shepherd-Barron, *Invariant theory for S5 and the rationality of M6*,
Compositio Mathematica **70** (1989), no.1,13–25。

* [原文 PDF](https://www.numdam.org/article/CM_1989__70_1_13_0.pdf)
* [書誌ページ](https://www.numdam.org/item/CM_1989__70_1_13_0/)

PDF は表紙を含む14ページ。原文の頁数と PDF の0-based index は次の対応で確認した。
本文 OCR に加え、印刷頁16,17,18,20,21,25を画像で読んだ。逐語引用はしない。

|箇所|印刷頁|PDF index|適用の確認点|
|---|---:|---:|---|
|Lemma8|16|4|character twist による射影商の同一性|
|Proposition9|17|5|5次元表現の商と quintic del Pezzo の symmetric square の対応|
|Theorem10|17–18|5–6|その symmetric-square 商の有理性|
|Proposition11|18–21|6–9|二つの opposite biform の **全 Aut(Q)** 商の有理性|
|Lemma14|20–21|8–9|torus-normalizer slice を quadric bundle に還元する証明|
|Theorem15|21|9|5次元 χ の射影商の有理性|
|character table|25|13|今回の [3,2] と原論文 χ′ の照合|

Proposition11 は A=P(H0(O_Q(1,2)))、B=P(H0(O_Q(2,1)))、
Q=P1×P1、G=Aut(Q) と定義している。G の非恒等成分による factor 交換も含む。
その証明は二段の slice と Lemma14 を通り、最後に Tsen の定理を用いる。
ここでは公刊された定理を用い、全 slice の係数を新しく再計算したとは主張しない。

## 今回の対象との一致

既受理の物理 Spin4 切片は厳密に

```
L=V(2,1)⊕V(1,2),  V(a,b)=Symᵃ(A)⊗Symᵇ(B)
```

である。追加 Spin5 成分は (U,V)↦(−V(y,x),U(y,x))。
各形式を非零定数倍で割ると、この符号は消え、原論文の **全** Aut(Q) 作用に一致する。
単独 biform の定理を直和へ移しているのでも、connected quotient から有限群の
成分を落としているのでもない。

したがって原論文 Proposition11 により、既存の四次元体

```
B0=C(P(V(2,1))×P(V(1,2)))^Aut(Q)
```

は C 上の有理関数体である。

S5 を通る別の照合もできる。既受理の projective-frame 還元では
B0=C(P(V5))^S5、V5 は [3,2] 型。
今回保存した character は原論文の χ′ に一致する。Theorem15 の χ と χ′ は
Lemma8 の sign twist で同じ射影商を持つため、この経路でも同じ結論を得る。
主たる適用は表現名の照合を要しない Proposition11 とする。

## 二倍率を戻した結論

既受理の `Q_FIRSTJET_FIVEPOINT_VECTOR_DESCENT_20261104_v63` は、中心 μ2² を
除いた二つの squared tautological line の直和が genuine Aut(Q)-linearized
rank2 vector bundle であることを使い、元の不変式体 F に対し

```
F=C(W16)^Sp4 = B0(s1,s2)
```

を証明している。この二変数を取り落として F=B0 とはしない。
Proposition11 の有理性と合わせると

```
C(W16)^Sp4 ≅ C(t1,t2,t3,t4,t5,t6).
```

従って W16 の Sp4 商は **C 上有理、次元6** である。
これは既存の stable rationality の上界だけより強い結論である。

## 結論しないこと

この適用は特定の六つの純超越生成元や、その有理的な逆写像を今回構成したものではない。
有限13回路表示から六つへ明示的に整理する作業も未実施である。
多項式不変式環が polynomial、全 Hilbert 級数、最小生成次数は結論しない。

NEW25 がこの商で定める divisor の具体式、その quotient/source の指数、次数24の式の存在、
残余曲線の正確な次数・種数、Poisson 成分の網羅や元の affine fiber の尽くし込みは別課題。
今回それらを解決したとは扱わない。

原 PDF 全体や画像をリポジトリへ複写しておらず、引用は書誌と命題・頁の参照のみ。
新しい β 点、高次数 probe、有限群不変式探索、GB は行っていない。

## 書誌記録と再現

```
@article{ShepherdBarron1989S5,
  author = {Shepherd-Barron, N. I.},
  title = {Invariant theory for S5 and the rationality of M6},
  journal = {Compositio Mathematica},
  volume = {70}, number = {1}, pages = {13--25}, year = {1989},
  url = {https://www.numdam.org/article/CM_1989__70_1_13_0.pdf}
}
```

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_firstjet_shepherd_barron_rationality_20261104_v63.py
```

primary-source 独立受理と既存の全 descent 独立受理を gate にする。
命題番号・頁・表現の一致と結論の範囲を同名 report に固定する。
