# 型7の残る五全源：実曲線・既存99写像・新23成分の独立監査

2026-11-07。PRIMARY E。**全項目 PASS**。
第一・第二・第七全源は既存 MPS 閉和集合、第六全源は既存99番目
`R4_GA9_DISTINCT_IRREDUCIBLE_QUINTIC_23` に含まれる。
第五全源は、新しい23次元の全 Poisson 成分に含まれ、その成分は prior99 と相異なる。
第三・第四源の二つの22次元成分の局所支持は別作者の対象であり、今回の再計算に含めない。
全体の番号・下界・台帳は統合担当が更新する。

## 1. 入力と独立性

独立 checker は `code/classification_nil22_long_curves_review.py`。
producer、SymPy、native CAS、private work/IFT/jet を import 又は実行していない。
既存の独立 Fraction 疎多項式 helper から全150係数と全350式を構成した。
公開 `CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_CURVE_RECOGNITION_20261107.json`
の literal curve が主入力。公開 ORBIT_MODELS のうち使うのは actual anchor、全源とその座標だけであり、
同じ report に収録された第三・第四の形式切片は今回の証明の依存でない。
全入力 hash と checker 全文は結果に保存する。

- 契約：`reports/CLASSIFICATION_NIL22_LONG_CURVES_REVIEW_JOB_20261107.json`
- 結果：`reports/CLASSIFICATION_NIL22_LONG_CURVES_REVIEW_20261107.json`
- gate PASS 後の全実行：2.271秒、最大 RSS 23,616KiB。

全源有限化には、独立再実行済み
`CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_REVIEW_V2_20261107.json`
の全七源・全零補集合・飽和次元≤22を依存として用いる。
prior98 の一般因子監査には独立確認済み
`CLASSIFICATION_GA9_DISTINCT_COMPONENT_REVIEW_20261107.json` を依存として用い、
その二十九の全 MPS 族因数恒等式を理由なく再実行してはいない。

## 2. actual curve と全源22軌道閉包

五曲線の全350 Jacobi は QQ[t] 上で零。
各 t=0 の全六括弧を、元の公開 nilpotent 表 row6 の全20係数へ canonical source 値を代入したものと
独立に一致させ、ORBIT_MODELS の integer scale も一致させた。
全 source parameter を零と見なして生成した別の点ではない。

各 anchor の全 GL5・Euler vertical・scalar 行列を組み、vertical rank24、全 rank47 と
保存47小行列式を再生した。従って actual projective orbit は22次元。
各正規源は前の全源定理の A4 又は A5 の主開部分で既約であり、その PGL 飽和も既約。
既証明の次元≤22と、この実軌道22を含むことから、全飽和閉包は anchor 軌道閉包そのものである。
この橋によって一点から全 parameter を覆う。選んだ切片だけから全支持を推論しない。

曲線の全 monomial に対し、重み σ と parameter 重み μ について

\[
 \operatorname{diag}(\lambda^\sigma)_*P(t)
 =\lambda P(t/\lambda^\mu)
\]

を検査した。μ は第一・第五・第六で1、第二・第七で2。
複素数体上で λ^μ=t を選べるので、全 t≠0 は実 projective orbit 一つとなる。
稠密な曲線部分の入る PGL 安定閉集合は t=0、さらにその全22軌道閉包を含む。

## 3. 第一・第二・第七の MPS 被覆

全曲線で h-Casimir を検査し、斉次発散ベクトルと全5×5 cofactor 行列から
特性多項式を直接 determinant 展開した。

|源|特性多項式|
|---|---|
|第一|u(t−u)(t+u)(2t−u)(2t+u)|
|第二|u(t−u²)(4t−u²)|
|第七|u(4t+u²)(16t+u²)|

t≠0 では零と四つの非零相異固有値、非零比は −2,−1,1,2。
h-Casimir を保つので、非零四次元 affine 部分を対角化し平行移動を除いて既済 semisimple 型12へ入る。
実 Pfaffian の t=0 は前の全源定理で既約と証明された quartic そのもの。
最高四次も非零であり、その斉次化も既約。四次 reducible locus は閉なので一般 t でも既約。
従って既済 `HYPERPLANE_QUARTIC_SEMISIMPLE_SUPPORT.md` の**全型12**支持定理が適用される。
有限個の既存 MPS 成分の和は閉かつ PGL 安定なので、曲線原点と第一・第二・第七の全源を含む。
特定の MPS 番号や一意性を仮定する必要はない。

## 4. 第五・第六の23次元全成分性

両者で構造的代表 t=1 を取り、primitive integer scale を全150係数で一致させた。
全150変数・全350式から直接組んだ D と、全 GL5・vertical・scalar 行列 M は

\[
 \operatorname{rank}D=102,\quad \operatorname{rank}M=48,\quad DM=0.
\]

保存された102・48角小行列式を正確な Fraction 消去で一致させた。
vertical rank は24であり、射影 actual orbit と全接空間はいずれも23次元。
raw Jacobi 空間でも実軌道/gauge/scalar の48次元下界と接次元48が一致するので、点は regular。
従って点の成分は唯一で滑らか、その23軌道閉包が全 Poisson 空間の既約成分になる。
これは IFT の完全切片・高次 jet を要求しない直接の成分性判定である。
前節の曲線と全22軌道閉包の橋により、対応する第五・第六の全源がその成分に含まれる。

## 5. 第六の既存99への全150写像

二次同次 tensor の canonical 代表を

\[
 P^{\rm can}=P-\frac15 E\wedge\operatorname{div}P
\]

と取った。発散の trace は零なので、この補正は発散を正確に零とし、同じ projective tensor を表す。
この操作は GL5-equivariant。全150係数で producer の代表と一致した。

公開 report の有理可逆 eigenbasis B と
`D0=diag(1,1,−1,1,−1/48)` に対し T=D0 B⁻¹ とする。
既存99の独立受理済み九括弧 tensor を K と書くと、直接

\[
 T_*P_6^{\rm can}=-\frac{125}{2}K
\]

を全150係数で確認した。さらに元の P6 についても

\[
 T_*P_6+\frac{125}{2}K
 =\frac15 E\wedge\operatorname{div}(T_*P_6)
\]

という literal Euler/scalar identity を全150係数で確認した。
因子 D だけの一致又は非零固有値比だけの認識ではない。
従って第六の非零曲線点は既存99の dense23軌道に入り、第六の全源閉包もその既存成分へ入る。

## 6. 第五の複素絶対既約 quintic と非正規性

canonical 代表に公開有理 eigenbasis を適用して全150係数を一致させ、primitive tensor が次の九括弧であることを確認した。

```
{h,y}=−6h²,       {h,z}=3hx,
{h,w}=−3hy+162hz−45xy,
{x,y}=9hx,        {x,z}=−2x²,
{x,w}=−5hz+2xy−243xz,
{y,z}=−9hz+3xy,
{y,w}=6hw−3y²+567yz,
{z,w}=2xw+7yz−243z².
```

省略した {h,x} は零。全 Jacobi は零、signed Euler-Pfaffian は厳密に75F、

\[
 F=h^3z^2+h^2wx^2+h^2xyz+27hx^2yz+3x^3y^2.
\]

w について一次、係数 h²x² に対し余項の h=0 制限は3x³y²、x=0 制限は h³z²。
従って C[h,x,y,z] 上で係数と余項の gcd=1。
Gauss により F は C 上で既約であり、QQ factorization からの推測ではない。

∂wF=h²x² なので singular support は h=0 又は x=0。
h=0,x≠0 では ∂yF=6x³y から y=0。
x=0,h≠0 では ∂zF=2h³z から z=0。
逆包含は全五偏微分の代入で確認した。従って

\[
 \operatorname{Sing}(F)_{\rm red}
 =\{h=x=0\}\cup\{h=y=0\}\cup\{x=z=0\}.
\]

三つの相異なる射影平面であり、既約三次元 hypersurface に codimension1 の特異支持があるので非正規。

## 7. 第五成分と全 prior99 の相違

第五点の canonical 代表について、全25 GL matrix と tensor scalar の150×26 stabilizer 行列を独立に作った。
rank24、保存した完全 kernel 二列の積は全150行で零、二列も独立である。
projective scalar matrix を除く stabilizer は正確に一次元。
公開生成元を eigenbasis で対角化すると固有値は
`(−512,−512,128,128,768)`、multiplicity は2,2,1である。

既存99の九括弧 canonical model に対しても同じ全150×26行列を作り、rank24、projective stabilizer 一次元を再生した。
その生成元は diagonal `(−50,−25,0,25,50)`、tensor character50であることを全括弧で検査した。
両成分は dense actual23軌道を持つ。もし同じ成分なら二つの稠密局所閉軌道はともに開で交わり、同じ軌道になる。
しかし唯一の projective stabilizer line の半単純生成元は共役・非零倍・scalar identity shift の下で
固有値 multiplicity を保ち、2,2,1 と1,1,1,1,1は一致しない。従って別成分。
独立の divisor 不変量としても、既存99の被約特異支持二平面に対し第五は三平面である。

残る prior98 は既済の独立全一般源監査をそのまま適用する。
六 rank2 項目は閉 rank 条件、三十六項目は全成分次元の差で区別される。
残る同次元 rank4 は二十九の MPS 成分と二十七の nonMPS 成分で、
各**全一般源**の実 quintic は因数分解又は非被約であることが受理済み。
次数五の reducible/zero cone は閉なので、その成分全体も同閉条件を満たす。
第五点の複素既約 F はその条件に反する。
これにより第五は全 prior99 と相異なる。
同時に審査される第三・第四の22次元成分とは23対22という全成分次元の差で区別される。

全 parameter の帰属、成分性、複素既約性、distinctness をここまで別々の矢印として確認した。
新成分の最終番号・claim ID・全体計数は root の統合に委ね、共有台帳・Git・原稿は変更していない。
