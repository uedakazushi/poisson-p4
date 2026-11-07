# projected Veroneseを保存する全Poisson源：全像次元22以下

2026-11-07、主義務G。滑らかなprojected Veronese S⊂P4をPoisson不変にする
全global Poisson bivectorの固定S源はaffine次元≤2、全PGL射影像は次元≤22である。
S上でbivectorが零という条件へ制限しない。一般rank二のrestrictionを全て保持する。
この結果は未知W優越成分（必要次元25–27）をこの全支持から排除するが、
個別E成分帰属や全源でrank四が不可能という強い結論は主張しない。

## 1. 完全なtangent条件とEuler商

`q=[a²−b²,b²−c²,ab,ac,bc]` を使う。
q_a,q_b,q_cはcone tangent空間の基底でEuler qを含む。
従ってprojective πがSへtangentであることは正確に

`q∧Π(q)=r₃(a,b,c) q_a∧q_b∧q_c`

と同値である。右辺のr₃は全十係数のcubicであり、その部分族を選ばない。
次数六の十式の全28係数から、ambient全150係数とr₃の十係数に対する280×160行列を作った。
全行列rank121、完全kernel39、Πへの射影kernelは零。
Euler∧全linear場のrank24がこのkernelに入ることを確認してからquotientを取った。
canonical divergence-zero sourceは正確に15次元である。

全350 projective Jacobiは15parameterの二次式で、rowspace rank54。
54個の独立な元JacobiをQQの非零scaleだけでprimitive整数生成元へ取り直した。
これらは全350式と同じ理想を生成する。
実anticanonical quinticはEulerと全4×4 Pfaffianの正しい外積で再構成した。
その係数はJacobiの二次rowspaceには入らないので、単純なrow membershipによるrank四排除は使わない。

## 2. 全源次元の有限証明

上の全54生成元について、既設Singularの有限liftstdで242polynomialを得た。
全10,079のQQ polynomial liftを保存し、それぞれ
`G_j=Σ_i L_ij J_i` を厳密に検証した。
各G_jのgrevlex leading monomialが生成するmonomial idealは、元Jacobi理想のinitial idealに含まれる。
従ってこのmonomial idealの次元が全源の次元上界である。
十五変数におけるminimum vertex coverの大きさは13であり、

`dim 固定S全Poisson源 ≤15−13=2`。

この証明にGの完全GB性、primary decomposition、nil指数又は有限jetの積分可能性は不要である。
零Pfaffian、rank二、全ての退化parameterと全Jacobi supportを含む。
producer実行は63.985秒。最初のruntimeで標準libraryを見つけなかった事実は別に保存し、
公開環境に既存のlibrary pathを指定しただけで、再インストールは行っていない。
QQ generatorのdenominator clearingは非零定数scaleなので理想を変えない。

## 3. 全PGL像とrank二二重曲面への橋

全smooth projected Veroneseは一つのPGL₅軌道である。
Sを保つPGL群はP2のautomorphismでprojection centerの非退化conicを保つもの、identity成分SO₃の三次元群である。
SがP4をspanするのでambient作用のkernelは零で、逆向きにもこのSO₃作用は埋込みへ延びる。
従ってPGL₅×P(固定源)の像へのfiberには三次元群の自由作用があり、

`dim 全PGL射影像 ≤24+(2−1)−3=22`。

点でのstabilizer又はorbitrankを全源へ外挿していない。
`IRREDUCIBLE_NONNORMAL_RANK2_CLEAN_SURFACE.md` は、Sが滑らか、一般ordinary node、
blowup normal-cone交叉にvertical成分が無いという全rank二二重曲面枝を、この唯一軌道へ還元する。
正規化X・導手のsmoothnessは使わず、exceptional residue/Gysin、次数≤6、整数NSとdegree四の初等的識別を使う。
両者を合わせ、この全subbranchの未知優越成分が排除された。
特異S、非node、vertical成分の補集合と、他のrank零曲面は別に残る。

## 4. 別作者算術

`code/classification_projected_veronese_tangent_review.py` は、既存の独立Fraction engineだけを再利用し、
producer、CAS、Singularをimportしない。
全280×160行、非零121minor、kernel39、Euler24、canonical商15、全350式、
全54生成元、実quintic、全10,079 lifts、全242leading monomialとcover13を再構成した。
PASS、11.508秒、最大RSS59,680KiB。

* `CLASSIFICATION_PROJECTED_VERONESE_TANGENT_20261107.json`：全線形源・Jacobi・実D。
* `CLASSIFICATION_PROJECTED_VERONESE_TANGENT_DIMENSION_20261107.json`：全liftとleading次元証明。
* `CLASSIFICATION_PROJECTED_VERONESE_TANGENT_REVIEW_20261107.json`：全独立再構成、reviewer bytesも保存。
* 契約：`jobs/projected_veronese_tangent{,_dimension,_review}.json`。

旧rank零producer/certificatesは変更していない。
全既約非正規quinticの尽くし又は成分数上界をこの一枝から主張しない。
