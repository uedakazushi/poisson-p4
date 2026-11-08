# 概念整理の研究ノート

2026-11-09 / v1。最初に RESEARCH_POLICY.md を読む。

## 現在の数学

原稿 papers/poisson_divisor_symmetries_ja.tex は独立に読める。主要な証明をここに重複しない。

- thm:abelian：P^(2n)の一般に非退化なPoisson構造で、実退化因子の被約既約因子数がrなら、次数関係Σa_i d_i=0の閉対数一形式からr−1次元の可換な大域Poisson対称性代数を得る。非被約・非正規・非SNCも許す。各因子のPoisson性、留数、完備性を用いる。
- cor:orbit：P4ではprojective orbit dimension≤25−r。等号や成分性は別証明。
- thm:fibre：固定した重複度付き反標準因子F上で、V_F=ker(Fdα−dF∧α)を先に求め、α²/2∈C×Fを課すと一般rank4のPoisson構造と正確に対応する。非被約基底上の全scheme同値は未主張。
- prop:modular：アフィン四空間でZ=div P=P#dlog(Pf P)。二因子型のmodular場は標準対称性の生成元に比例する。
- §4.3：二重二次因子で逆形式がd(β/Q)に還元でき、dβが非退化なら、L_ZΩ=ΩとZQ=cQへ帰着する。全閉形式の有界極primitive表示は本稿では未証明。

通常の証明を収録し今回読解した。新CAS実験ではない。基本構造原理の新規性は未確認。

## 概念的な対応表

| 幾何 | 強制される構造 | 残る本質 |
|---|---|---|
| 5個の座標超平面 | 4次元トーラス、交代行列 | 同一因子上の射影5パラメータ、orbit20+5=25という族次元 |
| H+四次、2H+三次 | 1次元の非零線形Poisson対称性 | 半単純/冪零の表現、固定因子の閉形式、境界包含 |
| H1+H2+三次 | 2次元の可換線形対称性 | 軌道≤22、残る真正パラメータと配置 |
| H+2Q | 対称性と有界極二形式 | primitive表示の一様補題、共形シンプレクティック線形代数 |
| 既約五次 | 因子由来対称性は零 | 正規化、導手、葉・residue |
| 一般rank2 | Pfaffian写像の不定locus | 飽和foliation、零因子、rank4への変形 |

この表は新しい完全分類や排他的な成分分割を意味しない。

## 次に実際の重複を削る

H+2Qの線形核表示を、dQのKoszul複体と係数次数を使う一つの極次数還元補題で証明する。
de Rham類が零という事実だけではβの次数≤1や極の上限は出ないので区別する。

二因子H+Q4とQ2+C3で同じ標準対称性の半単純/冪零/共鳴表現を用い、既存の別々の入口を統合する。
固定作用の不変テンソル空間を取り、実divisorとPfaffian条件を戻す。
未分類の可換Lie部分代数全体の問題へ不用意に拡大しない。

各成分を (実D, 因子次数/重複度, a_pi, generic stabilizer, V_Fの残自由度, 変形型) で比較する。
r−1は大域ベクトル場の次元であって一般軌道の階数ではない。
既約五次とrank2は別の構造定理として扱う。

## 文献と既存証明への接続

Pym, Constructions and classifications of projective Poisson varieties, Letters in Mathematical Physics 108 (2018), 573–632, DOI 10.1007/s11005-017-0984-5, arXiv:1701.08852。§5.1、Propositions5.14/5.16。本文 https://pmc.ncbi.nlm.nih.gov/articles/PMC5818673/ を確認。
M. Matviichuk, B. Pym, T. Schedler, A local Torelli theorem for log symplectic manifolds, arXiv:2010.08692v2 (2021-07-15)、§3.1、Example4.14、Theorem1.4。本文 https://arxiv.org/html/2010.08692v2 を確認。40成分を全空間の尽くしとして引用しない。

既存の適用先：research/classification/DOUBLE_QUADRIC_SUPPORT.md、HYPERPLANE_QUARTIC_COFACTOR_ENTRANCE.md、GA9_CONSTANT_COMPONENT.md、T37_SIX_TARGET_LEAF_OBSTRUCTIONS_SEPARATE_REVIEW.md。
今回必要部分を読んだ。全caseの再検証は未実施。新基本定理の証明はこれらの旧ノートを前提にしない。
主稿の全体再構成は未完。旧証明を本ノートの短い要約で置換しない。

## 組版

```sh
mkdir -p work/conceptual-build
lualatex -interaction=nonstopmode -halt-on-error -output-directory=work/conceptual-build papers/poisson_divisor_symmetries_ja.tex
lualatex -interaction=nonstopmode -halt-on-error -output-directory=work/conceptual-build papers/poisson_divisor_symmetries_ja.tex
```

今回2回成功、7頁。参照・引用の未解決とbox警告は零。全頁の縮小描画を確認した。
luatexjaのフォントスケール・addjfontfeaturesの警告は残る。新PDFは配布/Library用で、Gitにはこの再生成可能なTeXを保存する。
TeX SHA256 a316ad4a9a810b1b1853f0dba9656d7f8a18c0b01e98f5278c5b8f07938b075a。
PDF SHA256 d71312953165fb82289328cdd2ea81632cbdf9fe54b8b5750e5ce241ab8f2886。
この組版確認を数学全体の独立認証と呼ばない。
