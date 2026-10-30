# 一次文献と、この研究で用いる範囲

**v25-setup — 2026-10-30。今回のセットアップでは新しい文献調査・新規定理の主張を行っていない。**

この記録は v24 までの原稿と、2026-10-30 の評価で確認した出典の適用範囲をまとめたものである。文献の網羅的な新規性調査は未完了であり、「現在知られている全結果」「本研究が最初」といった優先権の結論は含まない。保存原稿の引用と、一次本文の独立確認を区別する。

## 評価で確認した中心文献

### CM：degree two の foliations

M. Corrêa and A. Muniz, *Holomorphic foliations of degree two and arbitrary dimension*, arXiv:2207.12880v5, 19 January 2026.

- 一次本文：[arXiv HTML v5](https://arxiv.org/html/2207.12880v5)。
- 本研究で使う箇所：本文 Theorem 4.1 の分岐分類、関連する低次数 foliation の記述。
- 適用範囲：v14 の非安定分岐の帰着を支える。stable branch についての文献の結論は、代数的可積分性と uniruled leaves の記述までである。
- 推論してはいけないこと：文献だけから本研究の rank2 の8成分・一般 rank2 の6成分の完全性を引用すること。v14 の安定 reflexive sheaf の議論は別に必要である。
- 引用上の注意：2026-10-30 の監査では導入と本文の記述に不整合、一般指数の印刷上のずれを認めた。実際に用いる本文の文言と仮定を保存して引用する。
- 内部証拠：[大域監査](../audits/2026-10-30/global_audit.txt)、[v14 rank2_global_eight.tex](../legacy/v14/sections/rank2_global_eight.tex)、[stable_rank2.tex](../legacy/v14/sections/stable_rank2.tex)。

### MPS-Torelli：SNC locus と交わる40成分

M. Matviichuk, B. Pym and T. Schedler, *A local Torelli theorem for log symplectic manifolds*, arXiv:2010.08692v2 (2021).

- 一次資料：[arXiv v2](https://arxiv.org/abs/2010.08692v2)、[PDF](https://arxiv.org/pdf/2010.08692v2)。
- 本研究で使う箇所：Theorem 1.4。v23 の同定では Table 2、Example 6.24 も参照している。
- 適用範囲：40という個数は、SNC log symplectic locus と交わる成分についての結果である。
- 推論してはいけないこと：全 Poisson 空間の成分数が40であること、すべての一般 rank4 成分がその40個に含まれること。
- 研究上の位置：保存草稿の追加11成分との比較の基準。追加成分の一般点が非 SNC であることだけでは、40成分との相違の証明にならない。閉包全体に及ぶ識別条件を使う。
- 新規性に関する留保：定理直後の「すべての log symplectic 構造がこれらの成分に属するか」という問いとの関係は公表候補Bの動機になる。ただし、以後の文献を網羅的に確認したとは主張しない。
- 内部証拠：[大域監査](../audits/2026-10-30/global_audit.txt)、[v8 count_and_frontier.tex](../legacy/v08/sections/count_and_frontier.tex)、[v23 ss_literature.tex](../legacy/v23/sections/ss_literature.tex)。

### Pym-Elliptic：purely elliptic な場合の分類

B. Pym, *Elliptic singularities on log symplectic manifolds and Feigin–Odesskii Poisson brackets*, Compositio Mathematica **153** (2017), 717–744, arXiv:1507.05668.

- 一次資料：[arXiv](https://arxiv.org/abs/1507.05668)、[DOI](https://doi.org/10.1112/S0010437X16008174)。
- 本研究で使う箇所：purely elliptic な $\mathbf P^4$ の分類、および v22 が引用する Theorem 4.5 と Table 1 の局所標準形。
- 適用範囲：modular residue の非消滅・横断性に関する仮定を満たす elliptic な領域。
- 推論してはいけないこと：「実際の退化因子が正規」という仮定だけで $q_{5,1}$ 型になることを、この文献の定理として引用すること。その強化は v14–v22 の追加議論が担う。
- 内部証拠：[幾何・方針監査](../audits/2026-10-30/geometry_strategy.txt)、[v22 setup.tex](../legacy/v22/sections/setup.tex)、[v20 inactive.tex](../legacy/v20/sections/inactive.tex)。

### Pym-Survey：分類問題の背景

B. Pym, *Constructions and classifications of projective Poisson varieties*, Letters in Mathematical Physics **108** (2018), 573–632, arXiv:1701.08852.

- 一次資料：[arXiv](https://arxiv.org/abs/1701.08852)、[公開本文](https://pmc.ncbi.nlm.nih.gov/articles/PMC5818673/)。
- 本研究で使う範囲：射影 Poisson 構造、退化因子、log symplectic 構造、既約成分による分類という問題設定の背景。
- Theorem 5.12 の条件：退化因子の特異 locus の余次元3、および Poisson tensor の零 locus 上での modular residue の非消滅。
- 推論してはいけないこと：この総説の定理を正規因子だけの分類や、2026年時点での完全な文献状況と読み替えること。
- 内部証拠：[幾何・方針監査](../audits/2026-10-30/geometry_strategy.txt)。

## 保存済みの全係数同定で使う文献

### MPS-Quantum：$T_{37}$ の同定

M. Matviichuk, B. Pym and T. Schedler, *Creating quantum projective spaces by deforming q-symmetric algebras*, arXiv:2411.10425v1 (2024).

- 一次資料：[arXiv HTML v1](https://arxiv.org/html/2411.10425v1)。
- 保存原稿で使う箇所：Example 5.6 の明示係数と Theorem 7.5 の Poisson 変形。
- 本研究での役割：v23 に収録された全係数比較により、SS と NN の境界を同じ既知24次元成分へ同定する。
- 範囲の注意：$T_{37}$ は本研究の内部番号であり、この論文の定理番号ではない。SS と NN を追加の二成分として数えない。
- 確認状態：係数同定の根拠は [v23 ss_literature.tex](../legacy/v23/sections/ss_literature.tex) と同版の `proofs/literature/` に保存されている。今回のセットアップで文献全体を改めて査読したわけではない。

## 原稿内の引用として保存し、統合稿で再確認する資料

次表の書誌情報と用途は保存原稿から引き継ぐ。2026-10-30 の評価で、以下の全本文・全適用条件を独立に再査読したという意味ではない。論文A・B・Cへ取り込む際は、指定する版・節・定理の原文を確認する。

| 引用 | 一次資料 | 保存原稿が指定する用途 | 原稿側の入口 |
| --- | --- | --- | --- |
| D. A. Buchsbaum and D. Eisenbud, *Algebra structures for finite free resolutions, and some structure theorems for ideals of codimension 3*, Amer. J. Math. **99** (1977), 447–485 | [DOI](https://doi.org/10.2307/2373926) | 高さ3の Pfaffian イデアルと Gorenstein 分解 | [v14 normal_quintics.tex](../legacy/v14/sections/normal_quintics.tex) |
| M. Gualtieri and B. Pym, *Poisson modules and degeneracy loci*, arXiv:1203.4293 | [arXiv](https://arxiv.org/abs/1203.4293) | 退化 locus と modular residue。v14 の参考文献は Theorem 7.9 を指定 | [v14 主原稿の文献一覧](../legacy/v14/p4_poisson_completeness_v14.tex) |
| T. A. Fisher, *Invariants for the elliptic normal quintic* | [著者公開PDF](https://www.dpmms.cam.ac.uk/~taf1000/papers/invenq.pdf) | §10、式(38)–(39) の明示例と v14 の厳密計算 | [v14 normal_quintics.tex](../legacy/v14/sections/normal_quintics.tex) |
| F. Loray, J. V. Pereira and F. Touzet, *Foliations with trivial canonical bundle on Fano 3-folds*, Math. Nachr. **286** (2013), 921–940 | [arXiv](https://arxiv.org/abs/1112.3871) | v14 の文献一覧は Theorem 6.2 を指定。別版での用途と混同しない | [v14 主原稿の文献一覧](../legacy/v14/p4_poisson_completeness_v14.tex) |
| J. V. Pereira and J. P. dos Santos, *Distributions with locally free tangent sheaf*, arXiv:2408.10057 | [arXiv](https://arxiv.org/abs/2408.10057) | locally free tangent sheaf に関する背景。reflexive な場合への無条件拡張はしない | [v14 主原稿の文献一覧](../legacy/v14/p4_poisson_completeness_v14.tex) |
| The Stacks Project, Tag 031S | [Tag 031S](https://stacks.math.columbia.edu/tag/031S) | 正規性の $R_1+S_2$ 判定 | [v22 主原稿の文献一覧](../legacy/v22/p4_poisson_completeness_v22.tex) |

これは保存原稿の全引用を再掲した表ではない。v17 の形式的線形化、v20 の割線幾何、その他の局所計算の引用は、それぞれの原稿の参考文献一覧から正確な書誌と適用箇所を取り込む。未確認の記憶で新しい定理番号や書誌情報を補わない。

## 次の引用監査で満たすこと

1. 各主定理について「既知の入力」「本研究の追加議論」「未監査の依存先」を分ける。
2. arXiv の版を固定して定理の仮定と結論を記録し、出版版で番号が異なる場合は対応を付ける。
3. rank2 の8/6尽くし、追加11成分、正規因子分類についてそれぞれ新規性を調べる。現状の参考文献表だけで完了扱いにしない。
4. 新規性の評価と数学的正しさの監査を分ける。既知結果と分かっても証明が誤りになるわけではなく、証明書が再現しても新規性が保証されるわけではない。
