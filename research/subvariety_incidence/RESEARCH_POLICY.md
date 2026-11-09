<!-- この研究固有の方針はルート RESEARCH_POLICY.md のGitHub正本・Library緊急避難原則に従う。 -->
# Poisson不変部分多様体とambient完全分類

2026-11-10 / v1。現在の研究モジュール。最初に本書を読み、`papers/poisson_subvariety_incidence_ja.tex` と `research/subvariety_incidence/NOTES.md` を参照する。

## 研究の位置付けと目標

`poisson-p4` の147成分完全分類を改めて数え直すのではなく、滑らかな超平面、二次・三次超曲面、4次del Pezzo曲面上のintrinsic Poisson構造と、周囲P4のPoisson構造との延長incidenceを研究し、rank2/rank4成分の幾何・次元・退化因子の特異性を説明する。既存の418頁の分類主稿は置換しない。新規性の文献調査は未完了。

## 今回の主結果と証明状態

- **滑らかな三次超曲面**：非零intrinsic Poisson構造のambient延長は一意でJacobian rank2型。差はH0(Ω²(2))=0を用いる。したがってrank4の実退化因子に滑らかな三次既約因子はない。既知のintrinsic分類（LPT/Mayanskiy）を使用。ambient一意性は既存研究にもあり、新規性は主張しない。
- **滑らかな次数4以上の超曲面**：非零ambient Poisson構造が不変に保つことはない。本文の消滅・制限議論を参照。
- **一般の滑らかなCI(2,2)曲面と非零Poisson切断の組**：全ambient延長スキームはJacobian rank2の被約一点。全20変数源、150係数式の二次部分の普遍Q(u,v)上の行関係と、非零minor 70875/88、4443337729552/7688671875を用いる。これは一般開の定理であり、全てのsmooth pairの定理ではない。検算の実行記録は同梱bundleのcertificatesにある。
- **rational成分の次元**：滑らかな三次＋2-planeのincidenceは40次元で一般ファイバーA4、36次元のRat(1,1,3)を説明。CI(2,2)＋一次切断のincidenceは30次元で一般ファイバーA2、28次元のRat(1,2,2)を説明。
- **単純な不変超平面**：横断的な大域Poisson対称性が存在することがrank4延長の必要十分条件。重複超平面への不存在は主張しない。
- **滑らかな二次超曲面**：既存研究から、intrinsic Aff型の二点延長と、同一restrictionのrank2/rank4 pencilを統合。新しい成分数は主張しない。

## 未解決の核心と次の具体的操作

1. **4次del Pezzoの真の延長例外locus**。一般開でのminor消失は計算pivotの失敗であり真の例外とは限らない。特別な反標準切断hと五退化memberを持つquadric pencilを層別し、全20変数Jacobi理想の他のpivotで覆えるか、実際にrank4解があるかを判定する。真の解はPfaffianと既存147成分への帰属を調べる。
2. Q3の10変数延長空間をintrinsic対称性の跳躍で層別し、rank2/rank4の同時延長の一般性を説明する。
3. P3の六intrinsic成分で横断的Poisson対称性を計算し、単純Hと重複Hを区別する。
4. 必要な一次文献を本文と照合し、新規性が未確認な主張をIntroductionで明示する。

## 再現・保存

LuaLaTeX原稿：`papers/poisson_subvariety_incidence_ja.tex`。
独立コード：`code/subvariety_incidence/verify_del_pezzo_rigidity.py`、`code/subvariety_incidence/check_incidence_examples.py`。
Python3 / SymPy 1.14.0。コードは単体実行可能で、結果の正確なJSONはGitHubのresearch/subvariety_incidence/certificates/に収録し、手動スクリプトで再生成・照合する。
原稿・研究方針・ノート・検証スクリプトと正確な証明書はGitHub mainで揃い、`bash scripts/verify_github_source.sh incidence` が11頁PDFを再生成する。Libraryの旧配布物は不要。

## 一般原則と保存許可

研究の重要性・正確性・新規性・実現可能性に応じて方針を変更し、変更理由と未解決の核心を現行方針に統合する。証明済み、条件付き、計算観察、未証明、既知と新規性未確認を明確に区別する。重要な先行研究は一次資料で確認し、文献情報を推測で補わない。論文は日本語LuaLaTeXで十分な証明を保持し、ノートは作業履歴ではなく現行知識と障害を整理する。研究再開は現行方針・ノート・論文だけで可能にする。**現在のユーザー指示によりGitHub mainが唯一の正本であり、Library保存はGitHubへの実push失敗時の緊急避難に限る。GitHub Actionsは使わない。** bundleには現行の必要ソース・PDF・コード・入力・証明書・再現手順のみを収め、旧版・ログ・cache・nested ZIP・重複を含めない。容量の固定上限より検証可能性を優先する。生成・検証・GitHub mainへの非force pushと先端読戻しを別々に報告し、既存Libraryファイルは明示指示なしに削除しない。
