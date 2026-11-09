# Poisson部分多様体の延長と147成分の幾何学的理解

2026-11-10 / v2。研究再開：ルート `RESEARCH_POLICY.md` → `REPRODUCE_GITHUB.md` → 本方針 → `NOTES.md` → `papers/poisson_subvariety_incidence_ja.tex`。

## 研究価値の最大化（上位原則）

数学的重要性、新規性、概念的理解、他の成分への波及、証明の正確性、実現可能性を総合して研究全体の価値を最大化する。現在の優先順位は固定のノルマではなく、新しい証拠・文献・障害によってより有望な課題が判明したら、元目標の解決点と残る問題を簡潔に整理して方針を更新する。進捗ログや旧方針の累積は不要。

## 現在の証明済み結果と限界

- **本文に一般証明**：滑らかなquadric `Q=(q=0)⊂P4` のstrict q-Casimir lift `P` に対し `β=i_(P∧P/2)Ω=f dq`（fは3次）、`D_P=2qf`、`v=div P∈so(q)`。一般rank4ならQ因子の重複度は1または2であり、`Q²|D_P ⇔ v|Q∧σ=0`。葉方向のmodular対称性が `2Q+H`、横断的対称性が `Q+C3` を区別する。さらにv=0ならrank≤2。主稿 `thm:quadric-modular-multiplicity`。
- **新しい系、有限計算と通常の証明の合成**：一般 `Rat(1,2)` 型ではstrict `so(q)` stabilizerが0、従って唯一のambient Poisson延長はJacobian型rank2。明示例の10×10 minor `27980515917168640` を正確計算で確認し、非空開へ拡張。特殊Ratの全点を分類した意味ではない。
- **既存結果の統合**：滑らかなcubicの延長一意性、一般CI(2,2) del Pezzoの被約一点延長、incidence次元36/28、単純Hの横断対称性、quadricのrank2/rank4の特殊例。仮定と証明は現行論文に保持。これらをすべて新規成果とは数えない。
- **文献**：Loray–Pereira–Touzet, Math. Nachr. 286 (2013) 921–940, DOI 10.1002/mana.201100354（Q3のintrinsic Rat/Log/Aff分類）；Mayanskiy, J. Algebra 424 (2015) 21–45, DOI 10.1016/j.jalgebra.2014.08.049。新定理の文献上の優先権は未確認。今回の正確な例示検算と一般証明・外部査読は区別する。

## 現在の研究順序

1. **滑らかなQ3の全延長ファイバー**。一般Ratは解決。次にLog(1,1,1)とAffの一般元・特殊元について `aut(Q,σ)⊂so(q)` を一次線形代数で求め、`P=s(σ)+K_T` のJacobiをその小さな対称性空間で解く。ファイバーの既約成分・非被約性、rank2/4、`Q+C3` / `2Q+H`、既存147成分への帰属を調べる。旧稿のLog/Affの一般ファイバー表は再検討対象であり、今回の共通定理だけで証明済みとは呼ばない。
2. **4次del Pezzoの真の例外集合**。一般延長剛性の消失minorと真のrank4延長を区別する。
3. **P3の6成分の横断Poisson対称性**を分類し、単純/重複Hを幾何学的に分ける。
4. それらの共通構造を147成分全体の幾何学的特徴づけに統合。
5. Feigin–Odesskii互換性、特異cubic/非正規因子、一般次元は価値が高まれば優先順位を再評価。

**次の具体的着手**：Log型一般元のstrict stabilizerと全二次Jacobiを計算し、modular重複度の定理に接続する。予想外の障害が出れば、その障害を説明する一般定理へ優先順位を変更してよい。

## 再現と保存の義務

唯一の正本はGitHub main。現行論文・研究ノート・方針・3つのSymPy producer・3つのexact JSON証明書・再生成手順を同一mainに置く。`bash scripts/verify_github_source.sh incidence`で手動検証。GitHub Actionsは使用しない。Libraryは実際のGitHub push失敗時の緊急避難のみ。既存Libraryファイルは削除しない。旧原稿、バックアップ、作業ログを現行bundleに積み重ねず、証明を軽量化のため削らない。
