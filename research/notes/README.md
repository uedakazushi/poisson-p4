# 研究ノートと評価資料の入口

**v25-setup、2026-10-30。新しい数学的定理を追加せず、v24までの研究を再開できるよう整理した状態。**

まずルートの [README](../../README.md) と実行方針を読む。このディレクトリの評価は、原稿・計算の現在地と継続方針を判断する資料であり、全先行版の独立査読や完全分類の証明に代わるものではない。

## 現状評価

| 資料 | 内容 |
| --- | --- |
| [p4_poisson_status_strategy_ja_20261030_v1.pdf](p4_poisson_status_strategy_ja_20261030_v1.pdf) | v24時点の成果、監査結果、未解決範囲、継続方針の評価 |
| [p4_poisson_status_strategy_ja_20261030_v1.tex](p4_poisson_status_strategy_ja_20261030_v1.tex) | 同評価文書の LuaLaTeX 原稿 |

この評価と v24 の `HANDOFF.txt` には役割の違いがある。v24 は最新の端点計算を詳しく記録しているが、全一般 rank4 成分を捉える v2 の $S_{61}$ は古い版に保存されている。研究全体の方針を v24 の局所課題だけから決めない。

## 2026-10-30 の監査記録

| 記録 | 確認した範囲 |
| --- | --- |
| [global_audit.txt](../../audits/2026-10-30/global_audit.txt) | rank2 の8成分・一般 rank2 の6成分、追加11成分の出典、下界57、$S_{61}$/$S_{78}$ の位置 |
| [local_audit.txt](../../audits/2026-10-30/local_audit.txt) | v24 の実際の成分所属、26次元共通族、接空間と二次障害、局所研究の限界 |
| [geometry_strategy.txt](../../audits/2026-10-30/geometry_strategy.txt) | 正規因子定理の入力鎖、五係数軌道族の意味、独立監査の優先箇所と大域戦略 |
| [v24_recheck.json](../../audits/2026-10-30/v24_recheck.json) | v24 の1298検証項目の再現結果 |
| [v14 verification_summary.json](../../audits/2026-10-30/v14_recheck/verification_summary.json) | v14 の10検証群の再現結果 |

監査メモ中の絶対 scratch パスや `review_work/` は当時の実行記録である。このリポジトリでの現在の入口は上のリンクを使う。チェックが PASS した範囲と、人間が確認すべき幾何学的網羅性・論理依存関係を区別する。

`local_audit.txt` にある重みを利用した有限 Kuranishi モデルの提案は、既存データからの継続方針である。v24で全モデルを構成・分解済みと読み替えない。`geometry_strategy.txt` の Borel 点の再計算は、v2 の既知の普遍点の独立再確認であり、新定理として数えない。

## 原資料へ戻る入口

| 目的 | 原資料 |
| --- | --- |
| 全一般 rank4 の代数モデル $S_{61}$ | [v2主原稿](../../legacy/v02/p4_poisson_completeness_v2.tex)、[v2 HANDOFF](../../legacy/v02/HANDOFF.txt) |
| rank2 の閉じた分類部分 | [v14主原稿](../../legacy/v14/p4_poisson_completeness_v14.tex)、[rank2_global_eight.tex](../../legacy/v14/sections/rank2_global_eight.tex) |
| 追加成分を含む下界57の整理 | [v8 count_and_frontier.tex](../../legacy/v08/sections/count_and_frontier.tex)、[公表用原稿の索引](../../papers/README.md) |
| 正規な実際の退化因子 | [v22 setup.tex](../../legacy/v22/sections/setup.tex)、[v22主原稿](../../legacy/v22/p4_poisson_completeness_v22.tex) とそこに明記されたv14–v21の入力 |
| SS/NNと既知 $T_{37}$ の同定 | [v23 ss_literature.tex](../../legacy/v23/sections/ss_literature.tex) |
| 五重超平面端点の局所変形 | [v24 HANDOFF](../../legacy/v24/HANDOFF.txt)、[v24主原稿](../../legacy/v24/p4_poisson_completeness_v24.tex) |
| 一次文献の適用範囲 | [REFERENCES.md](../REFERENCES.md) |
| 版・ファイル・主張を機械可読で追う | [SOURCE_INDEX.json](../../papers/SOURCE_INDEX.json) |

保存用の版別原本は `sources/v02.zip` から `sources/v24.zip`。展開後の `legacy/` は読解・再現用であり、過去の主張に訂正が必要な場合は原本を変更せず、新しいノートと現在の研究状況へ訂正を記録する。旧版の未解決リストを最新状態より優先しない。

## 今後のノート

新しいノートには、日付、扱う問題、入力と仮定、得た結果、証明または証明書、残る問題を記す。確立した結果と探索中の案を分け、終了時には現在の研究方針へ反映する。日本語のまとまった数学原稿は LuaLaTeX を使う。定期的な履歴は版管理に任せ、旧版全体を新しいノートへ再帰的に複製しない。
