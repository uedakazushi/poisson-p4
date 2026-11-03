# 正規な実際の退化因子の認識論文：完成・再現記録

2026-11-03。独立した日本語 LuaLaTeX 論文を 91 ページの PDF として保存した。

- 論文ソース：[normal_actual_divisor_recognition_ja.tex](../papers/normal_actual_divisor_recognition_ja.tex)
- PDF：[normal_actual_divisor_recognition_ja.pdf](../papers/normal_actual_divisor_recognition_ja.pdf)
- ビルド手順：[build_normal_actual_divisor_paper_20261103.py](../code/build_normal_actual_divisor_paper_20261103.py)
- 機械可読の検証記録：[NORMAL_ACTUAL_DIVISOR_PAPER_20261103.json](NORMAL_ACTUAL_DIVISOR_PAPER_20261103.json)
- 独立検証：[査読ノート](../research/notes/NORMAL_ACTUAL_DIVISOR_PAPER_20261103_v62_REVIEW.md)、[検証結果](NORMAL_MODEL_RREF_20261103_v62_REVIEW.json)、[検証コード](../code/normal_model_rref_review_20261103_v62.py)

## 証明範囲

標数零の代数閉体上で、実際の退化因子が正規である Poisson 構造の閉点分類と、同じ認識条件の相対的なスキーム論的記述を定理・証明の順にまとめた。非活動成分、直線を含む可約台、有理正規四次曲線の 28 族、E6/E8 の幾何、八つの有理木の有限被覆を含む。旧 TeX への再帰的な入力はなく、係数順序、局所化条件、例外分岐、証明書の再現位置を本文内に明記した。

相対版で正規性を仮定する対象は幾何学的ファイバーである。冪零元を持つ基底上の全空間の正規性、非正規な実際の退化因子の分類、Poisson 成分全体の尽くし上げ、歴史的な独創性は主張していない。Pym および Gualtieri–Pym の定理は、正規性・非零 residue・局所積構造等の仮定を対応させて引用した。

## 検証

LuaHBTeX 1.18.0（既設 TeX Live）で 3 回コンパイルした。210 個のラベルと 12 個の引用が解決し、LaTeX/hyperref の警告、重複ラベル、文字欠落、overfull box はゼロ。PyMuPDF により 91 ページすべての文字領域、空白ページ、未解決記号を検査し、表紙・相対認識証明・係数表・再現手順を目視でも確認した。

独立検証では 211_triple、22_node、31_external の整数係数制約を別に作り、QQ 上の RREF により印刷した自由列番号、150 個すべての Poisson 係数、residue 全成分の一致を確認した。円錐曲線およびねじれ三次曲線の conormal frame も両チャートで検証した。長い既存の有限代数証明書は変更せず、五つの依存監査のハッシュと被覆確認を再利用した。全旧セルの Gröbner 計算を今回再実行したという意味ではない。

## 再ビルド

リポジトリの既設 Cloud 環境で、未使用の出力ディレクトリ名を選んで実行する。下のコマンドは保存済み PDF と検証記録を変更しない。

```sh
source /workspace/cloud-setup/activate.sh
python code/build_normal_actual_divisor_paper_20261103.py \
  --output-directory work/normal-paper-rebuild-fresh \
  --pdf-destination work/normal-paper-rebuilt.pdf \
  --report work/normal-paper-rebuilt.json
```

出力ディレクトリは存在していてはならない。ソースは単独でビルドできる。旧証明書を再現する場合に限り、論文末尾の手順と `sources/v14.zip` から `sources/v22.zip` の展開内容を使う。公開済み Cloud 環境のセットアップを再実行する必要はない。

## 保存ハッシュ

- TeX SHA256: `198a858a5b27e5e708a5358cf3f47a3d6377de523598b7facce73a8ead486079`
- PDF SHA256: `0fe03f289e18e2dfa8759fb9bd1c21ae50e4c2aeaa846225f989f81a1189ec85`
- ビルドスクリプト SHA256: `52a04eb2aa2a0b7286ea31942c0f8170f57c6891a8b3f0b31c059db9056ee105`

本記録では main 原稿、研究台帳、コミットは変更していない。main への統合・通常の push は親タスクが行う。
