# 現行分類の証明・入力・証明書を再現する

対象は複素P4の射影Poissonスキームの被約成分であり、PGL5商は取らない。
主稿は `papers/poisson_p4_classification_ja.pdf`、共通機構の別稿は
`papers/poisson_divisor_symmetries_ja.pdf`。全147 IDの定義は主稿の
「全成分の定義と一般族の一対一カタログ」にある。
必要な別稿の正規実因子認識91頁と一次還元118頁も、現行PDFとTeXを収録する。
現在の達成・保存状況は `research/classification/CURRENT.md` のみを参照する。

## 数学の証明と今回の検証の区別

主定理は各入口の全支持証明、実Poisson変形、全周囲の成分性、有限閉像の比較を合成する。
過去の全実験や全外部定理を今回独立に再証明したという主張ではない。
原証明を再利用した範囲、今回読んだ一次資料、今回実行した算術は
`reports/*integration_validation.json` と各共通機構の `*validation.json` に記録した。
96種の正確な凍結入力は、記録されたSHA256と一致するGit又は埋め込み原文から回収した。
旧producerの異なるbyte版は実際の証明書の入力であり、旧配布ZIPや不要なバックアップではない。
`frozen_producer_index.json` は元パス、hash、回収元、現行保存先を対応させる。

`reports/required_proof_dependency_closure.json` は、最終全支持と各成分の
primary claimからの依存、本文入力、正確な係数・コード・証明書を閉じたintegrity索引である。
これは新しいclaim台帳や数学の自動認証ではない。
`classification_roots.json` に選んだ根を記録する。
旧カタログ全体を写したmanifestの二次研究項目、実行用jobの更新、旧文章のhash、
実行時間・resume記録、明示的に未採択とされた誤ったV1 proofは数学の前提に数えない。
各除外と文章・制御情報の版差は索引に残す。
係数・実行コードのhash不一致はこの扱いで免除せず、正確な凍結入力へ解決している。

必須のlegacyコード・データは個別ファイルとして同梱する。
`legacy/v03/code/poisson_p4_components_check.py` 等を旧ZIPから取り出す操作は不要である。
`work/rank3-mixed-support-obstruction.json` という旧入力名は、同一SHA256の
`reports/CLASSIFICATION_RANK3_MIXED_COMPONENT_20261107.json` に正確に対応する。
`.json.gz` は一つのJSONデータ、`.in.gz` は一つのnative算術入力であり、旧bundleではない。

## 現行の検証

リポジトリのrootから実行する。Python依存は `requirements.txt`。
Python cacheと生成物は `work/` に置く。

```bash
export PYTHONPYCACHEPREFIX="$PWD/work/python-cache"
python code/audit_required_dependencies.py --output work/dependency-closure.json
python code/check_global_structure.py --output work/global-structure.json
python code/check_marked_reconstruction.py --output work/marked-reconstruction.json
python code/check_quadric_residue_reconstruction.py --output work/quadric-residue.json
python code/check_conductor_conormal_reconstruction.py --output work/conductor-conormal.json
python code/check_normalization_residue_reconstruction.py --output work/normalization-residue.json
python code/check_rank_two_reconstruction.py --output work/rank-two.json
python code/check_component_catalogue.py --output work/component-catalogue.json
bash scripts/build_structural_manuscripts.sh
```

51族の標本評価は有理演算の復号と座標規約の検査である。
全パラメータJacobiと稠密性は主稿及び指定した一般源の証明による。
LuaLaTeXは二回組版し、box・未解決参照・重複参照を検査する。
既存のfont/scaling/hyperref警告は残っており、無警告との主張ではない。
TeX cache・インストール済み第三者ライブラリはbundleの対象外である。

## 凍結された一つの証明書を隔離して再実行する

次は今回実際に行った例である。作成先は存在しない新しいディレクトリを指定する。
原稿・原証明書・現行コードは上書きしない。

```bash
python scripts/materialize_frozen_proof_inputs.py \
  --report reports/CLASSIFICATION_RANK3_MIXED_COMPONENT_REVIEW_20261107.json \
  --destination work/replay-mixed
python work/replay-mixed/code/classification_rank3_mixed_component_review.py \
  --output "$PWD/work/replay-mixed/work/recomputed-review.json"
```

この再実行では全350 Jacobi、正確な周囲行列とminor、余核functional、二次障害−140、
実因子を確認した。六つのrank二成分の元行列も別に再実行し、旧証明書との一致を確認した。
これらの範囲を全147の独立検証へ拡張していない。
他のproducerは、対応するjobのCLI・環・反転元を確認し、同じ方法で入力を隔離する。
native計算を含むものにはMacaulay2/Singular/FLINT等、そのjobに明記された依存が必要である。

旧weight6独立reviewのfresh実装比較用scratchは、必要時に次で生成できる。
この一時出力を恒久的な未同梱入力として要求しない。

```bash
python code/q_global_weight6_two_axes_20261103_v63.py \
  --output work/old-plane-v62/weight6-review-fresh
python code/q_weight6_two_axes_review_20261103_v63.py --output work/weight6-review
```

## Library配布一式

Libraryの全ZIP volumeを同じフォルダに置く。
いずれかから `UNPACK.py` を取り出し、全volumeを指定する。

```bash
python UNPACK.py --parts poisson-p4-v72-part*.zip --destination poisson-p4-current
```

各volumeは同じmanifestを持ち、内容は重複しないcontent-addressed chunkである。
全chunkと復元された全ファイルをSHA256で検証する。
必要なディレクトリ参照リンクも同じ相対参照で復元し、論理パスと実体パスの一致を確認する。
旧ZIPの入れ子、ログ、cache、旧PDFは含めない。
論文PDFはLibraryに単独でも保存する。生成、CRC/hash確認、main保存、Library保存は別のreceiptで記録する。
