# 公開済み Cloud 環境での研究再開

2026年11月3日（Asia/Tokyo）。開始時の main は
`4c3f623fb698f9adc6077015b25e2325b11072af`。通常の fetch 後の
`origin/main` と一致した。開始時の未追跡ファイルは、セットアップで生成された
`reports/20261102T142937.117200Z-all/` のみで、追跡済みファイルの変更はなかった。

## 実行環境の再利用

`/workspace/cloud-setup/RESULTS.md` と `activate.sh` を読み、各 shell で
`source /workspace/cloud-setup/activate.sh` を使用した。セットアップ・依存導入・
全基準 suite は再実行していない。ローカル移行時の停止記録は再開元の履歴として扱った。

今回の軽量確認は以下の通り。

- `scripts/doctor.py --strict`：Python 3.12.14、SymPy 1.14.0、NumPy 2.3.5、
  mpmath 1.3.0、23版の展開、日本語 LuaLaTeX の存在を確認。
- `/workspace/cloud-setup/native-smoke`：FLINT 3.1.0、GMP 6.3.0、MPFR 4.2.2 で成功。
- `M2 --script /workspace/cloud-setup/smoke.m2`：商・有理因数分解・有限拡大体の確認に成功。
- `sha256sum -c research/handoff/20261102/SHA256SUMS`：全項目が一致。
- GitHub と文献サイトの結果・実行権限は [接続記録](CLOUD_CONNECTIVITY_20261103.md)。

## Macaulay2 と native 計算の互換性

Cloud は Macaulay2 1.24.11、引継ぎ元は1.26.05。今回実際に使用したのは、
重み付き多項式環、`gb`、イデアル包含、剰余、交叉、Hilbert 級数、有限商の基底である。
非平方 `a=0` の中央8例について QQ と F32003 上で被約支持を検証し、
保存済みの全 Hilbert 級数と一致した。全 package や未使用 script の互換性は主張しない。
新しい runner は実行ファイルを PATH から選び、Homebrew の絶対パスを使用しない。
保存済みコードと原資料は変更していない。

重み付き環で `degree(R/I)` をベクトル空間長と同一視する新 runner の初期 assert は
失敗した。これは今回の実装の誤りであり、旧結果の不一致ではない。
当該有限商の `basis` は70列、Hilbert 級数の係数和も70である一方、`degree` は
234047636766720を返した。長さの検証を `numColumns basis(R/I)` に修正し、
成功した最終出力と区別した。初期ログは `work/cloud-a0-central-modp-20261103/5.log`、
最小確認は `work/cloud-resume-20261103/degree-compatibility.m2` に残した。

既設の native 実行ファイル `a-zero-det` でも保存入力の531行386列の行加群を再現した。
FLINT 3.1 系での行演算522662回、行列式次数17、全18係数は引継ぎ記録と一致した。
正式な再現 runner は同じ既設 toolchain から実行ファイルを再構築し、入力と出力を保存する。

## 研究記録

研究の再現入口は `code/s61_binary_a_zero_replay_20261103.py`。
中央支持、例外値、横断長の上界と、周囲の冪零性をつなぐ幾何学的な証明は
別に区別して記録する。下界57、一般ファイバー次元2以下は維持し、一般五点の尽くしは
新しい検証の PASS だけから結論しない。

ChatGPT Library への ZIP 保存を行う汎用ツールは、この実行環境では見つからなかった。
Pages の添付用ツールは存在するが、指定 Page なしでの Library ZIP 保存とは異なる。
Library 保存は未実行であり、通常の GitHub main 保存を使用する。
