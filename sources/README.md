# 原資料の保存形式

2026-10-30 / v25-setup。

ユーザー提供の `p4_poisson_completeness_v24.zip` に含まれる全通常ファイルを、
v02〜v24の各版の直接の内容に分けて `v02.zip` ... `v24.zip` に保存した。
本文・PDF・コード・証明書・旧ログも含め、ファイルのバイト列は変更していない。
ZIPコンテナの日時・並び・圧縮は正規化したため、元ZIPと同一のZIPではない。

`manifest.json` に元ZIPのSHA-256、各再梱包ZIPのSHA-256、全ファイルのサイズと
SHA-256、元の長いパスに戻すためのprefixを収録する。
元ZIPのハッシュは出所の照合用であり、再梱包から元ZIPの圧縮メタデータの復元は行わない。

```sh
python3 scripts/bootstrap_sources.py
python3 scripts/bootstrap_sources.py --check
```

展開後は `legacy/v02` ... `legacy/v24` を短いパスとして使える。
例えば `legacy/v24/previous_v23` は `../v23` への相対リンクとなり、旧コードの
再帰パス参照も維持する。原本のmanifestを検証する旧runnerは、symlinkを辿る一時コピー
で実行する必要があるため、通常は `scripts/verify.py --suite integrity` を使用する。

展開済みの原本が変更されていた場合、bootstrapは上書きせず停止する。
研究上の修正は `code/` 等へ別ファイルとして置き、変更元と理由を記録する。
元HANDOFFの優先課題は歴史的な状態であり、現在の優先順位はルートの研究方針に従う。
