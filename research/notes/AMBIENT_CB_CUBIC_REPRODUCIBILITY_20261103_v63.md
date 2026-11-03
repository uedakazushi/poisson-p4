# CB cubic workerの入力cacheを公開証拠から復元する

2026-11-03。完成したworkerのbyte同一性を保つため、過去のsourceはignored work内のpoint-setup入力を読む。
公開済みgzipからその入力を戻す小さいhelperを追加する。凍結済み数学packetのcode/reportは変更しない。

新しいcheckoutで全物理QQ(b) workerを明示的に再実行する場合、まず

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/restore_CB_cubic_worker_inputs_20261103_v63.py

を実行する。既存cacheがあればseconds以外の全fieldが同じことを確認し、そのまま保持する。
異なる内容を上書きしない。続いてcubic packetを新しいdirectoryへコピーして、既存のpublisherへ`--recompute`を付ければよい。
generic setup workerが新しいgeneric inputを生成し、generic syzygy workerがそれを使う。

これは入力復元のみであり、新しいmathematical check、source setup全体のやり直し、GB runではない。
通常の短いcompleted-certificate監査には復元も不要である。
